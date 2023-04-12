from folding import generate_fold_constraint, build_map
from utility import multicore, avg, arg_sort, inv_list, flatten
from vienna  import*
from fitness import evaluate_rz_activity, score_affinity, score_thermobalance
from representation    import Segment
import numpy as np

def coarse_grain_pheno(bppms, map, er_seg_ids, constant_rz_seg_ids, truth_vector, bool_TT, bool_constant_rz):
    if bool_TT:
        num_segs   = max( [v[0] for v in map.values()]) + 1 
    else: #TO
        num_segs  = len(er_seg_ids) + 1
        er_start_seg_id   = min(er_seg_ids)


    num_states    = len(bppms)
    t_sppms = np.zeros(( num_states, num_segs, num_segs ))
    for state_id in range(num_states):
        for key, val in bppms[state_id].items(): #!!!! dict.items() is super slow if bppm is full!!!
            if None not in key: #Ignore unpaired bases
                seg1_id = map[key[0]][0]
                seg2_id = map[key[1]][0]
                if bool_TT:
                    if truth_vector[state_id] == 1:
                        #When the target output is 1, only consider interactions between the extension region segments.
                        if seg1_id in er_seg_ids and seg2_id in er_seg_ids:
                            t_sppms[state_id][seg1_id][seg2_id] += val
                    else:
                        t_sppms[state_id][seg1_id][seg2_id] += val


                    #Nullify entries corresponding to interactions involving at least one constant rz segment
                    #Currently, this only works for TT, but not TO.
                    if bool_constant_rz:
                        if seg1_id in constant_rz_seg_ids or seg2_id in constant_rz_seg_ids:
                            t_sppms[state_id][seg1_id][seg2_id] = 0

                else:
                    if seg1_id in er_seg_ids or seg2_id in er_seg_ids:
                        #always consider ER-ER interactions regardless of target output
                        if seg1_id in er_seg_ids and seg2_id in er_seg_ids:
                            t_sppms[state_id][seg1_id - er_start_seg_id + 1][seg2_id - er_start_seg_id + 1] += val
                        #only consider ER-Rz and Rz-ER interactions if the target output is 0
                        if truth_vector[state_id] == 0:
                            if seg1_id in er_seg_ids and seg2_id not in er_seg_ids:
                                t_sppms[state_id][seg1_id - er_start_seg_id + 1][0] += val #Using [0] to store Rz seg vals. This is different from paper which uses [num_segs - 1]
                            if seg1_id not in er_seg_ids and seg2_id in er_seg_ids:
                                t_sppms[state_id][0][seg2_id - er_start_seg_id + 1] += val


    return t_sppms #This is symmetric since in Vienna, we make the bppm dict symmetric



def pheno_perf_fits(input_data): #Operates on a sub-population of individuals
    tasks, common_data = input_data
    base_to_seg, name_to_bases, er_seg_ids, constant_rz_seg_ids, rz_seg_names, truth_vector, num_inputs, bool300, folding_temp, bool_TT, bool_constant_rz, off_cap = common_data
    num_states         = 2**num_inputs
    num_indss_batch     = len(tasks)


    logic_tasks          = flatten([task[0] for task in tasks])
    affinity_tasks   = flatten([task[1] for task in tasks])
    #print("LOGIC", logic_tasks)
    #print("affinity", affinity_tasks)
    logic_results        = Fold(logic_tasks, bool300, folding_temp, 1)
    #affinity_results = Fold(affinity_tasks, bool300, folding_temp, 1)
    affinity_results = Fold(affinity_tasks, bool300, folding_temp, 2)

    results            = list()   
    for ind_idx_batch in range(num_indss_batch):
                
        potential_fitnesses = dict()
        bppms              = list()
        mfe_dbs             = list()
        free_energies       = list()
        ens_divs            = list()
        on_scores          = list()
        off_scores         = list()
        excess_scores      = list()
        affinity_scores = list()
        
        #process logic results
        for state_id in range(num_states):
            task_id = ind_idx_batch*num_states + state_id
            (mfe_db, bppm, mfe_ad, freeEnergy, ensDiv) = logic_results[task_id]

            on_score, off_score, excess_score = evaluate_rz_activity(bppm, name_to_bases, off_cap)            

            if truth_vector[state_id] == 1:
                on_scores.append( on_score )
            else:
                off_scores.append( off_score )
                excess_scores.append(excess_score)
           
            mfe_dbs.append(mfe_db)           
            bppms.append(bppm)
            free_energies.append(freeEnergy)
            ens_divs.append(ensDiv)  

        #Calculate thermobalance
        thermobalance_score = score_thermobalance(free_energies)    

        #Calculate affinity
        for input_id in range(num_inputs): #There is 1 affinity task per input
            task_id = ind_idx_batch*num_inputs + input_id
            (mfe_db, bppm, mfe_ad, freeEnergy, ensDiv) = affinity_results[task_id]
            affinity_scores.append( score_affinity( input_id , name_to_bases, bppm, affinity_tasks[task_id][2], affinity_tasks[task_id][3] ) )


        potential_fitnesses['avgON']     = avg(on_scores)
        potential_fitnesses['avgOFF']    = avg(off_scores)
        potential_fitnesses['minON']     = min(on_scores)
        potential_fitnesses['minOFF']    = min(off_scores)
        potential_fitnesses['ON']        = avg([ potential_fitnesses['minON'], potential_fitnesses['avgON'] ])
        potential_fitnesses['OFF']       = avg([ potential_fitnesses['minOFF'], potential_fitnesses['avgOFF'] ])
        potential_fitnesses['switch']    = potential_fitnesses['ON'] + potential_fitnesses['OFF']
        potential_fitnesses['switchMin'] = sum([ potential_fitnesses['minON'], potential_fitnesses['minOFF'] ])
        potential_fitnesses['ensDiv']    = -avg(ens_divs)
        potential_fitnesses['excess']    = avg(excess_scores)
        potential_fitnesses['affinity'] = avg( [ avg(affinity_scores), min(affinity_scores) ] ) 
        potential_fitnesses['thermobalance'] = thermobalance_score 

        t_sppms    = coarse_grain_pheno(bppms, base_to_seg, er_seg_ids, constant_rz_seg_ids, truth_vector, bool_TT, bool_constant_rz)

        results.append( ( potential_fitnesses, t_sppms, mfe_dbs) )
    return results
