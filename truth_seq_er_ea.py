from ea import *
from representation import SRS, SCC
from utility import multicore, search_list_objs, non_dominated_sort
from semantics import *
from rz_templates import build_template
import csv

class TruthSeqEr(EA):
    def __init__(self, template_params, num_inputs, truth_vector, local_fit_names, 
                 viability_params, bool_fast_folding, folding_temp, bool_TT, bool_stringent,
                 neighborhood_size, num_processes):
        EA.__init__(self)
        self.truth_vector        = truth_vector
        self.num_inputs          = num_inputs
        self.local_fit_names     = local_fit_names
        self.bool_fast_folding   = bool_fast_folding
        self.folding_temp        = folding_temp
        self.bool_TT             = bool_TT
        self.bool_constant_rz    = 0 #Vestigital. TO DO: remove without breaking anything
        self.bool_stringent      = bool_stringent
        self.off_cap              = 0.8
        self.min_cap, self.max_cap, self.bool_knee, self.knee_cap, self.knee_gen = viability_params
        self.neighborhood_size = neighborhood_size
        self.num_processes = num_processes

        self.func_num            = bin_to_dec(truth_vector)
        self.func_name           = 'f-' + str(self.func_num) + '-' + str(num_inputs) + '-I' 
        self.num_states = int(2**self.num_inputs)


        #Build template
        self.template = build_template(*template_params)

      
    def create_individual(self):
        ind              = SRS()
        ind.sccs         = copy.deepcopy(self.template)      
        ind.truth_vector = copy.deepcopy(self.truth_vector)

        for scc in ind.sccs:
            scc.update_allowed() #This will also sample the bccs since the empty bcc will violate the allowed length constraint
    
        ind.update_mutation_maps(['bccs'])
        ind.update_mutation_probs()

        ind.coarse_pheno = [] #!!!This might be needed for resuming

        return ind

    def copy_ind(self, ind):
        newInd                     = SRS()
        newInd.sccs                = copy.deepcopy(ind.sccs)
        newInd.id                  = ind.id
        newInd.bcc_to_scc          = copy.deepcopy(ind.bcc_to_scc)
        newInd.mutation_probs       = copy.deepcopy(ind.mutation_probs)
        newInd.truth_vector        = copy.deepcopy(ind.truth_vector)
        newInd.input_states        = copy.deepcopy(ind.input_states)
        newInd.name_to_bases       = copy.deepcopy(ind.name_to_bases)
        newInd.potential_fitnesses  = copy.deepcopy(ind.potential_fitnesses)
        newInd.coarse_pheno         = copy.deepcopy(ind.coarse_pheno) #!!!This might be needed for resuming
        return newInd

    def generate_semantics_offspring(self):
        population_full = [search_list_objs(self.population, 'id', ind.id)[0] for ind in self.parents] #Because I removed the deep copies, the parents don't have any phenotype

        if self.bool_fast_folding == 0 or (self.bool_fast_folding == 1 and self.cur_gen == self.num_gens - 1):
            bool300 = 1
            pool    = population_full + self.offspring #We need to make sure ALL individuals in the final population go through bool 300 folding
        else:
            bool300 = 0
            pool    = self.offspring

        batch = list()
        for ind in pool:
            #Prepare tasks
            ind.GenerateSegments()
            ind.GenerateStrands()
            ind.update_base_maps()
            ind.BuildColorString()
            ind.generate_folding_tasks(self.num_inputs)
            #batch += ind.tasks
            batch.append(ind.tasks)

        rz_seg_names  = [seg.name for seg in ind.segments if seg.strand == 0] #!!Hardcoded for cis-acting
        common_data   = (ind.base_to_seg, ind.name_to_bases, ind.er_seg_ids, ind.constant_rz_seg_ids, rz_seg_names, self.truth_vector, self.num_inputs, bool300, self.folding_temp, self.bool_TT, self.bool_constant_rz, self.off_cap)

        results                   = multicore(batch, pheno_perf_fits, self.num_processes, common_data)
        #results                    = pheno_perf_fits([batch, common_data]) #Single core


        #Store results back in individuals and prepare distance calculation
        phenos = list()
        for i, ind in enumerate(pool):
            result = results[i]
            ind.potential_fitnesses = result[0]
            ind.coarse_pheno        = result[1]
            ind.mfe_dbs             = result[2]

        #Novelty
        #print("Starting novelty")        
        novelty_pool    = population_full + self.offspring
        num_inds         = len(novelty_pool)
        distance_matrix = self.calculate_distance_matrix(novelty_pool)  

        for i in range(num_inds):
            distances = list()
            for j in range(num_inds):
                distances.append(distance_matrix[i, j])                      

            neighbor_idxs = arg_sort(distances)[0:self.neighborhood_size]
            novelty_pool[i].potential_fitnesses['novelty'] = sum([distances[n_idx] for n_idx in neighbor_idxs])

            #Local competition
            dataset = list()
            dataset.append( [novelty_pool[i].potential_fitnesses[f_name] for f_name in self.local_fit_names] )
            for n_idx in neighbor_idxs:
                dataset.append( [novelty_pool[n_idx].potential_fitnesses[f_name] for f_name in self.local_fit_names] )
            ff_idxs = non_dominated_sort(dataset)
            for ff_id, ff_idx in enumerate(ff_idxs): #I'm assuming this is used in the line a bit below for 'rank'
                if 0 in ff_idx:
                    break
            novelty_pool[i].potential_fitnesses['rank']    = 1/float(1 + ff_id)

        #Viability nullification
        if self.bool_knee == 1:
            if self.cur_gen < self.knee_gen:
                self.cur_viability_threshold = self.min_cap + self.cur_gen * (self.knee_cap - self.min_cap) / float(self.knee_gen)
            else:
                self.cur_viability_threshold = self.knee_cap + (self.cur_gen - self.knee_gen) * (self.max_cap - self.knee_cap) / float(self.num_gens -1 - self.knee_gen)
        else:
            self.cur_viability_threshold = self.min_cap + self.cur_gen * (self.max_cap - self.min_cap) / float(self.num_gens -1)
        #print(self.cur_viability_threshold)     
        
        pool = self.population + self.offspring
        pf_names = list(ind.potential_fitnesses.keys())
        for ind in pool:
            for pfName in pf_names:
                if '_via' not in pfName: #Without this condition, we end up with name_via_via_via ... 
                    if ind.potential_fitnesses['switchMin'] < self.cur_viability_threshold:
                        ind.potential_fitnesses[pfName + "_via"] = -1000
                    else:
                        ind.potential_fitnesses[pfName + "_via"] = ind.potential_fitnesses[pfName]

        
    def calculate_distance_matrix(self, pool):
        batch   = list()
        num_inds = len(pool)
        for i in range(num_inds): 
            for j in range(num_inds):
                batch.append((pool[i].coarse_pheno, pool[j].coarse_pheno))
        results = multicore(batch, dot_distance_batch, self.num_processes, [])

        taskId               = 0
        distance_matrix = np.zeros(( num_inds, num_inds )) #Used later for MDS
        for i in range(num_inds):
            for j in range(num_inds):
                distance_matrix[i, j] = results[taskId]
                taskId += 1   

        return distance_matrix

    def post_process(self):
        #Determine which individuals are acceptable
        acceptable_pool = list()
        for ind in self.population:
            is_reliable = ind.potential_fitnesses['affinity'] >= 0.85 and ind.potential_fitnesses['thermobalance'] >= 6 and ind.potential_fitnesses['thermobalance'] <= 10
            #!!!!! if ind.potential_fitnesses['switchMin'] >= 0.90 and (is_reliable or not self.bool_stringent):
            if 1: #!!! Debugging
                ind.potential_fitnesses['acceptable'] = 1
                acceptable_pool.append(ind)
            else:
                ind.potential_fitnesses['acceptable'] = 0

        #Marginal diversity gain(MDG)
        #Handle unacceptable individuals
        for ind in self.population:
            if ind.potential_fitnesses['acceptable'] == 0:               
                ind.potential_fitnesses['mdg']        = -len(acceptable_pool) #This is effectively -inf

        #Handle acceptable individuals
        pool              = acceptable_pool
        num_inds           = len(pool)
        if num_inds == 0: #no acceptable individuals
            pass 
        elif num_inds == 1:
            pool[0].potential_fitnesses['mdg'] = 0
        else:
            distance_matrix   = self.calculate_distance_matrix(pool)      
            mdg_idxs          = list()
            furthest_pair_idx = np.argmax(distance_matrix)                          #This is the (single) index of the pair in a 2d matrix
            mdg_idxs += [furthest_pair_idx % num_inds, furthest_pair_idx // num_inds] #From this pair index, we get the index of each element of the pair
            while len(mdg_idxs) < num_inds:
                candidate_distances = list()
                unplaced_idxs       = [idx for idx in range(num_inds) if idx not in mdg_idxs]
                for unplaced_idx in unplaced_idxs:
                    marginal_distances = list()
                    for mdg_idx in mdg_idxs:
                        marginal_distances.append(distance_matrix[unplaced_idx, mdg_idx])
                    candidate_distances.append( min(marginal_distances)  )
                temp_idx = arg_sort(candidate_distances)[-1]
                mdg_idxs.append( unplaced_idxs[temp_idx] )

            pool[ mdg_idxs[0] ].potential_fitnesses['mdg'] = 0
            pool[ mdg_idxs[1] ].potential_fitnesses['mdg'] = 0
            for mdg_rank in range(2, num_inds):
                pool[ mdg_idxs[mdg_rank] ].potential_fitnesses['mdg'] = -mdg_rank

        #Export to CSV
        mdg_vals        = [ind.potential_fitnesses['mdg'] for ind in pool]
        ind_idxs_sorted = arg_sort(mdg_vals)[::-1]
        inds_ranked     = [pool[idx] for idx in ind_idxs_sorted]

        designs_folder = 'designs'
        if not os.path.exists(designs_folder):
            os.makedirs(designs_folder)

        header = ['Design ID', 'Gate sequence']
        header += ['Input ' + str(i+1) + ' sequence' for i in range(self.num_inputs)] 
        rows = list()
        for ind_idx, ind in enumerate(inds_ranked):
            rows.append([ind_idx] + ind.strands)

        with open('designs/ribogate_designs_' + self.func_name + '.csv','w') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(header)
            csv_out.writerows(rows)

        header = ['Design ID', 'State index', 'MFE secondary structure in dot bracket notation',
                  'Gate sequence', 'Folding constraint']
        rows = list()
        for ind_idx, ind in enumerate(inds_ranked):
            folding_constraints = [task[1] for task in ind.tasks[0]]
            for i in range(self.num_states):
                rows.append( (ind_idx, i, ind.mfe_dbs[i], ind.strands[0], folding_constraints[i]) )
    
        with open('designs/ribogate_structures_' + self.func_name + '.csv','w') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(header)
            csv_out.writerows(rows)

    def evaluate_fitness(self, ind):
        pass

    def mutate_parent(self, ind):
        for m in range(self.mutation_rate):
            #Randomly select a (global) index of a bcc. Use the bcc_to_scc map to get the scc containing that bcc and the bcc's relative position within that scc
            num_bccs_total       = len (ind.bcc_to_scc.keys() )
            bcc_id_all           = np.random.choice( range( num_bccs_total),  1, p=ind.mutation_probs)[0]
            (scc_id, bcc_id_scc) = ind.bcc_to_scc[bcc_id_all] 
            scc                  = ind.sccs[scc_id]
            bcc                  = scc.bccs[bcc_id_scc]
            options              = [i for i in scc.allowed[bcc_id_scc] if i != bcc] #Every bcc but the current one                                 
            scc.bccs[bcc_id_scc] = np.random.choice(options)       


