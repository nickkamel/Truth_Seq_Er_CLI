from utility import avg
def evaluate_rz_activity(bppm, name_to_bases, off_cap):  #This is a core disruption score
    #core
    rz_elem_scores     = list()

    core_score_abs = 0
    for core_seg_name in ['C1', 'C2', 'C3']:
        core_seg_score_abs = 0
        core_seg_base_ids  = name_to_bases[core_seg_name]
        for base_id in core_seg_base_ids:
            unpaired_prob = bppm[(base_id, None)]
            paired_prob   = 1 - unpaired_prob
            core_seg_score_abs += paired_prob
        core_score_abs += core_seg_score_abs
        if core_seg_name == 'C1':
            off_score_uncapped = core_seg_score_abs / len(core_seg_base_ids)
            if off_score_uncapped > off_cap:
                off_score = off_cap / off_cap #NK. Looking back at this, is there reason I don't just do off_score = 1?
            else:
                off_score = off_score_uncapped / off_cap
    #Here the off score only looks at C1 and will not reward disrupting it above a certain amount (this is why have have the cap). The reason for this is to make the ribozyme more switchable
 
    num_core_base_ids  = len(name_to_bases['C1'] + name_to_bases['C2'] + name_to_bases['C3'])
    core_score         = core_score_abs / num_core_base_ids
    rz_elem_scores.append(core_score)

    #stems
    excess_score = 0
    for i in range(1, 4):
        stem_name     = 's' + str(i)
        stem_base_ids = list(zip(name_to_bases[stem_name + 'A'], name_to_bases[stem_name + 'B'][::-1]))
        stem_score_  = 0
        stem_score_A = 0
        stem_score_B = 0
        wUp = 1 #unpaired weight
        wMp = 1 #paired weight
        for base_id, partner_base_id in stem_base_ids: 
            paired_prob            = bppm[(base_id, partner_base_id)]
            mispaired_unpaired_prob = 1 - paired_prob
            unpaired_prob1         = bppm[(base_id, None)]
            mispaired_prob1        = mispaired_unpaired_prob - unpaired_prob1
            unpaired_prob2         = bppm[(partner_base_id, None)]
            mispaired_prob2        = mispaired_unpaired_prob - unpaired_prob2
            stem_score_A          = stem_score_A + wUp*unpaired_prob1 + wMp*mispaired_prob1
            stem_score_B          = stem_score_B + wUp*unpaired_prob2 + wMp*mispaired_prob2

        stem_score_A = stem_score_A / float( len(stem_base_ids) )
        stem_score_B = stem_score_B / float( len(stem_base_ids) )
        #stem_score  = max(stem_score_A, stem_score_B)
        stem_score  = avg([stem_score_A, stem_score_B])
        rz_elem_scores.append(stem_score)

    #excess_score is used for the case where we want to stick closer to Penchosvky and not disrupt stems 1 and 3 in neither the on nor off state
    excess_score = -avg([rz_elem_scores[1], rz_elem_scores[3]]) #average s1 and s3 stem scores
    on_score     = -avg(rz_elem_scores)

    return on_score, off_score, excess_score

#We want to see if the input actually binds to its target OBS when it is concatenated on the same strand with it
def score_affinity(input_id, name_to_bases, bppm, rz_len, input_len):
    obs_ids   = name_to_bases['OBS' + str(input_id + 1)] #The first input has an id of 0 here, hence the +1.
    input_ids = range(rz_len + 1, rz_len + input_len + 1)

    score = 0
    for base_id, partner_base_id in zip(obs_ids, input_ids[::-1]): 
        paired_prob = bppm[(base_id, partner_base_id)]
        score += paired_prob
    return score / float(input_len)

def score_thermobalance(free_energies):
    if len(free_energies) == 8:
        #0 is the most stable because when we co-fold, we force the bases to be unbound which make it "less stable"
        #Hasse: 
        #0 -> 1,2,4
        #1 -> 3,5
        #2 -> 3,6
        #3 -> 7
        #4 -> 5,6
        #5 -> 7
        #6 -> 7
        energy_gaps = [free_energies[1] - free_energies[0], free_energies[2] - free_energies[0], free_energies[4] - free_energies[0],
                        free_energies[3] - free_energies[1], free_energies[5] - free_energies[1], 
                        free_energies[3] - free_energies[2], free_energies[6] - free_energies[2], 
                        free_energies[7] - free_energies[3],
                        free_energies[5] - free_energies[4], free_energies[6] - free_energies[4],
                        free_energies[7] - free_energies[5],
                        free_energies[7] - free_energies[6]] 
    elif len(free_energies) == 4:
        #0 -> 1, 2
        #1 -> 3
        #2 -> 3
        energy_gaps = [free_energies[1] - free_energies[0], free_energies[2] - free_energies[0],
                        free_energies[3] - free_energies[1], free_energies[3] - free_energies[2]]
    elif len(free_energies) == 2:
        #0 -> 1
        energy_gaps = [free_energies[1] - free_energies[0]]

    return avg(energy_gaps)
