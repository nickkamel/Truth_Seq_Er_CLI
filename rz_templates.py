from representation import SCC
from utility import *
from collections import defaultdict
import copy

def build_template(input_lengths, s_lengths, bool_negator, bool_linker_flex, bool_start_GG, input_seqs):
    num_inputs       = len(input_lengths)
    #seg_names        = ['startGG', 's1A', 'C1', 's2A', 'L0', 'OBS1', 'L1', 'OBS2', 'L2', 'OBS3', 'L3', 's2B', 'C2', 's3A', 's3H', 's3B', 'C3', 's1B']
    seg_names        = ['s1A', 'C1', 's2A', 'L0', 'OBS1', 'L1', 'OBS2', 'L2', 'OBS3', 'L3', 's2B', 'C2', 's3A', 's3H', 's3B', 'C3', 's1B']
    negator_len = 14
    linker_flex_len  = 3

    #This is needed to generate the correct segs_pos for scc
    if num_inputs < 3:
        seg_names = remove_from_list(seg_names, ['OBS3', 'L3'])
    if num_inputs < 2:
        seg_names = remove_from_list(seg_names, ['OBS2', 'L2'])

    if bool_linker_flex == 0:
        seg_names = remove_from_list(seg_names, ['L1', 'L2', 'L3'])
        if bool_negator == 0:
            seg_names = remove_from_list(seg_names, ['L0'])


    #if bool_start_GG == 0:
    #    seg_names = remove_from_list(seg_names, ['startGG'])
    segs_pos = defaultdict(int)
    seg_pos = 0
    for seg_name in seg_names:
        segs_pos[seg_name] = seg_pos
        seg_pos += 1

    if bool_negator == 0:
        l0_len = linker_flex_len
    else:
        l0_len = negator_len

    #This is to make code clean
    stem_constraints = dict()
    if bool_start_GG:
        #stem_constraints['s1'] = ['GG' + 'N'*(s_lengths[0]-2)     , 'N'*s_lengths[0]]
        stem_constraints['s1'] = ['GGH' + 'N'*(s_lengths[0]-3)     , 'N'*s_lengths[0]]
    else:
        stem_constraints['s1'] = ['N'*s_lengths[0]          , 'N'*s_lengths[0]] 
    stem_constraints['s2'] = ['G' + 'N'*(s_lengths[1]-1), 'N'*(s_lengths[1]-1) + 'C'] 
    stem_constraints['s3'] = ['A' + 'N'*(s_lengths[2]-1), 'N'*(s_lengths[2]-1) + 'U'] 

    obs_in_constraints = list()
    for input_idx in range(num_inputs):
        input_length = input_lengths[input_idx]
        if input_seqs[input_idx] == None:
            input_constraint = 'N'*input_length
        else:
            input_constraint = input_seqs[input_idx]
        obs_constraint = 'N'*input_length
        obs_in_constraints.append( [obs_constraint, input_constraint] )

    linker_sccs = list()
    if bool_negator or bool_linker_flex:
        linker_sccs.append( SCC(['L0']            , 1, l0_len          , ['N'*l0_len]          , (0,) , (segs_pos['L0'],)                 , ('#4363d8',)            ) )
    if bool_linker_flex:
        linker_sccs.append( SCC(['L1']            , 1, linker_flex_len , ['N'*linker_flex_len] , (0,) , (segs_pos['L1'],)                 , ('#911eb4',)            ) )
    if bool_linker_flex and num_inputs > 1:
        linker_sccs.append( SCC(['L2']            , 1, linker_flex_len , ['N'*linker_flex_len] , (0,) , (segs_pos['L2'],)                 , ('#f032e6',)            ) )
    if bool_linker_flex and num_inputs > 2:
        linker_sccs.append( SCC(['L3']            , 1, linker_flex_len , ['N'*linker_flex_len] , (0,) , (segs_pos['L3'],)                 , ('#fabebe',)            ) )

    obs_in_sccs = [         SCC(['OBS1', 'Input1'], 2, input_lengths[0], obs_in_constraints[0] , (0,1), (segs_pos['OBS1'], 0)             , ('#f58231', '#FFFFFF') )]
    if num_inputs > 1:
        obs_in_sccs.append( SCC(['OBS2', 'Input2'], 2, input_lengths[1], obs_in_constraints[1] , (0,2), (segs_pos['OBS2'], 0)             , ('#46f0f0', '#FFFFFF') ) )
    if num_inputs > 2:
        obs_in_sccs.append( SCC(['OBS3', 'Input3'], 2, input_lengths[2], obs_in_constraints[2] , (0,3), (segs_pos['OBS3'], 0)             , ('#bcf60c', '#FFFFFF') ) )

    stem_core_sccs = [
                            SCC(['s1A', 's1B']    , 2, s_lengths[0]    , stem_constraints['s1'], (0,0), (segs_pos['s1A'], segs_pos['s1B']), ('#e6194b', '#808000') ), 
                            SCC(['C1']            , 1, 7               , ['CUGANGA']           , (0,) , (segs_pos['C1'],)                 , ('#3cb44b',)            ),
                            SCC(['s2A', 's2B']    , 2, s_lengths[1]    , stem_constraints['s2'], (0,0), (segs_pos['s2A'], segs_pos['s2B']), ('#ffe119', '#008080') ), 
                            SCC(['C2']            , 1, 3               , ['GAA']               , (0,) , (segs_pos['C2'],)                 , ('#e6beff',)            ),
                            SCC(['s3A', 's3B']    , 2, s_lengths[2]    , stem_constraints['s3'], (0,0), (segs_pos['s3A'], segs_pos['s3B']), ('#000000', '#800000') ),
                            SCC(['s3H']           , 1, 4               , ['NNNN']              , (0,) , (segs_pos['s3H'],)                , ('#fffac8',)            ),
                            SCC(['C3']            , 1, 1               , ['H']                 , (0,) , (segs_pos['C3'],)                 , ('#aaffc3',)            )
    ]

    template = linker_sccs + obs_in_sccs + stem_core_sccs


    return template
