from truth_seq_er_ea import TruthSeqEr

import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Your algorithm description here')

    # Required arguments
    parser.add_argument('--num_inputs', type=int, required=True, help='Number of inputs')
    parser.add_argument('--truth_vector', type=str, required=True, help='Truth vector as a sequence of 0s and 1s')

    # Optional arguments
    parser.add_argument('--num_processes', type=int, default=min(12,os.cpu_count()), help='Number of processes (default: 12)')
    parser.add_argument('--novelty_search', type=lambda s: s.lower() in ['yes', 'true', '1'], default=True, help='Enable novelty search (default: yes)')
    parser.add_argument('--neighborhood_size', type=int, default=30, help='Novelty neighborhood size (default: 30)')
    parser.add_argument('--viability_nullification', type=lambda s: s.lower() in ['yes', 'true', '1'], default=True, help='Enable viability nullification (default: yes)')
    parser.add_argument('--viability_threshold', type=float, default=0.90, help='Viability threshold (default: 0.90)')
    parser.add_argument('--viability_knee_value', type=float, default=0.45, help='Viability knee value (default: 0.45)')
    parser.add_argument('--viability_knee_generation', type=int, default=50, help='Viability knee generation (default: 50)')
    parser.add_argument('--num_genss', type=int, default=200, help='Number of generations (default: 200)')
    parser.add_argument('--num_indss', type=int, default=300, help='Number of individuals (default: 300)')
    parser.add_argument('--stringent_mode', type=lambda s: s.lower() in ['yes', 'true', '1'], default=False, help='Enable stringent mode (default: no). If enabled, the final population is filtered for thermobalance and affinity')
    parser.add_argument('--folding_temperature', type=float, default=37.0, help='Folding temperature')
    parser.add_argument('--input_seqs', type=str, default=None, help='Input constraints as a sequence of characters from the set {N, A, U, G, C}')


    return parser.parse_args()

def validate_args(args):
    min_num_inputs = 1
    max_num_inputs = 3
    if args.num_inputs > max_num_inputs or args.num_inputs < min_num_inputs:
        raise ValueError("The number of inputs must be betwen 1 and 3")

    # Convert the truth vector to a list of 0s and 1s
    truth_vector = [int(x) for x in args.truth_vector if x in ['0', '1']]
    args.truth_vector = truth_vector


    # Check if the length of the truth vector is equal to 2 to the power of the number of inputs
    if len(truth_vector) != 2 ** args.num_inputs:
        raise ValueError('The length of the truth vector is not equal to 2 to the power of the number of inputs')



    if args.input_seqs != None:
        # Check if the input constraints are valid
        valid_chars = set('NAUGC')
        input_seqs = [x for x in args.input_seqs if x in valid_chars]

        # Check if the input_seqs are valid
        valid_chars = set('NAUGC')
        min_length = 14
        max_length = 24

        input_seqs = args.input_seqs.split(',')
        if len(input_seqs) != args.num_inputs:
            raise ValueError('Invalid number of input constraints provided')

        for constraint in input_seqs:
            if not all(c in valid_chars for c in constraint):
                raise ValueError(f'Invalid character in input constraint: {constraint}')

            if not (min_length <= len(constraint) <= max_length):
                raise ValueError(f'Input constraint {constraint} length is not within the allowed range ({min_length}, {max_length})')
    else:
        input_seqs = ['N' * 22] * args.num_inputs
    args.input_seqs = input_seqs



    return args

if __name__ == '__main__':
    args = parse_args()

    try:
        args = validate_args(args)
    except ValueError as e:
        print(f"Error: {e}")
        exit(1)


    stem_lengths      = [7,8,7]
    sel_fits = ['ON', 'OFF']
    if args.novelty_search:
        sel_fits.append('novelty')
    if args.viability_nullification:
        sel_fits = [fit + "_via" for fit in sel_fits]
    bool_fast_folding = 0
    bool_TT           = 1 #transparent rz, transparent er

    bool_negator = not args.truth_vector[-1] #If the last entry of truth vector is 0, then we want it to kill the rz so we need the logic linker
    bool_linker_flex  = 0
    input_lengths     = [len(input_seq) for input_seq in args.input_seqs]
    bool_start_GG     = 0
    template_params   = [input_lengths, stem_lengths, bool_negator, bool_linker_flex, bool_start_GG, args.input_seqs]    
    viability_params  = [-1.00, 0.90, 1, -0.05, 49]

    pSM, pSP, sSM, cP, bS = 0, [], 0, 0, 1
    mutation_rate      = 4
    mP                = [1, 0]
    local_fit_names   = ['switch'] #vestigial
    initial_pop       = None
    ea_params          = (args.num_indss, args.num_genss, pSM, pSP, sSM, cP, mutation_rate, mP, sel_fits, bS) 

    job = TruthSeqEr(template_params, args.num_inputs, args.truth_vector, local_fit_names, 
                     viability_params, bool_fast_folding, args.folding_temperature, bool_TT,
                     args.stringent_mode, args.neighborhood_size, args.num_processes)

    print("Generating ribogate designs implementing function with truth vector", args.truth_vector)

    job.run(*ea_params)