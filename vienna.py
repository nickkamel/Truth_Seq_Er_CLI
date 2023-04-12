import subprocess
from collections import defaultdict
import os

def parse_dot_bracket(dot_bracket):
    partners = dict()
    strands = dot_bracket.split('&')
    num_strands = len(strands)
    strandNames = ['A','B'] 

    pairs = list()
    waiting_to_be_paired = list()
    waiting_to_be_paired_pseudo = list()
    for strand_id,strand in enumerate(strands):
        for base_id,base_bracket_type in enumerate(strand):
            base_name = (strand_id,base_id+1)
            if base_bracket_type == "(":
                waiting_to_be_paired.append(base_name)
            elif base_bracket_type == "[":
                waiting_to_be_paired_pseudo.append(base_name)
            else:
                if base_bracket_type == ")":
                    partner = waiting_to_be_paired.pop()
                elif base_bracket_type == "]":
                    partner = waiting_to_be_paired_pseudo.pop()
                elif base_bracket_type == ".":
                    partner = (None, None)

                pairs.append((partner,base_name)) #non symmetric and out of order (for example, if we have (...), the () pair will be added after the dots
                partners[base_name] = partner
                if base_bracket_type != ".": #we don't want symmetry for the dummy (0,0) values
                    partners[partner] = base_name #symmetry
    #print partners
    return partners

def Fold(tasks, bool300, folding_temp, num_strands):
    input_data     = ''
    for task in tasks:
        seq        = task[0]
        constraint = task[1]
        input_data = input_data + seq + "\n" + constraint + "\n"
    input_data     = input_data[:-1] #Get rid of the last newline

    if bool300 == 1:
        str300   = '300'
    else:
        str300   = ''
    if num_strands == 1:
        app_name = "\RNAFold"
    else:
        app_name = "\RNACofold"

    script_dir = os.path.dirname(os.path.abspath(__file__))
    p1 = subprocess.Popen([script_dir + app_name + str300 + ".exe","-p","-C","--noPS","--bppmThreshold=1e-5", "--temp=" + str(folding_temp)], stdin=subprocess.PIPE,stdout=subprocess.PIPE)   
    
    input_data    = str.encode(input_data) #For python 3
    results       = p1.communicate(input=input_data)[0]
    #results_lines = results.split('\r\n')
    results_lines = results.decode().split('\r\n') #For python 3
    results_lines = results_lines[:-1] #Get rid of the last newline
    #print(results_lines)

    results_starts = list()
    for line_idx, line in enumerate(results_lines):
        if line == "Secondary structure":
            results_starts.append(line_idx)
    results_starts.append(len(results_lines)) #Add a virtual start id for what would be the next result
    #print(results_starts)

    num_results = len(results_starts) - 1

    processed_results = list()
    for i in range(num_results):
        result_lines = results_lines[results_starts[i]:results_starts[i+1]]
        #print(result_lines)

        mfe_db       = result_lines[1]
        mfe_ad       = parse_dot_bracket(mfe_db)
        mfe_energy   = float(result_lines[3])
        if num_strands == 1:
            ens_div      = float(result_lines[-1])

        for line_idx, line in enumerate(result_lines):
            if line == "BPPM start":
                bppm_start_line_idx = line_idx + 1
            if line == "BPPM end":
                bppm_end_line_idx   = line_idx - 1

        #print(bppm_start_line, bppm_end_line)
        bppm     = defaultdict(int)
        adj_dict   = dict()
        #num_bases = len(seq) #This will fail for variable length sequences like the case of affinity tasks with variable length inputs
        num_bases = len(tasks[i][0]) 
        for i in range(num_bases):
            adj_dict[i + 1] = list()
        for line_idx in range(bppm_start_line_idx, bppm_end_line_idx + 1):
            line = result_lines[line_idx]
            split = line.split('\t')
            b_id_1 = int(split[0])
            b_id_2 = int(split[1])
            prob = float(split[2]) ** 2 #The values given in the ps file are the sqrt of the probability
            bppm[(b_id_1, b_id_2)] = prob
            bppm[(b_id_2, b_id_1)] = prob
            adj_dict[b_id_1].append(b_id_2)
            adj_dict[b_id_2].append(b_id_1)
        for b_id_1 in adj_dict.keys():
            pairedSum = 0
            for b_id_2 in adj_dict[b_id_1]:
                pairedSum += bppm[(b_id_1, b_id_2)]
            unpaired_prob = 1 - pairedSum

            bppm[(b_id_1, None)]  = unpaired_prob
            bppm[(None, b_id_1 )] = unpaired_prob

        if num_strands == 1:
            processed_results.append((mfe_db, bppm, mfe_ad, mfe_energy, ens_div))
        else:
            processed_results.append((mfe_db, bppm, mfe_ad, mfe_energy, None))
    return processed_results



