from folding import generate_fold_constraint, build_map, allowed_lut
import numpy as np
from utility import*

class SRS():
    def __init__(self):
        self.sccs       = list()
        self.segments   = list() #This is always sorted
        self.strands    = list()

        #For constraints, core score. 
        #input: name of a segment
        #output: list of base (global) ids contained in that segment
        self.name_to_bases        = dict() 
        #For coarse graining. 
        #input:   global id of base 
        ##output: pair of (global) id of segment containing that base and position of base within that segment
        self.base_to_seg          = dict()  

        #bcc: base connected component
        #scc: segment connected component
        #bcc_to_scc: inputs a scalar integer representing the global id of a bcc. Outputs a pair of integer representing the scc that contains it and its local id within that scc
        self.bcc_to_scc           = dict()

        self.mutation_probs        = list()
        self.geno                 = ''
        self.folding_constraints  = list()
        self.potential_fitnesses   = dict()

        self.folded_strand_ids      = [0, 1] #This is used for building the color string. The ribozyme and the substrate, not the repeats. 

    def GenerateSegments(self):
        unordered_segments = list()
        for scc in self.sccs:
            if scc.size == 1:
                seq = ''.join(scc.bccs)
                seg = Segment(scc.names[0], scc.length, scc.segs_strand[0], scc.segs_pos[0], scc.segs_Colors[0], seq)
                unordered_segments.append(seg)
            elif scc.size == 2:
                seq_1 = list()
                seq_2 = list()
                for bcc in scc.bccs:
                   seq_1.append(bcc[0])
                   seq_2.append(bcc[1])
                seq_2  = seq_2[::-1]
                seq_1  = ''.join(seq_1)
                seq_2  = ''.join(seq_2)
                #print(scc.names, scc.segs_strand, scc.segs_pos, scc.segs_Colors)
                seg1  = Segment(scc.names[0], scc.length, scc.segs_strand[0], scc.segs_pos[0], scc.segs_Colors[0], seq_1)
                seg2  = Segment(scc.names[1], scc.length, scc.segs_strand[1], scc.segs_pos[1], scc.segs_Colors[1], seq_2)
                unordered_segments += [seg1, seg2]
                         
        #Sort segments
        self.segments = list()
        self.num_strands = max([seg.strand for seg in unordered_segments]) + 1
        for sd_id in range(self.num_strands):
            strand_segs = [seg for seg in unordered_segments if seg.strand == sd_id]
            sorted_idxs = arg_sort([seg.pos for seg in strand_segs])
            for idx in sorted_idxs:
                self.segments.append( strand_segs[idx] )

        #Extension seg ids
        self.er_seg_ids = list()
        #The most recent pool had a bug where the linker segments were not considered extension segments. This is because the template names were Lx.. but here they were interx ->  potential_extension_names = ['inter0', 'OBS1', 'inter1', 'OBS2', 'inter2', 'OBS3', 'inter3', 'OBS']
        potential_extension_names = ['L0', 'OBS1', 'L1', 'OBS2', 'L2', 'OBS3', 'L3', 'OBS'] #!!!Potential problem with this approach is that it gives a lot of weight to flexibility linkers which have little computationa; effect
        seg_names = [seg.name for seg in self.segments]
        for peName in potential_extension_names:
            if peName in seg_names:
                idx = inv_list(seg_names, peName)
                self.er_seg_ids.append(idx)

        #These are used if we want to stick closer to Penchovsky and only disrupt stem 2
        self.constant_rz_seg_ids = list()
        constant_rz_seg_names = ['startGG', 's1A', 's3A', 's3H', 's3B', 'C3', 's1B']
        seg_names = [seg.name for seg in self.segments]
        for rc_name in constant_rz_seg_names:
            if rc_name in seg_names:
                idx = inv_list(seg_names, rc_name)
                self.constant_rz_seg_ids.append(idx)

    def GenerateStrands(self):
        self.strands = list()
        for sd_id in range(self.num_strands):
            strand = ''
            strand_segs = [seg for seg in self.segments if seg.strand == sd_id]
            for seg in strand_segs:
                strand += seg.seq
            self.strands.append(strand)

    def BuildColorString(self):
        self.colorString = ''
        accum = 0 #color string is 0-indexed
        for sd_id in self.folded_strand_ids:
            strand_segs = [seg for seg in self.segments if seg.strand == sd_id]
            for seg in strand_segs:
                self.colorString += str(accum) + '-' + str(accum + seg.length - 1) + ':' + seg.color + ' ' #?? -1 because color string end is inclusive
                accum += seg.length
            accum += 2 #FORNA requires a gap of 2 between strands in order to display colors correctly

    def generate_folding_tasks(self, num_inputs):
        gate_len = len(self.strands[0])

        #Logic tasks
        logic_tasks = list()
        if num_inputs == 1:
            input_states = [(0, ), (1, )]
        elif num_inputs == 2:
            input_states = [(0, 0), (0, 1), (1,0), (1,1)]
        elif num_inputs == 3:
            input_states = [(0, 0, 0), (0, 0, 1), (0, 1,0), (0, 1,1), (1, 0, 0), (1, 0, 1), (1, 1,0), (1, 1,1)]
  
        for input_state in input_states:
            constrained_bases = list()
            for col_id, col in enumerate(input_state):
                if col == 1:
                    constrained_bases += self.name_to_bases['OBS' + str(col_id + 1)] #These are 1-indexed and right inclusive
            state_constraint = generate_fold_constraint( gate_len, constrained_bases ) 
            logic_tasks.append( (self.strands[0], state_constraint) )                    

        #Affinity tasks
        affinity_tasks = list()
        for i in range(num_inputs):
            input_len  = len(self.strands[i + 1])
            #constraint = generate_fold_constraint(gate_len, []) + generate_fold_constraint(input_len, [ ])
            #affinity_tasks.append( (self.strands[0]  + self.strands[i + 1], constraint, gate_len, input_len ) )

            constraint = '.'*len(self.strands[0]) + '&' + '.'*len(self.strands[i+1])
            affinity_tasks.append((self.strands[0] + '&' + self.strands[i+1], constraint, gate_len, input_len))

        self.tasks = [logic_tasks, affinity_tasks]
        self.input_states = input_states #Used for server

 

    def update_base_maps(self): #Requires segments to be formed
        lengths                                     = [seg.length for seg in self.segments]
        (ct_base_seg_all, unused, self.base_to_seg) = build_map(lengths, 1)
        self.name_to_bases                          = dict()
        names                                       = [seg.name for seg in self.segments]
        for i, name in enumerate(names):
            self.name_to_bases[name] = ct_base_seg_all[i]

    def update_mutation_maps(self, types):
        if 'bccs' in types:
            scc_lengths     = [scc.length for scc in self.sccs]
            self.bcc_to_scc = build_map(scc_lengths, 0)[2]

    def update_mutation_probs(self):
        num_bccs = len ( self.bcc_to_scc.keys() )
        weights  = list()
        for i in range(num_bccs):
            (scc_id, bcc_id_loc) = self.bcc_to_scc[i]
            scc                  = self.sccs[scc_id]
            cur_allowed           = scc.allowed[bcc_id_loc]
            if scc.size == 1: # unpaired
                weight = len(cur_allowed) - 1                # -1 because we care about the #of allowed bccs that a bcc can mutate into 
            else:
                weight = ( len(cur_allowed) - 1 ) / float(2) #/2 because each base pair contains 2 bases
            weights.append(weight)

        denom              = sum(weights)
        self.mutation_probs = [w/float(denom) for w in weights]
       

class SCC():
    def __init__(self, names, size, length, constraints_raw, strand_ids, segs_pos, colors):
        self.names         = names
        self.size          = size
        self.length        = length
        if size == 1:
            self.constraints = list(constraints_raw[0])
        else:
            self.constraints = list( zip(list(constraints_raw[0]), list(constraints_raw[1][::-1])) )


        self.segs_strand      = strand_ids
        self.segs_pos         = segs_pos
        self.segs_Colors      = colors #This is a tuple pair. The second elem is empty if there is only one strand

        self.mismatches      = list()
        self.allowed         = list()
        self.bccs            = list()


    def update_allowed(self): #(during initialization, size change, mismatch change)
        self.allowed = list()

        for i in range( self.length ):
            if self.size == 1:                
                cur_allowed = allowed_lut( self.constraints[i] )
            else:
                boolMismatch = i in self.mismatches
                cur_allowed = allowed_lut( (self.constraints[i][0], self.constraints[i][1], boolMismatch) )
            self.allowed.append( cur_allowed )

        #Handle length constraint violations. This is required during initialization
        num_bccs = len(self.bccs)
        new_len  = len(self.allowed)
        if num_bccs < new_len:
            for i in range(num_bccs, new_len):
                self.bccs.append( np.random.choice(self.allowed[i]) )
        if num_bccs > new_len:
            self.bccs = self.bccs[0:new_len]

        #Handle bcc constraint violations
        for i, bcc in enumerate(self.bccs):
            if bcc not in self.allowed[i]:
                self.bccs[i] = np.random.choice(self.allowed[i])
               

class Segment():
    def __init__(self, name, length, strand, pos, color, seq):
        self.name   = name
        self.length = length
        self.strand = strand
        self.pos    = pos
        self.color  = color
        self.seq    = seq
