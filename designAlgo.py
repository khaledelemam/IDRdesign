import math
import time
import random  
from helpers import *
from run_features import*
import pickle
from Bio import pairwise2
from Bio.Align import substitution_matrices



def random_sequence(length):  
    sample_string = 'ARNDCQEGHILKMFPSTWYV' # define the specific string   
    # random.seed(10)
    result = ''.join((random.choice(sample_string)) for x in range(length))  
    print("Randomly generated sequence is: ", result) 
    return result


def calc_all_dist(all_seq, arr_main, transformation, dist_dict):
    for s in all_seq:
        arr_rand = calc_feature_vector(s, transformation)
        dist = calculate_dist(arr_main, arr_rand)
        dist_dict[s] = dist
    return dist_dict


def changeResidue(random_seq, all_seqs, index):
    aa_string = 'ARNDCQEGHILKMFPSTWYV'
    aa_list = list(aa_string)
    aa_temp = aa_list
    for _ in range(len(aa_list)):
        if len(aa_temp)!=0:
            residueChange = aa_temp.pop(0)
            random_seq2 = random_seq[:index] + residueChange +random_seq[index+1:]
            all_seqs.append(random_seq2)
    return all_seqs

def calc_feature_vector(seq , transformation):
        seq_mol_feats = run_feats(seq)
        seq_df = parseFeaturesDict(seq_mol_feats)
        # feats= ["id=aliphatic m=[ALMIV]", "isoelectric_point", "my_kappa","KL_hydropathy"]
        # temp = seq_df.loc[:,feats]    
        # idx = [seq_df.columns.get_loc(i) for i in feats]
        # seq_standardized = transform(temp, transformation,idx)
        seq_standardized = transform(seq_df, transformation)
        arr_seq=  seq_standardized.to_numpy()  
        return arr_seq

def align(seq1, seq2):
    blosum62 = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.globaldx(seq1, seq2,blosum62)
  
    return max(alignments)
       

def design_sequences(target, seqs_number, query = ''):
    
    transformation = pickle.load(open('transformation.pkl', 'rb'))
    arr_main = calc_feature_vector(target, transformation)
    designed_seqs = []
    
    while seqs_number > 0:
        t0 = time.time()

        if query =='':
             random_seq = random_sequence(len(target))
        else:
            random_seq = query
            print("Query sequence is: ", random_seq)

        temp = random_seq

        arr_rand = calc_feature_vector(random_seq, transformation)
        dist = calculate_dist(arr_main, arr_rand)
        print("Initial distance: " + str(dist))

        precision = 0.000001 
        previous_step_size = 1
        curr_dis = math.inf
        iter = 0
     
        while previous_step_size > precision:
            
            dist_dict ={}
            all_seqs = []
            

            #change residue 20 times per position the move on to next postion
            for index in range(len(random_seq)):
                all_seqs = changeResidue(random_seq, all_seqs, index)

            dist_dict = calc_all_dist(all_seqs, arr_main, transformation, dist_dict)
            impr_seq = min(dist_dict, key = dist_dict.get)
            random_seq = impr_seq
            prev_dis = curr_dis
            curr_dis = dist_dict[impr_seq][0]
            previous_step_size = abs(curr_dis - prev_dis)

          
            iter +=1
        
        res = min(dist_dict.items(),key = lambda x:x[1])
        designed_seqs.append(res)
       
        alignment_before= align(target, temp)
        alignment_after = align(target,res[0])
      

        t1 = time.time()
        print("Iterations: ",iter)
        print("Time: ",t1-t0)
        print("Alignment before: " , alignment_before)
        print("Alignment after: " , alignment_after )
        print('\n')

        seqs_number = seqs_number -1
    return designed_seqs



# designed = design_sequences('MLFRNIEVGRQAAKLLTRTSSRLAWQSIGASRNISTIRQQIRKTQ', 1, 'DHKLYSVTNKGMADFLAPCPRMTVLEVAPRVQEILCEWCSCIYVY')
# print(designed)


#Generate random sequence from amino acid distribution

# aa = 'ARNDCQEGHILKMFPSTWYVU'
# keys = list(aa)
# aaFreq = {key: 0 for key in keys}
# seqID, seqs , lengths, dct, total= parseFASTA_generator('UP000005640_9606_SPOTD_7res_joined_min_30AA.fasta', aaFreq)
# prob = [(c/total) for c in dct.values()]
# new_count = npr.multinomial(45, prob)
# seq = generate_seq(new_count, keys)
