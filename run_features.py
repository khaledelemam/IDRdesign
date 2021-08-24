from multi_ConDens import *
from sequence_features import *
import sys,os
from helpers import*
import pickle


def run_feats (input_seq):
    """
    # Initialize data structures
    """
    ### these features are all functions of the amino acid frequencies    
    aafeats = [sf.net_charge,    sf.WF_complexity, sf.KL_hydropathy, sf.isoelectric_point, sf.frac_charged_residues, sf.ED_ratio, sf.RK_ratio]
    aafeatnames = ["net_charge",    "WF_complexity", "KL_hydropathy", "isoelectric_point", "FCR","ED_ratio", "RK_ratio"]

    lcs = [sf.net_charge, sf.KL_hydropathy, sf.frac_charged_residues]
    lcnames = ["net_charge","KL_hydropathy", "FCR"]

    ### these features have to be calculated from the whole sequence
    seqfeats = [sf.sequence_charge_decoration, sf.my_kappa, sf.my_omega]
    seqfeatnames = ["SCD", "my_kappa", "my_omega"]

    """
    # Read in data from files
    """
    ### repeats from a file
    if not os.path.isfile("repeats.pkl"):
        REP = read_repeats_from_file ('repeats.txt')
        pickle.dump(REP, open( "repeats.pkl", "wb" ))
    # print (str(len(list(REP.keys())))+' repeats read in '+'repeats.txt');        
    else:
        REP = pickle.load(open('repeats.pkl', 'rb'))

    ### motifs from a file
    if not os.path.isfile("u_fixed_w_feat.pkl"):
        MOT = read_motifs_from_file ('u_fixed_w_feat.txt')
        pickle.dump(MOT, open( "u_fixed_w_feat.pkl", "wb" ))
    else:
        MOT = pickle.load(open('u_fixed_w_feat.pkl', 'rb'))

    featnamesorted = sorted( aafeatnames ) + sorted( list(REP.keys()) + list(MOT.keys()) )

    mol_feats = {}

    aastring="ACDEFGHIKLMNPQRSTVWY"

    def calc_feats(sequences_to_use, names_to_use):
        #start new data structures for each file
        mol_feats = {}
        predvar = {}
        obs={}
        
        for m in list(MOT.keys()) + lcnames:
            predvar[m]=[]
        for f in featnamesorted:
            obs[f]=[]


        for x in range(len(sequences_to_use)):
            ugseq=sequences_to_use[x]  
            nseq=names_to_use[x]   
            mol_feats.setdefault(nseq,{})
            for r in list(REP.keys()):
                pat = re.compile(REP[r][0]+"{2,}")
                all_m=[] #IP added 1/11/2019
                repeatlength=0
                for m in pat.finditer(ugseq):
                    all_m.append(m)
                    repeatlength += (len(m.group(0)) - 1) #first residue in the repeat isn't counted 
                obs[r].append(float(repeatlength))
                mol_feats[nseq][r]=float(repeatlength)
            
            
            for m in list(MOT.keys()):
                pat = re.compile("".join(MOT[m]))
                mol_feats[nseq][m]=float(len(pat.findall(ugseq)))
                obs[m].append(float(len(pat.findall(ugseq))))
            
            
            obsaa = { aa : ugseq.count(aa) for aa in aastring }
            
            #IP 1/11/2019: iterates over a set of functions stored in a list aafeats
            #aafeats[i](obsaa) is a call to a function under the index i
            #dict obsaa is the argument passed to each function
            for i in range(0,len(aafeats)):
                mol_feats[nseq][aafeatnames[i]]=aafeats[i](obsaa)
                if aafeats[i](obsaa)==None:
                    print("SEQ: %s"%ugseq)
                obs[aafeatnames[i]].append(aafeats[i](obsaa))
            
            """
            IP 1/11/2019 -- Alan does not compute this for Z-scores (?)
            -- added here based on the list seqfeats
            """
            for i in range(0,len(seqfeats)):
                mol_feats[nseq][seqfeatnames[i]]=seqfeats[i](ugseq)
        return mol_feats

    isDirectory = os.path.isdir(input_seq)
    isFile = os.path.isfile(input_seq)

    if isDirectory:
        for mfa in os.listdir(input_seq): #mfa in sys.argv[3:]:

            output_file=input_seq+mfa[:-4]+"_FEAT.out.txt"
            ALN = me.Alignment()
            ALN.read_mfa(input_seq+"/"+mfa)          

            print(str(len(ALN.names))+' sequences read in '+mfa)
            sequences_to_use,names_to_use = seq_choosing_heuristic_tmp(ALN) ### this should do all the quality control 
            print(str(len(sequences_to_use))+' sequences passed filtering') 

            mol_feats = calc_feats(sequences_to_use, names_to_use)
            
            ## IF SPLITTING EACH OUTPUT BY UNIPROT ID, add SPLIT=True
            write_out_features(mol_feats,output_file,split_each=False)
            
            print("Done with all sequences")
            print("Written molecular features out to >> %s"%(output_file))
            print("Exiting cleanly...")
            
    elif isFile:
        ALN = me.Alignment()
        ALN.read_mfa(input_seq)          
        sequences_to_use,names_to_use = seq_choosing_heuristic_tmp(ALN) ### this should do all the quality control 
        print(str(len(sequences_to_use))+' sequences passed filtering') 
        mol_feats = calc_feats(sequences_to_use, names_to_use)

    else:

        ALN = me.Alignment()
        ALN.read_seq(input_seq)
        sequences_to_use,names_to_use = seq_choosing_heuristic_tmp(ALN)
        mol_feats = calc_feats(sequences_to_use, names_to_use)
     
    return mol_feats
