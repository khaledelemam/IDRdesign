import io
import sys,os
import math
import random
import re
import sequence_features as sf
from sequence_features import *
import mol_evol as me
from mol_evol import *
import basic_stats as stats
from basic_stats import *

#DEFAULTS
###the reference species in the file
REF_NUM=0
### sequence quality control heuristics
### minimum number of amino acids in the reference IDR
L_MIN = 1
### how much can the orthologous sequence differ in length from the reference
### L_FACTOR of two means it must be at least half as long and less than twice as long
L_FACTOR=3
### include this sequence if it is D_RATIO further away from the reference than it is from any other sequence
### make this very big to include all input sequences
D_RATIO=5.0
### how much evolutionary distance (in subs. per site) do we want? We need a lot for statistical power
### make this very big to inclue all input sequences
D_TOTAL=30.0
### based on a set of long human IDRs from Iva ###
### needed for F81 model ###
disofreq= {
'A' : 0.0788699046851705,
'C' : 0.0118527921135247,
'D' : 0.0461598457818217,
'E' : 0.0787539937748114,
'F' : 0.0196131563519694,
'G' : 0.078318168751861,
'H' : 0.0228797833856907,
'I' : 0.0242815326616339,
'K' : 0.0554383853661703,
'L' : 0.0729074474562964,
'M' : 0.0176086703421587,
'N' : 0.0303228093095522,
'P' : 0.10630318378937,
'Q' : 0.052863102517791,
'R' : 0.057842120067218,
'S' : 0.12233082931423,
'T' : 0.0587153155919235,
'V' : 0.0446360036803002,
'W' : 0.00627773490505093,
'Y' : 0.0140252201534557
}

# INPUT: instance of the class Alignment (containing all fasta sequences)
def seq_choosing_heuristic_tmp(aln): ##try to get dis_tot, but don't use "redundant" sequences
	seqlist = []
	seqnames = []
	for s in aln.names:
		ugseq=aln.seq[s].replace("-","")
		if ('X' in ugseq): 
			continue #no bad data
		if len(ugseq)==0:
			continue
		seqlist.append(ugseq) # store the actual sequence
		seqnames.append(s) # store the name associated with sequence
	return seqlist,seqnames	

def seq_choosing_heuristic(aln,ref,D_RAT,D_TOT,L_FACT): ##try to get dis_tot, but don't use "redundant" sequences
	sp_dis = {}
	refs=aln.seq[ref].replace("-","") ### reference has to be checked before this is run...
	for s in aln.names:
		if (s is ref): 
			continue
		ugseq=aln.seq[s].replace("-","")
		if ('X' in ugseq): 
			continue #no bad data
		if ( float(L_FACT*len(ugseq)) < float(len(refs)) ): 
			continue #this sequence is too short
		if ( float(len(ugseq)) > float(L_FACT*len(refs)) ): 
			continue #this sequence is too long
		dis = me.fastF81_aa_dis(aln.seq[s],aln.seq[ref],disofreq)
		if ((dis>=float(10)) or (dis<=float(0))): 
			continue #distance looks bad
		sp_dis[dis]=s  ### include this sequence
		
	totd=float(0)
	seqlist = []		
	for d in sorted(list(sp_dis.keys())):
		if (len(seqlist)==0):
			seqlist.append(sp_dis[d])
			totd=totd+d
		elif (totd<D_TOT):
			#eprint(d,sp_dis[d],totd)
			mindis=D_TOT
			for s in seqlist:
				###distance of this sequence to the previous ons added
				dis = me.fastF81_aa_dis(aln.seq[s],aln.seq[sp_dis[d]],disofreq)
				if (dis<mindis): mindis=dis
			#eprint("found something ",mindis," away")
			### include this sequence if it is D_RATIO further away from the reference than it is from any other sequence
			if (D_RAT*mindis > d): 
				seqlist.append(sp_dis[d])
				totd=totd+d
#			else :
#				eprint	(sp_dis[d]+" got filtered because"+str(mindis)+" is too small relative to "+str(d) )
#		else:
#			eprint	(sp_dis[d]+" got filtered because"+str(totd)+" is big enough " )		
	return seqlist			
					

def var_aa_lin_comb(f,aavars):
	aalist=list("ACDEFGHIKLMNPQRSTVWY")
	var=float(0)
	zeros = { aa:float(0) for aa in aalist }
	for aa in aalist:
		onehot = zeros.copy() ##python is crazy!!! if you do onehot = zeros, you modify zeros ?!?!?!?
		onehot[aa] = float(1) 
		var+= aavars[aa] * (f(onehot))**2
		#eprint(aa,onehot,aavars[aa],f(onehot),(f(onehot))**2)
	return var

#IP 1/11/2019
# writes to a file content of featdict
# featdict: nested dictionary of features computed from seq
#			first layer: 
#				keys: sequence ids (i.e. >ENS..x)
#				values: features
#			second layer:
#				keys: feature ids
#				values: feature values

def write_out_features(featdict,outpath,split_each=False):
	if split_each:
		keys=[key for key in featdict.keys()]
		first_entry=keys[0]
		output=outpath.replace(".txt",'').replace('.out','')
		#for entry in featdict[first_entry]:
		#	fout.write("\t%s"%(entry))	
		covered=[]
		for seqid,feats in featdict.items():
			uniprot=seqid.split("_")[0]
			if uniprot not in covered:
				covered.append(uniprot)
				fout=open(output+uniprot[1:]+".out.txt","a")
				for entry in featdict[first_entry]:
					fout.write("\t%s"%(entry))
			
		for seqid,feats in featdict.items():
			uniprot=seqid.split("_")[0]
			fout=open(output+uniprot[1:]+".out.txt","a")
			fout.write("\n%s"%(seqid))
			for featid,featval in feats.items():
				fout.write("\t%s"%(featval))
		
		
	else:
		fout=open(outpath,"w") # open a file on path
		#if not os.path.isfile(outfile):
		keys=[key for key in featdict.keys()]
		first_entry=keys[0]
		
		for entry in featdict[first_entry]:
			fout.write("\t%s"%(entry))	
		for seqid,feats in featdict.items():
			fout.write("\n%s"%(seqid))
			for featid,featval in feats.items():
				fout.write("\t%s"%(featval))
		
		fout.close()	
					
def eprint(*args):
	sys.stderr.write(str(args)+"\n")
		
def expected_matches_PBN (s,m,SM): #compute the expected number of matches to m given SM starting with sequence s
	if (m[-1] == '$'): #special case of C-terminal motifs that 
		newm = m[0:len(m)-1]
		news = s[len(s)-len(newm):]
		return expected_matches_PBN(news,newm,SM)
	#pat = re.compile("".join(m))
	aalist=list("ACDEFGHIKLMNPQRSTVWY")
	mean_pbn=float(0)
	var_pbn=float(0)
	for pos in range(0,len(s)-len(m)+1):
		pospr=float(1)
		for i in range(0,len(m)):
			thispr=float(0)
			if (m[i] == "."):
				thispr=float(1)
			elif ("^" in m[i]):
				#print(m[i])
				b=s[pos+i]
				for a in aalist:
					if (re.search(m[i],a)):
						thispr += SM[b][a]
			else: 
				b=s[pos+i]
				for a in m[i]:
					if (a in aalist):
						thispr += SM[b][a]
			pospr*=thispr
		mean_pbn += pospr
		var_pbn += (pospr - pospr*pospr)
	return( [mean_pbn,var_pbn] )
	
def expected_matches_PBN_multi (s,motifs,SM): #compute the expected number of matches to motifs given SM starting with sequence s
												#motifs is a dict with keys names of motifs and elements patterns by position
	aalist=list("ACDEFGHIKLMNPQRSTVWY")
	motlist=list(motifs.keys())
	mean_pbn= { m:float(0) for m in motlist }
	var_pbn= { m:float(0) for m in motlist }
	
	for pos in range(0,len(s)):
		pospr= { m:float(1) for m in motlist }
		for i in range(0,20):
			if (pos+i >= len(s)): continue
			b=s[pos+i]
			thispr= { m:float(0) for m in motlist }
			for a in aalist:	
				for m in motlist:
					if (i>=len(motifs[m])): continue
					if (motifs[m][i] == "."): continue
					if ("^" in motifs[m][i]):
						if (re.search(motifs[m][i],a)):
							thispr[m] += SM[b][a]
					else: 
						if (a in motifs[m][i]):
							thispr[m] += SM[b][a]
			for m in motlist: pospr[m] *= thispr[m]
		for m in motlist:
			mean_pbn[m] += pospr[m]
			var_pbn[m] += (pospr[m] - pospr[m]*pospr[m])	
	return( [mean_pbn,var_pbn] )	
		
		
def expected_repeat_residues (s,m,SM): #compute the expected number of residues in repeats
	
	aalist=list("ACDEFGHIKLMNPQRSTVWY")
	startpr = [float(0)] #repeat can't start at the first poisition
	for pos in range(1,len(s)-3):
		startpr.append(float(1))
		# ... one non-matching residue before the re[eat
		b=s[pos-1]
		thispr=float(0)
		for a in aalist:
			if (not a in m):
				thispr += SM[b][a]
		startpr[pos]*=thispr 
		for i in range(0,2): # ...and then two matching
			thispr=float(0)
			b=s[pos+i]
			for a in aalist:
				if (a in m):
					thispr += SM[b][a]		
			startpr[pos]*=thispr  
			
	for pos in range(len(s)-3,len(s)): startpr.append(float(0)) #repeat can't start at the last two poisitions
	
	repeatpr = [float(0)] #repeat can't start at the first poisition
	for pos in range(1,len(s)):
		thispr=float(0)
		for a in aalist:
			if (a in m):
				thispr += SM[s[pos]][a]
		repeatpr.append( startpr[pos-1] + repeatpr[pos-1]*thispr )
		##either the repeat starts now and this is the second residue, first residue isn't counted
		## or it continues from the previous posiiton
	#print(startpr)
	#print(repeatpr)			
	return sum(repeatpr)

def read_repeats_from_file (file):
	R = { }
	f = open(file, 'r') #open the file
	mname=""
	for line in f:
		line=line.rstrip("\n\r")
		if ('id=' in line):
			mname=line.replace("\t"," ")
			R[mname] = [] # add this header line to the list
		else:
			vals=line.split("\t")
			R[mname].append(vals[1]) ## these are fixed width motifs
	f.close()	
	return R

def read_motifs_from_file (file):		
	M = { }
	f = open(file, 'r') #open the file
	mname=""
	for line in f:
		line=line.rstrip("\n\r")
		if ('id=' in line):
			mname=line.replace("\t"," ")
			M[mname] = [] # add this header line to the list
		else:
			vals=line.split("\t")
			if ("^" in vals[1]):
				pstring="["
				for aa in "ACDEFGHIKLMNPQRSTVWY":
					if (not aa in vals[1]):
						pstring += aa
				pstring += "]"
				M[mname].append(pstring) ## 20% faster if we don't have to use regular expression matching
				#M[mname].append(vals[1])
			else:
				M[mname].append(vals[1]) ## these are fixed width motifs
	f.close()		
	return M		
	
def var_diff_LRTstat(x,v):
	#mx=mean(x)
	s2x=stdev(x)**2
	return float(len(x))*(math.log(v/s2x) + s2x/v - float(1))

def Poisson_sample(lam): # sample from either a single poisson or a vector of poissons with different means
	if (type(lam) is list):
		x = []
		for val in lam:
			x.append(Poisson_sample(val))
		return x	
	else:
		if (lam<float(15)): # if the mean is bigger than this, poisson is approximately normal
			rn=random.random()
			cumpr=float(0)
			pr=math.exp(-lam)
			for i in range (0,50):
				cumpr+= pr
				if (cumpr>=rn):
					return i
				pr = pr*lam/(i+1)	#trick to avoid calcluating factorials or exponents
				#print(str(i) + " "+str(pr)+" "+str(lam)+" "+str(cumpr))
			eprint("returning 51")	
			return 51
		else:
			return int(round( random.gauss(lam - 0.5, math.sqrt(lam)) ))	# continuity correction	

	
def expected_repeat_residues_F81 (s,m,d,freq): #compute the expected number
#for faster F81 calculations
	ss=float(0)
	for x in freq: ss+=freq[x]*freq[x]
	beta=float(1)/(float(1)-ss)
	
	en=float(0)
	aalist=list("ACDEFGHIKLMNPQRSTVWY")
	SM = { }
	for a in aalist: # precompute the substitution matrix
		SM[a] = { b:fastF81_aa(b,a,d,freq[b],beta) for b in aalist }
	
	startpr = [float(0)] #repeat can't start at the first poisition
	for pos in range(1,len(s)-3):
		startpr.append(float(1))
		# ... one non-matching residue before the re[eat
		b=s[pos-1]
		thispr=float(0)
		for a in aalist:
			if (not a in m):
				thispr += SM[b][a]
		startpr[pos]*=thispr 
		for i in range(0,2): # ...and then two matching
			thispr=float(0)
			b=s[pos+i]
			for a in aalist:
				if (a in m):
					thispr += SM[b][a]		
			startpr[pos]*=thispr  #probability that a repeat could start if one more matches
	repeatpr = [float(0)] #repeat can't start at the first poisition
	for pos in range(1,len(s)-3):
		thispr=float(0)
		for a in aalist:
			if (a in m):
				thispr += SM[s[pos]][a]
		repeatpr.append( startpr[pos-1] + repeatpr[pos-1]*thispr )
		##either the repeat starts now and this is the second residue, first residue isn't counted
		## or it continues from the previous posiiton
	#print(startpr)
	#print(repeatpr)			
	return sum(repeatpr)
	
def expected_matches_F81_PBN (s,m,d,freq): #compute the expected number of matches to m after d subs/site starting with sequence s
	#print (m)
	#pat = re.compile("".join(m))
	aastring="ACDEFGHIKLMNPQRSTVWY"
	aalist = list(aastring)
	
	#for faster F81 calculations
	ss=float(0)
	for x in freq: ss+=freq[x]*freq[x]
	beta=float(1)/(float(1)-ss)
	SM = { }
	for a in aalist: # precompute the substitution matrix
		SM[a] = { b:fastF81_aa(b,a,d,freq[b],beta) for b in aalist }
	
	mean_pbn=float(0)
	var_pbn=float(0)
	for pos in range(0,len(s)-len(m)):
		pospr=float(1)
		
		for i in range(0,len(m)):
			thispr=float(0)
			if (m[i] == "."):
				thispr=float(1)
			elif ("^" in m[i]):
				b=s[pos+i]
				for a in aalist:
					if (re.search(m[i],a)):
						thispr += SM[b][a]
			else: 
				b=s[pos+i]
				for a in aalist:
					if (a in m[i]):
						thispr += SM[b][a]
			pospr*=thispr
		mean_pbn += pospr
		var_pbn += (pospr - pospr*pospr)
	return( [mean_pbn,var_pbn] )
			
			