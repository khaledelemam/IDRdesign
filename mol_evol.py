import math
import random

class Alignment:
	def __init__(self):
		self.names = [] #list of names of sequences in the alignment
		self.seq = {} #the actual sequence indexed by their names
		self.info = "" #any other information 	
	
	def read_json(self,file): #read from a json file	
		f = open(file, 'r') #open the file
		lines=[]
		foundseq=0
		for line in f:
			line=line.rstrip("\n\r")
			#print (lines)
			if ('}' in line) :
				#print("lines:",lines)
				seqname=""
				for l in lines:
					if ('species' in l):
						vals=l.split(": ")
						#print(vals)
						seqname=vals[-1].replace("\"","").replace(",","")
						if (not seqname in self.names): self.names.append(seqname) # add this name to the list if it isn't there
						##print ("found"+seqname)
						##(apparently json can have more than one alignment)
				for l in lines:	
					if ('seq\":' in l):
						vals=l.split(": ")
						if (not seqname in self.seq): ##this==the first bit of sequence for this species
							self.seq[seqname] = vals[-1].replace("\"","").replace(",","")
						else:
							self.seq[seqname] += vals[-1].replace("\"","").replace(",","")	
				lines=[]
				foundseq=0
			elif ('{' in line) :
				foundseq=1
			elif(foundseq==1) :
				lines.append(line);	
		f.close()		
					
				
	def read_mfa(self,file): #read from a multiple fasta file
		f = open(file, 'r') #open the file
		for line in f:
			line=line.rstrip("\n\r")
			if ('>' in line):
				self.names.append(line) # add this header line to the list
			else:
				seqname=self.names[-1] # this==the most recently found header line
				if (not seqname in self.seq):
					self.seq[seqname] = line
				else:
					self.seq[seqname] += line		
		f.close()
	# def read_seqList(self,seq_list):
	# 	for seq in seq_list:
			
		
	def read_seq(self, seq):
		self.names.append("")
		seqname=self.names[-1]
		self.seq[seqname] = seq

		
	def print_mfa(self): #print fasta format to the screen
		for s in self.names:
			if ('>' in s):
				 print(s)
			else:
				print ('>' + s)
			print(self.seq[s])
			
	def print_maf(self): #print maf format to the screen
						#assumes the score line has been stored as the info for the alignment
		if ('a score' in self.info): 
			print(self.info)
		else:
			print('a score= '+self.info)
		for s in self.names:
			print ("s "),;print(s+"\t"),;print(self.seq[s]+"\t\t\t")
		print ('\n')
	
	def sub(self,beg,end,slist): #return a sub alignment
		#print('trying to create an alignment from ',beg,'to',end,' for ',len(slist), ' of ',len(self.names))
		SA=Alignment()
		for s in slist:
		#	print(s)
		#	print (len(self.seq[s]),end) 
			if ((len(self.seq[s])>=end) and (s in self.names)): 
				SA.seq[s]=self.seq[s][beg:end]
				SA.names.append(s)
		#print (len(SA.names))		
		##remove any columns that are entirely gaps	
		#SA.print_mfa()	
		toremove=[]
		for i in range(0,len(SA.seq[SA.names[0]])):
			if (SA.isallgapped(i)==1): toremove.append(i)
		for pos in reversed(toremove):
			for s in SA.names:	SA.seq[s] = SA.seq[s][:pos] + SA.seq[s][pos+1:]
		SA.info = self.info + "_" + str(beg) + "_" + str(beg+end)	
		#print(toremove)	
		#SA.print_mfa()			
		return SA

	def isungapped(self,pos): #==this position ungapped
		for s in self.names:
			if (len(self.seq[s])<=pos): return 0
			if (self.seq[s][pos]=='-'): return 0
		return 1
				   
	def isallgapped(self,pos): #==this position all gaps
		for s in self.names:
			if (len(self.seq[s])<=pos): return 0
			if (not self.seq[s][pos]=='-'): return 0
		return 1		   	
		

def P_JC69_aa(A,B,d): # probability of changing from B to A in evolutionary distnace d
	if (A==B):
		return 0.05*(1 + float(19)*math.exp(-float(20)*d/float(19)))
	else:
		return 0.05*(1 - math.exp(-float(20)*d/float(19)))	

def JC69_nt_dis(s1,s2): #distance between two aligned sequences
	nid=0
	ndiff=0
	ratio = float(3)/float(4)
	for pos in range(0,len(s1)):
		if ('-' in s1[pos]+s2[pos]):
			pass
		elif (s1[pos]==s2[pos]):
			nid+=1
		else:
			ndiff+=1
	if (nid==0):
		return float(10)
	elif (ndiff==0):
		return float(0)
	else:				
		frac = float(ndiff)/float(nid+ndiff)
	
		if (frac>=ratio):
			return float(10)
		else:
			return -ratio*math.log(float(1) - frac/ratio)

def JC69_aa_dis(s1,s2): #distance between two aligned sequences
	nid=0
	ndiff=0
	ratio = float(19)/float(20)
	for pos in range(0,len(s1)):
		if ('-' in s1[pos]+s2[pos]):
			pass
		elif (s1[pos]==s2[pos]):
			nid+=1
		else:
			ndiff+=1
	if (nid==0):
		return float(10)
	elif (ndiff==0):
		return float(0)
	else:				
		frac = float(ndiff)/float(nid+ndiff)
	
		if (frac>=ratio):
			return float(10)
		else:
			return -ratio*math.log(float(1) - frac/ratio)

def P_F81_aa(A,B,d,freq):
	ss=float(0)
	for x in freq:
		#print(x)
		ss+=freq[x]*freq[x]
	beta=float(1)/(float(1)-ss)
	#BETA - constant : 1.075 -> ss=0.07
	if (A==B):
		return math.exp(-beta*d) + freq[A]*(1 - math.exp(-beta*d))
	else:
		return freq[A]*(1 - math.exp(-beta*d))

def fastF81_aa(A,B,d,fA,beta):
	if (A==B):
		return math.exp(-beta*d) + fA*(1 - math.exp(-beta*d))
	else:
		return fA*(1 - math.exp(-beta*d))		

def fastF81_aa_dis(s1,s2,freq): #can't find a trusted reference for this, but seems to be right... same as JC, but beta=3/4
	ss=float(0) 
	for x in freq.values(): 
		ss += x*x 
	beta=float(1)/(float(1)-ss)
	nid=0
	ndiff=0
	alnpairs = zip(list(s1),list(s2))
	for aa1, aa2 in alnpairs:
		if ((aa1=='-') or (aa2=='-')): 
			pass
		elif (aa1==aa2):
			nid+=1
		else:
			ndiff+=1
	if (nid==0):
		return float(10)
	elif (ndiff==0):
		return float(0)
	frac = float(ndiff)/float(nid+ndiff)	
	if (frac*beta>float(1)):
		return float(10)	
	else:				
		frac = float(ndiff)/float(nid+ndiff)
		return -math.log(float(1) - frac*beta)/beta
		
def F81_aa_dis(s1,s2,freq): #distance between two aligned sequences, assuming s1 evolved from s2
	initdis=JC69_aa_dis(s1,s2) ## use JC69 as an initial guess
	#brute force 1D ML estimation
	N=150 #how many values to try
	
	ss=float(0) ###for faster f81 calculations
	for x in freq.values(): ss += x*x 
	beta=float(1)/(float(1)-ss)
	#print (ss)
	#print (sum([ val*val for val in freq.values()]))
	alnpairs = zip(list(s1),list(s2))
	
	ilogL=float(0)		
	for aa1, aa2 in alnpairs:
		if ((aa1=='-') or (aa2=='-')): continue
		#ilogL += math.log(P_F81_aa(s1[pos],s2[pos],initdis,freq))
		ilogL += math.log(fastF81_aa(aa1,aa2,initdis,freq[aa1],beta))
	#print(ilogL)
	
	bestd=initdis
	bestL=ilogL
	for i in range(0,N):
		d=initdis + 0.01*initdis*(float(N)/float(2) - float(i))
		logL=0
		for aa1, aa2 in alnpairs:
			if ((aa1=='-') or (aa2=='-')): continue		
			logL += math.log(fastF81_aa(aa1,aa2,d,freq[aa1],beta))	
		#print (d,i,bestd),;print("\t"),;
		#print (logL-ilogL)	
		if(logL>bestL):
			bestL=logL
			bestd=d
			#if(bestd>float(9)):
			#	eprint(str(i)+": logL="+str(logL)+" F81d="+str(d))
	return bestd	

def sim_multi_prot_given_anc(seq,d,freq,n): ##simulate n sequences from seq
											##still could use an object approach to avoid 
											##passing back n sequences
	SM = { }
	aalist = list("ACDEFGHIKLMNPQRSTVWY")
	#for faster F81 calculations
	ss=float(0)
	for x in freq: ss+=freq[x]*freq[x]
	beta=float(1)/(float(1)-ss)
	for a in aalist: # precompute the substitution matrix
		SM[a] = [fastF81_aa(b,a,d,freq[b],beta) for b in aalist]
	newseq= []
	seqlist=list(seq)
	for i in range(0,n):
		newseq.append( "".join( [ aalist[random_int(SM[b])] for b in seqlist ] ) )
	### try the dreaded 2D list construction	not faster in this case ... 
	#RES = [ "".join([aalist[random_int(SM[b])] for b in seqlist ])  for i in range(0,n) ]
	return(newseq)
				

def random_int(pr):
	cumpr=float(0)
	rn = random.random()
	i=int(0);
	for p in pr:
		cumpr += p
		i+=1
		if (cumpr>=rn): return (i-1)
	#print(cumpr)				
						
def sim_prot_given_anc(seq,d,freq):
	newseq=[]
	aalist = list("ACDEFGHIKLMNPQRSTVWY")
	#for faster F81 calculations
	ss=float(0)
	for x in freq:
		#print(x)
		ss+=freq[x]*freq[x]
	beta=float(1)/(float(1)-ss)
	
	for b in list(seq):
		thispr=float(0)
		rn=random.random()
		for a in aalist:
			thispr += fastF81_aa(a,b,d,freq[a],beta)
			if (thispr>=rn):
				newseq.append(a)
				break;
	return("".join(newseq))			
		