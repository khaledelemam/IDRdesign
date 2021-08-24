import math

def log (x): ##vectorized log function
	if (type(x) is list):
		return [log(val) for val in x]
	elif (x>0):
		return math.log(float(x))
	else:
		return ""	

def mean (x):
	tot =float(0)
	for val in x: tot += float(val)
	return tot/float(len(x))

def stdev (x):
	tot =float(0)
	m=mean(x)
	for val in x: tot += float( val - m )**2
	if ((tot>=float(0)) and (len(x)>1)): 
		return math.sqrt( tot/float(len(x)-1) )
	else:
		return float(0)	

def t_stat (x, y): #Welch's t-statistic
	totxx =float(0)
	totx=float(0)
	for val in x:
		totxx += float(val)**2
		totx+=float(val)
	mx = totx/float(len(x))			
	totyy=float(0)
	toty=float(0)	
	for val in y:
		totyy += float(val)**2	
		toty+=float(val)
	my = toty/float(len(y))	
	var= (totxx/float(len(x)) - mx*mx)/float(len(x)) + (totyy/float(len(y)) - my*my)/float(len(y))
	if (var>=float(0)):	
		return (mx - my)/math.sqrt(var)
	else:
		return 0
			
def paired_t_stat (x, y): #paired t-statistic
	if (len(x) != len(y)):
		print("problem with paired t-stat for ",len(x),len(y))
		return 0
		
	totxx =float(0)
	totx=float(0)
	for val1, val2 in zip(x,y):
		diff = float(val1)-float(val2)
		totxx += diff*diff
		totx += diff
	mx = totx/float(len(x))			
	var= (totxx/float(len(x)) - mx*mx)/float(len(x))
	if (var>float(0)):	
		return mx/math.sqrt(var)
	else:
		return 0

def correl (x,y): #pearson's
	if (len(x) != len(y)):
		print("problem with corr for ",len(x),len(y))
		return 0
	totxx =float(0)
	totx=float(0)
	totyy=float(0)
	toty=float(0)	
	totxy=float(0)	
	for valx,valy in zip(x,y):
		totxx += float(valx)**2
		totx+=float(valx)
		totyy += float(valy)**2	
		toty+=float(valy)
		totxy+=float(valy)*float(valx)
	n=float(len(x))
	mx = totx/n			
	my = toty/n	
	if ((totxx < n*mx*mx) or (totyy < n*my*my)): return 0
	return (totxy - n*mx*my)/( math.sqrt( totxx - n*mx*mx )	* math.sqrt( totyy - n*my*my )	)			