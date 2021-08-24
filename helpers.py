import math
import pandas as pd
import random
import csv



def parseFeaturesDict(featdict):
    IDs=[key for key in featdict.keys()]
    first_entry=IDs[0]
    featureNames = [entry for entry in featdict[first_entry]]
    values = []
    all_values =[]
    output = []
    for seqid,feats in featdict.items():
        values = []
        for featid,featval in feats.items():
            values.append(featval)
        all_values.append(values)
    df = pd.DataFrame(all_values, columns=featureNames,index = IDs)
    return df
  


def eucdis(v1, v2):
    '''
    Takes in 2 feature vectors
    Returns the euclidean distance

    '''
    dist = [(a - b)**2 for a, b in zip(v1, v2)]
    dist = math.sqrt(sum(dist))
    return dist


def z_score(df):
    '''
    Takes in a dataframe of features and returns a list of 2 lists
    List 1: normalized dataframe
    List 2: mean and std dev used for tranformation
    '''
    # copy the dataframe
    # df_std = df.copy()
    transformation = []
    # apply the z-score method
    for column in df.columns:
        # print( df_std[column].std() )
        transformation.append((df[column].mean(), df[column].std()))
        if df[column].std() != 0: #don't divide by 0
            df[column] = (df[column] - df[column].mean()) / df[column].std()
 
    return df, transformation


def transform(df, transformation, idx =[]):
    '''
    Takes in dataframe of one sequence and tranformation vector (mean,std dev)
    Returns normalized dataframe
    '''
    # copy the dataframe
    # df_std = df.copy()

    ## search for sepcific features
    if idx != []:
        # print("HERE")
        i= idx.pop(0)
        for column in df.columns:
            # print(i)
            if transformation[i][1] != 0: #don't divide by 0
                df[column] = (df[column] - transformation[i][0]) / transformation[i][1]
            if len(idx) != 0:
                i = idx.pop(0)
            else:
                break

    else:
        i= 0
        for column in df.columns:
            if transformation[i][1] != 0: #don't divide by 0
                df[column] = (df[column] - transformation[i][0]) / transformation[i][1]
            i+=1
    return df


def calculate_dist(arr_db, arr_query):
    '''
    Takes in 2 feature vectors
    Returns list with euclidean distances
    '''
    all_dist = []
    for i in range(len(arr_db)):
        dist = eucdis(arr_query[0],arr_db[i])
        all_dist.append(dist)
    return all_dist


def parseFASTA_generator(file, dct):
 
    seqID =[]
    seqs=[]
    lengths = []
    total = 0
    with open(file) as f:
        line = f.readline()
        seqID.append(line.strip())
        # print(line)
        count = 1
        while line:
            line = f.readline()
            if line.strip():  # ignore blank lines
                if (count % 2 == 0): #ID
                    seqID.append(line.strip())
                else: # seq
                    seq = line.strip()
                    # print(seq)
                    length = len(seq)
                    lengths.append(length)
                    # print(list(seq))
                    for i in seq:
                        if i != "X":
                            dct[i] += 1
                            total+=1
                    # print(length)
                    seqs.append(seq)
            count +=1
        return seqID, seqs, lengths, dct, total

def generate_seq(new_count, keys):

    randStr = ''
    for count,letter in zip(new_count, keys):
        for i in range(count):
            randStr += letter
    randLst = list(randStr)
    random.shuffle(randLst)
    res = ''.join(randLst)
    return res


def parseFASTA(file):
    '''
    Takes in a fasta file
    Returns  2 lists
    List 1: protein seq ID
    List 2: protein sequences

    '''
    seqID =[]
    seqs=[]
    with open(file) as f:
        line = f.readline()
        seqID.append(line.strip())
        count = 1
        while line:
            line = f.readline()
            if line.strip():  # ignore blank lines
                if (count % 2 == 0): #ID
                    seqID.append(line.strip())
                else: # seq
                    seqs.append(line.strip())
            count +=1
        return seqID, seqs

def write_csv(dists):
    with open('design_stats.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Iterations", "Distance"])
        for i in range(len(dists)):
            writer.writerow([i+1, dists[i]])




