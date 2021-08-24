from scipy.stats import shapiro, wilcoxon, fisher_exact, ttest_rel
import numpy as np
from helpers import parseFASTA, calculate_dist
from designAlgo import calc_feature_vector, align
import pickle
import csv

### Cox15

#cox15 dist
cox15_random_dists =[3.077619,3.131759,15.93054,2.837519,4.348167,4.243072,2.312528,3.781911,5.047213,3.352399,3.921691,2.526132,3.353345,4.404619,4.201034,2.935876,3.893772,3.353119,4.26132,5.033628]
cox15_designed_dists =[0.737343,0.03157,0.745498,0.208408,0.637559,0.096407,0.412454,0.280466,0.674779,0.130959,0.359936,0.102049,0.209018,0.220069,0.093402,0.672966,0.083446,0.307186,0.692461,0.119536]

#cox15 align
cox15_random_alignPerc = [35.56,31.11,31.11,35.56,31.11,31.11,31.11,28.89,26.67,28.89,31.11,35.56,28.89,33.33,28.89,26.67,31.11,28.89,31.11,28.89]
cox15_designed_alignPerc =[42.22,44.44,35.56,40,46.67,35.56,37.78,35.56,42.22,35.56,44.44,46.67,35.56,40,40,35.56,40,40,40,40]
cox15_homolog_alignPerc = [82.22222222,73.3333333,80,44.4444444,33.3333333,28.8888889,46.6666667,48.8888889,48.8888889,40,42.2222222,40,44.4444444,44.4444444,40,46.6666667,46.6666667,44.4444444]


cox15_random_align62 = [87,81,85,82,78,88,72,72,66,78,79,84,83,85,69,79,85,78,77,76]
cox15_designed_align62 = [99,115,92,89,103,97,90,88,95,91,104,106,90,95,95,97,95,99,93,92]
cox15_homolog_align62 = [184,175,184,109,82,65,108,119,116,101,100,99,106,100,95,111,120,108]

#cox15 mitofates
cox15_mitofates = [[0,11],[20,9]]

#cox15 solubility
cox15_random_solubility = [2.12,0.3,0.51,1.32,-0.057,-2.2,1.32,0.13,-0.65,1.51,1.02,1.57,0.3,0.056,0.4,-1.13,1.29,1.31,0.27,1.62]
cox15_designed_solubility = [1.43,1.49,1.56,1.31,1.91,1.16,1.88,1.24,1.93,1.29,1.75,1.39,1.02,1.91,1.02,1.56,1.59,1.78,1.84,1.81]

#homolog lengths = [45, 45, 45, 45, 55, 40, 34, 48, 58, 58, 49, 50, 51, 48, 48, 47, 49, 49, 49]


#### Ded1-N

#ded1 align
ded1_random_align = [32.95,28.41,29.55,27.27,27.27,31.82,27.27,25,31.82,26.14,26.14,29.55,30.68,23.86,23.86,22.73,27.27,28.41,27.27,27.27]
ded1_designed_align = [45.45,39.77,43.18,45.45,51.14,45.45,45.45,50,43.18,47.73,40.91,47.73,36.36,44.32,38.64,48.86,42.05,38.64,50,44.3]
ded1_homolog_align = [94.3181818,90.9090909,90.9090909,57.9545455,61.3636364,57.9545455,67.0454545,53.4090909,53.4090909,64.7727273,55.6818182,57.9545455,59.0909091,53.4090909,53.4090909,54.5454545,56.8181818,56.8181818]

ded1_random_align62 = [174,158,133,157,148,160,174,157,149,173,140,144,166,156,139,134,136,144,155,155]
ded1_designed_align62 = [240,212,231,223,235,271,245,231,264,218,248,210,245,198,233,210,251,219,206,267]
ded1_homolog_align62 = [467,451,452,305,322,299,346,281,278,335,295,296,312,289,293,300,297,304]


#ded1 dist
ded1_random_dist = [13.83147,13.71706,14.77798,13.83787,15.40684,14.85746,15.73366,13.26091,16.02031,12.38556,14.55626,13.41441,14.05659,15.48634,14.66737,13.43559,13.71289,14.49181,12.18131,14.22911]
ded1_designed_dist =[0.113410321,0.443221708,2.209749543,0.348226247,0.466345987,0.344720351,0.13188713,0.796642876,2.155668503,2.195832734,0.224755706,0.105549912,0.138645605,2.158112294,0.739023872,0.508345711,2.150114397,2.150114397,0.117750441,2.155264813]
ded1_homolog_dist = [3.404236016,4.340661194,3.193576078,18.59639289,9.347133959,12.03798025,12.41796835,16.85220668,7.096840542,7.248142466,11.604065,11.90501795,12.4448682,11.87828835,11.20290797,6.860104643,12.61284546,13.70622864]

#ded1 pscore
ded1_random_pscore= [1.14,1,0.93,0.94,0.8,1.3,0.81,0.74,1.07,1.03,0.95,1.14,0.73,0.85,0.91,0.99,0.85,0.71,0.99,0.97]
ded1_designed_pscore = [8.35,5.22,6.43,6.37,6.48,8.46,7.93,5.75,5.76,6.49,5.99,5.19,5.63,6.75,6.58,5.34,7.16,6.26,6.71,8.65]

#ded1 solubility
ded1_random_solubility =[0.19,1.41,-0.13,0.87,-0.62,0.85,0.35,1.38,0.02,0.7,-0.07,-0.4,-0.23,-0.6,-1.45,1.44,-0.66,0.77,0.32,-0.76]
ded1_designed_solubility= [2.62,2.65,2.61,2.59,2.52,2.53,2.55,2.57,2.54,2.5,2.37,2.6,2.56,2.29,2.51,2.5,2.46,2.48,2.45,2.51]

#homolog lengths [88, 92, 92, 95, 91, 84, 95, 91, 81, 101, 106, 79, 86, 105, 92, 98, 97, 95, 96]



def homolog_dist(out,homolog_fasta,sequence):
    transformation = pickle.load(open('transformation.pkl', 'rb'))
    arr_main = calc_feature_vector(sequence, transformation)
    seqID ,seqs = parseFASTA(homolog_fasta)

    with open(out, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Distance"])
            for s in seqs:
                arr_rand = calc_feature_vector(s,transformation)
                dis = calculate_dist(arr_main, arr_rand)
                writer.writerow(dis)


def homolog_align(out, homolog_fasta, sequence):
    seqID ,seqs = parseFASTA(homolog_fasta)
    with open(out, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Alignment"])
            for s in seqs:
                a = align(s,sequence)
                writer.writerow([a.score])

def fishertest(data):
    #fisher exact t-test for MitoFates
    # table:  Random  Designed
    # Yes      0        11
    # No       20       9

    s, pvalue = fisher_exact(data)
    print(pvalue)

    if pvalue < 0.05:
        print("Significant difference")
    else:
        print("Non-significant difference")




def calc_stats(v1, v2):

    differences = []
    for i,j in zip(v1, v2):
        differences.append(abs(i-j))

    data = np.asarray(differences)

    #test if differences are normally distibuted
    stat, p = shapiro(data)
    print(p)

    
    if p > 0.05: #Gaussian
        #preform paired t-test
        s, pvalue = ttest_rel(v1, v2)
        print("Gaussian, ttest")

    else: #not Gaussian
        #preform Wilcoxon signed-rank test
        s, pvalue = wilcoxon(v1, v2)
        print("Not Gaussian,wilcox") 
    print(pvalue)
    if pvalue < 0.05:
        print("Significant difference")
    else:
        print("Non-significant difference")



# calc_stats(cox15_designed_align62, cox15_random_align62)
# fishertest(cox15_mitofates)


    


#ded1 design dist, ded1 random dist = 3.58e-21 -> significant, ttest
#ded1 design align, ded1 random align  =  2.42e-11 -> significant
#ded1 design align, ded1 homolog align =  7.63 e-6 -> significant, ttest

#ded1 design align62, ded1 random align62  =  4.84e-12 -> significant, ttest
#ded1 design align62, ded1 homolog align62 =  7.63 e-6 -> significant, wilcox

#ded1 design pscore ded1 random pscore = 8.05e-16 -> significant, ttest
#ded1 design solubility, ded1 random solubility = 3.73e-11 -> siginifcant, ttest


#cox15 design dist, cox15 random dist = 1.91e-21 -> significant, wilcox
#cox15 design align, cox15 random align  =  4.82e-11 -> significant
#cox15 design align, cox15 homolog align =  0.035 -> significant, wilcox

#cox15 design align62, cox15 random align62  =  9.61e-9 -> significant, ttest
#cox15 design align62, cox15 homolog align62 =  0.023 -> significant, wilcox


#cox15 design mitofates, cox15 random mitofates = 0.00015 -> significant, fisher
#cox15 design solubility, cox15 random solubility = 0.0003 -> siginifcant, wilcox


#79,49,31,40 -- 83,48,36,44 (h,r,d) cox15
#63,28,44 -- 67,31,48 (h,r,d) ded1