import numpy as np
#import matplotlib.pyplot as plt
#from sys import argv

f = 'YY1_interestingHets.txt'



# split the file into chr1-22
#chr_snppos_dic = {}
#for i in xrange(1,22):
#    chr_snppos_dic[i] = zip(snppos[chromosome == i], winning_parent[chromosome == i])

def ASE_Concordance(c, d, chromosome, snppos,winning_parent ):  #c is the chrom number
    """
    output is a vector (list) [1, -1, ...]  1 is discordant, -1 is concordant
    Discordant: For Asym snppos, within d bp, there are both paternal and maternal Asym.
    """
    #c = 1
    #print c
    c_snppos = snppos[chromosome == c]  #the snps on that chrom
    #print c_snppos.shape
    c_winning_parent = winning_parent[chromosome == c]
    c_concordance = np.full(c_snppos.shape[0], 0, int)  #set initial to "NA"
    for i in range(c_snppos.shape[0]): 
        #print i
        s = (c_snppos - c_snppos[i] <= d) & (c_snppos - c_snppos[i] >= 0)
        #print c_winning_parent[s]
        if c_winning_parent[s].shape[0] > 1: # have more than one ASE snp in 1Kb (d=1000)
            c_concordance[s] = -1  # Concordance cluster
    
    for i in range(c_snppos.shape[0]): 
        #print i
        s = (c_snppos - c_snppos[i] <= d) & (c_snppos - c_snppos[i] >= 0)
        #print c_winning_parent[s]
        if np.unique(c_winning_parent[s]).shape[0] > 1:  # have more than one parent in 1000bp (d=1000) ex:PPMM, PMPPP,
                c_concordance[s] = 1 # Discordance cluster
    return list(c_concordance)


def TF_ASB_Discordance(f):
    chromosome = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #1-22
    snppos = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
    winning_parent =  np.loadtxt(f, dtype=str ,delimiter='\t', usecols=[13], skiprows=1)  #M or P
    Total_c_concordance = []
    for chrom in xrange(1,23):
        Total_c_concordance += ASE_Concordance(chrom, 1000,chromosome, snppos,winning_parent )
    
    c = np.array(Total_c_concordance)
    
    output_list=[]
    tf = f.split('_')[0]
    with open ('TF_ASB_Discordant_FrominterestingHets_d1K/'+tf+'_ASB_Discordant_FrominterestingHets_d1K.txt', 'w') as out:
        output_list.append('\t'.join(['chrm', 'snppos' ,'TF.winning_parent' , 'TF.ASB.Discordance']))
        for i in xrange(snppos.shape[0]):
            output_list.append('\t'.join([str(chromosome[i]),str(snppos[i]), winning_parent[i], str(Total_c_concordance[i])]))
        out.write('\n'.join(output_list))

import os
import glob
file_names = glob.glob('*_interestingHets.txt')
file_names.sort()
for f in file_names:
    print f
    TF_ASB_Discordance(f)