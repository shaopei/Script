import numpy as np
import matplotlib.pyplot as plt
from sys import argv

#f = argv[1]
#f = '/Users/shaopei/Desktop/part_D/Danko_lab_work/Transcrition_directionality_project/TD_workstation/alleleseq_Groseq_against_Diploid.2015.feb5.snpIndelSV/interestingHets_noWeird.txt'
f = 'interestingHets_noWeird.txt'
#f = 'counts.txt'


chromosome = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #1-22
snppos = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
winning_parent =  np.loadtxt(f, dtype=str ,delimiter='\t', usecols=[13], skiprows=1)  #M or P

# split the file into chr1-22
#chr_snppos_dic = {}
#for i in xrange(1,22):
#    chr_snppos_dic[i] = zip(snppos[chromosome == i], winning_parent[chromosome == i])

def ASE_Concordance(c, d):  #c is the chrom number
    """
    output is a vector (list) [1, -1, ...]  1 is discordant, -1 is concordant
    Discordant: For Asym snppos, within d bp, there are both paternal and maternal Asym.
    """
    #c = 1
    print c
    c_snppos = snppos[chromosome == c]  #the snps on that chrom
    print c_snppos.shape
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

Total_c_concordance = []
for chrom in xrange(1,23):
    Total_c_concordance += ASE_Concordance(chrom, 1000)

c = np.array(Total_c_concordance)

output_list=[]
with open ('Groseq_Asym_Discordant_FrominterestingHets_noWeird_d1K.txt', 'w') as out:
    output_list.append('\t'.join(['chrm', 'snppos' ,'Groseq.winning_parent' , 'Groseq.ASE.Discordance']))
    for i in xrange(snppos.shape[0]):
        output_list.append('\t'.join([str(chromosome[i]),str(snppos[i]), winning_parent[i], str(Total_c_concordance[i])]))
    out.write('\n'.join(output_list))

discordance_map ={-1:"C", 1:"D", 0:"N"}
output_list=[]
with open ('Groseq_Asym_Discordance_FrominterestingHets_noWeird_d1K.bed', 'w') as out:
    output_list.append('\t'.join(['#chrom', 'chromStart' ,'chromEnd','Groseq.ASE.Discordance']))
    for i in xrange(snppos.shape[0]):
        output_list.append('\t'.join(['chr'+str(chromosome[i]),str(snppos[i]-1),str(snppos[i]), discordance_map[Total_c_concordance[i]]]))
    out.write('\n'.join(output_list))


output_list=[]
with open ('Groseq_winning_parent_FromCounts_index0.bed', 'w') as out:
    output_list.append('\t'.join(['#chrom', 'chromStart' ,'chromEnd','Groseq.winning_parent']))
    for i in xrange(snppos.shape[0]):
        output_list.append('\t'.join(['chr'+str(chromosome[i]),str(snppos[i]-1),str(snppos[i]), winning_parent[i]]))
    out.write('\n'.join(output_list))


    
output_list=[]
with open ('/workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/snp.call.all', 'U') as snp:
    with open ('snp.call.all.bed', 'w') as out:
        output_list.append('\t'.join(['#chrom', 'chromStart' ,'chromEnd','genotype']))
        for l in snp.readlines():
            ll = l.strip().split('\t')
            output_list.append('\t'.join(['chr'+ll[0],str(int(ll[1])-1), ll[1], ll[5]]))
        out.write('\n'.join(output_list))



