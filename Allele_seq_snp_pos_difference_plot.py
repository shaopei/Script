import numpy as np
import matplotlib.pyplot as plt

f = '/Users/shaopei/Desktop/part_D/Danko_lab_work/Transcrition_directionality_project/TD_workstation/alleleseq_Groseq_against_Diploid.2015.feb5.snpIndelSV/interestingHets_noWeird.txt'

chromosome = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=[0], skiprows=1) #1-22
snppos = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
winning_parent =  np.loadtxt(f, dtype=str ,delimiter='\t', usecols=[13], skiprows=1)  #M or P

snppos_plus1 = np.zeros(snppos.shape)
snppos_plus1 [1:] = snppos[:-1]
delta_snppos = snppos - snppos_plus1


#plt.hist(delta_snppos[chromosome== '1'])
plt.hist(delta_snppos[abs(delta_snppos) <= 10000])
plt.title(f.split('/')[-1]+'_delta_snppos')
plt.savefig(f.split('/')[-1]+'_delta_snppos_less_than10K.png')
#plt.show()

plt.hist(delta_snppos[abs(delta_snppos) <= 1000])
plt.title(f.split('/')[-1]+'_delta_snppos')
plt.savefig(f.split('/')[-1]+'_delta_snppos_less_than1K.png')



#from sys import argv
#f = argv[1]
#
#for i in range(10,23):
#    print  'sed -i \'s/'+`i`+'_paternal/chr'+`i`+'/g\' '+f
#
#for i in range(1,10):
#    print  'sed -i \'s/'+`i`+'_paternal/chr'+`i`+'/g\' '+f
