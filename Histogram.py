import numpy as np
f = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.counts_TSS_nonpair_Present_PeakCallFiltered.bed'

chromosome = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
chromStart = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
chromEnd = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
Con_Dis = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)
Tss_plus = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[4], skiprows=1)
Tss_minus = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[5], skiprows=1)

regionLength = chromEnd - chromStart
Tss_count = Tss_plus + Tss_minus

import matplotlib.pyplot as plt

plt.hist(Tss_count,bins=range(20))
plt.show()








plt.hist(regionLength[np.logical_and(Con_Dis =='C' ,  regionLength>100)], bins=range(1,5000,100))
plt.title('regionLength >100, C')
plt.xlabel('regionLength')
plt.ylabel('freq')
plt.savefig('regionLength_C.hist.png')
plt.hist(regionLength[np.logical_and(Con_Dis =='D' ,  regionLength>100)], bins=range(1,5000,100))
plt.title('regionLength >100, D')
plt.xlabel('regionLength')
plt.ylabel('freq')
plt.savefig('regionLength_D.hist.png')


