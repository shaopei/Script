import numpy as np





i_f = 'interestingHets_noWeird.txt'
i_chromosome = np.loadtxt(i_f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #1-22
i_snppos = np.loadtxt(i_f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
i_winning_parent =  np.loadtxt(i_f, dtype=str ,delimiter='\t', usecols=[13], skiprows=1)  #M or P

with open ('Groseq_interestingHets_noWeird_merge1K.bed', 'U') as merged_input:
    with open ('Groseq_interestingHets_noWeird_merge1K_Discordance.bed', 'w') as out:
        out.write(merged_input.readline().strip())
        out.write('\tGroseq.Discordance\n')
        for l in merged_input.readlines():
            out.write(l.strip())
            m_chr, m_chromStart, m_chromEnd = l.strip().split('\t')
            c = int(m_chr.split('chr')[-1])
            c_snppos = i_snppos[i_chromosome == c]
            c_winning_parent = i_winning_parent[i_chromosome == c]
            m_winning_parent = c_winning_parent[np.logical_and((c_snppos - 1 - int(m_chromStart) >= 0 ),(int(m_chromEnd) - c_snppos >= 0))]
            if m_winning_parent.shape[0] <1:
                print "ERROR, there is NO snp in this region"
            elif m_winning_parent.shape[0] == 1:
                out.write('\t'+m_winning_parent[0]+'\n')
            elif len(set(m_winning_parent)) > 1:
                out.write('\tD\n')
            else:
                out.write('\tC\n')



###
#merged_f = 'Groseq_interestingHets_noWeird_merge1K.bed'
merged_f = 'Groseq_interestingHets_noWeird_merge1K_Discordance.bed'
chromosome = np.loadtxt(merged_f, dtype=str ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
chromStart = np.loadtxt(merged_f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
chromEnd = np.loadtxt(merged_f, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
Con_Dis = np.loadtxt(merged_f, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)
regionLength = chromEnd - chromStart

import matplotlib.pyplot as plt
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
plt.show()



