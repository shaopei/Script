import numpy as np
######test if TF ASB are enriched in Groseq Discordance regions
f = 'd.SymCals.all.txt'
#f = 'd.SymCals.USF2.txt'
f_int = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.counts_100plus.bed'
f_end = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF_FishersExactTest_100plus.txt'
TF_index = xrange(7,125)
TF_index_3 = xrange(4,122)

#data = [l.strip().split('\t') for l in open(f, "U").readlines()]
chromosome = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
snppos = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
TF = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=TF_index, skiprows=1)



with open ('Groseq_interestingHets_noWeird_merge1K_Discordance.bed', 'U') as merged_input:
    with open (f_int, 'w') as out:
        out.write(merged_input.readline().strip())
        for i in open(f,'U').readline().strip().split('\t')[7:125]:
            out.write('\t')
            out.write('.'.join(i.strip('"').split('.')[1:]))
        out.write('\n')
        for l in merged_input.readlines():
            m_chr, m_chromStart, m_chromEnd = l.strip().split('\t')[0:3]
            if int(m_chromEnd)- int(m_chromStart) > 100:
                out.write(l.strip())
                c = int(m_chr.split('chr')[-1])
                c_snppos = snppos[chromosome == c]
                c_sym = TF[chromosome == c]
                m_sym = c_sym[np.logical_and((c_snppos - 1 - int(m_chromStart) >= 0 ),(int(m_chromEnd) - c_snppos >= 0))]
                for i in xrange(m_sym.shape[1]):
                    temp = set(m_sym[:,i])
                    if '"Asym"' in temp:
                        out.write('\tAsym')
                    elif '"Sym"' in temp:
                        out.write('\tSym')
                    else:
                        out.write('\tNA')
                out.write('\n')


#f_end = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF_FishersExactTest.txt'
#TF_index_3 = xrange(4,122)

#from fisher import pvalue
import scipy.stats as stats
Merged_TF_SCls  = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=TF_index_3, skiprows=1)
Merged_Concordance= np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)
TF_names = open(f_int, "U").readline().strip().split('\t')[4:122]


D_Asym = (Merged_TF_SCls[Merged_Concordance == 'D'] =='Asym').sum(axis=0)
D_Sym = (Merged_TF_SCls[Merged_Concordance == 'D'] =='Sym').sum(axis=0)
C_Asym = (Merged_TF_SCls[Merged_Concordance == 'C'] =='Asym').sum(axis=0)
C_Sym = (Merged_TF_SCls[Merged_Concordance == 'C'] =='Sym').sum(axis=0)
C_NA = (Merged_TF_SCls[Merged_Concordance == 'C'] =='NA').sum(axis=0)
D_NA = (Merged_TF_SCls[Merged_Concordance == 'D'] =='NA').sum(axis=0)

with open (f_end, 'w') as out:
    #out.write('\t'.join(['TF_names', 'p.left_tail', 'p.right_tail', 'p.two_tail','s_oddsratio', 's_pvalue','C_Asym', 'C_Sym', 'D_Asym', 'D_Sym']))
    out.write('\t'.join(['TF_names', 'oddsratio', 'pvalue_twotail','D_Asym', 'D_Sym', 'C_Asym', 'C_Sym','NA_oddsratio', 'NA_pvalue_twotail','C_NA', 'D_NA']))
    out.write('\n')
    for tf in range(len(TF_names)):
        #p = pvalue(C_Sym[tf], C_Asym[tf], D_Sym[tf], D_Asym[tf])  
        s_oddsratio, s_pvalue = stats.fisher_exact([[D_Asym[tf], D_Sym[tf]],[C_Asym[tf], C_Sym[tf]]])
        na_oddsratio, na_pvalue = stats.fisher_exact([[D_Asym[tf], D_NA[tf]],[C_Asym[tf], C_NA[tf]]])
        #print TF_names[tf], p.left_tail, p.right_tail, p.two_tail
        #out.write('\t'.join([TF_names[tf], str(p.left_tail), str(p.right_tail), str(p.two_tail),str(s_oddsratio), str(s_pvalue),str(C_Sym[tf]), str(C_Asym[tf]), str(D_Sym[tf]), str(D_Asym[tf])]))
        out.write('\t'.join([TF_names[tf],str(s_oddsratio), str(s_pvalue),str(D_Asym[tf]), str(D_Sym[tf]), str(C_Asym[tf]), str(C_Sym[tf]), str(na_oddsratio), str(na_pvalue),str(C_NA[tf]), str(D_NA[tf])]))
        out.write('\n')

######test if TF ASB are enriched in Groseq ASE region?
import numpy as np
import scipy.stats as stats
f = 'd.SymCals.all.txt'
f_out = 'Groseq.Asym_TF.Asym_FishersExactTest.txt'
TF_SCls = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=xrange(7,125), skiprows=1)
Groseq_Scls = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=[58], skiprows=1)
TF_names = open(f, "U").readline().strip().split('\t')[7:125]


GAsym_TFAsym = (TF_SCls[Groseq_Scls == '"Asym"'] =='"Asym"').sum(axis=0)
GAsym_TFSym = (TF_SCls[Groseq_Scls == '"Asym"'] =='"Sym"').sum(axis=0)
GSym_TFAsym = (TF_SCls[Groseq_Scls== '"Sym"'] =='"Asym"').sum(axis=0)
GSym_TFSym = (TF_SCls[Groseq_Scls == '"Sym"'] =='"Sym"').sum(axis=0)

import scipy.stats as stats
with open (f_out, 'w') as out:
    #out.write('\t'.join(['TF_names', 'p.left_tail', 'p.right_tail', 'p.two_tail','s_oddsratio', 's_pvalue','C_Asym', 'C_Sym', 'D_Asym', 'D_Sym']))
    out.write('\t'.join(['TF_names', 'oddsratio', 'pvalue_twotail','GAsym_TFAsym','GSym_TFAsym','GAsym_TFSym','GSym_TFSym']))
    out.write('\n')
    for tf in range(len(TF_names)):
        s_oddsratio, s_pvalue = stats.fisher_exact([[GAsym_TFAsym[tf],GSym_TFAsym[tf]],[GAsym_TFSym[tf],GSym_TFSym[tf]]])
        out.write('\t'.join([TF_names[tf],str(s_oddsratio), str(s_pvalue),str(GAsym_TFAsym[tf]),str(GSym_TFAsym[tf]),str(GAsym_TFSym[tf]),str(GSym_TFSym[tf])]))
        out.write('\n')
        


######test if TF ASB are enriched in Groseq Discordance snps (NOT regions)
import numpy as np
import scipy.stats as stats
#f = 'd.SymCals.all.with.Groseq.ASE.Discordance.txt'
f = 'd.SymCals.all_Groseq_Asym_Discordant_FrominterestingHets_noWeird_d1K.txt'
f_out = 'Groseq.Discordance.snp_TF.Asym_FishersExactTest.txt'
TF_SCls = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=xrange(7,125), skiprows=1)
Groseq_Scls = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=[244], skiprows=1)
TF_names = open(f, "U").readline().strip().split('\t')[7:125]


G_D_TFAsym = (TF_SCls[Groseq_Scls == '1'] =='"Asym"').sum(axis=0)
G_D_TFSym = (TF_SCls[Groseq_Scls == '1'] =='"Sym"').sum(axis=0)
G_C_TFAsym = (TF_SCls[Groseq_Scls== '-1'] =='"Asym"').sum(axis=0)
G_C_TFSym = (TF_SCls[Groseq_Scls == '-1'] =='"Sym"').sum(axis=0)

import scipy.stats as stats
with open (f_out, 'w') as out:
    #out.write('\t'.join(['TF_names', 'p.left_tail', 'p.right_tail', 'p.two_tail','s_oddsratio', 's_pvalue','C_Asym', 'C_Sym', 'D_Asym', 'D_Sym']))
    out.write('\t'.join(['TF_names', 'oddsratio', 'pvalue_twotail','G_D_TFAsym','G_D_TFSym','G_C_TFAsym','G_C_TFSym']))
    out.write('\n')
    for tf in range(len(TF_names)):
        s_oddsratio, s_pvalue = stats.fisher_exact([[G_D_TFAsym[tf],G_D_TFSym[tf]],[G_C_TFAsym[tf],G_C_TFSym[tf]]])
        out.write('\t'.join([TF_names[tf],str(s_oddsratio), str(s_pvalue),str(G_D_TFAsym[tf]),str(G_D_TFSym[tf]),str(G_C_TFAsym[tf]),str(G_C_TFSym[tf])]))
        out.write('\n')


######test if TSS pairs are enriched in Groseq Discordance regions
import numpy as np

f_plus = 'tss_paired_gm12878_plus.bed'
Tss_plus_chromosome = np.loadtxt(f_plus, dtype=str ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
Tss_plus_chromStart_End= np.loadtxt(f_plus, dtype=int ,delimiter='\t', usecols=[1,2], skiprows=1)
f_minus = 'tss_paired_gm12878_minus.bed'
Tss_minus_chromosome = np.loadtxt(f_minus, dtype=str ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
Tss_minus_chromStart_End= np.loadtxt(f_minus, dtype=int ,delimiter='\t', usecols=[1,2], skiprows=1)


f_int = 'TSS_pairs_in_Groseq_merge1K_regions.txt'
with open ('Groseq_interestingHets_noWeird_merge1K_Discordance.bed', 'U') as merged_input:
    with open (f_int, 'w') as out:
        out.write(merged_input.readline().strip()) ##chrom  chromStart      chromEnd        Groseq.Discordance
        out.write('\tTss_paired_plus\tTss_paired_minus\n')
        for l in merged_input.readlines():
            m_chr, m_chromStart, m_chromEnd = l.strip().split('\t')[0:3]
            m_chromStart, m_chromEnd = int(m_chromStart), int(m_chromEnd)
            out.write(l.strip())
            c_Tss_plus_chromStart_End = Tss_plus_chromStart_End[Tss_plus_chromosome == m_chr]
            c_Tss_minus_chromStart_End = Tss_minus_chromStart_End[Tss_minus_chromosome == m_chr]
            m_Tss_plus = c_Tss_plus_chromStart_End[np.logical_and((c_Tss_plus_chromStart_End[:,0] <= m_chromEnd),(m_chromStart <= c_Tss_plus_chromStart_End[:,1]))]
            print m_Tss_plus.shape[0]
            m_Tss_minus = c_Tss_minus_chromStart_End[np.logical_and((c_Tss_minus_chromStart_End[:,0] <= m_chromEnd),(m_chromStart <= c_Tss_minus_chromStart_End[:,1]))]
            out.write('\t')
            out.write('\t'.join([str(m_Tss_plus.shape[0]), str(m_Tss_minus.shape[0])]))
            out.write('\n')

chromStart = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
chromEnd = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
regionLength = chromEnd - chromStart
Groseq_Concordance = np.loadtxt(f_int , dtype=str ,delimiter='\t', usecols=[3], skiprows=1) #C,D,P,M
Tss_plus_count = np.loadtxt(f_int , dtype=int ,delimiter='\t', usecols=[4], skiprows=1)
Tss_minus_count = np.loadtxt(f_int , dtype=int ,delimiter='\t', usecols=[5], skiprows=1)
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



plt.hist(Tss_plus_count[Groseq_Concordance=='D'], bins=range(10), alpha = 0.5, color = 'blue')
plt.hist(Tss_plus_count[Groseq_Concordance=='C'], bins=range(10), alpha = 0.5, color = 'yellow')
plt.title('TSS counts per region')
plt.xlabel('number of TSS')
plt.ylabel('freq')
plt.savefig('Tss_count_inConNDis.png')
plt.close()

plt.hist(regionLength[Groseq_Concordance=='C'], bins=range(1,5000,100), alpha = 0.2)
plt.hist(regionLength[Groseq_Concordance=='D'], bins=range(1,5000,100), alpha = 0.2)
plt.hist(regionLength[np.logical_and(Groseq_Concordance=='C', Tss_plus_count>0)], bins=range(1,5000,100), alpha = 0.2)
plt.hist(regionLength[np.logical_and(Groseq_Concordance=='D', Tss_plus_count>0)], bins=range(1,5000,100), alpha = 0.2)
#plt.title('Groseq_Concordance==D, Tss_plus_count>0')
plt.xlabel('regionLength')
plt.ylabel('freq')
plt.savefig('regionLength_TSS.hist.png')
plt.close()

##########test if TSS (NOT pairs) are enriched in Groseq Discordance regions
import numpy as np
f_tss = 'tss_all_gm12878.bed'
Tss_chromosome = np.loadtxt(f_tss, dtype=str ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
Tss_chromStart_End= np.loadtxt(f_tss, dtype=int ,delimiter='\t', usecols=[1,2], skiprows=1)
Tss_strand = np.loadtxt(f_tss, dtype=str ,delimiter='\t', usecols=[5], skiprows=1)

f_int = 'TSS_NOTpairs_in_Groseq_merge1K_regions.txt'
with open ('Groseq_interestingHets_noWeird_merge1K_Discordance.bed', 'U') as merged_input:
    with open (f_int, 'w') as out:
        out.write(merged_input.readline().strip()) ##chrom  chromStart      chromEnd        Groseq.Discordance
        out.write('\tTss_plus\tTss_minus\n')
        for l in merged_input.readlines():
            m_chr, m_chromStart, m_chromEnd = l.strip().split('\t')[0:3]
            m_chromStart, m_chromEnd = int(m_chromStart), int(m_chromEnd)
            out.write(l.strip())
            c_Tss_chromStart_End = Tss_chromStart_End[Tss_chromosome == m_chr]
            c_Tss_strand = Tss_strand[Tss_chromosome == m_chr]
            m_Tss_plus = c_Tss_chromStart_End[np.logical_and(c_Tss_strand =='+',np.logical_and((c_Tss_chromStart_End[:,0] <= m_chromEnd),(m_chromStart <= c_Tss_chromStart_End[:,1])))]
            #print m_Tss_plus.shape[0]
            m_Tss_minus = c_Tss_chromStart_End[np.logical_and(c_Tss_strand =='-',np.logical_and((c_Tss_chromStart_End[:,0] <= m_chromEnd),(m_chromStart <= c_Tss_chromStart_End[:,1])))]
            out.write('\t')
            out.write('\t'.join([str(m_Tss_plus.shape[0]), str(m_Tss_minus.shape[0])]))
            out.write('\n')

chromStart = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
chromEnd = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
regionLength = chromEnd - chromStart
Groseq_Concordance = np.loadtxt(f_int , dtype=str ,delimiter='\t', usecols=[3], skiprows=1) #C,D,P,M
Tss_plus_count = np.loadtxt(f_int , dtype=int ,delimiter='\t', usecols=[4], skiprows=1)
Tss_minus_count = np.loadtxt(f_int , dtype=int ,delimiter='\t', usecols=[5], skiprows=1)
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



plt.hist(Tss_plus_count[Groseq_Concordance=='D'], bins=range(15), alpha = 0.5, color = 'blue')
plt.hist(Tss_plus_count[Groseq_Concordance=='C'], bins=range(15), alpha = 0.5, color = 'yellow')
plt.title('TSS counts per region')
plt.xlabel('number of TSS_all (NOT pair)')
plt.ylabel('freq')
plt.savefig('Tss.all_count_inConNDis.png')
plt.close()

plt.hist(regionLength[Groseq_Concordance=='C'], bins=range(1,5000,100), alpha = 0.2)
plt.hist(regionLength[Groseq_Concordance=='D'], bins=range(1,5000,100), alpha = 0.2)
plt.hist(regionLength[np.logical_and(Groseq_Concordance=='C', Tss_plus_count>0)], bins=range(1,5000,100), alpha = 0.2)
plt.hist(regionLength[np.logical_and(Groseq_Concordance=='D', Tss_plus_count>0)], bins=range(1,5000,100), alpha = 0.2)
#plt.title('Groseq_Concordance==D, Tss_plus_count>0')
plt.xlabel('regionLength')
plt.ylabel('freq')
plt.savefig('regionLength_TSS.hist.png')
plt.close()

######test if TF ASB are enriched in Groseq Discordance regions with at least one TSS (non_pair or pair depends on input 'TSS_NOTpairs_in_Groseq_merge1K_regions.txt')

f = 'd.SymCals.all.txt'
#f = 'd.SymCals.USF2.txt'
f_int = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.counts_TSS_nonpair_Present.bed'
f_end = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF_FishersExactTest_TSS_nonpair_Present.txt'
TF_index = xrange(7,125)
TF_index_3 = xrange(6,124)

#data = [l.strip().split('\t') for l in open(f, "U").readlines()]
chromosome = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
snppos = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
TF = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=TF_index, skiprows=1)



with open ('TSS_NOTpairs_in_Groseq_merge1K_regions.txt', 'U') as merged_input:
    with open (f_int, 'w') as out:
        out.write(merged_input.readline().strip())
        for i in open(f,'U').readline().strip().split('\t')[7:125]:
            out.write('\t')
            out.write('.'.join(i.strip('"').split('.')[1:]))
        out.write('\n')
        for l in merged_input.readlines():
            m_chr, m_chromStart, m_chromEnd,_, TSS_plus_count, TSS_minus_count = l.strip().split('\t')[0:6]
            if int(TSS_plus_count) + int(TSS_minus_count ) > 0:
                out.write(l.strip())
                c = int(m_chr.split('chr')[-1])
                c_snppos = snppos[chromosome == c]
                c_sym = TF[chromosome == c]
                m_sym = c_sym[np.logical_and((c_snppos - 1 - int(m_chromStart) >= 0 ),(int(m_chromEnd) - c_snppos >= 0))]
                for i in xrange(m_sym.shape[1]):
                    temp = set(m_sym[:,i])
                    if '"Asym"' in temp:
                        out.write('\tAsym')
                    elif '"Sym"' in temp:
                        out.write('\tSym')
                    else:
                        out.write('\tNA')
                out.write('\n')

import scipy.stats as stats
Merged_TF_SCls  = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=TF_index_3, skiprows=1)
Merged_Concordance= np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)
TF_names = open(f_int, "U").readline().strip().split('\t')[6:124]


D_Asym = (Merged_TF_SCls[Merged_Concordance == 'D'] =='Asym').sum(axis=0)
D_Sym = (Merged_TF_SCls[Merged_Concordance == 'D'] =='Sym').sum(axis=0)
C_Asym = (Merged_TF_SCls[Merged_Concordance == 'C'] =='Asym').sum(axis=0)
C_Sym = (Merged_TF_SCls[Merged_Concordance == 'C'] =='Sym').sum(axis=0)
C_NA = (Merged_TF_SCls[Merged_Concordance == 'C'] =='NA').sum(axis=0)
D_NA = (Merged_TF_SCls[Merged_Concordance == 'D'] =='NA').sum(axis=0)

with open (f_end, 'w') as out:
    #out.write('\t'.join(['TF_names', 'p.left_tail', 'p.right_tail', 'p.two_tail','s_oddsratio', 's_pvalue','C_Asym', 'C_Sym', 'D_Asym', 'D_Sym']))
    out.write('\t'.join(['TF_names', 'oddsratio', 'pvalue_twotail','D_Asym', 'D_Sym', 'C_Asym', 'C_Sym','NA_oddsratio', 'NA_pvalue_twotail','C_NA', 'D_NA']))
    out.write('\n')
    for tf in range(len(TF_names)):
        #p = pvalue(C_Sym[tf], C_Asym[tf], D_Sym[tf], D_Asym[tf])  
        s_oddsratio, s_pvalue = stats.fisher_exact([[D_Asym[tf], D_Sym[tf]],[C_Asym[tf], C_Sym[tf]]])
        na_oddsratio, na_pvalue = stats.fisher_exact([[D_Asym[tf], D_NA[tf]],[C_Asym[tf], C_NA[tf]]])
        #print TF_names[tf], p.left_tail, p.right_tail, p.two_tail
        #out.write('\t'.join([TF_names[tf], str(p.left_tail), str(p.right_tail), str(p.two_tail),str(s_oddsratio), str(s_pvalue),str(C_Sym[tf]), str(C_Asym[tf]), str(D_Sym[tf]), str(D_Asym[tf])]))
        out.write('\t'.join([TF_names[tf],str(s_oddsratio), str(s_pvalue),str(D_Asym[tf]), str(D_Sym[tf]), str(C_Asym[tf]), str(C_Sym[tf]), str(na_oddsratio), str(na_pvalue),str(C_NA[tf]), str(D_NA[tf])]))
        out.write('\n')

######test if TSS pairs are enriched in Groseq Discordance regions
import numpy as np

f_plus = 'tss_SS_gm12878_plus.bed'
f_minus = 'tss_SS_gm12878_minus.bed'
f_int = 'TSS_SSpairs_in_Groseq_merge1K_regions.txt'

Tss_plus_chromosome = np.loadtxt(f_plus, dtype=str ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
Tss_plus_chromStart_End= np.loadtxt(f_plus, dtype=int ,delimiter='\t', usecols=[1,2], skiprows=1)
Tss_minus_chromosome = np.loadtxt(f_minus, dtype=str ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
Tss_minus_chromStart_End= np.loadtxt(f_minus, dtype=int ,delimiter='\t', usecols=[1,2], skiprows=1)


with open ('Groseq_interestingHets_noWeird_merge1K_Discordance.bed', 'U') as merged_input:
    with open (f_int, 'w') as out:
        out.write(merged_input.readline().strip()) ##chrom  chromStart      chromEnd        Groseq.Discordance
        out.write('\tTss_paired_plus\tTss_paired_minus\n')
        for l in merged_input.readlines():
            m_chr, m_chromStart, m_chromEnd = l.strip().split('\t')[0:3]
            m_chromStart, m_chromEnd = int(m_chromStart), int(m_chromEnd)
            out.write(l.strip())
            c_Tss_plus_chromStart_End = Tss_plus_chromStart_End[Tss_plus_chromosome == m_chr]
            c_Tss_minus_chromStart_End = Tss_minus_chromStart_End[Tss_minus_chromosome == m_chr]
            m_Tss_plus = c_Tss_plus_chromStart_End[np.logical_and((c_Tss_plus_chromStart_End[:,0] <= m_chromEnd),(m_chromStart <= c_Tss_plus_chromStart_End[:,1]))]
            print m_Tss_plus.shape[0]
            m_Tss_minus = c_Tss_minus_chromStart_End[np.logical_and((c_Tss_minus_chromStart_End[:,0] <= m_chromEnd),(m_chromStart <= c_Tss_minus_chromStart_End[:,1]))]
            out.write('\t')
            out.write('\t'.join([str(m_Tss_plus.shape[0]), str(m_Tss_minus.shape[0])]))
            out.write('\n')

chromStart = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
chromEnd = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
regionLength = chromEnd - chromStart
Groseq_Concordance = np.loadtxt(f_int , dtype=str ,delimiter='\t', usecols=[3], skiprows=1) #C,D,P,M
Tss_plus_count = np.loadtxt(f_int , dtype=int ,delimiter='\t', usecols=[4], skiprows=1)
Tss_minus_count = np.loadtxt(f_int , dtype=int ,delimiter='\t', usecols=[5], skiprows=1)
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



plt.hist(Tss_plus_count[Groseq_Concordance=='D'], bins=range(10), alpha = 0.5, color = 'blue')
plt.hist(Tss_plus_count[Groseq_Concordance=='C'], bins=range(10), alpha = 0.5, color = 'yellow')
plt.title('TSS SS counts per region')
plt.xlabel('number of SS TSS')
plt.ylabel('freq')
plt.savefig('Tss_SS_inConNDis.png')
plt.close()

(Tss_plus_count>0).sum()
(Tss_minus_count>0).sum()

(Tss_plus_count[Groseq_Concordance=='C']>0).sum()
(Tss_plus_count[Groseq_Concordance=='D']>0).sum()
(Tss_minus_count[Groseq_Concordance=='C']>0).sum()
(Tss_minus_count[Groseq_Concordance=='D']>0).sum()

######test if TF ASB are enriched in Groseq Discordance regions with at least one US/SU,UU,SS TSS 

f = 'd.SymCals.all.txt'
#f = 'd.SymCals.USF2.txt'

TF_index = xrange(7,125)
TF_index_3 = xrange(6,124)

#data = [l.strip().split('\t') for l in open(f, "U").readlines()]
chromosome = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
snppos = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
TF = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=TF_index, skiprows=1)

f_int0 = 'TSS_UUpairs_in_Groseq_merge1K_regions.txt'
f_int = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.counts_TSS_UUpairs_Present.bed'
f_end = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF_FishersExactTest_TSS_UUpairs_Present.txt'

with open (f_int0 , 'U') as merged_input:
    with open (f_int, 'w') as out:
        out.write(merged_input.readline().strip())
        for i in open(f,'U').readline().strip().split('\t')[7:125]:
            out.write('\t')
            out.write('.'.join(i.strip('"').split('.')[1:]))
        out.write('\n')
        for l in merged_input.readlines():
            m_chr, m_chromStart, m_chromEnd,_, TSS_plus_count, TSS_minus_count = l.strip().split('\t')[0:6]
            if int(TSS_plus_count) + int(TSS_minus_count ) > 0:
                out.write(l.strip())
                c = int(m_chr.split('chr')[-1])
                c_snppos = snppos[chromosome == c]
                c_sym = TF[chromosome == c]
                m_sym = c_sym[np.logical_and((c_snppos - 1 - int(m_chromStart) >= 0 ),(int(m_chromEnd) - c_snppos >= 0))]
                for i in xrange(m_sym.shape[1]):
                    temp = set(m_sym[:,i])
                    if '"Asym"' in temp:
                        out.write('\tAsym')
                    elif '"Sym"' in temp:
                        out.write('\tSym')
                    else:
                        out.write('\tNA')
                out.write('\n')

import scipy.stats as stats
Merged_TF_SCls  = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=TF_index_3, skiprows=1)
Merged_Concordance= np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)
TF_names = open(f_int, "U").readline().strip().split('\t')[6:124]


D_Asym = (Merged_TF_SCls[Merged_Concordance == 'D'] =='Asym').sum(axis=0)
D_Sym = (Merged_TF_SCls[Merged_Concordance == 'D'] =='Sym').sum(axis=0)
C_Asym = (Merged_TF_SCls[Merged_Concordance == 'C'] =='Asym').sum(axis=0)
C_Sym = (Merged_TF_SCls[Merged_Concordance == 'C'] =='Sym').sum(axis=0)
C_NA = (Merged_TF_SCls[Merged_Concordance == 'C'] =='NA').sum(axis=0)
D_NA = (Merged_TF_SCls[Merged_Concordance == 'D'] =='NA').sum(axis=0)

with open (f_end, 'w') as out:
    #out.write('\t'.join(['TF_names', 'p.left_tail', 'p.right_tail', 'p.two_tail','s_oddsratio', 's_pvalue','C_Asym', 'C_Sym', 'D_Asym', 'D_Sym']))
    out.write('\t'.join(['TF_names', 'oddsratio', 'pvalue_twotail','D_Asym', 'D_Sym', 'C_Asym', 'C_Sym','NA_oddsratio', 'NA_pvalue_twotail','C_NA', 'D_NA']))
    out.write('\n')
    for tf in range(len(TF_names)):
        #p = pvalue(C_Sym[tf], C_Asym[tf], D_Sym[tf], D_Asym[tf])  
        s_oddsratio, s_pvalue = stats.fisher_exact([[D_Asym[tf], D_Sym[tf]],[C_Asym[tf], C_Sym[tf]]])
        na_oddsratio, na_pvalue = stats.fisher_exact([[D_Asym[tf], D_NA[tf]],[C_Asym[tf], C_NA[tf]]])
        #print TF_names[tf], p.left_tail, p.right_tail, p.two_tail
        #out.write('\t'.join([TF_names[tf], str(p.left_tail), str(p.right_tail), str(p.two_tail),str(s_oddsratio), str(s_pvalue),str(C_Sym[tf]), str(C_Asym[tf]), str(D_Sym[tf]), str(D_Asym[tf])]))
        out.write('\t'.join([TF_names[tf],str(s_oddsratio), str(s_pvalue),str(D_Asym[tf]), str(D_Sym[tf]), str(C_Asym[tf]), str(C_Sym[tf]), str(na_oddsratio), str(na_pvalue),str(C_NA[tf]), str(D_NA[tf])]))
        out.write('\n')


####
######test if TF DISCORDANT ASB are enriched in Groseq Discordance regions with at least one TSS

f = 'TF.Discordance.call.batch1.txt'
#TSS_f = 'TSS_pairs_in_Groseq_merge1K_regions.txt'
#f_int = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.Discordance_TSS_pair_Present.bed'
#f_end = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.Discordance_FishersExactTest_TSS_pair_Present.txt'

TSS_f = 'TSS_NOTpairs_in_Groseq_merge1K_regions.txt'
f_int = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.Discordance_TSS_NOTpair_Present.bed'
f_end = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.Discordance_FishersExactTest_TSS_NOTpair_Present.txt'


TF_index = xrange(2,120)
TF_index_3 = xrange(6,124)

#data = [l.strip().split('\t') for l in open(f, "U").readlines()]
chromosome = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
snppos = np.loadtxt(f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
TF = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=TF_index, skiprows=1)



with open (TSS_f, 'U') as merged_input:
    with open (f_int, 'w') as out:
        out.write(merged_input.readline().strip())
        for i in open(f,'U').readline().strip().split('\t')[2:120]:
            out.write('\t')
            out.write('.'.join(i.strip('"').split('.')[3:]))
        out.write('\n')
        for l in merged_input.readlines():
            m_chr, m_chromStart, m_chromEnd,_, TSS_plus_count, TSS_minus_count = l.strip().split('\t')[0:6]
            if int(TSS_plus_count) + int(TSS_minus_count ) > 0:
                out.write(l.strip())
                c = int(m_chr.split('chr')[-1])
                c_snppos = snppos[chromosome == c]
                c_sym = TF[chromosome == c]
                m_sym = c_sym[np.logical_and((c_snppos - 1 - int(m_chromStart) >= 0 ),(int(m_chromEnd) - c_snppos >= 0))]
                for i in xrange(m_sym.shape[1]):
                    temp = set(m_sym[:,i])
                    print temp
                    if '1' in temp:
                        out.write('\tD')
                    elif '-1' in temp:
                        out.write('\tC')
                    else:
                        out.write('\tNA')
                out.write('\n')

import scipy.stats as stats
Merged_TF_SCls  = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=TF_index_3, skiprows=1)
Merged_Concordance= np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)
TF_names = open(f_int, "U").readline().strip().split('\t')[6:124]


D_Asym = (Merged_TF_SCls[Merged_Concordance == 'D'] =='D').sum(axis=0)
D_Sym = (Merged_TF_SCls[Merged_Concordance == 'D'] =='C').sum(axis=0)
C_Asym = (Merged_TF_SCls[Merged_Concordance == 'C'] =='D').sum(axis=0)
C_Sym = (Merged_TF_SCls[Merged_Concordance == 'C'] =='C').sum(axis=0)
C_NA = (Merged_TF_SCls[Merged_Concordance == 'C'] =='NA').sum(axis=0)
D_NA = (Merged_TF_SCls[Merged_Concordance == 'D'] =='NA').sum(axis=0)

with open (f_end, 'w') as out:
    #out.write('\t'.join(['TF_names', 'p.left_tail', 'p.right_tail', 'p.two_tail','s_oddsratio', 's_pvalue','C_Asym', 'C_Sym', 'D_Asym', 'D_Sym']))
    out.write('\t'.join(['TF_names', 'oddsratio', 'pvalue_twotail','Groseq.D_TF.ASB.D', 'Groseq.D_TF.ASB.C', 'C_TF.D', 'C_TF.C','NA_oddsratio', 'NA_pvalue_twotail','C_NA', 'D_NA']))
    out.write('\n')
    for tf in range(len(TF_names)):
        #p = pvalue(C_Sym[tf], C_Asym[tf], D_Sym[tf], D_Asym[tf])  
        s_oddsratio, s_pvalue = stats.fisher_exact([[D_Asym[tf], D_Sym[tf]],[C_Asym[tf], C_Sym[tf]]])
        na_oddsratio, na_pvalue = stats.fisher_exact([[D_Asym[tf], D_NA[tf]],[C_Asym[tf], C_NA[tf]]])
        #print TF_names[tf], p.left_tail, p.right_tail, p.two_tail
        #out.write('\t'.join([TF_names[tf], str(p.left_tail), str(p.right_tail), str(p.two_tail),str(s_oddsratio), str(s_pvalue),str(C_Sym[tf]), str(C_Asym[tf]), str(D_Sym[tf]), str(D_Asym[tf])]))
        out.write('\t'.join([TF_names[tf],str(s_oddsratio), str(s_pvalue),str(D_Asym[tf]), str(D_Sym[tf]), str(C_Asym[tf]), str(C_Sym[tf]), str(na_oddsratio), str(na_pvalue),str(C_NA[tf]), str(D_NA[tf])]))
        out.write('\n')



