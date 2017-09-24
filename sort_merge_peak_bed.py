from sys import argv
#from operator import itemgetter

f = argv[1]#'ENCFF002CGO.bed'
l = [l.strip().split('\t') for l in open(f,'U').readlines()]

for ll in l:
    ll[0] = ll[0].split('chr')[-1]

l.sort()
#sorted_list = sorted(l, key=itemgetter(0,1,2))

doprint=0
p_chrom, p_chromStart, p_chromEnd = l[0][0:3]
for ll in l:
    chrom, chromStart, chromEnd = ll[0:3]
    if (p_chrom == chrom) and (int(p_chromEnd) >= int(chromStart)):
        if (int(p_chromEnd) <= int(chromEnd)):
            p_chromEnd = chromEnd
    else:
        if (doprint == 1):
            print p_chrom+'\t'+p_chromStart+'\t'+p_chromEnd
        p_chrom= chrom
        p_chromStart= chromStart
        p_chromEnd= chromEnd
    doprint = 1

print p_chrom+"\t"+p_chromStart+"\t"+chromEnd


######

import numpy as np
Binding = np.loadtxt('counts.txt', dtype=float ,delimiter='\t', usecols=[16], skiprows=1)


###
##
import os
import glob
file_names = glob.glob('*.bed')
template ="python /workdir/sc2457/AlleleSeq_pipeline_v1.2a/CombineSnpCounts.py 5 /workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/snp.call.all %s /workdir/sc2457/alleleseq.gersteinlab.org/NA12878_diploid_dec16.2012.alleleseq.input/rd_4327183snps_na12878_hg19.txt Groseq.count.files/%s.count.txt counts.log /workdir/sc2457/Groseq.gm12878_Allele_seq_result/SRR1552485.scythed.sickled.cnt"
run_sh_lines=[]
for f in file_names:
    new_file = f.split('.')[0]
    run_sh_lines.append(template %(f, new_file))


out = open('CombineSnpCounts.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()


for i in range(10):
    out = open('CombineSnpCounts_'+str(i)+'.sh',"w")
    out.write("\n".join(run_sh_lines[32*i:32*(i+1)]))
    out.close()
##


import os
import glob
file_names = glob.glob('CombineSnpCounts_*.sh')
template ='nohup sh %s &> %s.log &'

run_sh_lines=[]
for f in file_names:
    new_file = f.split('.')[0]
    run_sh_lines.append(template %(f, new_file))

out = open('CombineSnpCounts.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()



####determine if a snp are inlcuding in any regions in a bed file
import multiprocessing
import pickle
import numpy as np
#snp_f = '/workdir/sc2457/Groseq.gm12878_Allele_seq_result/interestingHets_noWeird.txt'
#snp_f ='/workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/snp.call.all'
snp_f = '/workdir/sc2457/TF_count_batch1N2/d.SymCals.all.txt'
snp_chrom = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
snp_snppos = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
#snp_in_region =  np.zeros(snp_snppos.shape, dtype=np.int)

snp_in_region_list=[]

#import os
#import glob
#file_names = glob.glob('*.bed')

file_names = [l.strip() for l in open('Larger_peak_file_list.txt','U').readlines()]
#file_names = ['ENCFF001USG.bed']
#for f in file_names:
#    title.append(f)
#    snp_in_region =  np.zeros(snp_snppos.shape, dtype=np.int)
#    with open(f+'.bed','U') as peak_bed:
#        regions =[ l.strip().split('\t') for l in peak_bed.readlines()]
#        for r in regions:
#            bed_chrom, bed_chromStart, bed_chromEnd = r[0:3]
#            try:
#                bed_chrom = int(bed_chrom.split('chr')[-1])
#                bed_chromStart = int (bed_chromStart)
#                bed_chromEnd = int(bed_chromEnd)
#                snp_in_region[np.logical_and(snp_chrom == bed_chrom,
#                               np.logical_and((snp_snppos - 1 - int(bed_chromStart) >= 0 ),(int(bed_chromEnd) - snp_snppos >= 0)))] = 1
#            except ValueError:
#                pass
#                
#        print f, snp_in_region.sum()
#        snp_in_region_list.append(snp_in_region)



def get_sno_in_region_vector(f):
    snp_in_region =  np.zeros(snp_snppos.shape, dtype=np.int)
    with open(f+'.bed','U') as peak_bed:
        regions =[ l.strip().split('\t') for l in peak_bed.readlines()]
        for r in regions:
            bed_chrom, bed_chromStart, bed_chromEnd = r[0:3]
            try:
                bed_chrom = int(bed_chrom.split('chr')[-1])
                bed_chromStart = int (bed_chromStart)
                bed_chromEnd = int(bed_chromEnd)
                snp_in_region[np.logical_and(snp_chrom == bed_chrom,
                               np.logical_and((snp_snppos - 1 - int(bed_chromStart) >= 0 ),(int(bed_chromEnd) - snp_snppos >= 0)))] = 1
            except ValueError:
                pass
                
        print f, snp_in_region.sum()
        return f, snp_in_region

file_names = [l.strip() for l in open('Larger_peak_file_list.txt','U').readlines()]
#get_sno_in_region_vector(file_names[1])

pool = multiprocessing.Pool(processes=32)
#pool_output = pool.map(subject_occurrence_per_file_memory_demanding , f_list )
pool_output = pool.map(get_sno_in_region_vector, file_names)
    # pool_output looks like [(f, snp_in_region), (f, snp_in_region),...]
pool.close() # no more tasks
pool.join()

SRF = get_sno_in_region_vector('ENCFF002CHW')
STAT1 = get_sno_in_region_vector('ENCFF001VFL')

with open('larger.peak.call.file_summary_with_TF.SymCals.all',"w") as pickle_out:
    pickle.dump(pool_output, pickle_out )

a = np.concatenate((snp_chrom.reshape((2024311,1)), snp_snppos.reshape((2024311,1))),axis=1)
title = ['chrom', 'snppos']

for p in pool_output:
    fp, snp_presence_p = p
    title.append(fp)
    a = np.concatenate((a, snp_presence_p.reshape((2024311,1))),axis=1)

title.append('ENCFF002CHW')
title.append('ENCFF001VFL')
a = np.concatenate((a, SRF[1].reshape(2024311,1)),axis=1)
a = np.concatenate((a, STAT1[1].reshape(2024311,1)),axis=1)

#tf_peak_dic = {}
#for i in range(len(title)-2):
#    tf_peak_dic[title[i+2]] = a[:,i+2]

#np.savetxt('d.SymCals.all_withTFPeakCall.txt', a, fmt='%d', delimiter='\t')
#np.savetxt('d.SymCals.all_withTFPeakCall.header.txt', np.array(title).reshape(1,104), fmt='%s', delimiter='\t')

b = np.concatenate((np.array(title).reshape(1,106), a),axis=0)
#np.savetxt('d.SymCals.all_withTFPeakCall.txt', b, fmt='%s', delimiter='\t')


TF_names = [tf.strip('\"') for tf in open(snp_f,'U').readline().strip().split('\t')[7:125]]
TF = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=xrange(7,125), skiprows=0)
TF[0,:] = np.array(TF_names)


combined_matrix = np.concatenate((b[:,0:2],TF),axis=1)
combined_matrix = np.concatenate((c,b[:,2:]),axis=1)

np.savetxt('d.SymCals.all_withTFPeakCall.txt', combined_matrix, fmt='%s', delimiter='\t')

combined_matrix = c 




######test if TF ASB are enriched in Groseq Discordance regions with at least one TSS (non_pair or pair depends on input 'TSS_NOTpairs_in_Groseq_merge1K_regions.txt')

snp_f='d.SymCals.all.YY1_all.txt'
TF_names = [tf.strip('\"') for tf in open(snp_f,'U').readline().strip().split('\t')[7:126]]
#TF = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=xrange(7,126), skiprows=0)
#TF[0,:] = np.array(TF_names)

chromosome = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
snppos = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
TF_i = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=xrange(7,126), skiprows=1)
#peak_calls = b[1:,2:]
peak_calls = np.loadtxt('d.SymCals.all_withTFPeakCall.txt', dtype=str ,delimiter='\t', usecols=xrange(120,224), skiprows=1)
peak_calls_file_names = open('d.SymCals.all_withTFPeakCall.txt','U').readline().strip().split('\t')[120:224]

f_merge = '/workdir/sc2457/TF_count_batch1N2/TSS_NOTpairs_in_Groseq_merge1K_regions.txt'
f_int = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.counts_TSS_nonpair_Present_PeakCallFiltered.bed'
f_end = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF_FishersExactTest_TSS_nonpair_Present_PeakCallFiltered.txt'

f_merge = '/workdir/sc2457/TF_count_batch1N2/TSS_pairs_in_Groseq_merge1K_regions.txt'
f_int = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.counts_TSS_pair_Present_PeakCallFiltered.bed'
f_end = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF_FishersExactTest_TSS_pair_Present_PeakCallFiltered.txt'



with open (f_merge, 'U') as merged_input:
    with open (f_int, 'w') as out:
        out.write(merged_input.readline().strip())
        out.write('\t')
        out.write('\t'.join(TF_names))
        out.write('\t')
        out.write('\t'.join(peak_calls_file_names))
        out.write('\n')
        for l in merged_input.readlines():
            m_chr, m_chromStart, m_chromEnd,_, TSS_plus_count, TSS_minus_count = l.strip().split('\t')[0:6]
            if int(TSS_plus_count) + int(TSS_minus_count ) > 0:
                out.write(l.strip())
                c = int(m_chr.split('chr')[-1])
                c_snppos = snppos[chromosome == c]
                c_sym = TF_i[chromosome == c]
                c_peak_call = peak_calls[chromosome == c]
                m_sym = c_sym[np.logical_and((c_snppos - 1 - int(m_chromStart) >= 0 ),(int(m_chromEnd) - c_snppos >= 0))]
                m_peak_call = c_peak_call[np.logical_and((c_snppos - 1 - int(m_chromStart) >= 0 ),(int(m_chromEnd) - c_snppos >= 0))]
                for i in xrange(m_sym.shape[1]):
                    temp = set(m_sym[:,i])
                    if '"Asym"' in temp:
                        out.write('\tAsym')
                    elif '"Sym"' in temp:
                        out.write('\tSym')
                    else:
                        out.write('\tNA')
                for i in xrange(m_peak_call.shape[1]):
                    temp = set(m_peak_call[:,i])
                    #print temp 
                    if '1' in temp:
                        out.write('\t1')
                    else:
                        out.write('\t0')
                out.write('\n')

import numpy as np 
import scipy.stats as stats
#f_int = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.counts_TSS_nonpair_Present_PeakCallFiltered.bed'
Merged_TF_SCls  = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=xrange(6,125), skiprows=1)
Merged_TF_PeakCls  = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=xrange(125,229), skiprows=1)
Merged_Concordance= np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)
TF_names = ['.'.join(tf.split('.')[1:]) for tf in open(f_int, "U").readline().strip().split('\t')[6:125]]
Peak_file_names = open(f_int, "U").readline().strip().split('\t')[125:]

tf_peak_file_dic={}
for l in open('Larger_peak_file_list_TF.dic', 'U').readlines()[1:]:
    peak_file_name, tagert_label = l.strip().split('\t')[0:2]
    tf_peak_file_dic[tagert_label] = peak_file_name

TF_names_withPeakFile = []
for tf in TF_names:
    new_tf = tf.split('.')[0]
    try :
        #tf_peak_file_dic[tf]
        TF_names_withPeakFile.append((tf,tf_peak_file_dic[new_tf]))
    except KeyError:
        print 'KeyError: ', tf

TF_names_nparray={}
for i in range(len(TF_names)):
    TF_names_nparray[TF_names[i]] = Merged_TF_SCls[:,i]

TF_peak_nparray={}
for i in range(len(Peak_file_names)):
    TF_peak_nparray[Peak_file_names[i]] = Merged_TF_PeakCls[:,i]

TF_names_withPeakFile.sort()
tf_test_result_combine=[]
for tf,tf_peak_file in TF_names_withPeakFile:
    #print tf,tf_peak_file
    tf_SymCals = TF_names_nparray[tf]
    tf_peak = TF_peak_nparray[tf_peak_file]
    test_list = [('D','Asym'), ('D', 'Sym'), ('C', 'Asym'), ('C','Sym'), ('D','NA'),('C','NA')]
    test_result=[]
    for i in range(len(test_list)):
        test_result.append((tf_SymCals[np.logical_and(tf_peak =='1',Merged_Concordance==test_list[i][0])] == test_list[i][1]).sum())
    tf_test_result_combine.append((tf, tf_peak_file, test_result)) #D_Asym, D_Sym, C_Asym, C_Sym, D_NA, C_NA
    #D_Asym = (tf_SymCals[np.logical_and(tf_peak =='1',Merged_Concordance=='D')] == 'Asym').sum()


with open (f_end, 'w') as out:
    #out.write('\t'.join(['TF_names', 'p.left_tail', 'p.right_tail', 'p.two_tail','s_oddsratio', 's_pvalue','C_Asym', 'C_Sym', 'D_Asym', 'D_Sym']))
    out.write('\t'.join(['TF_names', 'oddsratio', 'pvalue_twotail','D_Asym', 'D_Sym', 'C_Asym', 'C_Sym','D_NA', 'C_NA']))
    out.write('\n')
    for tf, tf_peak_file, tf_DC_count in tf_test_result_combine:
        #p = pvalue(C_Sym[tf], C_Asym[tf], D_Sym[tf], D_Asym[tf])  
        s_oddsratio, s_pvalue = stats.fisher_exact([[tf_DC_count[0], tf_DC_count[1]],[tf_DC_count[2], tf_DC_count[3]]])
        print tf, s_oddsratio, s_pvalue
        out.write('\t'.join([tf, str(s_oddsratio), str(s_pvalue)]))
        for i in tf_DC_count:
            out.write('\t%d' % i)
        out.write('\n')

