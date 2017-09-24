import numpy as np
#snp_f = '/workdir/sc2457/Groseq.gm12878_Allele_seq_result/interestingHets_noWeird.txt'
#snp_chrom = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
#snp_snppos = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
#winning_parent = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=[13], skiprows=1)
##m_snp = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=[7], skiprows=1)
##p_snp = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=[8], skiprows=1)
##ACGT_reads_count = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[9,10,11,12], skiprows=1)

def read_file_to_nparray (snp_f):
    temp = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[0,1], skiprows=1) #chr1-chr22
    winning_parent = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=[13], skiprows=1)
    #mat_allele =  np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=[7], skiprows=1)
    #pat_allele =  np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=[8], skiprows=1)
    cACGT = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=[7,8,9,10,11,12], skiprows=1) #chr1-chr22
    return temp[:,0], temp[:,1], winning_parent,cACGT #chr, snppos, winning_parent

plus_groseq_fp = '/workdir/sc2457/Groseq.gm12878_Allele_seq_result/StranSpecificASE/interestingHets_noW_plus.txt'
minus_groseq_fp = '/workdir/sc2457/Groseq.gm12878_Allele_seq_result/StranSpecificASE/interestingHets_noW_minus.txt'

plus_chrom, plus_snppos, plus_winning_parent,plus_cACGT = read_file_to_nparray(plus_groseq_fp)
minus_chrom, minus_snppos, minus_winning_parent, minus_cACGT = read_file_to_nparray(minus_groseq_fp)


tss_plus_f = '/workdir/sc2457/TF_count_batch1N2/tss_paired_gm12878_plus.bed'
tss_minus_f = '/workdir/sc2457/TF_count_batch1N2/tss_paired_gm12878_minus.bed'


tss_plus = [l.strip().split('\t') for l in open(tss_plus_f, 'U').readlines()]
tss_minus = [l.strip().split('\t') for l in open(tss_minus_f, 'U').readlines()]


def if_snp_in_a_region_NOT_strandSpecific (rChrom, rStart, rEnd,sChrom=snp_chrom, snppos=snp_snppos, winning_parent=winning_parent):
    """
    sChrom and snppos are np.array
    rChrom, rStart, rEnd are integer
    """
    
    # keep snppos and winning_parent in the same chrom as rChrom
    c_snppos = snppos[sChrom == rChrom]
    c_winning_parent = winning_parent[sChrom == rChrom]
    # the ones in the region
    in_ones = np.logical_and((c_snppos - 1 - rStart >= 0 ),(rEnd - c_snppos >= 0))
    return in_ones, c_snppos[in_ones], c_winning_parent[in_ones]


def if_snp_in_a_region(rChrom, rStart, rEnd, strand):
    """
    sChrom and snppos are np.array
    stand = 'plus' or 'minus'
    rChrom, rStart, rEnd are integer
    """
    assert strand in ['plus', 'minus']
    if strand == 'plus':
        sChrom,snppos, winning_parent, cACGT = plus_chrom, plus_snppos, plus_winning_parent, plus_cACGT
    else:
        sChrom,snppos, winning_parent, cACGT  = minus_chrom, minus_snppos, minus_winning_parent, minus_cACGT
        
    # keep snppos and winning_parent in the same chrom as rChrom
    c_snppos = snppos[sChrom == rChrom]
    c_winning_parent = winning_parent[sChrom == rChrom]
    c_cACGT = cACGT[sChrom == rChrom]
    # the ones in the region
    in_ones = np.logical_and((c_snppos - 1 - rStart >= 0 ),(rEnd - c_snppos >= 0))
    return in_ones, c_snppos[in_ones], c_winning_parent[in_ones], c_cACGT[in_ones]


def plus_Tss_three_stage_check(pChrom, pStart, pEnd, d_s0, d_s1, d_l, d_after_rStart=20, s0_count=0, s1_count=0, l_count=0):
    pw_use,pcACGT_use = None,None
    pt, ps, pw, pcACGT = if_snp_in_a_region(pChrom, pStart+d_after_rStart , pEnd + d_s0, 'plus')
    if pt.sum() > 0:
        pw_use = pw[-1]  #use last one
        ps_use = ps[-1]
        pcACGT_use = pcACGT[-1]
        s0_count += 1
        #if len(set(pw)) > 1:
        #print 'strict0',s0_count , pw,'pw_use:', pw_use, pChrom, pStart, pEnd
    else:
        pt, ps, pw ,pcACGT = if_snp_in_a_region(pChrom, pEnd + d_s0 , pEnd + d_s1, 'plus')
        if pt.sum() > 0:
            pw_use = pw[0]  #use first one
            pcACGT_use = pcACGT[0]
            s1_count += 1
        #    print 'strict1',s1_count , pw,'pw_use:', pw_use, pChrom, pStart, pEnd
        else:
            pt, ps, pw,pcACGT  = if_snp_in_a_region(pChrom, pEnd + d_s1 , pEnd + d_l, 'plus')
            if pt.sum() > 0:
                pw_use = pw[0]  #use first one
                pcACGT_use = pcACGT[0]
                l_count += 1
        #        print 'loose',l_count , pw,'pw_use:', pw_use, pChrom, pStart, pEnd
    return pt, ps, pw, pw_use, pcACGT_use


def minus_Tss_three_stage_check(mChrom, mStart, mEnd, d_s0, d_s1, d_l, d_after_rStart=20, s0_count=0, s1_count=0, l_count=0):
    mw_use, mcACGT_use = None,None
    mt, ms, mw,mcACGT = if_snp_in_a_region(mChrom, mStart - d_s0 , mEnd - d_after_rStart, 'minus' )
    if mt.sum() > 0:
        mw_use = mw[0]  #use 1st one
        mcACGT_use = mcACGT[0]
        s0_count += 1
        #if len(set(pw)) > 1:
        #print 'strict0',s0_count , mw,'mw_use:', mw_use, mChrom, mStart, mEnd
    else:
        mt, ms, mw,mcACGT = if_snp_in_a_region(mChrom, mStart - d_s1 , mStart - d_s0, 'minus')
        if mt.sum() > 0:
            mw_use = mw[-1]  #use last one
            mcACGT_use = mcACGT[-1]
            s1_count += 1
        #    print 'strict1',s1_count , mw,'mw_use:', mw_use, mChrom, mStart, mEnd
        else:
            mt, ms, mw,mcACGT = if_snp_in_a_region(mChrom, mStart - d_l , mStart - d_s1, 'minus')
            if mt.sum() > 0:
                mw_use = mw[-1]  #use last one
                mcACGT_use = mcACGT[-1]
                l_count += 1
        #        print 'loose',l_count , mw,'mw_use:', mw_use, mChrom, mStart, mEnd
    return mt, ms, mw, mw_use, mcACGT_use

ACGT_map={'A':2, 'C':3, 'G':4, 'T':5}
def directionality(pcACGT_use, mcACGT_use):
    mat_count_plus = int(pcACGT_use[ACGT_map[pcACGT_use[0]]])
    mat_count_minus =int(mcACGT_use[ACGT_map[mcACGT_use[0]]])
    pat_count_plus = int(pcACGT_use[ACGT_map[pcACGT_use[1]]])
    pat_count_minus =int(mcACGT_use[ACGT_map[mcACGT_use[1]]])
    print mat_count_plus, mat_count_minus, pat_count_plus, pat_count_minus
    
    
    
    
    





def get_Tss_dl_regions_with_groseq_ASB(d_l, out_fp, d_region):
    with open(out_fp, 'w') as out:
        out.write('\t'.join(['#chrom','chromStart','chromEnd','Groseq.Discordance','Tss_plus_winParent','Tss_minus_winParent'])
            +'\n')
        d_s0 = 50
        d_s1 = 100
        #d_l = 500
        d_after_rStart = 20
        count, s0_count, s1_count, l_count = 0,0,0,0
        Con_count, Dis_count = 0, 0
        for i in xrange(len(tss_plus)):
            pChrom, pStart, pEnd = tss_plus[i][0:3]
            try:
                pChrom = int(pChrom.split('chr')[-1])
                pStart, pEnd = int(pStart), int(pEnd)
                # check different regions of Tss
                pt, ps, pw, pw_use,pcACGT_use = plus_Tss_three_stage_check(pChrom, pStart, pEnd, d_s0, d_s1, d_l,d_after_rStart)
                mChrom, mStart, mEnd = tss_minus[i][0:3]
                mChrom = int(mChrom.split('chr')[-1])
                assert mChrom == pChrom
                mStart, mEnd = int(mStart) , int(mEnd)
                mt, ms, mw, mw_use,mcACGT_use = minus_Tss_three_stage_check(mChrom, mStart, mEnd, d_s0, d_s1, d_l,d_after_rStart)
                
                if min(pt.sum(),mt.sum()) >= 1:
                    count += 1
                    if pw_use == mw_use:
                        #print pcACGT_use, mcACGT_use
                        Con_count +=1
                        out.write('\t'.join(['chr'+str(pChrom), str(min(pStart,mStart) - d_region), str(max(pEnd, mEnd) + d_region),'Concordant', pw_use, mw_use]+list(pcACGT_use)+list(mcACGT_use))
                                  +'\n')
                    elif pw_use != mw_use:
                        directionality(pcACGT_use, mcACGT_use)
                        Dis_count +=1
                        out.write('\t'.join(['chr'+str(pChrom), str(min(pStart,mStart) - d_region), str(max(pEnd, mEnd) + d_region),'Discordant', pw_use, mw_use])
                                  +'\n')
                        #print 'D:', str(pChrom)+':'+str(pStart)+'-'+str(pEnd), pw_use, mw_use
                    #print 'chr'+str(pChrom), str(min(pStart,mStart) - d_l), str(max(pEnd, mEnd) + d_l), pw_use, mw_use
                    #print count, pChrom, str(pStart)+'-'+str(pEnd), pw, mw
            except ValueError:
                pass
    print '\t'.join(['d='+str(d_l), 'Con', str(Con_count), 'Dis', str(Dis_count)])


get_Tss_dl_regions_with_groseq_ASB(150, 'GM12878_GroSeq_d150_dRegion250_withStrandSpecific.txt', 250) 



#####test if TF ASB are enriched in Groseq Discordance regions with at least one TSS (non_pair or pair depends on input 'TSS_NOTpairs_in_Groseq_merge1K_regions.txt')

snp_f='/local/workdir/sc2457/ENCODE/fastq/PeakCallBedFiles/d.SymCals.all.YY1_all_W.txt'
TF_names_i = [tf.strip('\"') for tf in open(snp_f,'U').readline().strip().split('\t')[7:130]]
chromosome = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
snppos = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
TF_i = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=xrange(7,130), skiprows=1)
peak_calls = np.loadtxt('/local/workdir/sc2457/ENCODE/fastq/PeakCallBedFiles/d.SymCals.all_withTFPeakCall.txt', dtype=str ,delimiter='\t', usecols=xrange(120,225), skiprows=1)
peak_calls_file_names = open('/local/workdir/sc2457/ENCODE/fastq/PeakCallBedFiles/d.SymCals.all_withTFPeakCall.txt','U').readline().strip().split('\t')[120:225]

####
snp_f='/workdir/sc2457/ENCODE/reDoTF/reDoTF_counts_file/SRF/d.SymCals.all.txt'
TF_names_i = [tf.strip('\"') for tf in open(snp_f,'U').readline().strip().split('\t')[7:169]]
TF_i = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=xrange(7,169), skiprows=1)
TF_i_old = np.copy(TF_i)
#chromosome, snppos = snp_chrom, snp_snppos
chromosome = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[0], skiprows=1) #chr1-chr22
snppos = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)

peak_calls = np.loadtxt('/workdir/sc2457/ENCODE/reDoTF/ENCODE_bed/snps_with_in_PeakCallbed.txt', dtype=int ,delimiter='\t', usecols=xrange(2,440), skiprows=1)
peak_calls_file_names = open('/workdir/sc2457/ENCODE/reDoTF/ENCODE_bed/snps_with_in_PeakCallbed.txt','U').readline().strip().split('\t')[2:440]


#f_merge, f_int, f_end = 'f_merge', 'f_int', 'f_end'
#get_Tss_dl_regions_with_groseq_ASB(150, f_merge, 250)


import scipy.stats as stats

#def get_TF_FishersExactTest_TSS_PeakCallFiltered (f_merge, f_int, f_end):

metadata = [l.strip().split('\t') for l in open('/workdir/sc2457/ENCODE/reDoTF/ENCODE_bed/all.bed.metadata.tsv', 'U').readlines()]
ExperimentAcc_BedFileAcc_dic = {}
BedFileAcc_ExperimentAcc_dic = {}
Target_BedFileAcc_dic = {}
for l in metadata[1:]:
    experiment_accession = l[3]
    file_format = l[1]
    file_accession = l[0]
    target = l[15].split('-')[0]
    if file_format.count('bed') > 0:
        if experiment_accession not in ExperimentAcc_BedFileAcc_dic:
            ExperimentAcc_BedFileAcc_dic[experiment_accession] = []
        ExperimentAcc_BedFileAcc_dic[experiment_accession].append(file_accession)
        BedFileAcc_ExperimentAcc_dic[file_accession] = experiment_accession
        if target not in Target_BedFileAcc_dic:
            Target_BedFileAcc_dic[target] = []
        Target_BedFileAcc_dic[target].append(file_accession)


TFnamesExpAcc_withPeakFileFromSameTarget = []
TFnamesExpAcc_withPeakFileFromSameExpAcc = []
for tf in TF_names_i:
    target, ExpAcc = tf.split('.')[1:]
    try :
        #tf_peak_file_dic[tf]
        TFnamesExpAcc_withPeakFileFromSameTarget.append((tf,Target_BedFileAcc_dic[target]))
        TFnamesExpAcc_withPeakFileFromSameExpAcc.append((tf, ExperimentAcc_BedFileAcc_dic[ExpAcc]))
    except KeyError:
        print 'KeyError: ', tf

TF_names_nparray_snp={}
for i in range(len(TF_names_i)):
    TF_names_nparray_snp[TF_names_i[i]] = TF_i[:,i]

TF_peak_nparray={}
for i in range(len(peak_calls_file_names)):
    TF_peak_nparray[peak_calls_file_names[i]] = peak_calls[:,i]#.astype(int)

def keep_snps_within_peak_files():
    for tf, tf_peak_files in TFnamesExpAcc_withPeakFileFromSameTarget:
        tf_peak = np.copy(TF_peak_nparray[tf_peak_files[0]])
        for tf_p_f in tf_peak_files[1:]:
            #print TF_peak_nparray[tf_p_f]
            tf_peak = tf_peak + TF_peak_nparray[tf_p_f]
        TF_names_nparray_snp[tf][tf_peak == 0] = 'NotInPeak'
        if tf.count('SRF') > 0:
            print tf, tf_peak

def use_All_snps():
    TF_i = TF_i_old
    #TF_i = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=xrange(7,165), skiprows=1)

keep_snps_within_peak_files()
with open (f_merge, 'U') as merged_input:
    with open (f_int, 'w') as out:
        out.write(merged_input.readline().strip())
        out.write('\t')
        out.write('\t'.join(TF_names_i))
        out.write('\t')
        out.write('\t'.join(peak_calls_file_names))
        out.write('\n')
        for l in merged_input.readlines():
            m_chr, m_chromStart, m_chromEnd,_, TSS_plus_count, TSS_minus_count = l.strip().split('\t')[0:6]
            #if int(TSS_plus_count) + int(TSS_minus_count ) > 0:
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
                if 1 in temp:
                    out.write('\t1')
                else:
                    out.write('\t0')
            out.write('\n')








Merged_TF_bed = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=xrange(0,3), skiprows=1)
Merged_TF_SCls  = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=xrange(6,168), skiprows=1)
Merged_TF_PeakCls  = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=xrange(168,606), skiprows=1)
Merged_Concordance= np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)
TF_names = ['.'.join(tf.split('.')[1:]) for tf in open(f_int, "U").readline().strip().split('\t')[6:168]]
#TF_names = ['.'.join(tf.split('.')[1:]) for tf in TF_names_i ]
Peak_file_names = open(f_int, "U").readline().strip().split('\t')[168:]


metadata = [l.strip().split('\t') for l in open('/workdir/sc2457/ENCODE/reDoTF/ENCODE_bed/all.bed.metadata.tsv', 'U').readlines()]
ExperimentAcc_BedFileAcc_dic = {}
BedFileAcc_ExperimentAcc_dic = {}
Target_BedFileAcc_dic = {}
for l in metadata[1:]:
    experiment_accession = l[3]
    file_format = l[1]
    file_accession = l[0]
    target = l[15].split('-')[0]
    if file_format.count('bed') > 0:
        if experiment_accession not in ExperimentAcc_BedFileAcc_dic:
            ExperimentAcc_BedFileAcc_dic[experiment_accession] = []
        ExperimentAcc_BedFileAcc_dic[experiment_accession].append(file_accession)
        BedFileAcc_ExperimentAcc_dic[file_accession] = experiment_accession
        if target not in Target_BedFileAcc_dic:
            Target_BedFileAcc_dic[target] = []
        Target_BedFileAcc_dic[target].append(file_accession)


TFnamesExpAcc_withPeakFileFromSameTarget = []
TFnamesExpAcc_withPeakFileFromSameExpAcc = []
for tf in TF_names:
    target, ExpAcc = tf.split('.')
    try :
        #tf_peak_file_dic[tf]
        TFnamesExpAcc_withPeakFileFromSameTarget.append((tf,Target_BedFileAcc_dic[target]))
        TFnamesExpAcc_withPeakFileFromSameExpAcc.append((tf, ExperimentAcc_BedFileAcc_dic[ExpAcc]))
    except KeyError:
        print 'KeyError: ', tf

TF_names_nparray={}
for i in range(len(TF_names)):
    TF_names_nparray[TF_names[i]] = Merged_TF_SCls[:,i]

TF_peak_nparray={}
for i in range(len(Peak_file_names)):
    TF_peak_nparray[Peak_file_names[i]] = Merged_TF_PeakCls[:,i]

TF_names_withPeakFile = TFnamesExpAcc_withPeakFileFromSameTarget
TF_names_withPeakFile.sort()
tf_test_result_combine=[]

for tf,tf_peak_files in TF_names_withPeakFile:
#for tf in TF_names:
    #print tf,tf_peak_file
    tf_SymCals = TF_names_nparray[tf]
    tf_peak = TF_peak_nparray[tf_peak_files[0]]
    for tf_p_f in tf_peak_files[1:]:
        #print TF_peak_nparray[tf_p_f][8]
        tf_peak = tf_peak + TF_peak_nparray[tf_p_f]
    test_list = [('D','Asym'), ('D', 'Sym'), ('C', 'Asym'), ('C','Sym'), ('D','NA'),('C','NA')]
    test_result=[]
    if tf.count('SRF')>0 or tf.count('MYB')>0:
        print tf, 'D.Asym'
        D_Asym_bed= Merged_TF_bed[np.logical_and(tf_SymCals== 'Asym', np.logical_and(tf_peak >= 1,Merged_Concordance=='D'))]
        for chrom,chromStart,chromEnd  in D_Asym_bed:
            print chrom+':'+chromStart+'-'+chromEnd
        print tf, 'C.Asym'
        C_Asym_bed= Merged_TF_bed[np.logical_and(tf_SymCals== 'Asym', np.logical_and(tf_peak >= 1,Merged_Concordance=='C'))]
        for chrom,chromStart,chromEnd  in C_Asym_bed:
            print chrom+':'+chromStart+'-'+chromEnd
    for i in range(len(test_list)):
        test_result.append((tf_SymCals[np.logical_and(tf_peak>0, Merged_Concordance==test_list[i][0])] == test_list[i][1]).sum())
    tf_test_result_combine.append((tf, tf_peak_files, test_result)) #D_Asym, D_Sym, C_Asym, C_Sym, D_NA, C_NA
    #D_Asym = (tf_SymCals[np.logical_and(tf_peak =='1',Merged_Concordance=='D')] == 'Asym').sum()





with open (f_end, 'w') as out:
    #out.write('\t'.join(['TF_names', 'p.left_tail', 'p.right_tail', 'p.two_tail','s_oddsratio', 's_pvalue','C_Asym', 'C_Sym', 'D_Asym', 'D_Sym']))
    out.write('\t'.join(['TF_names', 'oddsratio', 'pvalue_twotail','D_Asym', 'D_Sym', 'C_Asym', 'C_Sym','D_NA', 'C_NA']))
    out.write('\n')
    for tf, tf_peak_files, tf_DC_count in tf_test_result_combine:
        #p = pvalue(C_Sym[tf], C_Asym[tf], D_Sym[tf], D_Asym[tf])  
        s_oddsratio, s_pvalue = stats.fisher_exact([[tf_DC_count[0], tf_DC_count[1]],[tf_DC_count[2], tf_DC_count[3]]])
        print tf, s_oddsratio, s_pvalue
        out.write('\t'.join([tf, str(s_oddsratio), str(s_pvalue)]))
        for i in tf_DC_count:
            out.write('\t%d' % i)
        out.write('\n')



#f_region = '/workdir/sc2457/TF_count_batch1N2/TSS_NOTpairs_in_Groseq_merge1K_regions.txt'
#f_int = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF.counts_TSS_nonpair_Present_PeakCallFiltered.bed'
#f_end = 'Groseq.interestingHets.noWeird.merge1K.Discordance_TF_FishersExactTest_TSS_nonpair_Present_PeakCallFiltered.txt'


def run(d):
    f_region = 'TssPair.d'+str(d)+'_dRegion250_withStrandSpecific.GroseqDiscordanceInterestingHetsNoWeird.bed'
    f_int = 'TssPair.d'+str(d)+'_dRegion250_withStrandSpecific.GroseqDiscordanceInterestingHetsNoWeird_TF.counts.PeakCallFiltered.bed'
    f_end = 'TssPair.d'+str(d)+'_dRegion250_withStrandSpecific.GroseqDiscordanceInterestingHetsNoWeird_TF.counts.PeakCallFiltered_FishersExactTest.txt'
    
    get_Tss_dl_regions_with_groseq_ASB(d, f_region, 250)  #fix the TF binding region to Tss +- 250
    get_TF_FishersExactTest_TSS_PeakCallFiltered (f_region, f_int, f_end)

d_l_list=[100,120,150,180,200]
for i in d_l_list:
    run(i)


data_dic={}
with open('Merged_TssPair_dRegion250_withStrandSpecific.GroseqDiscordanceInterestingHetsNoWeird_TF.counts.PeakCallFiltered_FishersExactTest.txt', 'w') as output:
    for d in d_l_list:
        f_end = 'TssPair.d'+str(d)+'_dRegion250_withStrandSpecific.GroseqDiscordanceInterestingHetsNoWeird_TF.counts.PeakCallFiltered_FishersExactTest.txt'
        data_dic[d] = [l.strip().split('\t') for l in open(f_end, 'U').readlines()]
            
    for d in d_l_list:
        output.write('\t'.join(data_dic[d][0][0:2] + [str(d)+'_pvalue_twotail']+data_dic[d][0][3:]))
        output.write('\t')
    output.write('\n')
    for i in range(1,len(data_dic[100])):
        for d in d_l_list:
            output.write('\t'.join(data_dic[d][i]))
            output.write('\t')
        output.write('\n')
    


