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
    cACGT = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=[7,8,9,10,11,12], skiprows=1) #mat_all pat_all cA cC cG cT
    return temp[:,0], temp[:,1], winning_parent,cACGT #chr, snppos, winning_parent

#plus_groseq_fp = '/workdir/sc2457/Groseq.gm12878_Allele_seq_result/StranSpecificASE/fdr0.2/interestingHets_noW_plus.txt'
#minus_groseq_fp = '/workdir/sc2457/Groseq.gm12878_Allele_seq_result/StranSpecificASE/fdr0.2/interestingHets_noW_minus.txt'
plus_groseq_fp = 'interestingHets_noW_plus.txt'
minus_groseq_fp = 'interestingHets_noW_minus.txt'


plus_chrom, plus_snppos, plus_winning_parent,plus_cACGT = read_file_to_nparray(plus_groseq_fp)
minus_chrom, minus_snppos, minus_winning_parent, minus_cACGT = read_file_to_nparray(minus_groseq_fp)


tss_plus_f = '/workdir/sc2457/TF_count_batch1N2/tss_paired_gm12878_plus.bed'
tss_minus_f = '/workdir/sc2457/TF_count_batch1N2/tss_paired_gm12878_minus.bed'


tss_plus = [l.strip().split('\t') for l in open(tss_plus_f, 'U').readlines()]
tss_minus = [l.strip().split('\t') for l in open(tss_minus_f, 'U').readlines()]

#
#def if_snp_in_a_region_NOT_strandSpecific (rChrom, rStart, rEnd,sChrom=snp_chrom, snppos=snp_snppos, winning_parent=winning_parent):
#    """
#    sChrom and snppos are np.array
#    rChrom, rStart, rEnd are integer
#    """
#    
#    # keep snppos and winning_parent in the same chrom as rChrom
#    c_snppos = snppos[sChrom == rChrom]
#    c_winning_parent = winning_parent[sChrom == rChrom]
#    # the ones in the region
#    in_ones = np.logical_and((c_snppos - 1 - rStart >= 0 ),(rEnd - c_snppos >= 0))
#    return in_ones, c_snppos[in_ones], c_winning_parent[in_ones]


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
    #print pcACGT_use
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


 #mat_all pat_all cA cC cG cT
ACGT_map={'A':2, 'C':3, 'G':4, 'T':5}
def directionality(pcACGT_use, mcACGT_use):
    mat_count_plus = pcACGT_use[ACGT_map[pcACGT_use[0]]] 
    mat_count_minus =mcACGT_use[ACGT_map[mcACGT_use[0]]]
    pat_count_plus = pcACGT_use[ACGT_map[pcACGT_use[1]]]
    pat_count_minus =mcACGT_use[ACGT_map[mcACGT_use[1]]]
    return [mat_count_plus, pat_count_plus,mat_count_minus,pat_count_minus]
    
    
    
    
    





def get_Tss_dl_regions_with_groseq_ASB(d_l, out_fp, d_region):
    with open(out_fp, 'w') as out:
        out.write('\t'.join(['chrom','chromStart','chromEnd','Groseq.Discordance',
                             'Tss_plus_winParent','Tss_minus_winParent',
                             'mat_count_plus', 'pat_count_plus','mat_count_minus',  'pat_count_minus' ,
                             'plus_Mat', 'plus_Pat',
                             'Tss_plus_A','Tss_plus_C','Tss_plus_G','Tss_plus_T',
                             'minus_Mat', 'minus_Pat',
                             'Tss_minus_A','Tss_minus_C','Tss_minus_G','Tss_minus_T'])
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
                        out.write('\t'.join(['chr'+str(pChrom), str(min(pStart,mStart) - d_region), str(max(pEnd, mEnd) + d_region),'Concordant', pw_use, mw_use]+directionality(pcACGT_use, mcACGT_use)+list(pcACGT_use)+list(mcACGT_use))
                                  +'\n')
                    elif pw_use != mw_use:
                        #print directionality(pcACGT_use, mcACGT_use)
                        Dis_count +=1
                        out.write('\t'.join(['chr'+str(pChrom), str(min(pStart,mStart) - d_region), str(max(pEnd, mEnd) + d_region),'Discordant', pw_use, mw_use]+directionality(pcACGT_use, mcACGT_use)+list(pcACGT_use)+list(mcACGT_use))
                                  +'\n')
                        #print 'D:', str(pChrom)+':'+str(pStart)+'-'+str(pEnd), pw_use, mw_use
                    #print 'chr'+str(pChrom), str(min(pStart,mStart) - d_l), str(max(pEnd, mEnd) + d_l), pw_use, mw_use
                    #print count, pChrom, str(pStart)+'-'+str(pEnd), pw, mw
            except ValueError:
                pass
    print '\t'.join(['d='+str(d_l), 'Con', str(Con_count), 'Dis', str(Dis_count)])


def compare_two_CD_region_list():
    file1_fp='/workdir/sc2457/Groseq.gm12878_Allele_seq_result/StranSpecificASE/BetaBinfdr0.2/GM12878_GroSeq_d150_dRegion250_withStrandSpecific.txt'
    file2_fp='/workdir/sc2457/Groseq.gm12878_Allele_seq_result/StranSpecificASE/fdr0.1/GM12878_GroSeq_d150_dRegion250_withStrandSpecific_fdr0.1.txt'
    
    file1_regions_C_set=set()
    file2_regions_C_set=set()
    overlap_regions_C_set=set()
    file1_regions_D_set=set()
    file2_regions_D_set=set()
    overlap_regions_D_set=set()
    
    
    with open(file1_fp, 'U') as f1:
        for l in f1:
            ll = l.strip().split('\t')
            if ll[3] == 'Concordant':
                file1_regions_C_set.add((ll[0], ll[1], ll[2], ll[3]))
            elif ll[3] == 'Discordant':
                file1_regions_D_set.add((ll[0], ll[1], ll[2], ll[3]))
    
    
    with open(file2_fp, 'U') as f2:
        for l in f2:
            ll = l.strip().split('\t')
            if ll[3] == 'Concordant':
                file2_regions_C_set.add((ll[0], ll[1], ll[2], ll[3]))
                if (ll[0], ll[1], ll[2], ll[3]) in file1_regions_C_set:
                    overlap_regions_C_set.add((ll[0], ll[1], ll[2], ll[3]))
            elif ll[3] == 'Discordant':
                file2_regions_D_set.add((ll[0], ll[1], ll[2], ll[3]))
                if (ll[0], ll[1], ll[2], ll[3]) in file1_regions_D_set:
                    overlap_regions_D_set.add((ll[0], ll[1], ll[2], ll[3]))
    
    len(file1_regions_D_set)
    len(file2_regions_D_set)
    len(overlap_regions_C_set)
    len(overlap_regions_D_set)
##








if __name__ == '__main__':
    get_Tss_dl_regions_with_groseq_ASB(150, 'GM12878_GroSeq_d150_dRegion250_withStrandSpecific.txt', 250) 


