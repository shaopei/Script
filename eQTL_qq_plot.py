# use TSS_d_regions_CandD_directionality_index.py to generate a file containing the Concordant and Discordant regions
# get_Tss_dl_regions_with_groseq_ASB(150, 'GM12878_GroSeq_d150_dRegion250_withStrandSpecific.txt', 250)

#divide eQTL into three groups:
#ones in Concordant region
#ones in Discordant region
#ones near TSS but not belong to the two groups above
#the rest

import numpy as np
import scipy.stats as stats


eQTL_f = '/workdir/sc2457/eQTL/All_eQTL.gff.v3'
eQTL_chrom = np.loadtxt(eQTL_f, dtype=str ,delimiter='\t', usecols=[0], skiprows=0) #converters={0: lambda x: int(x.strip('chr'))})
#set(eQTL_chrom)
#set(['chr13', 'chr12', 'chr11', 'chr10', 'chr17', 'chr16', 'chr15', 'chr14', 'chr19', 'chr18', 'chr6_qbl_hap2', 'chr22', 'chr20', 'chr21', 'chr6_cox_hap1', 'chr7', 'chr6', 'chr5', 'chr4', 'chr3', 'chr2', 'chr1', 'chr9', 'chr8'])

eQTL_paper = np.loadtxt(eQTL_f, dtype=str ,delimiter='\t', usecols=[1], skiprows=0) 
eQTL_snppos = np.loadtxt(eQTL_f, dtype=int ,delimiter='\t', usecols=[3], skiprows=0)
eQTL_Q_vlaue = np.loadtxt(eQTL_f, dtype=float ,delimiter='\t', usecols=[5], skiprows=0)

eQTL_f2 = 'eQTL_Qvalue_cutoff_hapmap3_cis_hg19.txt'
#eQTL_f2 = 'eQTL_Qvalue_cutoff_hapmap3_trans_hg19.txt'
eQTL_chrom = np.loadtxt(eQTL_f2, dtype=str ,delimiter='\t', usecols=[1], skiprows=1)
eQTL_paper = np.loadtxt(eQTL_f2, dtype=str ,delimiter='\t', usecols=[4], skiprows=1)  #GeneSymbol
eQTL_snppos = np.loadtxt(eQTL_f2, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
eQTL_Q_vlaue = np.loadtxt(eQTL_f2, dtype=float ,delimiter='\t', usecols=[5], skiprows=1)




#-1 * np.log10()


directional_region_fp = '/workdir/sc2457/eQTL/GM12878_GroSeq_d150_dRegion250_withStrandSpecific.txt'
d_chrom = np.loadtxt(directional_region_fp, dtype=str ,delimiter='\t', usecols=[0], skiprows=1)
d_chromStart = np.loadtxt(directional_region_fp, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
d_chromEnd = np.loadtxt(directional_region_fp, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
d_Con_Discordance = np.loadtxt(directional_region_fp, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)

tss_plus_f = '/workdir/sc2457/TF_count_batch1N2/tss_paired_gm12878_plus.bed'
tss_minus_f = '/workdir/sc2457/TF_count_batch1N2/tss_paired_gm12878_minus.bed'
tss_plus = [l.strip().split('\t') for l in open(tss_plus_f, 'U').readlines()]
tss_minus = [l.strip().split('\t') for l in open(tss_minus_f, 'U').readlines()]

def if_snp_in_a_region_NOT_strandSpecific (rChrom, rStart, rEnd, sChrom=eQTL_chrom, snppos=eQTL_snppos):
    """
    sChrom and snppos are np.array
    rChrom, rStart, rEnd are integer
    """
    # keep snppos and winning_parent in the same chrom as rChrom
    c_snppos = snppos[sChrom == rChrom]
    # the ones in the region
    in_ones = np.logical_and((c_snppos - 1 - rStart >= 0 ),(rEnd - c_snppos >= 0))
    return in_ones, c_snppos[in_ones], eQTL_paper[sChrom == rChrom][in_ones], eQTL_Q_vlaue[sChrom == rChrom][in_ones]

# find eQTL in Concordant region

def snps_in_regions(x_regions):
    count = 0
    for i in xrange((len(d_chrom[x_regions]))):
        _, snp, paper, qValue = if_snp_in_a_region_NOT_strandSpecific(d_chrom[x_regions][i], d_chromStart[x_regions][i], d_chromEnd[x_regions][i])
        if len(snp)>0:
            count +=1
            print d_chrom[x_regions][i], d_chromStart[x_regions][i], d_chromEnd[x_regions][i], abs(d_chromStart[x_regions][i] - d_chromEnd[x_regions][i]), len(snp), snp, paper, qValue
    return count

snps_in_regions(d_Con_Discordance=='Concordant')
snps_in_regions(d_Con_Discordance=='Discordant')

#C_regions = d_Con_Discordance=='Concordant'
#count = 0
#for i in xrange((len(d_chrom[C_regions]))):
#    _, snp, paper, qValue = if_snp_in_a_region_NOT_strandSpecific(d_chrom[C_regions][i], d_chromStart[C_regions][i], d_chromEnd[C_regions][i])
#    if len(snp)>0:
#        count +=1
#        print d_chrom[C_regions][i], d_chromStart[C_regions][i], d_chromEnd[C_regions][i], abs(d_chromStart[C_regions][i] - d_chromEnd[C_regions][i]), len(snp), snp, paper, qValue
#print count


Tss_regions = []
d_region = 250
for i in xrange(len(tss_plus)):
    pChrom, pStart, pEnd = tss_plus[i][0:3]
    try:
        pChrom = int(pChrom.split('chr')[-1])
        pStart, pEnd = int(pStart), int(pEnd)
        mChrom, mStart, mEnd = tss_minus[i][0:3]
        mChrom = int(mChrom.split('chr')[-1])
        assert mChrom == pChrom
        mStart, mEnd = int(mStart) , int(mEnd)
        #print '\t'.join(['chr'+str(pChrom), str(min(pStart,mStart) - d_region), str(max(pEnd, mEnd) + d_region)])
        Tss_regions.append(['chr'+str(pChrom), min(pStart,mStart) - d_region, max(pEnd, mEnd) + d_region])
    except ValueError:
        pass


count = 0
for i in xrange((len(Tss_regions))):
    rChrom, rStart, rEnd = Tss_regions[i]
    _, snp, paper, qValue = if_snp_in_a_region_NOT_strandSpecific(rChrom, rStart, rEnd)
    if len(snp)>0:
        count +=1
        #print rChrom, rStart, rEnd, rEnd - rStart ,len(snp),snp , paper, qValue

C_oddsratio, C_pvalue = stats.fisher_exact([[3,65],[423,21505]])
D_oddsratio, D_pvalue = stats.fisher_exact([[5,52],[423,21505]])
CtoD_oddsratio, CtoD_pvalue = stats.fisher_exact([[3,65],[5,52]])

#cis
C_oddsratio, C_pvalue = stats.fisher_exact([[19,68-19],[1330,21928-1330]])
D_oddsratio, D_pvalue = stats.fisher_exact([[15,57-15],[1330,21928-1330]])
CtoD_oddsratio, CtoD_pvalue = stats.fisher_exact([[19,68-19],[15,57-15]])

#trans
C_oddsratio, C_pvalue = stats.fisher_exact([[3,68-3],[251,21928-251]])
D_oddsratio, D_pvalue = stats.fisher_exact([[2,57-2],[251,21928-251]])
CtoD_oddsratio, CtoD_pvalue = stats.fisher_exact([[3,68-3],[2,57-2]])

 #'_with_eOTL', '_without_eQTL', 'C_Asym', 'C_Sym'





###map snp ID to chromosome location
import glob
file_names = glob.glob('bed_chr_*.bed')

snpID_chr_map={}
for bed_file in file_names:
    with open(bed_file, 'U') as bed_f:
        data = [l.strip().split('\t') for l in bed_f.readlines()[1:]]
        for l in data:
            if l[3] in snpID_chr_map:
                print 'duplicate', l
            snpID_chr_map[l[3]] = (l[0], l[2])

eQTL_all = '/workdir/sc2457/eQTL/seeQTL/eqtl_all.txt'
eQTL_all_withChrom = '/workdir/sc2457/eQTL/seeQTL/eqtl_all_withChromLocation.txt'
with open(eQTL_all, 'U') as e:
    with open(eQTL_all_withChrom, 'w') as out:
        for l in e:
        #for l in e.readlines()[1:]:
            ll = l.strip().split('\t')
            if ll[0] in snpID_chr_map:
                out.write('\t'.join(snpID_chr_map[ll[0]]))
            else:
                out.write('\t'.join(('','')))
            out.write('\t')
            out.write(l)
            



###filter vcf file and keep snp overlape with GM12878

for i in range(1,23):
    print ('nohup cat All_20151104.vcf | awk \'BEGIN{OFS="\t"} ($1 == %s) {print $1,$2,$3,$4,$5}\' > chr%s_20151104.vcf &' % (i,i))


gmSNP_fp='/workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/snp.call.' #all'

#dbSNP_vcf_fp='/workdir/sc2457/dbSNP/'#All_20151104.vcf'  #GRCh37.p13
#dbSNP_vcf_chrom = np.loadtxt(dbSNP_vcf_fp, dtype=str ,delimiter='\t', usecols=[0], skiprows=56) 
#dbSNP_vcf_chromPos = np.loadtxt(dbSNP_vcf_fp, dtype=int ,delimiter='\t', usecols=[1], skiprows=56) 
#dbSNP_vcf_snpID_REF_ALT = np.loadtxt(dbSNP_vcf_fp, dtype=str ,delimiter='\t', usecols=[2:5], skiprows=56)


SNP_in_GM=[]
for i in xrange(1,23):
    dbSNP_vcf_fp='/workdir/sc2457/dbSNP/'+'chr'+str(i)+'_20151104.vcf'#All_20151104.vcf'  #GRCh37.p13
    print dbSNP_vcf_fp
    dbSNP_vcf_chromPos = np.loadtxt(dbSNP_vcf_fp, dtype=int ,delimiter='\t', usecols=[1], skiprows=56) 
    dbSNP_vcf_snpID_REF_ALT = np.loadtxt(dbSNP_vcf_fp, dtype=str ,delimiter='\t', usecols=[2,3,4], skiprows=56)
    print gmSNP_fp+str(i)
    gmSNP_data = np.loadtxt(gmSNP_fp+str(i), dtype=str ,delimiter='\t', usecols=[0,1,2,3,4,5], skiprows=0)
    with open(gmSNP_fp+str(i)+'_withSNP_ID', 'w') as out:
        out.write('\t'.join(['Chrom', 'ChromPos','snpID','REF','ALT', 'mat','pat']))
        out.write('\n')
        for j in xrange(gmSNP_data.shape[0]):
            #print j
            _,pos,ref,_,_,genotype = gmSNP_data[j]
            m, d = genotype[0], genotype[1]
            #print i, pos, ref, genotype, m, d
            d_pos = dbSNP_vcf_chromPos == int(pos)
            if d_pos.any:
                d_snpID_REF_ALT = dbSNP_vcf_snpID_REF_ALT[d_pos]
                in_ones = np.logical_and((d_snpID_REF_ALT[:,1] == ref), np.logical_or(d_snpID_REF_ALT[:,2] == m, d_snpID_REF_ALT[:,2] == d))
                #print i, dbSNP_vcf_chromPos[d_pos][in_ones], d_snpID_REF_ALT[in_ones]
                if len(d_snpID_REF_ALT[in_ones]) ==1:
                    out.write('\t'.join([str(i)]+[str(dbSNP_vcf_chromPos[d_pos][in_ones][0])]+list(d_snpID_REF_ALT[in_ones][0])+[m,d]))
                    out.write('\n')
                    SNP_in_GM.append(d_snpID_REF_ALT[in_ones][0,0])
                else:
                    for k in xrange(len(d_snpID_REF_ALT[in_ones])):
                        out.write('\t'.join([str(i)]+[str(dbSNP_vcf_chromPos[d_pos][in_ones][k])]+list(d_snpID_REF_ALT[in_ones][k])+[m,d]))
                        out.write('\n')
                        SNP_in_GM.append(d_snpID_REF_ALT[in_ones][i,0])
            
        
        
        




gmSNP_fp='/workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/snp.call.' #all'
import multiprocessing
def annotate_SNPID(i):
    #SNP_in_GM_seperate=[]
    dbSNP_vcf_fp='/workdir/sc2457/dbSNP/'+'chr'+str(i)+'_20151104.vcf'#All_20151104.vcf'  #GRCh37.p13
    print dbSNP_vcf_fp
    dbSNP_vcf_chromPos = np.loadtxt(dbSNP_vcf_fp, dtype=int ,delimiter='\t', usecols=[1], skiprows=56) 
    dbSNP_vcf_snpID_REF_ALT = np.loadtxt(dbSNP_vcf_fp, dtype=str ,delimiter='\t', usecols=[2,3,4], skiprows=56)
    print gmSNP_fp+str(i)
    #gmSNP_data = np.loadtxt(gmSNP_fp+str(i), dtype=str ,delimiter='\t', usecols=[0,1,2,3,4,5], skiprows=0)
    with open(gmSNP_fp+str(i)+'_withSnpID', 'w') as out:
        out.write('\t'.join(['Chrom', 'ChromPos','snpID','REF','ALT', 'mat','pat']))
        out.write('\n')
        for l in open(gmSNP_fp+str(i), 'U'):
        #for j in xrange(gmSNP_data.shape[0]):
            #print j
            _,pos,ref,_,_,genotype,_ = l.strip().split('\t') #gmSNP_data[j]
            m, d = genotype[0], genotype[1]
            #print i, pos, ref, genotype, m, d
            d_pos = dbSNP_vcf_chromPos == int(pos)
            if d_pos.any:
                d_snpID_REF_ALT = dbSNP_vcf_snpID_REF_ALT[d_pos]
                in_ones = np.logical_and((d_snpID_REF_ALT[:,1] == ref), np.logical_or(d_snpID_REF_ALT[:,2] == m, d_snpID_REF_ALT[:,2] == d))
                #print i, dbSNP_vcf_chromPos[d_pos][in_ones], d_snpID_REF_ALT[in_ones]
                if len(d_snpID_REF_ALT[in_ones]) ==1:
                    out.write('\t'.join([str(i)]+[str(dbSNP_vcf_chromPos[d_pos][in_ones][0])]+list(d_snpID_REF_ALT[in_ones][0])+[m,d]))
                    out.write('\n')
                    SNP_in_GM_seperate.append(d_snpID_REF_ALT[in_ones][0,0])
                elif len(d_snpID_REF_ALT[in_ones]) >1:
                    for k in xrange(len(d_snpID_REF_ALT[in_ones])):
                        out.write('\t'.join([str(i)]+[str(dbSNP_vcf_chromPos[d_pos][in_ones][k])]+list(d_snpID_REF_ALT[in_ones][k])+[m,d]))
                        out.write('\n')
                        SNP_in_GM_seperate.append(d_snpID_REF_ALT[in_ones][k,0])
    #return SNP_in_GM_seperate

pool = multiprocessing.Pool(processes=21)
#pool_output = pool.map(subject_occurrence_per_file_memory_demanding , f_list )
pool_output = pool.map(annotate_SNPID, range(1,22))
# pool_output looks like [(f1, q1,s1), (f2, q2,s2),...]
pool.close() # no more tasks
pool.join()



### this takes too long, use sort, join in linux instead ###
#eQTL_all = '/workdir/sc2457/eQTL/seeQTL/eqtl_all.txt'
##eQTL_all_readin = np.loadtxt(eQTL_all, dtype=str ,delimiter='\t', usecols=range(0,5), skiprows=1)
#gmSNP_fp = '/workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/withSNP_ID/'
#import numpy as np
#import multiprocessing
#def annotate_SNPID_eQTL(i):
#    print gmSNP_fp+'snp.call.'+str(i)+'_withSnpID'
#    gmSNP_fp_readin = np.loadtxt(gmSNP_fp+'snp.call.'+str(i)+'_withSnpID', dtype=str ,delimiter='\t', usecols=range(0,7), skiprows=1)
#    with open(gmSNP_fp+'snp.call.'+str(i)+'_withSNP.id_eQTL.pValue', 'w') as out:
#        out.write('\t'.join(['Chrom', 'ChromPos','snpID','REF','ALT', 'mat','pat','gene','t-stat','p-value','FDR']))
#        out.write('\n')
#        for l in open(eQTL_all, 'U'):
#            ll = l.strip().split('\t')
#            ins = gmSNP_fp_readin[:,2] == ll[0]  #ll[0] is snpID
#            if ins.any():
#                assert len(gmSNP_fp_readin[ins]) == 1,list(gmSNP_fp_readin[ins][0])+ll
#                #print gmSNP_fp_readin[ins], ll
#                out.write('\t'.join(list(gmSNP_fp_readin[ins][0])+ll[1:]))
#                out.write('\n')
#    print gmSNP_fp+'snp.call.'+str(i)+'_withSnpID', 'Done'
#
#import multiprocessing
#pool = multiprocessing.Pool(processes=12)
#pool_output = pool.map(annotate_SNPID_eQTL, range(1,10))
#pool.close() # no more tasks
#pool.join()
### this takes too long, use sort, join in linux instead ###



####

##annotate study according line number##
eQTL_all = '/workdir/sc2457/eQTL/seeQTL/eqtl_all.txt'
data = [l.strip().split('\t') for l in open(eQTL_all, 'U').readlines()[1:]]


l_data=open('/workdir/sc2457/eQTL/seeQTL/eqtl_all_study_line_number.txt', 'U').readlines()
line_index=[]
for l in l_data:
    print l
    ll = l.strip().split(':')
    line_number = ll[0]
    study_name = ll[1].split('eqtl/eQTL_results_')[-1].split('.txt')[0]
    print study_name, line_number
    line_index.append((study_name, line_number))

with open('/workdir/sc2457/eQTL/seeQTL/eqtl_all_study.txt', 'w') as out:
    for i in xrange(len(line_index)-1):
        start_line = int(line_index[i][1])-1
        end_line = int(line_index[i+1][1])-2
        print line_index[i][0], start_line, end_line
        for j in xrange(start_line, end_line):
            if len(data[j]) == 5:
                out.write('\t'.join(data[j]+[line_index[i][0]]))
                out.write('\n')
##annotate study according line number##


    
## for eacg snpID, keep the line with lowerst p-value
with open('/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers.txt', 'U') as input_f:
    with open('/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_minPvalue.txt', 'w') as out:
        snpID, p_value = 'rs1000000', 0.343128287974117
        the_line_with_lowest_pValue = ''
        for l in input_f:
            ll=l.strip().split('\t')
            if snpID == ll[2]:
                if float(ll[10]) < p_value:
                    the_line_with_lowest_pValue = l
                    p_value = float(ll[10])
            else:
                print snpID, p_value
                #print the_line_with_lowest_pValue
                out.write(the_line_with_lowest_pValue)
                snpID = ll[2]
                p_value = float(ll[10])
                the_line_with_lowest_pValue = l



## for each snpID, keep the line with lowerst FDR
with open('/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll.txt', 'U') as input_f:
    with open('/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_minFDR.txt', 'w') as out:
        snpID, fdr = 'rs1000000', 0.986078756632218
        the_line_with_lowest_FDR = ''
        for l in input_f:
            ll=l.strip().split('\t')
            if snpID == ll[2]:
                if float(ll[11]) < fdr:
                    the_line_with_lowest_FDR = l
                    fdr = float(ll[11])
            else:
                print snpID, fdr
                #print the_line_with_lowest_pValue
                out.write(the_line_with_lowest_FDR)
                snpID = ll[2]
                fdr = float(ll[11])
                the_line_with_lowest_FDR = l


###find snps within Concordant, Dicordant, or OtherTssPair regions


import numpy as np
#eQTL_f = '/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_minPvalue_ChromSorted.txt'
#eQTL_f = '/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_minPvalue_ChromSorted.txt'
#eQTL_f = '/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers.txt'
eQTL_f='/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_geneWithin1K.txt'  #only use the snp gene pairs within 1K
eQTL_chrom = np.loadtxt(eQTL_f, dtype=int ,delimiter='\t', usecols=[0], skiprows=0)
eQTL_paper = np.loadtxt(eQTL_f, dtype=str ,delimiter='\t', usecols=[12], skiprows=0) 
eQTL_snppos = np.loadtxt(eQTL_f, dtype=int ,delimiter='\t', usecols=[1], skiprows=0)
eQTL_snpID = np.loadtxt(eQTL_f, dtype=str ,delimiter='\t', usecols=[2], skiprows=0)
eQTL_p_vlaue = np.loadtxt(eQTL_f, dtype=str ,delimiter='\t', usecols=[10], skiprows=0)
eQTL_FDR = np.loadtxt(eQTL_f, dtype=str ,delimiter='\t', usecols=[11], skiprows=0)

directional_region_fp = '/workdir/sc2457/eQTL/GM12878_GroSeq_d150_dRegion250_withStrandSpecific.txt'
#directional_region_fp = '/workdir/sc2457/Groseq.gm12878_Allele_seq_result/StranSpecificASE/fdr0.2/GM12878_GroSeq_d150_dRegion250_withStrandSpecific_fdr0.2.txt'
d_chrom = np.loadtxt(directional_region_fp, dtype=str ,delimiter='\t', usecols=[0], skiprows=1)
d_chromStart = np.loadtxt(directional_region_fp, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
d_chromEnd = np.loadtxt(directional_region_fp, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
d_Con_Discordance = np.loadtxt(directional_region_fp, dtype=str ,delimiter='\t', usecols=[3], skiprows=1)

tss_plus_f = '/workdir/sc2457/TF_count_batch1N2/tss_paired_gm12878_plus.bed'
tss_minus_f = '/workdir/sc2457/TF_count_batch1N2/tss_paired_gm12878_minus.bed'
tss_plus = [l.strip().split('\t') for l in open(tss_plus_f, 'U').readlines()]
tss_minus = [l.strip().split('\t') for l in open(tss_minus_f, 'U').readlines()]

def if_snp_in_a_region_NOT_strandSpecific (rChrom, rStart, rEnd, sChrom=eQTL_chrom, snppos=eQTL_snppos):
    """
    sChrom and snppos are np.array
    rChrom, rStart, rEnd are integer
    """
    rChrom = int(rChrom.split('chr')[-1])
    c_snppos = snppos[sChrom == rChrom]
    # the ones in the region
    in_ones = np.logical_and((c_snppos - 1 - rStart >= 0 ),(rEnd - c_snppos >= 0))
    return in_ones, c_snppos[in_ones], eQTL_paper[sChrom == rChrom][in_ones], eQTL_p_vlaue[sChrom == rChrom][in_ones], eQTL_FDR[sChrom == rChrom][in_ones], eQTL_snpID[sChrom == rChrom][in_ones]


# find eQTL in Concordant region

def snps_in_regions_prep(x_regions):
    snpID_list=[]
    count = 0
    for i in xrange((len(d_chrom[x_regions]))):
        _, snp, paper, pValue, fdr, snpid = if_snp_in_a_region_NOT_strandSpecific(d_chrom[x_regions][i], d_chromStart[x_regions][i], d_chromEnd[x_regions][i])
        if len(snp)>0:
            count +=1
            for k in xrange(len(snp)):
                print d_chrom[x_regions][i], d_chromStart[x_regions][i], d_chromEnd[x_regions][i], abs(d_chromStart[x_regions][i] - d_chromEnd[x_regions][i]), len(snp), snp[k], paper[k], pValue[k], fdr[k],snpid[k]
                snpID_list.append(snpid[k])
    return count, snpID_list

C_count, C_snpID_list = snps_in_regions_prep(d_Con_Discordance=='Concordant')
D_count, D_snpID_list = snps_in_regions_prep(d_Con_Discordance=='Discordant')
#C_count =31
# D_count = 28

CD_overlap_snpID_list=[]
for s in C_snpID_list:
    if s in D_snpID_list:
        print s
        CD_overlap_snpID_list.append(s)


# remove SNPs in both Concordant and Discordant regions
def snps_in_regions(x_regions, out_fp, list_to_exclude):
    snpID_list=[]
    with open(out_fp, 'w') as out:
        out.write('\t'.join(['chrom','chromStart','chromEnd', 'region_length', 'SnpNo', 'SnpLoc', 'study', 'p_value', 'FDR', 'SnpID']))
        out.write('\n')
        count = 0
        for i in xrange((len(d_chrom[x_regions]))):
            _, snp, paper, pValue, fdr, snpid = if_snp_in_a_region_NOT_strandSpecific(d_chrom[x_regions][i], d_chromStart[x_regions][i], d_chromEnd[x_regions][i])
            if len(snp)>0:
                count +=1
                for k in xrange(len(snp)):
                    if snpid[k] not in list_to_exclude:
                        out.write('\t'.join([str(d_chrom[x_regions][i]), str(d_chromStart[x_regions][i]), str(d_chromEnd[x_regions][i]), str(abs(d_chromStart[x_regions][i] - d_chromEnd[x_regions][i])), str(len(snp)), str(snp[k]), str(paper[k]), pValue[k], fdr[k],snpid[k] ]))
                        out.write('\n')
                        print d_chrom[x_regions][i], d_chromStart[x_regions][i], d_chromEnd[x_regions][i], abs(d_chromStart[x_regions][i] - d_chromEnd[x_regions][i]), len(snp), snp[k], paper[k], pValue[k], fdr[k],snpid[k]
                        snpID_list.append(snpid[k])
    return count, snpID_list

# the output Concordant and Discordant file has no overlap in SNPs
C_count, C_snpID_list = snps_in_regions(d_Con_Discordance=='Concordant', eQTL_f+'_Concordant', CD_overlap_snpID_list)
D_count, D_snpID_list =snps_in_regions(d_Con_Discordance=='Discordant', eQTL_f+'_Discordant', CD_overlap_snpID_list)

CD_overlap_snpID_list_after=[]
for s in C_snpID_list:
    if s in D_snpID_list:
        print s
        CD_overlap_snpID_list_after.append(s) # len(CD_overlap_snpID_list_after) should be zero

Tss_regions = []
d_region = 250
for i in xrange(len(tss_plus)):
    pChrom, pStart, pEnd = tss_plus[i][0:3]
    try:
        pChrom = int(pChrom.split('chr')[-1])
        pStart, pEnd = int(pStart), int(pEnd)
        mChrom, mStart, mEnd = tss_minus[i][0:3]
        mChrom = int(mChrom.split('chr')[-1])
        assert mChrom == pChrom
        mStart, mEnd = int(mStart) , int(mEnd)
        #print '\t'.join(['chr'+str(pChrom), str(min(pStart,mStart) - d_region), str(max(pEnd, mEnd) + d_region)])
        Tss_regions.append(['chr'+str(pChrom), min(pStart,mStart) - d_region, max(pEnd, mEnd) + d_region])
    except ValueError:
        pass
#len(Tss_regions) = 21928

# remove Concordant and Discordant regions
for di in xrange((len(d_chrom))):
    Tss_regions.remove([d_chrom[di],d_chromStart[di], d_chromEnd[di]])

count = 0
with open(eQTL_f+'_otherTssPairs', 'w') as out:
    out.write('\t'.join(['chrom','chromStart','chromEnd', 'region_length', 'SnpNo', 'SnpLoc', 'study', 'p_value','FDR','SnpID']))
    out.write('\n')
    for i in xrange((len(Tss_regions))):
        rChrom, rStart, rEnd = Tss_regions[i]
        _, snp, paper, pValue, fdr,snpid  = if_snp_in_a_region_NOT_strandSpecific(rChrom, rStart, rEnd)
        if len(snp)>0:
            count +=1
            for k in xrange(len(snp)):
                out.write('\t'.join([rChrom, str(rStart), str(rEnd), str(abs(rStart - rEnd )), str(len(snp)), str(snp[k]), str(paper[k]), pValue[k], fdr[k], snpid[k] ]))
                out.write('\n')
                print rChrom, rStart, rEnd, rEnd - rStart ,len(snp),snp , paper, pValue, fdr, snpid 


# each region, only keep one SNP with lowest p-value
def keep_the_SNP_with_smallest_p_value(input_fp):#, output_fp):
    with open(input_fp, 'U') as input_f:
        with open(input_fp+'_SmallestPvalue', 'w') as out:
            out.write('\t'.join(['chrom','chromStart','chromEnd', 'region_length', 'SnpNo', 'SnpLoc', 'study', 'p_value','FDR','SnpID']))
            out.write('\n')
            chrom, chromStart, chromEnd = input_f.readline().strip().split('\t')[0:3]
            p_value=1
            the_line_with_lowest_pValue = ''
            for l in input_f:
                ll=l.strip().split('\t')
                if ll[0] ==chrom and ll[1] == chromStart and ll[2] == chromEnd:
                    if float(ll[7]) < p_value:
                        the_line_with_lowest_pValue = l
                        p_value = float(ll[7])
                else:
                    out.write(the_line_with_lowest_pValue)
                    chrom, chromStart, chromEnd = ll[0:3]
                    p_value = float(ll[7])
                    the_line_with_lowest_pValue = l
                    


name_list=['_Concordant', '_Discordant','_otherTssPairs']
for name in name_list:
    keep_the_SNP_with_smallest_p_value(eQTL_f+name)



### keep the line with snpID in Concordant region
import numpy as np

def keep_the_SNP_in_the_set(input_fp, output_fp, set_name):
    with open(input_fp, 'U') as input_f:
        with open(output_fp, 'w') as out:
            for l in input_f:
                ll=l.strip().split('\t')
                if ll[2] in set_name:
                    out.write(l)

Concordant_Snp_ID = np.loadtxt('/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_minPvalue_ChromSorted.txt_Concordant', dtype=str ,delimiter='\t', usecols=[9], skiprows=1)
keep_the_SNP_in_the_set('eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers.txt', 'eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_Concordant.txt', set(Concordant_Snp_ID))


Discordant_Snp_ID = np.loadtxt('/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_minPvalue_ChromSorted.txt_Discordant', dtype=str ,delimiter='\t', usecols=[9], skiprows=1)
keep_the_SNP_in_the_set('eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers.txt', 'eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_Discordant.txt', set(Discordant_Snp_ID))

otherTssPairs_Snp_ID = np.loadtxt('/workdir/sc2457/eQTL/seeQTL/eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_minPvalue_ChromSorted.txt_otherTssPairs', dtype=str ,delimiter='\t', usecols=[9], skiprows=1)
keep_the_SNP_in_the_set('eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers.txt', 'eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_otherTssPairs.txt', set(otherTssPairs_Snp_ID))

EntrezGeneID_GencodeGeneID_map ={}
with open('/workdir/sc2457/eQTL/seeQTL/Homo_sapiens.gene_info') as input_f:
    for l in input_f:
        ll = l.strip().split('\t')
        dbXrefs = ll[5].split('|')
        print dbXrefs
        for d in dbXrefs:
            if d.split(':')[0] == 'Ensembl':
                EntrezGeneID_GencodeGeneID_map[ll[1]] = d.split(':')[1]
                print d.split(':')

GencodeGeneID_ChromLocation_map={}
with open('/workdir/sc2457/eQTL/seeQTL/gencode.v19.annotation.gtf_genes') as input_f:
    for l in input_f:
        print l
        ll = l.strip().split(' ')
        GencodeGeneID_ChromLocation_map[ll[4].split('.')[0]] = ll[0:3]

snp_gene_map={}
with open('eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers.txt', 'U') as input_f:
    with open('eqtl_all_cis_study_inGM12878SnpCallAll_RemoveMyers_geneWithin1K.txt', 'w') as out:
        for l in input_f:
            ll = l.strip().split('\t')
            if ll[8] in EntrezGeneID_GencodeGeneID_map: #ll[8] is the EntrezGeneID
                if EntrezGeneID_GencodeGeneID_map[ll[8]] in GencodeGeneID_ChromLocation_map:
                    chromosome, chromoStart, chromoEnd = GencodeGeneID_ChromLocation_map[EntrezGeneID_GencodeGeneID_map[ll[8]]]
                    if chromosome.split('chr')[-1] == ll[0]:  #ll[0] is chromosome number eg. 1,2,..,22
                        if abs(int(ll[1]) - int(chromoStart)) < 1000:
                            out.write(l)
                            if ll[2] not in snp_gene_map:
                                snp_gene_map[ll[2]] = [(ll[8],abs(int(ll[1]) - int(chromoStart)))]
                            elif (ll[8],abs(int(ll[1]) - int(chromoStart))) not in snp_gene_map[ll[2]]:
                                snp_gene_map[ll[2]].append((ll[8],abs(int(ll[1]) - int(chromoStart))))

for p in snp_gene_map:
    if len(snp_gene_map[p]) > 1:
        print p, snp_gene_map[p]


