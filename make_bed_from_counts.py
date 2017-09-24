import numpy as np

def switch_bed_col (f, length, target_col_index, output_fp):
    temp = np.loadtxt(f, dtype=str ,delimiter='\t', usecols=xrange(length), skiprows=0) #chr1-chr22
    strand = np.copy(temp[:,5].reshape(temp.shape[0], 1))
    strand [strand == '+'] = 'plus'
    strand [strand == '-'] = 'minus'
    a = np.append(temp, strand, axis=1)
    np.savetxt(output_fp, a[:,target_col_index], fmt='%s', delimiter='\t')


#switch_bed_col ('tss_all_gm12878.bed', 6, [0,1,2,6,4,5,3], 'tss_all_gm12878_strand.bed')


def read_file_to_nparray (snp_f):
    temp = np.loadtxt(snp_f, dtype=int ,delimiter='\t', usecols=[0,1], skiprows=1) #chr1-chr22
    winning_parent_SymCls = np.loadtxt(snp_f, dtype=str ,delimiter='\t', usecols=[13, 14], skiprows=1)
    return temp[:,0], temp[:,1], winning_parent_SymCls #chr, snppos, winning_parent


def bed_from_counts_file(f, output_fp):
    chrm, chrmEnd, winning_parent_SymCls = read_file_to_nparray (f)
    chrmStart = chrmEnd -1
    winP = winning_parent_SymCls [:,0]
    winP[winning_parent_SymCls [:,1] == 'Sym'] = 'Sym'
    np.savetxt(output_fp, np.column_stack((chrm, chrmStart,chrmEnd, winP)), fmt='%s', delimiter='\t')

#bed_from_counts_file('interestingHets_noW_minus.txt', 'Groseq_interestingHets_noW_minus.bed' )
bed_from_counts_file('SRF_interestingHets.txt', 'SRF_batch1N2_interestingHets.bed')

bed_from_counts_file('SRF.ENCSR000BGE_counts.txt', 'SRF.ENCSR000BGE_counts.bed')
bed_from_counts_file('SRF.ENCSR000BMI_counts.txt', 'SRF.ENCSR000BMI_counts.bed')