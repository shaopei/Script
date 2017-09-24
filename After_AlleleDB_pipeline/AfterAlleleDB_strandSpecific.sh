Min_count=$1
echo "Min_count         =$1" 
Max_Pvalus=1
plus_input_file=interestingHets_plus.txt
minus_input_file=interestingHets_minus.txt
PL=/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline

tss_plus_f=${PL}/tss_paired_gm12878_plus.bed
tss_minus_f=${PL}/tss_paired_gm12878_minus.bed



R --vanilla --slave --args $(pwd) ${plus_input_file} ${Min_count} ${Max_Pvalus} < ${PL}/filter_counts_file.R 
R --vanilla --slave --args $(pwd) ${minus_input_file} ${Min_count} ${Max_Pvalus} < ${PL}/filter_counts_file.R 

python ${PL}/TSS_d_regions_CandD_directionality_index_20170612.py ${plus_input_file:0:-4}_MinCount${Min_count}_MaxPvalue${Max_Pvalus}.txt \
${minus_input_file:0:-4}_MinCount${Min_count}_MaxPvalue${Max_Pvalus}.txt ${tss_plus_f} ${tss_minus_f} \
GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount${Min_count}_MaxPvalue${Max_Pvalus}.txt

#python /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/TSS_d_regions_CandD_directionality_index_20170612.py interestingHets_plusNoW.txt interestingHets_minusNoW.txt /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/tss_paired_gm12878_plus.bed /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/tss_paired_gm12878_minus.bed test.txt

