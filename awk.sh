awk 'BEGIN{OFS="\t"} ($15 <= 0.05){print $0}'

awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}'
awk 'BEGIN{OFS="\t"} ($16 <= 0.05){print $0}'

awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}'

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6=="+"?"-":"+"}' 

awk 'BEGIN{OFS="\t"} ($4 != $5) {print $0}'


awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}'

cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == $i) {print $1,$2,$3,$4,$5}' > chr$s_20151104.vcf

nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 2) {print $1,$2,$3,$4,$5}' > chr2_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 3) {print $1,$2,$3,$4,$5}' > chr3_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 4) {print $1,$2,$3,$4,$5}' > chr4_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 5) {print $1,$2,$3,$4,$5}' > chr5_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 6) {print $1,$2,$3,$4,$5}' > chr6_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 7) {print $1,$2,$3,$4,$5}' > chr7_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 8) {print $1,$2,$3,$4,$5}' > chr8_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 9) {print $1,$2,$3,$4,$5}' > chr9_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 10) {print $1,$2,$3,$4,$5}' > chr10_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 11) {print $1,$2,$3,$4,$5}' > chr11_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 12) {print $1,$2,$3,$4,$5}' > chr12_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 13) {print $1,$2,$3,$4,$5}' > chr13_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 14) {print $1,$2,$3,$4,$5}' > chr14_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 15) {print $1,$2,$3,$4,$5}' > chr15_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 16) {print $1,$2,$3,$4,$5}' > chr16_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 17) {print $1,$2,$3,$4,$5}' > chr17_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 18) {print $1,$2,$3,$4,$5}' > chr18_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 19) {print $1,$2,$3,$4,$5}' > chr19_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 20) {print $1,$2,$3,$4,$5}' > chr20_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 21) {print $1,$2,$3,$4,$5}' > chr21_20151104.vcf &
nohup cat All_20151104.vcf | awk 'BEGIN{OFS="\t"} ($1 == 22) {print $1,$2,$3,$4,$5}' > chr22_20151104.vcf &



awk '{print NR,$0}' eqtl_all.txt > eqtl_all_lineNumber.txt


zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz | awk 'BEGIN{OFS="\t"} (NR >103 && $18 >0.5 && $18<1.5) {print "chr"$1,$2-1, $2,$14,$15,$16}'


zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz | awk 'BEGIN{OFS="\t"} (NR <104){print $0}; (NR >103){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$25}' \
| awk 'BEGIN{OFS="\t"} (NR <105){print $0}; (NR >104 && $7=="PASS"){print $0}'

 |head -n 120


conditional-expression ? action1 : action2 ;


zcat SRR4041365_rcR1.fastq.gz| awk '(NR%2 ==0){print $0}; (NR%2==1){print rev($0)}'
