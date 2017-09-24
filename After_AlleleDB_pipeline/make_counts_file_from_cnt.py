fastq_metadata_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/metadata_fastq.tsv' # for downloaded fastq
output_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_Target_counts/'
Min_count = '5'
PL = '/workdir/sc2457/alleleDB/alleledb_pipeline/'
BASE='/workdir/sc2457/'
SNPS =BASE+'SNP/1000genome_vol1.ftp.release.20130502/snp.call.all'
CNVS= BASE+'alleleseq.gersteinlab.org/NA12878_diploid_dec16.2012.alleleseq.input/rd_4327183snps_na12878_hg19.txt'
BNDS= 'hits.bed'


#make dictionary from fastq_metadata
fastq_metadata = [l.strip().split('\t') for l in open(fastq_metadata_fp, 'U').readlines()]
ExperimentAcc_FastqFileAcc_dic = {}
FastqFileAcc_ExperimentAcc_dic = {}
Target_FastqFileAcc_dic = {}
FastqFileAcc_Target_dic = {}
ExperimentAcc_Target_dic={}
Target_ExperimentAcc_dic={}
for l in fastq_metadata[1:]:
    experiment_accession = l[3]
    file_accession = l[0]+'.cnt'
    target = l[16].split('-')[0]
    if experiment_accession not in ExperimentAcc_FastqFileAcc_dic:
        ExperimentAcc_FastqFileAcc_dic[experiment_accession] = []
    ExperimentAcc_FastqFileAcc_dic[experiment_accession].append(file_accession)
    FastqFileAcc_ExperimentAcc_dic[file_accession] = experiment_accession
    if target not in Target_FastqFileAcc_dic:
        Target_FastqFileAcc_dic[target] = []
    Target_FastqFileAcc_dic[target].append(file_accession)
    FastqFileAcc_Target_dic[file_accession] = target
    ExperimentAcc_Target_dic[experiment_accession] = target
    if target not in Target_ExperimentAcc_dic:
        Target_ExperimentAcc_dic[target] =set()
        Target_ExperimentAcc_dic[target].add(experiment_accession)
#
#with open('make_ExpAcc_counts.sh', 'w') as out:
#    for exp in ExperimentAcc_FastqFileAcc_dic:
#        out.write(' '.join(['python', PL+'CombineSnpCounts.py', Min_count, SNPS, BNDS, CNVS, output_fp+exp+'_counts.txt',output_fp+exp+'_counts.log']+ExperimentAcc_FastqFileAcc_dic[exp])+'\n')


with open('make_Target_counts.sh', 'w') as out:
    for target in Target_FastqFileAcc_dic:
        out.write(' '.join(['python', PL+'CombineSnpCounts.py', Min_count, SNPS, BNDS, CNVS, output_fp+target+'_counts.txt',output_fp+target+'_counts.log']+Target_FastqFileAcc_dic[target])+'\n')


import glob
file_names = glob.glob('make_Target_counts_s*')
for f in file_names:
    print 'nohup sh '+f+' &> '+f+'.log&'



