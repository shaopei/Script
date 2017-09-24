## counts all Asym and Sym in each region
## don't use those one


import scipy.stats as stats

Con_Asym = 'Concordant.Asym.counts.txt'
Con_Sym = 'Concordant.Sym.counts.txt'
Dis_Asym = 'Discordant.Asym.counts.txt'
Dis_Sym = 'Discordant.Sym.counts.txt'
fastq_metadata_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/metadata_fastq.tsv' # for downloaded fastq
fastq_metadata_fp ='/workdir/sc2457/ENCODE/reDoTF/metadata.tsv' # for old 


#compile counts
file_list = [Dis_Asym, Dis_Sym,Con_Asym, Con_Sym]
FastqFileAcc_Asym_count={}
for f in file_list:
    with open(f, 'U') as count_f:
        for l in count_f:
            fastq_file_acc = l.split(':')[0].split('.')[0]
            #print fastq_file_acc
            count = l.strip().split(':')[-1]
            if fastq_file_acc not in FastqFileAcc_Asym_count:
                FastqFileAcc_Asym_count[fastq_file_acc]=[]
            FastqFileAcc_Asym_count[fastq_file_acc].append(int(count))

for f in FastqFileAcc_Asym_count:
    assert len(FastqFileAcc_Asym_count[f])==4, f+' length incorrect:'+ str(len(FastqFileAcc_Asym_count[f]))
    

#map fastq to target
fastq_metadata = [l.strip().split('\t') for l in open(fastq_metadata_fp, 'U').readlines()]
FastqFileAcc_ExperimentAcc_dic = {}
FastqFileAcc_Target_dic = {}
ExperimentAcc_Target_dic = {}
for l in fastq_metadata[1:]:
    experiment_accession = l[3]
    file_accession = l[0]
    #target = l[16].split('-')[0] # for new meta
    target = l[15].split('-')[0]  #for old
    #print target
    FastqFileAcc_ExperimentAcc_dic[file_accession] = experiment_accession
    FastqFileAcc_Target_dic[file_accession] = target
    ExperimentAcc_Target_dic[experiment_accession] = target

# perform Fisher's excat test
with open('Fisher_excat_test.txt' ,'w') as out:
    out.write('\t'.join(str(s) for s in ['FastqFileAcc','Target', 'ExperimentAcc',
                                                 's_oddsratio', 's_pvalue', 'Dis_Asym', 'Dis_Sym','Con_Asym', 'Con_Sym'])+'\n')
    for f in FastqFileAcc_Asym_count:
        Dis_Asym, Dis_Sym,Con_Asym, Con_Sym = FastqFileAcc_Asym_count[f]
        Dis_Sym = Dis_Sym - 1  #grep Sym get 1 from the header
        Con_Sym = Con_Sym - 1  #grep Sym get 1 from the header
        s_oddsratio, s_pvalue = stats.fisher_exact([[Dis_Asym, Dis_Sym],[Con_Asym, Con_Sym]])
        try:
            out.write('\t'.join(str(s) for s in [f,FastqFileAcc_Target_dic[f], FastqFileAcc_ExperimentAcc_dic[f],
                                                 s_oddsratio, s_pvalue, Dis_Asym, Dis_Sym,Con_Asym, Con_Sym])+'\n')
        except KeyError:
            out.write('\t'.join(str(s) for s in [f,ExperimentAcc_Target_dic[f], f, s_oddsratio, s_pvalue, Dis_Asym, Dis_Sym,Con_Asym, Con_Sym])+'\n')



