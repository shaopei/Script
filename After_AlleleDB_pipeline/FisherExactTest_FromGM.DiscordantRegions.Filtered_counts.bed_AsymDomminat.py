# need to be under the folder like 6_CNV.Peak.BlackList.Filtered_counts_bed_in_GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1
fastq_metadata_fp = '/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/metadata_fastq.tsv' # for downloaded fastq
#fastq_metadata_fp ='/workdir/sc2457/ENCODE/reDoTF/metadata.tsv' # for old 


import glob
import scipy.stats as stats

file_names = glob.glob('*cordantRegions.Filtered_counts.bed')

#for f in file_names:
#    with open(f, 'U') as f_in:
#        with open(f[0:-4]+'_Asym.dominant.combine.txt', 'w') as f_out:
#            SymCals = []
#            Region = []
#            for l in f_in.readlines()[1:]:
#                ll = l.strip().split('\t')
#                if ll[6:9] != Region and len(Region)>0:
#                    if SymCals.count('Asym') >0:
#                        print Region, 'Asym'
#                        f_out.write('\t'.join(Region+['Asym'])+'\n')
#                    elif SymCals.count('Sym') >0:
#                        print Region, 'Sym'
#                        f_out.write('\t'.join(Region+['Sym'])+'\n')
#                    SymCals = []
#                SymCals.append(ll[4])
#                Region = ll[6:9]
#            #last line
#            if SymCals.count('Asym') >0:
#                print Region, 'Asym'
#                f_out.write('\t'.join(Region+['Asym'])+'\n')
#            elif SymCals.count('Sym') >0:
#                print Region, 'Sym'
#                f_out.write('\t'.join(Region+['Sym'])+'\n')

#examine each Concordant or Docordant region, if there is one Asym --> Asym, if only Sym --> Sym
FastqFileAcc_Symcals_combined={}
for f in file_names:
    fastq_acc = f.split('.')[0]
    Con_Dis = f.split('.')[-3]
    if fastq_acc not in FastqFileAcc_Symcals_combined:
        FastqFileAcc_Symcals_combined[fastq_acc]={Con_Dis:[]}
    else:
        FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis]=[]
    with open(f, 'U') as f_in:
        SymCals = []
        Region = []
        for l in f_in.readlines()[1:]:
            ll = l.strip().split('\t')
            if ll[6:9] != Region and len(Region)>0:
                if SymCals.count('Asym') >0:
                    print Region, 'Asym'
                    FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis].append('Asym')
                elif SymCals.count('Sym') >0:
                    print Region, 'Sym'
                    FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis].append('Sym')
                SymCals = []
            SymCals.append(ll[4])
            Region = ll[6:9]
        #last line
        if SymCals.count('Asym') >0:
            print Region, 'Asym'
            FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis].append('Asym')
        elif SymCals.count('Sym') >0:
            print Region, 'Sym'
            FastqFileAcc_Symcals_combined[fastq_acc][Con_Dis].append('Sym')

#compile Asym and Sym counts
FastqFileAcc_Asym_count={}
for f in FastqFileAcc_Symcals_combined:
    Dis_Asym = FastqFileAcc_Symcals_combined[f]['DiscordantRegions'].count('Asym')
    Dis_Sym = FastqFileAcc_Symcals_combined[f]['DiscordantRegions'].count('Sym')
    Con_Asym = FastqFileAcc_Symcals_combined[f]['ConcordantRegions'].count('Asym')
    Con_Sym  = FastqFileAcc_Symcals_combined[f]['ConcordantRegions'].count('Sym')
    FastqFileAcc_Asym_count[f] =  [Dis_Asym, Dis_Sym,Con_Asym, Con_Sym]


#map fastq to target
fastq_metadata = [l.strip().split('\t') for l in open(fastq_metadata_fp, 'U').readlines()]
FastqFileAcc_ExperimentAcc_dic = {}
FastqFileAcc_Target_dic = {}
ExperimentAcc_Target_dic = {}
for l in fastq_metadata[1:]:
    experiment_accession = l[3]
    file_accession = l[0]
    target = l[16].split('-')[0] # for new meta
    #target = l[15].split('-')[0]  #for old
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
        s_oddsratio, s_pvalue = stats.fisher_exact([[Dis_Asym, Dis_Sym],[Con_Asym, Con_Sym]])
        try:
            out.write('\t'.join(str(s) for s in [f,FastqFileAcc_Target_dic[f], FastqFileAcc_ExperimentAcc_dic[f],
                                                 s_oddsratio, s_pvalue, Dis_Asym, Dis_Sym,Con_Asym, Con_Sym])+'\n')
        except KeyError:
            out.write('\t'.join(str(s) for s in [f,ExperimentAcc_Target_dic[f], f, s_oddsratio, s_pvalue, Dis_Asym, Dis_Sym,Con_Asym, Con_Sym])+'\n')




