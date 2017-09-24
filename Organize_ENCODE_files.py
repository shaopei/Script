

metadata = [l.strip().split('\t') for l in open('metadata.tsv', 'U').readlines()]

Experiment_File_dic = {}
File_Experiment_dic = {}
Experiment_acc_target_dic = {}
for l in metadata[1:]:
    Experiment_accession = l[3]
    Experiment_acc_target_dic[Experiment_accession] = l[15].split('-')[0]
    File_Experiment_dic[l[0]] = Experiment_accession
    file_accession = l[0]
    if Experiment_accession not in Experiment_File_dic:
        Experiment_File_dic[Experiment_accession] =[]
    Experiment_File_dic[Experiment_accession].append(file_accession)





import os
import glob
file_names = glob.glob('*.fastq.gz')
#file_names = glob.glob('*.cnt')

run_sh_lines=[]
#mkdir_temp = 'mkdir %s_%s'
#cp_temp = 'cp PIPELINE_TF_2015OldSV.mk %s_%s/.'
#for exp_acc in Experiment_File_dic.keys():
#    run_sh_lines.append(mkdir_temp %(Experiment_acc_target_dic[exp_acc], exp_acc))
#    run_sh_lines.append(cp_temp % (Experiment_acc_target_dic[exp_acc], exp_acc))

mv_temp = 'mv %s %s_%s/.'
for f in file_names:
    run_sh_lines.append(mv_temp % (f, Experiment_acc_target_dic[File_Experiment_dic[f.split('.')[0]]], File_Experiment_dic[f.split('.')[0]] ))

with open('run_1-4.sh', 'w') as out:
    out.write('\n'.join(run_sh_lines))


import os
import glob
temp1 = 'cd /workdir/sc2457/ENCODE/reDoTF/%s\nmake -f PIPELINE_TF_2015OldSV.mk > PIPELINE_TF_2015OldSV.log'
run_sh_lines=[]


#for folder in ['CEBPZ_ENCSR347NOB', 'PAX5_ENCSR000BHD', 'IRF4_ENCSR000BGY/']:
#    run_sh_lines.append(temp1 % (folder))
#
#with open('run_2.sh', 'w') as out:
#    out.write('\n'.join(run_sh_lines))


folder_names = glob.glob('*_ENCSR*')
for folder in folder_names:
    if folder not in ['CEBPZ_ENCSR347NOB', 'PAX5_ENCSR000BHD', 'IRF4_ENCSR000BGY/']:
        run_sh_lines.append(temp1 % (folder))

with open('run_3.sh', 'w') as out:
    out.write('\n'.join(run_sh_lines))
    
    
import os
import glob
folder_names = glob.glob('*_ENCSR*')

temp1 = 'cp %s/counts.txt reDoTF_counts_file/%s.%s_counts.txt'
run_sh_lines=[]
for folder in folder_names:
    a,b = folder.split('_')
    run_sh_lines.append(temp1 % (folder,a,b))

with open('cp_counts.sh', 'w') as out:
    out.write('\n'.join(run_sh_lines))




####
data =[l.strip().split('\t') for l in open('file_acc.txt', 'U').readlines()[1:]]
temp1 = 'cp %s_%s.fastq.gz file_Acc/%s.fastq.gz'
temp2 = 'cp %s_%s.cnt file_Acc/%s.cnt'

run_sh_lines=[]
for l in data:
    file_accession, target, replicate = l[0:3]
    run_sh_lines.append(temp1 %(target,replicate, file_accession))
    run_sh_lines.append(temp2 %(target,replicate, file_accession))

with open('cp.sh', 'w') as out:
    out.write('\n'.join(run_sh_lines))


import os
import glob
folder_names = glob.glob('*.gz')
with open('downloaded_file_names', 'w') as out:
    out.write('\n'.join(folder_names))


NotToInclude =[l.strip() for l in open('NotToDownloadFilename.txt', 'U').readlines()]
with open('files2.txt', 'U') as input_f:
    with open('files3.txt', 'w') as out:
        for l in input_f.readlines():
            file_name = l.strip().split('/')[-1]
            if file_name not in NotToInclude:
                out.write(l)




metadata = [l.strip().split('\t') for l in open('/workdir/sc2457/ENCODE/reDoTF/metadata.tsv', 'U').readlines()]
file_data = [l.strip().split() for l in open('file_size.txt', 'U').readlines()[1:-2]]
File_size_dic = {}
for l in metadata[1:]:
    File_size_dic[l[0]] = l[36]

for l in file_data:
    size = l[4]
    file_accession = l[-1].split('.')[0]
    #print file_accession, size
    try:
        #if File_size_dic[file_accession] == size:
            #print 'correct'
        if File_size_dic[file_accession] != size:
            print file_accession, size, File_size_dic[file_accession]
    except KeyError:
        print 'KeyError: '+file_accession


#####
WaitForDownloadExpAcc = []
metadata = [l.strip().split('\t') for l in open('metadata.tsv', 'U').readlines()]
wait_file_names = [l.strip().split('/')[-1].split('.') for l in open('files3.txt', 'U').readlines()]

Experiment_File_dic = {}
File_Experiment_dic = {}
Experiment_acc_target_dic = {}
for l in metadata[1:]:
    Experiment_accession = l[3]
    Experiment_acc_target_dic[Experiment_accession] = l[15].split('-')[0]
    File_Experiment_dic[l[0]] = Experiment_accession
    file_accession = l[0]
    if Experiment_accession not in Experiment_File_dic:
        Experiment_File_dic[Experiment_accession] =[]
    Experiment_File_dic[Experiment_accession].append(file_accession)


for f in wait_file_names:
    if f[1] == 'fastq':
        WaitForDownloadExpAcc.append(File_Experiment_dic[f[0]])

WaitForDownloadExpAcc= set(WaitForDownloadExpAcc)

import os
import glob
temp1 = 'cd /workdir/sc2457/ENCODE/reDoTF/%s\npwd\npwd >> /workdir/sc2457/ENCODE/reDoTF/whereIam \nmake -f PIPELINE_TF_2015OldSV.mk > PIPELINE_TF_2015OldSV.log'
run_sh_lines=[]
wait_sh_lines=[]

folder_names = glob.glob('*_ENCSR*')
for folder in folder_names:
    exp_acc = folder.split('_')[-1]
    if exp_acc not in WaitForDownloadExpAcc:
        run_sh_lines.append(temp1 % (folder))
    elif exp_acc in WaitForDownloadExpAcc:
        wait_sh_lines.append(temp1 % (folder))
    else:
        print 'ERROR!!'

with open('run_alleleseq.sh', 'w') as out:
    out.write('\n'.join(run_sh_lines))

with open('run_alleleseq_WAIT.sh', 'w') as out:
    out.write('\n'.join(wait_sh_lines))

