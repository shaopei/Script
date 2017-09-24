import os
import glob

file_names = glob.glob('*.fastq.gz')
# mkdir TF
# mv TF_*.fastq.gz TF/.
# cd TF/
# make -f PIPELINE_TF.mk > PIPELINE_TF.log
# cd ../TF
# make -f PIPELINE_TF.mk > PIPELINE_TF.log

template = "mkdir %s\nmv %s_*.fastq.gz %s/."

run_sh_lines=[]
for f in file_names:
    ff=f[:-11]
    run_sh_lines.append(template %(ff, ff, ff))

out = open('mv.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()



##
import os
import glob
file_names = glob.glob('*')
#template = "rm %s/PIPELINE_TF.mk"
template = "cp PIPELINE_TF_2015OldSV.mk %s/."
run_sh_lines=[]
for f in file_names:
    run_sh_lines.append(template %(f))

out = open('cp_PIPELINE.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()

##
import os
import glob
file_names = glob.glob('*')
#template = "cd /workdir/sc2457/ENCODE/fastq/batch1/%s\npwd\nmake -f PIPELINE_TF_2015OldSV.mk > PIPELINE_%s_2015OldSV.log 2>&1"
template = "cd /workdir/sc2457/ENCODE/fastq/PE/batch1/%s\npwd\nmake -f PIPELINE_TF_2015OldSV.mk > PIPELINE_%s_2015OldSV.log 2>&1"
run_sh_lines=[]
for f in file_names:
    run_sh_lines.append(template %(f, f))

out = open('run_PE_PIPELINE.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()

##
import os
import glob
file_names = glob.glob('*_p1.fastq.gz')
template = "nohup zcat %s_p1.fastq.gz %s_p2.fastq.gz > ../%s.fastq &"
run_sh_lines=[]
for f in file_names:
    ff = f[:-12]
    run_sh_lines.append(template %(ff,ff,ff))

out = open('combine_PE_reads.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()

##
import os
import glob
file_names = glob.glob('*')
template = "mv %s/counts.txt SE_count_batch1/%s_counts.txt"
template2 = "mv %s/counts.log SE_count_batch1/%s_counts.log"
template3 = "mv %s/FDR.txt SE_count_batch1/%s_FDR.txt"
template4 = "mv %s/interestingHets.txt SE_count_batch1/%s_interestingHets.txt"

run_sh_lines=[]
for f in file_names:
    run_sh_lines.append(template %(f, f))
    run_sh_lines.append(template2 %(f, f))
    run_sh_lines.append(template3 %(f, f))
    run_sh_lines.append(template4 %(f, f))

out = open('mv_batch1_result.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()

##
import os
import glob
file_names = glob.glob('*')
template = "mv %s/counts.txt PE_count_batch1/%s_counts.txt"
template2 = "mv %s/counts.log PE_count_batch1/%s_counts.log"
template3 = "mv %s/FDR.txt PE_count_batch1/%s_FDR.txt"
template4 = "mv %s/interestingHets.txt PE_count_batch1/%s_interestingHets.txt"
template5 = "mv %s_counts.txt %s.PE_counts.txt"

run_sh_lines=[]
for f in file_names:
    #run_sh_lines.append(template %(f, f))
    #run_sh_lines.append(template2 %(f, f))
    #run_sh_lines.append(template3 %(f, f))
    ff = f.split('_')[0]
    run_sh_lines.append(template5 %(ff, ff))

out = open('mv_batch1_result.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()

##
import os
import glob

file_names = glob.glob('*.fastq.gz')
file_names.sort()
with open('batch2_list.txt','w') as out:
    out.write('\n'.join(file_names))
# mkdir TF
# mv TF_*.fastq.gz TF/.
# cd TF/
# make -f PIPELINE_TF.mk > PIPELINE_TF.log
# cd ../TF
# make -f PIPELINE_TF.mk > PIPELINE_TF.log

template = "mv %s %s/."
template2 = "cd /workdir/sc2457/ENCODE/fastq/batch1N2/%s\npwd\nmake -f PIPELINE_TF_2015OldSV.mk > PIPELINE_%s_batch2_2015OldSV.log 2>&1"

run_sh_lines=[]
TF = set()
for f in file_names:
    ff=f[:-11]
    TF.add(ff)
    run_sh_lines.append(template %(f, ff))

TF=list(TF)
TF.sort()
for tf in TF:
    run_sh_lines.append(template2 %(tf, tf))

out = open('run_batch2.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()
##
import os
import glob
file_names = glob.glob('*')
template = "cp %s/counts.txt SE_count_batch2/%s_counts.txt"

run_sh_lines=[]
for f in file_names:
    run_sh_lines.append(template %(f, f))

out = open('cp_batch2_result.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()

####

#cd ATF2
#mkdir ATF2_2
#mv ATF2_2.cnt ATF2_2.fastq.gz ATF2_2/.
#cp PIPELINE_TF_2015OldSV.mk ATF2_2/.
#cd ATF2_2/ 
#make -f PIPELINE_TF_2015OldSV.mk > PIPELINE_ATF2_2_2015OldSV.log
template = "cd /workdir/sc2457/ENCODE/fastq/batch1N2/%s\nmkdir %s\nmv %s.cnt %s.fastq.gz %s/.\ncp PIPELINE_TF_2015OldSV.mk %s/.\ncd %s/\nnohup make -f PIPELINE_TF_2015OldSV.mk &> PIPELINE_%s_2015OldSV.log &"

file_list = [l.strip() for l in open('batch2_list.txt','U').readlines()]
run_sh_lines=[]
for f in file_list:
    folder = f.split('_')[0]
    new_file = f.split('.')[0]
    run_sh_lines.append(template%(folder,new_file, new_file, new_file, new_file, new_file, new_file, new_file))

out = open('seperate_batch2_result.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()


#
import os
import glob
file_names = glob.glob('*')
template = "mv %s_counts.txt %s-1_counts.txt"

run_sh_lines=[]
for f in file_names:
    new_file = f.split('_')[0]
    run_sh_lines.append(template %(new_file, new_file))

out = open('mv.sh',"w")
out.write("\n".join(run_sh_lines))
out.close()


