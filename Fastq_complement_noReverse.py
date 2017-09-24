import sys
import gzip

in_fastq = sys.argv[1]
out_fastq = sys.argv[2]

def Comp(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	return "".join([seq_dict[base] for base in seq])

with open(in_fastq) as f:
    with gzip.open(out_fastq ,"wb") as out:
        NR = 1
        for l in f:
            if NR%4 == 2:
                out.write(Comp(l.strip()))
                out.write('\n')
            else:
                out.write(l)
            NR +=1
        
    
