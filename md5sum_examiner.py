from sys import argv

assert len(argv) ==3, 'python md5sum_examiner.py metadata.tsv ../ENCODE_chipseq_md5sum'

meta_fp = argv[1]  # 'metadata.tsv'
file_need_to_check = argv[2] #  '../ENCODE_chipseq_md5sum'

file_md5sum_dic={}
with open(meta_fp,'U') as meta_f:
    head = meta_f.readline().strip().split('\t')
    m_index = head.index('md5sum')
    for l in meta_f:
        ll = l.strip().split('\t')
        #print ll[0], ll[m_index]
        file_md5sum_dic[ll[0]] = ll[m_index]

with open(file_need_to_check, 'U') as check_f:
    for l in check_f:
        md5sum, file_path = l.strip().split('  ')
        File_accession = file_path.split('/')[-1].split('.')[0]
        if file_md5sum_dic[File_accession] != md5sum:
            print File_accession, md5sum, file_md5sum_dic[File_accession]
    print 'End'