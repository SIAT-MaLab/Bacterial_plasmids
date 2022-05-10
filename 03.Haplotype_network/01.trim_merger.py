#python trim_merger.py overlap_ends.fa merger_overlap_len out_trim.fa 
from Bio import SeqIO
import sys
args = sys.argv
merge_len = {}
with open(args[2]) as merger_list:
    for line in merger_list.readlines():
        merge_len[line.split(' ')[0]] = int(line.split(' ')[1])
with open(args[3],'w') as f_out:
    with open(args[1]) as f_in:
         for record in SeqIO.parse(f_in, "fasta"):
                merge_len2 = merge_len[record.name]
                res_name = '>' + record.name+'\n'
                res_seq = str(record.seq)[merge_len2:]+'\n'
                f_out.write(res_name+res_seq)

