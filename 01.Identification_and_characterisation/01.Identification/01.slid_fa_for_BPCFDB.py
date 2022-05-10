import random
from Bio import SeqIO 
import sys
import math
args = sys.argv
fh=open(args[1])
slid_len = int(args[2])
fout = open(args[3],'w')

def subseq(all_len,rec_name,m):
    point = (m-1) * slid_len
    sub_seq = str(rec.seq)[point:point+slid_len]
    outname(rec_name,m)
    fout.write(sub_seq+'\n') 
def outname(rec_name,i):
    name_part = '%04d' % i        
    res_name = ">" + rec_name + "_sub" + name_part + "\n"
    fout.write(res_name)

for rec in SeqIO.parse(fh,'fasta'):
    all_len = len(str(rec.seq))
    n = math.ceil(all_len / slid_len)
    if(n==1):
        outname(rec.name,0)
        fout.write(str(rec.seq)+'\n')
    else:
        for m in range(1,n):            
            subseq(all_len,rec.name,m)
        sub_seq = str(rec.seq)[all_len-slid_len:]
        outname(rec.name,n)
        fout.write(sub_seq+'\n')
fh.close()
fout.close()
