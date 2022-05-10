#This script tends sto gnerate plasmids with uneven length distribution. It is because many NCBI plasmid have length less than 200 kbp. 
#Therefore, we run this script for two times and combine the results: 1. All plasmids/chromsomes are as "fh"; 2. All plasmids with length > 10k as "fh".
#After generating fragments within 2k-200kbp, we add random sampling for even length distribution of plasmid/chromosome fragments. 
import random
from Bio import SeqIO 
import sys
import math
args = sys.argv
fh=open(args[1])
fout = open(args[2],'w')

def subseq(all_len,cut_min_len,rec_name,m):
    sub_len = random.randint(cut_min_len,200000)
    if sub_len > all_len:
        sub_len = random.randint(cut_min_len,all_len)
    point = random.randint(sub_len,all_len)
    sub_seq = str(rec.seq)[point-sub_len:point]
    outname(rec_name,m)
    fout.write(sub_seq+'\n') 
def outname(rec_name,i):
    name_part = '%04d' % i        
    res_name = ">" + rec_name + "_sub" + name_part + "\n"
    #res_name = ">" + rec_name + "_sub2_" + name_part + "\n"
    fout.write(res_name)

for rec in SeqIO.parse(fh,'fasta'):
    all_len = len(str(rec.seq))
    if(all_len<2000):
        outname(rec.name,0)
        fout.write(str(rec.seq)+'\n')
    else:
        n = math.ceil(all_len / 200000)
        for m in range(1,n+1):            
            subseq(all_len,2000,rec.name,m)
        #for m in range(1,5*n+1):
            #subseq(all_len,10000,rec.name,m)
fh.close()
fout.close()
