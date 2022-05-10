# Usage: python rearrange.py clstr_417.fa out.fa 
# This script can't rearrange the sequences in occasions that the searching region spreads on two ends of one sequence or have mutations. 
# We manually check and rearrange the sequences that are overlooked by this script. 
from Bio import SeqIO
import sys
args = sys.argv
merge_len = {}
with open(args[2],'w') as f_out:
    with open(args[1]) as f_in:
         for record in SeqIO.parse(f_in, "fasta"):
            #clstr_417
            #a = record.seq.find('ATCTGTTTGCTTTTAAATTCA')
            #clstr_599
            a = record.seq.find('ACTAGAAAGTTTCGTGGTAAATGCTTTTTGACTTGTTATATGTAGGTTTT')
            if a != -1:
                seq = record.seq
            if a ==-1:
                seq = record.seq.reverse_complement()
                #a = record.seq.find('ATCTGTTTGCTTTTAAATTCA')
                a = seq.find('ACTAGAAAGTTTCGTGGTAAATGCTTTTTGACTTGTTATATGTAGGTTTT')
            if a != -1:
                res_name = '>' + record.name+'\n'
                res_seq = str(seq)[a:]+str(seq)[:a]+'\n'
                f_out.write(res_name+res_seq)
                
