import os 
import logging 
import argparse 

def parse_args():
    parser = argparse.ArgumentParser(description="Identification of plasmids from metagenomic/genomic contigs")
    parser.add_argument('-s', required = True, help='SRR_ID')
    parser.add_argument('-c', required = True, help='contig_ID')
    parser.add_argument('-l', required = True, type = int, help='The overlap length of two ends of the contig')
    parser.add_argument('-t', required = True, help='The directionary of trimed fasta file')
    parser.add_argument('-o', required = True, help='Output directionary')
    return parser.parse_args()

def create_subfa(contigID, overlap_len):
    outfile = contigID + "_sub.fa"
    fout = open(outfile,'w')
    content = open("../"+contigID).read().split('\n')
    name = content[0]
    seq = content[1]
    seq_sub1 = seq[:300]
    seq_sub2 = seq[len(seq)-overlap_len-300:len(seq)-overlap_len]
    seq = seq_sub2 + seq_sub1 
    fout.write(name+"\n"+seq+"\n")
    fout.close()
    
def run_cmd(cmd):
    logging.debug(cmd)
    res = os.system(cmd)
    if res != 0:
        logging.error("Error: "+ cmd)
        exit(1)

def seq_blast(contigID, trimD, SRR_ID,overlap_len):
    db_file_pre = trimD + '/' + SRR_ID
    cmd1 = "makeblastdb -dbtype nucl  -in " + db_file_pre + "_1.trim.fa  -input_type fasta -out " + db_file_pre  + "_1.trim.fa_db"
    run_cmd(cmd1)
    cmd2 = "makeblastdb -dbtype nucl  -in " + db_file_pre + "_2.trim.fa  -input_type fasta -out " + db_file_pre  + "_2.trim.fa_db"
    run_cmd(cmd2)
    cmd3 = "blastn  -query " + contigID + "_sub.fa  -db " + db_file_pre + "_1.trim.fa_db -out read1.out -evalue 1e-10 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\""
    run_cmd(cmd3)
    cmd4 = "blastn  -query " + contigID + "_sub.fa  -db " + db_file_pre + "_2.trim.fa_db -out read2.out -evalue 1e-10 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\""
    run_cmd(cmd4)
    
def main():
    args = parse_args()
    SRR_ID = args.s
    contigID = args.c
    overlap_len = args.l
    trimD = os.path.abspath(args.t)
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s: %(message)s", datefmt = '%Y-%m-%d %H:%M:%S %a')
    
    outdir = args.o
    if not os.path.exists(outdir):
        logging.info("Making output dir: "+ outdir)
        os.makedirs(outdir)
    else:
        logging.info("Output dir exists: "+ outdir)
    os.chdir(outdir)
    logging.info("Create subseq of input contig...")
    create_subfa(contigID, overlap_len)
    seq_blast(contigID, trimD, SRR_ID, overlap_len)
    
if __name__ == "__main__":
    main()
