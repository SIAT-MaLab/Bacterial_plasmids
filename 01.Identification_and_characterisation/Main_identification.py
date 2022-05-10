# Software requirment: Diamond/prodigal/mash
# Database requirment: NCBI chromosome/plasmid fragments; PGMs and CGMs
from Bio import SeqIO
import os 
import sys 
import logging 
import argparse
import pandas as pd
import subprocess
import random 

def parse_args():
    parser = argparse.ArgumentParser(description="Identification of plasmids from metagenomic/genomic contigs")
    parser.add_argument('-f', required = True, help='Input fasta file')
    parser.add_argument('-o', required = True, help='Output directory')
    parser.add_argument('-t', type = int, default = 20, help='Number of threads, default: 20')
    return parser.parse_args()

def run_mash(args_f, name, mash_db, threads):
    command = 'mash sketch -i -k 17 '+args_f+' -o '+name + ' 2>/dev/null'
    logging.debug(command)
    res = os.system(command)
    if res != 0:
        logging.error("Mash-sketch run failed")
        exit(1)
    command = 'mash dist '+name+'.msh '+mash_db+' -d 0.15 -p ' + str(threads) + ' > '+ name + '_mash_dist'
    logging.debug(command)
    res = os.system(command)
    if res != 0:
        logging.error("Mash-dist run failed")
        exit(1)

def parse_mash_result(args_f, name, mash_ref_type):
    command = "grep '^>' "+args_f+" | awk '{print $1}' > "+ name +"_query_name ; sed -i 's/^>//' " + name +"_query_name"
    logging.debug(command)
    res = os.system(command)
    if res != 0:
        logging.error("Getting query name failed")
        exit(1)
        
    query_name = pd.read_csv(name +"_query_name",sep=' ',header=None)
    query_name.columns = ['contig']
    dist = pd.read_csv(name + '_mash_dist',sep='\t',header=None)
    dist = dist.iloc[:,[0,1,2]]
    dist.columns = ['contig','contig2','dist']
    ref_type = pd.read_csv(mash_ref_type,sep=' ',header=None)
    ref_type.columns = ['contig2','type']
    
    dist = dist.groupby('contig')
    top5 = dist.apply(lambda x: x.nsmallest(5,'dist',keep='all')).reset_index(drop=True)
    top5_refInfo = pd.merge(top5,ref_type,on='contig2',how='left').drop_duplicates()
    
    type_count = top5_refInfo.groupby(['contig','type']).size().reset_index(name='counts')
    res = type_count.pivot(index='contig',columns='type',values='counts').fillna(0)
    res.columns.name = None
    res = res.reset_index()
    res['predict'] = res.apply(lambda x: 'plasmid' if x['plasmid'] > x['chromosome'] else 'chromosome',axis=1)
    res = pd.merge(query_name,res,on='contig',how='left').drop_duplicates()
    res.to_csv(name + '_mash_result',index=False)
    return(res)

def run_prodigal(args_f, name, unknown_seq, threads):
    dic_allSeq = {} 
    for rec in SeqIO.parse(args_f,'fasta'):
        dic_allSeq[rec.name]= str(rec.seq)
    random.shuffle(unknown_seq)
    split_items = [unknown_seq[i::threads] for i in range(threads)]
    i = 0 
    i = 0 
    ps = []
    out_protein_file = [] 
    in_fasta_part = [] 
    log_file = [] 
    for seqList in split_items:
        i = i + 1
        fout = open('unknown_seq_part'+str(i)+'.fa','w')
        for seq_name in seqList:
            seq = dic_allSeq.get(seq_name,'wrong')
            out_name = '>' + str(seq_name)+'\n'
            fout.write(out_name)
            fout.write(seq+'\n')
        fout.close()
        in_fasta_part.append('unknown_seq_part'+str(i)+'.fa')
        command = 'prodigal -p meta -d /dev/null -o /dev/null -s /dev/null -a ' + 'unknown_seq_part'+str(i)+'_protein.faa -i '+ 'unknown_seq_part'+str(i)+'.fa  2> prodigal' + str(i) + '_log'
        log_file.append('prodigal'+str(i)+'_log')
        logging.debug(command)
        p = subprocess.Popen(command,shell=True)
        ps.append(p)
        out_protein_file.append('unknown_seq_part'+str(i)+'_protein.faa')
    for p in ps:
        p.wait()
    command = 'cat ' + (' ').join(out_protein_file) + ' > ' + name + '_protein.faa'
    os.system(command)
    command = 'cat ' + (' ').join(in_fasta_part) + ' > ' + name + '_prodigal_input.fa'
    os.system(command)
    command = 'cat ' + (' ').join(log_file) + ' > ' + name + '_prodigal.log'
    os.system(command)
    command = 'rm ' + (' ').join(in_fasta_part)
    os.system(command)
    command = 'rm ' + (' ').join(out_protein_file)
    os.system(command)
    command = 'rm ' + (' ').join(log_file)
    os.system(command)

def run_diamond(name, pfam_db, threads):
    command = 'diamond blastp --db ' + pfam_db +' --query ' +  name + '_protein.faa --out ' + name + '_diamond_out.tab --outfmt 6 --max-target-seqs 1 --evalue 1e-3 -p' + str(threads) + ' 1> '+name + '_diamond.log'
    res = os.system(command)
    if res != 0:
        logging.error("Diamond run failed")
        exit(1)

def parse_diamond_result(name,pfam_protein,pfam_type):
    input = pd.read_csv(name + '_diamond_out.tab',sep='\t',header=None)
    input = input.iloc[:,0:2]
    input.columns = ['query','gene']
    
    cluster = pd.merge(input,pfam_protein,on='gene',how='left').drop_duplicates()
    cluster['contig'] = [i.rsplit("_", 1)[0] for i in cluster['query']]
    contig_cluster = cluster[['contig','cluster']].drop_duplicates()
    tmp = contig_cluster['cluster'].isin(pfam_type['cluster'])
    cluster_use = contig_cluster[tmp]
    cluster_use = pd.merge(cluster_use,pfam_type,on='cluster',how='left').drop_duplicates()
    
    type_count = cluster_use.groupby(['contig','type']).size().reset_index(name='counts')
    res = type_count.pivot(index='contig',columns='type',values='counts').fillna(0)
    res['predict'] = res.apply(lambda x: 'plasmid' if x['plasmid'] >0 and x['chromosome'] == 0 else 'chromosome',axis=1)
    res['contig'] = res.index
    diamond_res = res[['contig','predict']]
    return(diamond_res)

def main():
    args = parse_args()
    dirname = os.path.split(os.path.abspath(__file__))[0]
    mash_db = dirname + '/reference/repre_slid.fna.msh'
    pfam_db = dirname + '/reference/Pfam-A.dmnd'
    mash_ref_type = dirname + '/reference/BPCFDB_ref_type'
    base = os.path.basename(args.f)
    name_file = os.path.splitext(base)[0]
    threads = args.t
    
    outdir = args.o
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    name = os.path.join(outdir, name_file)
    
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s: %(message)s", datefmt = '%Y-%m-%d %H:%M:%S %a')
    
    logging.info('*' * 15 + 'BEGIN' + '*' * 15)
    logging.info('*' * 15 + 'First stage' + '*' * 15)
    logging.info('Use Mash to identify PLSs against BPCFDB.')
    run_mash(args.f,name,mash_db,threads)
    logging.info('Parse mash result...')
    mash_res = parse_mash_result(args.f,name,mash_ref_type)
    
    logging.info('*' * 15 + 'Second stage' + '*' * 15)
    logging.info('Identify PLSs with PGMs but no CGMs.')
    logging.info('Unknown contigs from mash result are analysed in this stage.\n\
    Plasmids identified here are considered potentially novel.')
    
    # Get unknown contigs from mash result and predict genes
    tmp = mash_res['predict'].isna()
    unknown_mash = list(mash_res[tmp]['contig'])
    run_prodigal(args.f,name,unknown_mash,threads)
    
    # Run diamond to compare genes with pfam fasta database
    run_diamond(name, pfam_db, threads)
    pfam_protein = pd.read_csv(dirname + '/reference/CGM_PGM_protein',sep=' ',header=None)
    pfam_protein.columns = ['gene','cluster']
    pfam_type = pd.read_csv(dirname + '/reference/CGM_PGM_type')
    diamond_res = parse_diamond_result(name,pfam_protein,pfam_type)
    
    logging.info('Combine result from mash and plasmid-like gene analysis')
    mash_predict = mash_res[ ~ tmp][['contig','predict']]
    mash_predict['source'] = 'mash'
    
    diamond_predict = pd.merge(mash_res[tmp]['contig'],diamond_res,on='contig',how='left').drop_duplicates()
    diamond_predict = diamond_predict.fillna('chromosome')
    diamond_predict['source'] = 'Pfam_gene'
    
    result = pd.concat([mash_predict,diamond_predict])
    result.to_csv(name + '_result.csv',index=False)
    logging.info('Result is wrote in ' + name + '_result.csv')
    logging.info('*' * 15 +'END'+ '*' * 15)

if __name__ == "__main__":
    main()
