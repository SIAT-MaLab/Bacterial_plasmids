from Bio import SeqIO
import pandas as pd 
import networkx as nx 
import json 
import os 
import logging 
import argparse 
import pyfastg

def parse_args():
    parser = argparse.ArgumentParser(description="Identification of plasmids from metagenomic/genomic contigs")
    parser.add_argument('-s', required = True, help='SRR_ID')
    parser.add_argument('-c', required = True, help='contig_ID')
    parser.add_argument('-l', required = True, type = int, help='The length of contig')
    parser.add_argument('-cr', required = True, type = float, help='The cov_ref of contig')
    parser.add_argument('-p', required = True, help='The file containing nodes list of path')
    parser.add_argument('-q', required = True, help='The directionary of file containing sequences of begining and ending nodes')
    parser.add_argument('-t', required = True, help='The directionary of trimed fasta file')
    parser.add_argument('-a', required = True, help='The assembly directionary')
    parser.add_argument('-o', required = True, help='Output directionary')
    return parser.parse_args()

def extract_node_name(name_long):
    name1 = name_long[:-1]
    name2 = name1.split(':')[0]
    rc = False
    if name2.endswith("'"):
        rc = True
    name = name2.split('_')[1]
    if rc:
        name += "-"
    else:
        name += "+"
    return name

def extract_path(ID_list):
    name1 = open(ID_list).read().split('\n')[:-1][0]
    name2 = name1.split(':')[1]
    path = json.loads(name2)
    return path
    
    
def extract_node_fa(dic_use, path, outfile, ori_file):
    fout = open(outfile,'w')
    nodes = list(set(path[1:-1]))
    for rec2 in nodes:
        seq = dic_use.get(rec2,'wrong')
        res_name = '>' + rec2+'\n'
        fout.write(res_name)
        fout.write(seq+'\n')
    fout.flush()
    fout.close()
    cmd = "cat " + ori_file + " >> " + outfile
    run_cmd(cmd)

def run_cmd(cmd):
    logging.debug(cmd)
    res = os.system(cmd)
    if res != 0:
        logging.error("Error: "+ cmd)
        exit(0)

def read_blast(node_fa_file, SRR_ID, trimD):
    cmd1 = "makeblastdb  -dbtype nucl  -in " + node_fa_file + "  -input_type fasta -out " + node_fa_file + "_db"
    run_cmd(cmd1)
    file_pre = trimD + '/' + SRR_ID
    cmd2 = "blastn  -query " + file_pre + "_1.trim.fa  -db " + node_fa_file + "_db -out read1.out -evalue 1e-10 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\""
    run_cmd(cmd2)
    cmd3 = "blastn  -query " + file_pre + "_2.trim.fa  -db " + node_fa_file + "_db -out read2.out -evalue 1e-10 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\""
    run_cmd(cmd3)


def pair_info(file1, file2):
    cmd1 = "awk \'$4/$13>0.9 {if($9<$10){print $1,$2,\"+\"}else{print $1,$2,\"-\"}}\' " + file1 + " > tmp1"
    run_cmd(cmd1)
    cmd2 = "awk \'$4/$13>0.9 {if($9<$10){print $1,$2,\"+\"}else{print $1,$2,\"-\"}}\' " + file2 + " > tmp2"
    run_cmd(cmd2)
    
    df1 = pd.read_csv('tmp1', sep=' ', header=None)
    header1 = ['ID', 'read1', 'ori1']
    df1.columns = header1
    df2 = pd.read_csv('tmp2', sep=' ', header=None)
    header2 = ['ID', 'read2', 'ori2']
    df2.columns = header2 
    df_merge = pd.merge(df1,df2)
    
    df_merge['read1_sub'] = df_merge['read1'].str[0:-1]
    df_merge['read2_sub'] = df_merge['read2'].str[0:-1]
    tmp1 = df_merge[df_merge['read1_sub'] != df_merge['read2_sub']]
    tmp2 = tmp1[tmp1['ori1'] != tmp1['ori2']]
    
    res1_df = tmp2[tmp2['ori1'] == "+"]
    res2_df = tmp2[tmp2['ori1'] == "-"]
    return res1_df, res2_df

def path_make(g, G, contigID, path, cov_ref):
    path_new = nx.shortest_path(g, contigID + "_2", contigID + "_1")
    path_new[0] = path[0]
    path_new[-1] = path[-1]
    a = []
    for i in path_new[1:-1]:
        covX, MC = cov_M(G, path, i, cov_ref)
        a.append(MC)
    res = 1 in a 
    return path_new, res
    
class NotPassMC(Exception):
    pass

def pair2link(g, path, contigID, contig_len, fastg_file, cov_ref):
    G = pyfastg.parse_fastg(fastg_file)
    path_new, path_flag = path_make(g, G, contigID, path, cov_ref)
    if path_flag:
        raise NotPassMC
    G_edges = list(G.edges)
    base_len = contig_len
    for i in range(0,len(path_new)-1):
        node1 = path_new[i]
        node2 = path_new[i+1]
        if((node1, node2) not in G_edges):
            if(i==0):
                base_len = base_len + 160 
            else:
                G.nodes[node1]['length'] = G.nodes[node1]['length'] + 160 
    for i in path_new[1:len(path_new)-1]:
        all_len = base_len + G.nodes[i]['length']
    complete = contig_len / all_len 
    return complete 

def cov_M(g, path, node, cov_ref):
    nodes_A = list(g[node])
    MC = 0
    if(set(nodes_A).issubset(set(path))): 
        if(node in nodes_A):
            cov_fold = max(int(g.nodes[node]['cov'] / cov_ref + 0.2), 1)
        else:
            cov_fold = 1
    else:
        cov_left = g.nodes[node]['cov'] - cov_BA(g, node, path)
        cov_fold = max(int(cov_left / cov_ref + 0.2), 1)
        if(cov_fold > 1):
            node_M = node[:len(node)-1]
            node_1 = node_M + "+"
            node_2 = node_M + "-"
            if(node_1 in path and node_2 in path and cov_fold == 2):
                cov_fold = 1
            else:
                cov_fold = 1 
                MC = 1 
    return cov_fold, MC
    
def cov_BA(g, node, path):
    adj_B = list(g.in_edges(node))
    cov_B = 0
    for i in adj_B:
        if i[0] not in path:
            cov_B = cov_B + g.nodes[i[0]]['cov']
    adj_A = list(g[node])
    cov_A = 0
    for i in adj_A:
        if i not in path:
            cov_A = cov_A + g.nodes[i]['cov']
    cov_BA = min(cov_B, cov_A)
    return cov_BA
    
def main():
    args = parse_args()
    base = args.s
    SRR_ID = base
    contigID = args.c
    contig_len = args.l
    fastg_file = os.path.abspath(args.a) + "/" + args.s + "/assembly_graph.fastg"
    trimD = os.path.abspath(args.t)
    nodes_file = os.path.abspath(args.p)
    ori_file = os.path.abspath(args.q) + "/" + contigID + "_sub.fa"
    cov_ref = args.cr
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s: %(message)s", datefmt = '%Y-%m-%d %H:%M:%S %a')
    
    outdir = args.o
    if not os.path.exists(outdir):
        logging.info("Making output dir: "+ outdir)
        os.makedirs(outdir)
    else:
        logging.info("Output dir exists: "+ outdir)
    os.chdir(outdir)
    logging.info('01.Parsing the fastg file...')
    
    fh=open(fastg_file)
    dic_use = {} 
    for rec in SeqIO.parse(fh,'fasta'):
        name = extract_node_name(rec.name)
        dic_use[name]= str(rec.seq)
    fh.close()
    logging.info('Parsing: done')

    logging.info("02.Extracting the sequences of focused nodes in assembly graph...")
    node_fa_file = base + "_nodes.fa"
    path = extract_path(nodes_file)
    extract_node_fa(dic_use, path, node_fa_file, ori_file)
    logging.info("Extracting the sequences: done!")
    
    logging.info("03.Blast the sequences of nodes against the reads...")
    read_blast(node_fa_file, base, trimD)
    logging.info("Blast: done!")
    logging.info("04.Paring the blast result...")
    #path = extract_path(nodes_file)
    df1, df2 = pair_info("read1.out", "read2.out")
    try:
        g = nx.from_pandas_edgelist(df1,'read1','read2',create_using=nx.DiGraph())
        complete = pair2link(g, path, contigID, contig_len, fastg_file, cov_ref)
        logging.info("Done. Completeness: " + SRR_ID + " " + contigID + " " + str(complete))
    except:
        try:
            g = nx.from_pandas_edgelist(df2,'read2','read1',create_using=nx.DiGraph())
            complete = pair2link(g, path, contigID, contig_len, fastg_file, cov_ref)
            logging.info("Done. Completeness(reverse): " + SRR_ID + " " + contigID + " " + str(complete))
        except nx.NetworkXNoPath:
            logging.info('Done. NetworkXNoPath: ' + SRR_ID + " " + contigID )
        except nx.NodeNotFound:
            logging.info('Done. NodeNotFound in new path: ' + SRR_ID + " " + contigID )
        except NotPassMC:
            logging.info('Done. NotPassMC: ' + SRR_ID + " " + contigID )
    
if __name__ == "__main__":
    main()

