# Example: python complete.py SRR9205994 2979009+ 3031990+ 4383 VVJT01000021.1 21.133999
# python complete.py fastg_file begin_node end_node PLS_len contig_ID cov_ref
import pyfastg 
import networkx as nx 
import json
import sys 
args = sys.argv
g = pyfastg.parse_fastg("03.assembly/" + args[1] + "/assembly_graph.fastg")

def completeness(node1, node2, begin, PLS_len=args[4], contig_ID=args[5],SRR_ID=args[1], cov_ref = args[6]):
    path = nx.shortest_path(g, node1, node2, weight='length')
    if(node1 != begin):
        path = [begin] + path
    cov_ref = float(cov_ref)
    len_out, MC_flag = sum_len(path, cov_ref)
    complete = int(PLS_len) / (len_out + int(PLS_len))
    MC_out = "SingleCopy:" + json.dumps(path) if MC_flag==0 else "CopyError:"+json.dumps(path)
    print(SRR_ID, contig_ID, complete, MC_out)

def sum_len(path, cov_ref):
    len_all = 0
    if(len(path)==2):
        MC_flag = 0
    else:
        MC_flag = 0
        for i in path[1:len(path) - 1]:
            covX, MC = cov_M(path, i, cov_ref)
            len_all = len_all + g.nodes[i]['length'] * covX
            if(MC==1):
                MC_flag =1
    return len_all, MC_flag
    
def cov_M(path, node, cov_ref):
    nodes_A = list(g[node])
    MC = 0
    if(set(nodes_A).issubset(set(path))): 
        if(node in nodes_A):
            cov_fold = max(int(g.nodes[node]['cov'] / cov_ref + 0.2), 1)
        else:
            cov_fold = 1
    else:
        cov_left = g.nodes[node]['cov'] - cov_BA(node, path)
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
    
def cov_BA(node, path):
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
    

try: 
    if(args[2]!=args[3]):
        completeness(args[2],args[3],args[2])
    else:
        end_node = list(g[args[2]])
        for node in end_node: 
            completeness(node,args[3],args[2])
except nx.NetworkXNoPath:
    print(args[1],args[5],'- NetworkXNoPath')
#except CopyError as e:
#    print(args[1],args[5],'- UnsolvedCopy:'+str(e))
