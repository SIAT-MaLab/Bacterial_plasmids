import numpy as np 
import pandas as pd 
import logging 
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch
file = 'all_refined_plas.dist' 
iniCols = ['genome1','genome2','dist','p','kmers']
uCols = ['genome1','genome2','dist']
Mdb = pd.read_csv(file, names=iniCols, usecols=uCols, sep='\t')

def cluster_hierarchical(db, linkage_cutoff):
    names = list(db.columns)
    arr =  np.asarray(db)
    try:
        arr = ssd.squareform(arr)
    except:
        logging.error("The database passed in is not symmetrical!")
        sys.exit()
    linkage = sch.linkage(arr, method= 'complete')
    fclust = sch.fcluster(linkage,linkage_cutoff, criterion='distance')
    Cdb = _gen_cdb_from_fclust(fclust,names)
    return Cdb

def _gen_cdb_from_fclust(fclust,names):
    Table={'cluster':[],'genome':[]}
    for i, c in enumerate(fclust):
        Table['cluster'].append(c)
        Table['genome'].append(names[i])
    return pd.DataFrame(Table)

linkage_db = Mdb.pivot("genome1","genome2","dist")
Cdb = cluster_hierarchical(linkage_db, linkage_cutoff= 0.01)
Cdb.to_csv('mash_cluster99',index=False)
#Cdb = cluster_hierarchical(linkage_db, linkage_cutoff= 0.001)
#Cdb.to_csv('mash_cluster999',index=False)