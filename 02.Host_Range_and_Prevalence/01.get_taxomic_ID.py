from Bio import Entrez
from pandas import DataFrame,Series 
import pandas as pd 
import json
Entrez.email = '***********'

f1 = open('genome_assembly_ID')
f2=open('genome_metadata','a')
for line in f1.readlines():
    term = line.strip()
    handle = Entrez.esearch(db="assembly",term=term)
    record = Entrez.read(handle)
    for id in record['IdList']:
        handle2 = Entrez.esummary(db="assembly", id=id)
        record2 = Entrez.read(handle2,validate=False)
        res = json.dumps(record2['DocumentSummarySet']['DocumentSummary'][0])
        f2.write(res+"\n")
f1.close()
f2.close()

keys = ['Synonym','Taxid']
final_res=[]
with open('genome_metadata') as f4:
    for line in f4.readlines():
        line = line.strip()
        use = json.loads(line)
        use2 = [use.get(key) for key in keys]
        res_line = [use2[0].get('Genbank'),use2[1]]
        final_res.append(res_line)
final_res_df = DataFrame(final_res)
final_res_df.columns=['genome','Taxonomic_ID']
final_res_df.to_csv("genome_taxonomic_ID",sep="\t",index=False)