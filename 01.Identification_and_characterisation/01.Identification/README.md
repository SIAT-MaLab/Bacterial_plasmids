## Identification
### 1. Constructing BPCFDB 
+ Sliding NCBI plasmids/chromosomes into 200 kbp non-overlapping fragments -- 01.slid_fa_for_BPCFDB.py 
+ Clustering fragments into clusters by MCL and choosing representatives to make BPCFDB
> mash skecth -i slid.fa  
> mash dist slid.fa.msh slid.fa.msh > slid\_fa.dist  
> cut -f 1,2,3 slid\_fa.dist  
> mcl slid_fa.dist > slid\_fa.mcl
### 2. Choosing CGMs and PGMs
+ Randomly splitting NCBI plasmids/chromosomes into 2 kbp - 200 kbp fragments -- 02.rand_2To200k.py
+ Computing Pfam gene frequency on plasmid/chromosome fragments by Diamond -- 02.pfam_ratio.xlsx
+ Choosing CGMs and PGMs according to gene frequency 
### 3. Testing the workflow 
+ Mimic contigs -- 02.rand\_2To200k.py & sel_contig.R --> 02.mimic_contig.fa.7z


 


