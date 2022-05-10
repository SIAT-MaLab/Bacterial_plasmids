## Haplotype network
#### Selecting clstr599 and clstr417 as examples to construct haplotype network.  
+ Selecting member contigs with circular form 
+ Removing overlapping region at two ends of each circular contig -- 01.trim\_merger.py
+ Arranging contigs with same beginning point and same direction to avoid the influence of circular permutated sequences on sequence alignment -- 02.rearrannge.py
+ Aligning sequences with MUSCLE 
+ Removing member contigs with >1 bp gap in the alignment comparnig with the member contigs with median length
> clstr417/isolates\_4139align.fa  
> clstr417/metagenomes\_4139align.fa  
> clstr599/isolates\_align.fa  
> clstr599/metagenomes\_align.fa  
> Note: There was a reverse complement mutation (10bp) in the intergenic region between rep gene and chm2bp gene. This mutation was modified manually as single base mutation in the alignment file.  

+ Constructing haplotype network 
> clstr417/isolates.nex  
> clstr417/metagenomes.nex  
> clstr599/isolates.nex  
> clstr599/metagenomes.nex  
> Note: Sigle base gap in the *nex file were modified as SNP manually so that the information of those gap were used to construct haplotype network.   
