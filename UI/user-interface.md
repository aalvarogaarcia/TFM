

<a name="6"></a>
#### 6. Detect hits using BLASTn with NCIB Influenza Virus Database:
<a name="6.1"></a>  
##### 6.1 Download and build NCBI Influenza Virus database:
```Bash
# 1. Download latest NCBI Influenza DB sequences and metadata
wget https://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz

# 2. Gunzip NCBI FLU FASTA
# Replace FASTA headers >gi|{gi}|gb|{accession}|{description} with >{accession}|{description} for easier parsing and processing
zcat influenza.fna.gz | sed -E 's/^>gi\|[0-9]+\|gb\|(\w+)\|/>/' > influenza.fna

# 3. Make BLASTDB:
# In case you installed BLAST through conda:
conda activate blast

makeblastdb -in influenza.fna

mkdir blast_db

mv influenza.fna.* blast_db/

conda deactivate
```


