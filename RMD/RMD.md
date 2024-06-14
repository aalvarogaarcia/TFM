
<a name="7"></a> 
#### 7. Align non-host reads to a reference genome:
<a name="7.1"></a>  
##### 7.1 Select the reference based on BLASTn results:
```Bash

strain=$( cat sample.strain.txt )

# Path to multiple Influenza strains reference genomes separated by segments:
infA_h1n1_dir=/path/to/reference/Influenza_A_virus_H1N1_California/separated_segments
infA_h3n2_dir=/path/to/referenceInfluenza_A_virus_H3N2_Wisconsin/separated_segments
infB_dir=/path/to/referenceInfluenza_B_virus_Yamagata/separated_segments

# Select reference folder:
if [[ ${strain} == "A-H1N1" ]]; then
  reference_dir=${infA_h1n1_dir}
elif [[ ${strain} == "A-H3N2" ]]; then
  reference_dir=${infA_h3n2_dir}
elif [[ ${strain} == "B" ]]; then
  reference_dir=${infB_dir}
else
  echo "LOG: There is no FASTA file for strain ${strain}."
  exit 1
fi
```
