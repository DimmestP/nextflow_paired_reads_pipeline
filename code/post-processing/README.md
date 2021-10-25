# post-processing 

### find_polyA_site.sh
Bash script for extracting all reads mapped to the 3'UTR of the constructs of interest (and the normalising genes PGK1, RPS3, TSA1)

```
# how to run script
# enable script to be executed
chmod +x find_polyA_site.sh

# run script whilst passing the path to the folder with the
# find_polyA_site nextflow pipeline output and the 
# experiment name
./find_polyA_site.sh /homes/wallacelab/datastore/wallace_rna/data/2021/10-Oct/Sam/chimera_quantseq_pipeline_output QuantSeq
```
