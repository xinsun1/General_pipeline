# Pipeline for aDNA

Things matters for aDNA ***Poor data quality***:
- deamination, C to T/U mutation at the flank region of some reads
- short DNA fragment
- low sequencing depth
- different sequencing library method
- etc. (to add)


This pipeline include the following analysis
1. Mapping
  1.1 Paleomix mapping
  1.2 Quality control  
3. Genotyping
  1. Genotype likelihood
  2. Pseudohaplid calling
4. Population structure analysis
  1. PCA with projection
  2. ADMIXTURE with genotype likelihood
5. Phylogeny
  1. Treemix
6. Gene flow test - admixtools
7. 
Before doing this, remember to perform data QC first.

## 1. Genotype likelihood
test




## 2. pseudohaploid calling

