# Pipeline for aDNA - genotyope likelihood calling and pseudohaploid calling

**Before doing this, remember to perform data QC first.**

Things matters for aDNA ***Poor data quality***:
- deamination, C to T/U mutation
- low sequencing depth
- short DNA fragment
- others (to add)

Reasons why GL or pseudohaploid
- both will reduce the 



## 1. Genotype likelihood
test




## 2. pseudohaploid calling
Two types of pseudohaploid 
- random base
- consensus base

### 2.1 haploid calling
An example code to call consensus base using ANGSD
``` bash

```

### 2.2 data filtering
Things consider to filter:
- only keep SNPs in high quality(modern/high depth) samples
- MAF
- missingness




