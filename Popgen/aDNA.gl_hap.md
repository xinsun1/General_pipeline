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


``` bash

# choose SNP ID for filter
awk '{print $2}' rand_n87_raw.b.tv.se.mac2.mis20.bim > rand_n87_raw.b.tv.se.mac2.mis20.snpID
# extract selected SNPs
plink2 --bfile rand_n87_raw.b --chr-set 29 --extract rand_n87_raw.b.tv.se.mac2.mis20.snpID --make-bed --out rand_n87_raw.b.tv.se.mac2.mis20.all

plink2 --bfile rand_n87_raw.b.tv.se.mac2.mis20.all --chr-set 29 --geno 0.8 --maf 0.05 --make-bed --indep-pairwise 10kb 1 0.5 --out rand_n87_raw.b.tv.se.mac2.mis20.all.maf05.mis80.ld

plink2 --bfile rand_n87_raw.b.tv.se.mac2.mis20.all.maf05.mis80.ld --chr-set 29 --extract rand_n87_raw.b.tv.se.mac2.mis20.all.maf05.mis80.ld.prune.in --make-bed --out rand_n87_raw.b.tv.se.mac2.mis20.all.maf05.mis80.ld.pruned

# MORE

```


Add extra sample (outgrup)
`run.add_outgroup.sh`
```bash
#!/bin/bash

# BATCH="rand_q30.plink2.tv.maf05.mis25.ld_pruned"
# OUT="L._culpaeus"
# OUT_FA=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/11-angsd/2-Dog/1-hap_fa/AndeanFox/AndeanFox.minD1.q20.cons.fa.gz

BATCH=$1
OUT=$2
OUT_FA=$3


#### 0. set env

#### 1.generate loci list
awk '{split($1,a,":"); print a[1]":"a[2]"-"a[2]}' ${BATCH}.snp > ${BATCH}.list_snp.samtools

#### 2. get gt
samtools faidx -r ${BATCH}.list_snp.samtools $OUT_FA | grep -v ">"  > ${BATCH}.${OUT}.gt

#### 3. generate eig gt
paste ${BATCH}.snp ${BATCH}.${OUT}.gt | awk '{if($7=="N"){$8=9}else{if($7==$5){$8=2}else{if($7==$6){$8=0}else{$8=9}}}; print $8}' > ${BATCH}.${OUT}.gt.012

#### 4. paste
paste -d "" ${BATCH}.geno ${BATCH}.${OUT}.gt.012 > ${BATCH}.${OUT}.geno

ln -s ${BATCH}.snp ${BATCH}.${OUT}.snp
cp ${BATCH}.indi ${BATCH}.${OUT}.indi
echo "${OUT}" | awk '{print $1,"U",$1}' OFS='\t' >> ${BATCH}.${OUT}.indi
```

