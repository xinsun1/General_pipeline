# Quality control process

QC for aDNA including the following process
- check mapdamage result
- check error rate with `ANGSD`
- check relatedness (not necessary but recommended)
- be aware of the genotype missing rate/sequencing depth for each sample
- etc

## 1. check mapdamage result
The mapdamage plot is given by `paleomix`
Check these two files `Fragmisincorporation_plot.pdf` and `Length_plot.pdf` to get an impression of you data quality.

## 2. error rate estimation with ANGSD
We use a "error free" sample to compare the relative error rate for each sample with [ANGSD](http://www.popgen.dk/angsd/index.php/Error_estimation)

An example code
``` bash
#!/bin/bash

WDIR=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/d-Bos
REF=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/d-Bos/0.ref/ARS-UCD1.2/ARS-UCD1.2.fasta
BAM_LIST=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/d-Bos/1.qc.sample_meta/list.round1.bam

CHR_LIST=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/d-Bos/list_chr.angsd.ARD-UCD1.2.chr1_29

# load angsd before use
module purge
module load tools
module load htslib/1.8 angsd/0.931


# gen anc
cd $WDIR/0.angsd_fa
BAM=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/d-Bos/0.transfer/Others/buffalo01.sort.dedup_realign_q25.bam
ID=buffalo01

# angsd -nThreads 4 -i $BAM -dofasta 2 -doCounts 1 -minQ 30 -minmapq 20 \
#                   -setminDepthInd 3 -remove_bads 1 -uniqueOnly 1 \
#                   -ref $REF -out $ID.minD3.q20.f2.cons

ANC=$WDIR/0.angsd_fa/$ID.minD3.q20.f2.cons.fa.gz

module load samtools/1.14
samtools faidx $ANC

# err estimate
##### run per chr
cd $WDIR/1.qc.err
for i in $(less $CHR_LIST)
do
        angsd -doAncError 1 -anc $ANC -ref $REF -out err_chr.$i -bam $BAM_LIST -minMapQ 30 -minQ 20 -remove_bads 1 -uniqueOnly 1 -checkBamHeaders 0 -r $i &
done
wait
```



