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
#### run per chr
cd $WDIR/1.qc.err
for i in $(less $CHR_LIST)
do
        angsd -doAncError 1 -anc $ANC -ref $REF -out err_chr.$i -bam $BAM_LIST -minMapQ 30 -minQ 20 -remove_bads 1 -uniqueOnly 1 -checkBamHeaders 0 -r $i &
done
wait
```

Merge results for illustration
`estError.R` can be downloaded [here](https://github.com/ANGSD/angsd/tree/master/R)
``` bash
for i in {1..29}; do echo "err_chr.${i}:.ancErrorChr"; done > list.filename
awk 'FNR!=1 {for(i=1;i<=NF;i++){a[FNR-1][i]+=$i};tf=NF;tl=FNR-1;next}END{for(i=1;i<=tl;i++){for(j=1;j<=tf;j++){printf a[i][j]"\t";};print ""}}' $(< list.filename) > err_all.ancError.all

# estError.R from ANGSD install directory
Rscript ../../program/angsd-master/R/estError.R file=err_all.ancError.all
```

Illustration with R/python
Check the ANGSD manual

Example code for R
``` R
#### error rate ####
# 'plot angsd err estimation result'

## 0. set env
setwd('~/Documents/Projects/Bos/1.qc/')

library(tidyverse)

raw = read.table(file='./errorEst.txt', stringsAsFactors = F, header = F, fill = T )

ind = read.table('list.round1')

meta_sample=read_csv(file="../0.meta/Bos_dataset - meta_information.csv",
                     skip_empty_rows = TRUE)
meta_sample=meta_sample[! is.na(meta_sample$ID),] # remove empty rows
rownames(meta_sample) = meta_sample$ID


n=85


mut = raw[2:(n+1), 2:ncol(raw)]
mut = as.data.frame(sapply(mut, as.numeric))
mut = cbind(mut, ind)
colnames(mut) = c(raw[1,1: (ncol(raw) -1)], 'id')

p = mut %>% pivot_longer(!id, names_to = "mut_type", values_to = "mut_rate") 
# %>% subset(! mut_type  %in%  c("C -> T", "G -> A", "A -> G", "T -> C")) # transversion only

overall = data.frame(a=raw[(n+2):nrow(raw),1], id=ind[,1]) %>% 
    separate(a, c(NA, "mut_rate", NA), sep = " ", remove = T) %>%
    mutate(mut_rate=as.numeric(mut_rate))

p_theme = theme(panel.background = element_rect(fill = NA, colour = "black", linetype = 1),
                panel.grid.major = element_line(color = NA),
                panel.grid.minor = element_line(color = NA),
                axis.title.y = element_blank(), axis.title.x = element_blank(),
                axis.text.x = element_text(size=14,angle = 90, hjust=1, vjust = 0.2),
                legend.position = "none"
)


p %>%
    # add pop, dp group
    mutate(pop = meta_sample[id,]$data_group_bam_folder,
           dp = meta_sample[id,]$`DP_chr1-29`) %>%
    
    mutate(dp1=case_when(dp>=1 ~ ">= 1x",
                         TRUE ~ "<1x")) %>%
    
    ggplot() + geom_col(aes(x=id, y=mut_rate, fill=dp1)) +
    facet_grid(mut_type ~ pop , scales = "free_x", space = "free_x") +
    scale_fill_manual(values=rep(c('#e41a1c','#377eb8'),100)) + 
    theme(panel.background = element_rect(fill = NA, colour = "black", linetype = 1),
          panel.grid.major = element_line(color = NA),
          panel.grid.minor = element_line(color = NA),
          axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=10,angle = 90, hjust=1, vjust = 0.2),
          legend.position = "bottom"
    )
ggsave('err_mut_fixed_16x9.png',width = 16, height = 9, device = "png", dpi = 500)


p %>%
    # add pop, dp group
    mutate(pop = meta_sample[id,]$data_group_bam_folder,
           dp = meta_sample[id,]$`DP_chr1-29`) %>%
    
    mutate(dp1=case_when(dp>=1 ~ ">= 1x",
                         TRUE ~ "<1x")) %>%
    
    ggplot() +
    geom_col(aes(x=id, y=mut_rate, fill=mut_type),position=position_dodge()) +
    facet_grid(dp1 ~ pop , scales = "free_x", space = "free_x") +
    theme(panel.background = element_rect(fill = NA, colour = "black", linetype = 1),
          panel.grid.major = element_line(color = NA),
          panel.grid.minor = element_line(color = NA),
          axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=10, angle = 90, hjust=1, vjust = 0.2),
          legend.position = "right")
ggsave('err_all_free_16x9.png',width = 16, height = 9, device = "png", dpi = 500)


ggplot(overall) + geom_col(aes(x=id, y=mut_rate), fill="blue") + 
    theme(panel.background = element_rect(fill = NA, colour = "black", linetype = 1),
          panel.grid.major = element_line(color = "grey", linetype=2),
          panel.grid.minor = element_line(color = NA),
          axis.title.y = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=13, angle = 90, hjust=1, vjust = 0.2),
          legend.position = "right")

ggsave('err_all_free_16x16.png',width = 16, height = 16, device = "png", dpi = 300)
ggsave('err_overall_withS_16x9.png',width = 16, height = 9, device = "png", dpi = 300)

```

# 3. relatedness with ngsRelate
Sometimes, it is useful/necessary to remove related individuals in the dataset.
Imaging one of the population is just a family. 
It is essentially about the **accurate estimation of allele frequency in the population sampled**.

An example code
``` bash
#!/bin/bash

# set WDIR
WDIR=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/d-Bos

MAF_FILE=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/d-Bos/2.gl/gl_tv_mafF_misF.mafs.gz
B_FILE=/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/d-Bos/2.gl/gl_tv_mafF_misF.beagle.gz


# load module
module load tools
module load htslib/1.8


# prep freq
cd $WDIR/1.qc.rel
# if -ref cut f6
zcat $MAF_FILE | cut -f6 |sed 1d > freq

/home/projects/ku-cbd/data/norwegian_wolves/SNP_calling/program/ngsRelate/ngsRelate -f freq -G $B_FILE -O rel_maf05_misF -n 85 -p 40
```




