# Selection analysis summary, annotation and visualization with R

## 0. related resources online

Some good tuturial I found available online:
- ~~[theta fst annotation and visualization](https://jyanglab.com/AGRO-932/chapters/a1.2-lab/w5lab.html#1)~~
- [Great bioconductor tutorial](https://uclouvain-cbio.github.io/BSS2019/03_CoreApproachesInBioconductor.html)
- [Bioconductor annotation packages](https://bioconductor.org/packages/3.16/data/annotation/)
- [Bioconductor annotation pipelines](https://www.bioconductor.org/help/course-materials/2015/UseBioconductorFeb2015/A01.5_Annotation.html)
- [Tutorial for GO in R](https://yulab-smu.top/biomedical-knowledge-mining-book/index.html)


## 1. batch data processing

R functions to process PBS result
``` R
```
Gene annotation
``` R
library(tidyverse)
library(GenomicRanges)
library(ggrepel)

BiocManager::install("AnnotationHub", lib="/projects/mjolnir1/people/gnr216/a-software/R-lib/4.2")
BiocManager::install("TxDb.Cfamiliaris.UCSC.canFam3.refGene", lib="/projects/mjolnir1/people/gnr216/a-software/R-lib/4.2")
BiocManager::install("org.Cf.eg.db", lib="/projects/mjolnir1/people/gnr216/a-software/R-lib/4.2", force=TRUE)

.libPaths("/projects/mjolnir1/people/gnr216/a-software/R-lib/4.2")



```

## 2. visualization


## 3. summary
