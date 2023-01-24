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

``` r
# 1. load modules
library(tidyverse)
library(ggrepel)

# 2. load run batch
library(googlesheets4)
gs4_deauth()

meta_pbs = read_sheet(
    "https://docs.google.com/spreadsheets/d/19HddRvbTJVOEIEGF-SPiHdLgtlVRwOSJ8CuZWFTg3Ms/edit#gid=445479383",
    sheet="PBS")

meta_pbs_g = meta_pbs %>%
    dplyr::select(c("G1", "G2", "G3")) %>%
    filter(!is.na(G1)) # remove blank line

meta_pbs_pop = meta_pbs %>%
    dplyr::select(c("G_PBS", "pop", "sub",	"Pop.sub",	"F")) %>%
    filter(!is.na(G_PBS)) # remove blank line

# tidyverse didn't work, use loop
meta_pbs_g_run = c()
for(i in 1:nrow(meta_pbs_g)){
    tmp = meta_pbs_g %>%
        filter(row_number() == i) %>%
        mutate(
            rG1=list(meta_pbs_pop %>% filter(G_PBS == G1) %>% pull(Pop.sub)),
            rG2=meta_pbs_pop %>% filter(G_PBS == G2, `F`=="F") %>% pull(Pop.sub),
            rG3=meta_pbs_pop %>% filter(G_PBS == G3, `F`=="F") %>% pull(Pop.sub),
            ) %>%
        unnest(rG1)
    meta_pbs_g_run = rbind(meta_pbs_g_run, tmp)
}

# 3. load functions
run_pbs = function(p1, p2, p3, wdir){
    # p1(target), p2, p3:   name of the pop.sub
    # wdir:     subdir for fst results
    
    # 1. read fst results
    fst_dir = wdir
    list_fst_file = list.files(path = fst_dir, pattern = "windowed.weir.fst")
    # double \\ to escape
    fst_p1p2 = read_tsv(
        paste0(fst_dir,
               list_fst_file[
                   c(grep(paste0(p1,"__",p2, "\\."),list_fst_file),
                     grep(paste0(p2,"__",p1, "\\."),list_fst_file)
                     )]),
        col_names = T)
    fst_p1p3 = read_tsv(
        paste0(fst_dir,
               list_fst_file[
                   c(grep(paste0(p1,"__",p3, "\\."),list_fst_file),
                     grep(paste0(p3,"__",p1, "\\."),list_fst_file)
                   )]),
        col_names = T)
    fst_p2p3 = read_tsv(
        paste0(fst_dir,
               list_fst_file[
                   c(grep(paste0(p2,"__",p3, "\\."),list_fst_file),
                     grep(paste0(p3,"__",p2, "\\."),list_fst_file)
                   )]),
        col_names = T)
    
    ## 2. get pbs
    fst_p2p3 = fst_p2p3 %>% unite("name", CHROM:BIN_END, remove = FALSE) %>% 
        mutate(fst_new=ifelse(MEAN_FST < 0, 0, MEAN_FST)) %>% mutate(fst_mod=-log(1-fst_new)) %>%
        dplyr::select(-c("N_VARIANTS","WEIGHTED_FST","MEAN_FST","fst_new")) %>% tibble::column_to_rownames('name')
    
    fst_p1p3 = fst_p1p3 %>% unite("name", CHROM:BIN_END, remove = FALSE) %>% 
        mutate(fst_new=ifelse(MEAN_FST < 0, 0, MEAN_FST)) %>% mutate(fst_mod=-log(1-fst_new)) %>%
        dplyr::select(-c("N_VARIANTS","WEIGHTED_FST","MEAN_FST","fst_new")) %>% tibble::column_to_rownames('name')
    
    fst_p1p2 = fst_p1p2 %>% unite("name", CHROM:BIN_END, remove = FALSE) %>% 
        mutate(fst_new=ifelse(MEAN_FST < 0, 0, MEAN_FST)) %>% mutate(fst_mod=-log(1-fst_new)) %>%
        dplyr::select(-c("N_VARIANTS","WEIGHTED_FST","MEAN_FST","fst_new")) %>% tibble::column_to_rownames('name')
    # merge results
    shared_names = rownames(fst_p2p3)[rownames(fst_p2p3) %in% rownames(fst_p1p3) & rownames(fst_p2p3) %in% rownames(fst_p1p2)]
    merged = cbind(fst_p1p2[shared_names,], fst_p1p3[shared_names,4], fst_p2p3[shared_names,4])
    colnames(merged) = c("CHROM","BIN_START","BIN_END","fst_12","fst_13","fst_23")
    merged = merged %>% 
        mutate(pbs=(fst_13 + fst_12 - fst_23)/2) %>%
        subset(! CHROM == 'chrX') %>%
        mutate(CHR=as.numeric(substr(CHROM,4,5))) %>%
        mutate(POS=BIN_END/1000)
    
    
    don = merged %>% 
        
        # Compute chromosome size
        group_by(CHR) %>% 
        summarise(chr_len=max(POS)) %>% 
        
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        dplyr::select(-chr_len) %>%
        
        # Add this info to the initial dataset
        left_join(merged, ., by=c("CHR"="CHR")) %>%
        
        # Add a cumulative position of each SNP
        arrange(CHR, POS) %>%
        mutate(BPcum=POS+tot)
    return(don)
}

run_pbs_batch = function(meta_pbs_g_run, fst_dir, pdir){
    meta_pbs_g_run_g1 = meta_pbs_g_run %>% dplyr::select(G1) %>% distinct() %>% pull
    for(i in 1:length(meta_pbs_g_run_g1)){
        # 1. run batch per G1
        tmp_run = meta_pbs_g_run %>%
            filter(G1 == meta_pbs_g_run_g1[i]) %>%
            unite("batch", c(rG1,rG2,rG3), sep="-", remove = FALSE)
        tmp_don_s = c()
        for(n in 1:nrow(tmp_run)){
            # run pbs per batch
            print(name_batch)
            tmp_pbs = run_pbs(
                p1=tmp_run[n,]$rG1,
                p2=tmp_run[n,]$rG2,
                p3=tmp_run[n,]$rG3,
                wdir=fst_dir)
            # add batch information
            name_batch = tmp_run[n,]$batch
            tmp_pbs = tmp_pbs %>% 
                mutate(batch=name_batch)
            tmp_don_s = rbind(tmp_don_s, tmp_pbs)
        }
        # 2. plot per G1
        # 2.1 get plot center for each chr
        axisdf = tmp_don_s %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
        # 2.2 precompute quantile for each plot
        cutoff=0.9995 # set cutoff
        tmp_quantile = tmp_don_s %>%
            group_by(batch) %>%
            summarize(threshold = quantile(pbs, probs = cutoff)) 
        # 2.2 plot
        
        ggplot(tmp_don_s, aes(x=BPcum, y=pbs)) +
            
            # Show all points
            geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
            scale_color_manual(values = rep(c('#e41a1c', '#377eb8'), 22 )) +
            
            # custom X axis:
            scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, expand = expansion(mult = 0.01)) +
            # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
            
            geom_hline(data= tmp_quantile, linetype=2,
                       aes(yintercept=threshold)) +
            
            # one batch per line
            facet_grid(batch~., scales = "free") +
            
            # Custom the theme:
            theme(panel.background = element_rect(fill = NA, colour = "black", linetype = 1),
                  panel.grid.major = element_line(color = NA),
                  panel.grid.minor = element_line(color = NA),
                  legend.position = "none",
                  axis.title = element_text(size=20, face = "bold"),
                  axis.text = element_text(size = 16)) +
            xlab("Chromosome") + ylab("PBS")
        p_width=12
        p_height=min(3*(nrow(tmp_run) +1), 30) # max height: 30
        p_name= paste0(pdir,"pbs.", meta_pbs_g_run_g1[i], ".", p_width,"x", p_height, ".png")
        ggsave(p_name, width = p_width, height = p_height, device = "png", dpi = 500)
    }
}
run_pbs_anno_batch = function(meta_pbs_g_run, fst_dir, pdir, tx_db, org_db,
                              cutoff=0.9995, flank=FALSE, gff3=FALSE){
    # run pbs anno batch function:   run pbs with gene annotation for all combinations given in wdir
    # be careful with conflicts bewtween tidyverse and others
    # input:
    #   meta_pbs_g_run:     list of metainfor for run
    #   fst_dir:            fst result dir, e.g. './fst/'
    #   pdir:               pbs plot dir, e.g. './fst_plot/'
    #   tx_db:              tx_db loaded
    #   org_db:             org_db loaded
    #   flank:              if flank db is used, default=FALSE
    #   gff3:               if gff3 txdb is used, default=FALSE
    # output:             pbs plot
    # return:             df of outlier genes
    anno_gene_all = c()
    meta_pbs_g_run_g1 = meta_pbs_g_run %>% dplyr::select(G1) %>% distinct() %>% pull
    for(i in 1:length(meta_pbs_g_run_g1)){
        # 1. run batch per G1
        tmp_run = meta_pbs_g_run %>%
            filter(G1 == meta_pbs_g_run_g1[i]) %>%
            unite("batch", c(rG1,rG2,rG3), sep="-", remove = FALSE)
        tmp_don_s = c()
        tmp_anno_s = c()
        for(n in 1:nrow(tmp_run)){
            # run pbs per batch
            
            tmp_pbs = run_pbs(
                p1=tmp_run[n,]$rG1,
                p2=tmp_run[n,]$rG2,
                p3=tmp_run[n,]$rG3,
                wdir=fst_dir)
            # add batch information
            name_batch = tmp_run[n,]$batch
            print(name_batch)
            tmp_pbs = tmp_pbs %>% 
                mutate(batch=name_batch)
            tmp_don_s = rbind(tmp_don_s, tmp_pbs)
            # get annotation
            if(gff3){
                # gff3 used
                tmp_anno_full = run_anno_gff3(tmp_pbs, tx_db, org_db, cutoff)
            }else if(isFALSE(flank)){
                # no flank
                tmp_anno_full = run_anno_noflank(tmp_pbs, tx_db, org_db, cutoff)
            }else{
                # flank db
                tmp_anno_full = run_anno(tmp_pbs, tx_db, org_db, cutoff)
            }
            
            tmp_anno_uniq = tmp_anno_full %>%
                group_by(gene) %>%
                filter(pbs == max(pbs)) %>%
                distinct(gene, .keep_all = TRUE) %>%
                mutate(batch=name_batch)
            tmp_anno_s = rbind(tmp_anno_s, tmp_anno_uniq)
            tmp_anno_full = tmp_anno_full %>%
                mutate(batch=name_batch)
            anno_gene_all = rbind(anno_gene_all, tmp_anno_full)
        }
        # 2. plot per G1
        # 2.1 get plot center for each chr
        axisdf = tmp_don_s %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
        # 2.2 precompute quantile for each plot
        tmp_quantile = tmp_don_s %>%
            group_by(batch) %>%
            summarize(threshold = quantile(pbs, probs = cutoff)) 
        # 2.2 plot
        tmp_don_s %>%
            # remove pbs <0.1
            filter(pbs>=0.1) %>%
            ggplot(aes(x=BPcum, y=pbs)) +
            
            # Show all points
            geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
            scale_color_manual(values = rep(c('#e41a1c', '#377eb8'), 22 )) +
            
            # custom X axis:
            scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, expand = expansion(mult = 0.01)) +
            # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
            
            geom_hline(data= tmp_quantile, linetype=2,
                       aes(yintercept=threshold)) +
            
            # one batch per line
            facet_grid(batch~., scales = "free") +
            
            # plot gene 
            geom_label_repel(data=tmp_anno_s, aes(label=gene),
                             size=4, nudge_y = 0.05) +
            
            # Custom the theme:
            theme(panel.background = element_rect(fill = NA, colour = "black", linetype = 1),
                  panel.grid.major = element_line(color = NA),
                  panel.grid.minor = element_line(color = NA),
                  legend.position = "none",
                  axis.title = element_text(size=20, face = "bold"),
                  axis.text = element_text(size = 16)) +
            xlab("Chromosome") + ylab("PBS")
        p_width=12
        p_height=min(3*(nrow(tmp_run) +1), 30) # max height: 30
        p_name= paste0(pdir,"pbs.anno.", meta_pbs_g_run_g1[i],".", cutoff, ".", p_width,"x", p_height, ".png")
        ggsave(p_name, width = p_width, height = p_height, device = "png", dpi = 500)
    }
    return(anno_gene_all)
}



run_anno = function(res_pbs, tx_db, org_db, cutoff=0.9995){
    # 1. get outlier regions
    threshold = quantile(res_pbs$pbs, probs = cutoff)
    pbs_top =  res_pbs %>% filter(pbs >= threshold)
    # 2. grange overlap, input cf3, db
    
    pbs_top_gr = GRanges(pbs_top$CHROM,
                         IRanges(start=pbs_top$BIN_START, end=pbs_top$BIN_END))
    pbs_top_gr$pbs = pbs_top$pbs
    pbs_top_gr$BPcum = pbs_top$BPcum
    
    pbs_top_gr_gene = findOverlaps(pbs_top_gr, tx_db)
    pbs_top_gr_gene_id = AnnotationDbi::select(
        org_db,
        keys=tx_db[pbs_top_gr_gene@to,]@ranges@NAMES,
        columns = c("SYMBOL", "GENENAME")) %>%
        dplyr::select(SYMBOL)
    pbs_top_gene_df = data.frame(
        as(pbs_top_gr[pbs_top_gr_gene@from,], "data.frame"),
        gene=pbs_top_gr_gene_id$SYMBOL)
        # group_by(gene) %>%
        # filter(pbs == max(pbs))
    return(pbs_top_gene_df)
}
run_anno_noflank = function(res_pbs, tx_db, org_db, cutoff=0.9995){
    # 1. get outlier regions
    threshold = quantile(res_pbs$pbs, probs = cutoff)
    pbs_top =  res_pbs %>% filter(pbs >= threshold)
    # 2. grange overlap, input cf3, db
    pbs_top_gr = GRanges(pbs_top$CHROM,
                         IRanges(start=pbs_top$BIN_START, end=pbs_top$BIN_END))
    pbs_top_gr$pbs = pbs_top$pbs
    pbs_top_gr$BPcum = pbs_top$BPcum
    
    pbs_top_gr_gene = findOverlaps(pbs_top_gr, genes(tx_db))
    pbs_top_gr_gene_id = select(
        org_db,
        keys=genes(tx_db)[pbs_top_gr_gene@to,]$gene_id,
        columns = c("SYMBOL", "GENENAME")) %>%
        dplyr::select(SYMBOL)
    pbs_top_gene_df = data.frame(
        as(pbs_top_gr[pbs_top_gr_gene@from,], "data.frame"),
        gene=pbs_top_gr_gene_id$SYMBOL)
    # group_by(gene) %>%
    # filter(pbs == max(pbs))
    return(pbs_top_gene_df)
}
run_anno_gff3 = function(res_pbs, tx_db, org_db, cutoff=0.9995){
    # run_anno:     get gene annotation for pbs result with flank regions for genes
    # input:    
    #   res_pbs:    pbs_result by run_pbs
    #   tx_db:      loaded gff3 tx_db for region annotation
    #   org_db:     loaded org_db for gene id annotation, not used
    # return:       data.frame with annotation information
    
    # 1. get outlier regions
    threshold = quantile(res_pbs$pbs, probs = cutoff)
    pbs_top =  res_pbs %>% filter(pbs >= threshold)
    # 2. grange overlap, input cf3, db
    
    pbs_top_gr = GRanges(str_sub(pbs_top$CHROM, start=4L),
                         IRanges(start=pbs_top$BIN_START, end=pbs_top$BIN_END))
    pbs_top_gr$pbs = pbs_top$pbs
    pbs_top_gr$BPcum = pbs_top$BPcum
    
    pbs_top_gr_gene = findOverlaps(pbs_top_gr, tx_db)
    pbs_top_gene_df = data.frame(
        as(pbs_top_gr[pbs_top_gr_gene@from,], "data.frame"),
        gene=tx_db[pbs_top_gr_gene@to,]$Name,
        gene_id=tx_db[pbs_top_gr_gene@to,]$gene_id) %>%
        mutate(gene=case_when(is.na(gene) ~ gene_id,
                               TRUE ~ gene))
    # group_by(gene) %>%
    # filter(pbs == max(pbs))
    return(pbs_top_gene_df)
}



# 4. load dbs
.libPaths("/projects/mjolnir1/people/gnr216/a-software/R-lib/4.2")
library(GenomicRanges)
#library('TxDb.Cfamiliaris.UCSC.canFam3.refGene')
#cf3 = TxDb.Cfamiliaris.UCSC.canFam3.refGene

# add flank 1kb for cf3
# cf3_up = flank(genes(cf3), width=1000 ,start = TRUE, both = FALSE)
# cf3_down = flank(genes(cf3), width=1000 ,start = FALSE, both = FALSE)
# cf3_f1kb = c(genes(cf3), cf3_up, cf3_down)
# cf3_f1kb=trim(unlist(range(split(cf3_f1kb, ~gene_id))))

txdb_new = import.gff3("/projects/mjolnir1/people/gnr216/r.ref/wolf/Canis_lupus_familiaris.CanFam3.1.104.gff3")
txdb_new_gene = txdb_new[txdb_new$type=="gene", ]
txdb_new_gene_up = flank(txdb_new_gene, width=1000 ,start = TRUE, both = FALSE)
txdb_new_gene_down = flank(txdb_new_gene, width=1000 ,start = FALSE, both = FALSE)
# NB! this is a stupid way, just keep all three in one grange object. Shoulbe be merged. 
txdb_new_gene_f1kb = unlist(as(list(txdb_new_gene,txdb_new_gene_up, txdb_new_gene_down), "GRangesList"))


library("org.Cf.eg.db")
org_db = org.Cf.eg.db

# 5. run pbs
test_anno = run_pbs_anno_batch(
    meta_pbs_g_run,
    fst_dir = './fst/',
    pdir='./fst_plot/',
    cf3_f1kb,
    org_db,
    cutoff=0.9995,
    flank = TRUE
)



```

## 2. visualization


## 3. summary
