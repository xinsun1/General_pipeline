
## 1. pca plot function in R
GL result plot (with/without projection)
``` r
#### gl pca ####
### pcangsd result


# function to generate PCA results from PCAngsd result
# input:
#   cov_file:       covariant matrix from pcangsd result
#   plot_legend:    plot legend df
#   sample_meta:    sample meta information df, rownames=id
#   list_full:      sample id list for full cov (ordered by gl files)
#   list_pro:       sample id list for projection, default=FALSE
#   list_se:        sample id list for pca
#   list_id_plot:   sample id list for plot 
#   plot_group:     resulst plot group in meta, default_var = 'pop'
#   pro_dim:        dimensiion to project, default=c(1,2)
# output:
#   ggplot object without p_theme
plot_pcangsd = function(cov_file, plot_legend, sample_meta,
                        list_full, list_se,
                        plot_group='pop',
                        list_id_plot=FALSE,list_pro=FALSE, pro_dim=c(1,2)){
    
    # 1. set env
    library(tidyverse)
    library(vegan)
    library(ggrepel)
    
    # 2. read files
    cov_full = as.matrix(read.table(cov_file))
    list_se_idx = match(list_se, list_full)

    ## 3.1 perform pca
    cov_se = cov_full[list_se_idx, list_se_idx]
    se_de = eigen(cov_se) # decomposite covariance matrix
    evec_se = se_de$vectors[,pro_dim]
    eval_se = se_de$values
    var_explained = se_de$values/sum(se_de$values)
    evec_all = cbind(list_se_idx, evec_se) # sample_index, PCdim1, PCdim2
    
    ## 3.2 perfrom projection
    if(! isFALSE(list_pro)){
        list_pro_idx = match(list_pro, list_full)
        # project one at a time
        for(i in list_pro_idx){
            # same as step 3
            list_tmp = c(list_se_idx, i)
            cov_pro = cov_full[list_tmp, list_tmp]
            pro_de = eigen(cov_pro)
            evec_pro = pro_de$vectors[,pro_dim]
            
            # procrustes and project the new sample
            evec_pro_se = evec_pro[-nrow(evec_pro),] # remove the sample for projection
            tmp_projection = procrustes(evec_se, evec_pro_se) # generate procrustes parameters
            tmp_evec = predict(tmp_projection, t(evec_pro[nrow(evec_pro), pro_dim]))
            evec_all = rbind(evec_all, c(i, tmp_evec[1,]))
        }
    }
    colnames(evec_all) = c('sample_index', 'PCdim1', 'PCdim2')
    evec_all = as.data.frame(evec_all)
    evec_all = evec_all %>%
        # sort by sample_index
        arrange(sample_index) %>%  
        mutate(
            # add id
            id = list_full[sample_index]) %>%
        # add plot group
        mutate(
            pop = sample_meta[id,] %>% pull(as.name(plot_group))
            # pop = sample_meta[id,]$pop
            )
    
    # 4. plot result
    p = evec_all %>%
        ggplot(aes(x=PCdim1, y=PCdim2)) +
        geom_point(aes(fill= pop, shape=pop), size =3, alpha = 0.8) 
    
    if(! isFALSE(list_id_plot)){
        # if plot id to the result
        p = p +
            geom_text_repel(
                data = ~subset(., list_id_plot %in% evec_all$id),
                aes(label=id), 
                max.overlaps = Inf, # show overlapped points
                box.padding = 0.5, 
                seed=88)
    }
    
    # add legend and title
    p = p + 
        scale_fill_manual(breaks = plot_legend$pop,
                          values = plot_legend$color) +
        scale_shape_manual(breaks = plot_legend$pop,
                           values = plot_legend$shape) +
        xlab(str_c('PC',pro_dim[1],': ', sprintf("%1.2f%%", 100*var_explained[pro_dim[1]]))) +
        ylab(str_c('PC',pro_dim[2],': ', sprintf("%1.2f%%", 100*var_explained[pro_dim[2]])))
    
    # 5. return result
    return(p)
}

#### procrustes projection with PCangsd ####
# 1. set env
library(tidyverse)

# 2. prepare files
wdir='~/Documents/Projects/Tiger_ancient_historical/South_China_tiger/Progress_project/population_structure/gl_pca/'
setwd(wdir)

list_pro = c("RUSA21", "M2")
list_pro = FALSE
list_full = sample_meta %>%
    group_by(GL_order) %>%
    arrange() %>%
    pull(id)

list_se = list_full[! list_full %in% list_pro]

# 3. get output
p_pca = plot_pcangsd('./gl.tv.maf05.e6.i1k.cov', meta_legend, sample_meta,
             list_full, list_se, 
             list_id_plot = list_full, pro_dim = c(1,2)) + p_theme_pca
p_pca12 =
    plot_pcangsd(
    './gl.tv.maf05.e6.i1k.cov', meta_legend, sample_meta,
    list_full, list_se, 
    list_id_plot = list_full, pro_dim = c(1,2)) +
    p_theme_pca +
    theme(
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.background = element_rect(color="grey"),
        aspect.ratio = 1
        ) +
    guides(fill=guide_legend(ncol=1))

p_pca13 =
    plot_pcangsd(
        './gl.tv.maf05.e6.i1k.cov', meta_legend, sample_meta,
        list_full, list_se, 
        list_id_plot = list_full, pro_dim = c(1,3)) +
    p_theme_pca +
    theme(
        legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.background = element_rect(color="grey"),
        aspect.ratio = 1
    ) +
    guides(fill=guide_legend(ncol=1))

# 3.1 merege multiple plots
library(cowplot)
plot_grid(p_pca12, p_pca13, labels = NULL )

ggsave('pca_gl.grid.pca12-13.12x6.png',width = 12, height = 6, device = "png", dpi = 500)

```
