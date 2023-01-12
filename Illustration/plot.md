





# 1. relatedness plot
 with `lattice`
Example code
``` R 
#### mikkel data ####

## 0. set env
setwd("~/Documents/Projects/wolf project/9.inbreed/mikkel_rel/")
library(tidyverse)
library(RColorBrewer)

## 1.read raw data 
raw_rel = read_delim(file = './collected (1).txt',
                     delim = "\t")
raw_king = read_delim(file = './King.txt',
                     delim = "\t")


# 2. plot_king
rel_mod = raw_king %>%
    # get id
    separate(id, c("s1", "s2"), sep="-") %>%
    
    pivot_wider(names_from = s2, values_from = king) %>%
    column_to_rownames(var = 's1')


rel_mod = raw_rel %>%
    select(a, b, KING) %>%
    # add id name
    mutate(a=list_id[a+1,]$id,
           b=list_id[b+1,]$id) %>%
    
    pivot_wider(names_from = b, values_from = KING) %>%
    column_to_rownames(var = 'a') 



library("lattice")
color_list= c("#01415B", "#005148", "#019587", "#A6BC09", "#CCEA8D")
color_list= c("#363432", "#196774", "#90A19D" ,"#F0941F",  "#EF6024")


library(Cairo)
CairoPNG(filename = "rel_mikkel.png", width = 2000, height = 2000, dpi=300)
levelplot(as.matrix(rel_mod), scales=list(x=list(rot=90)),
          at=c(-1.55,0.0442,0.0884,0.177,0.354,0.5), 
          col.regions=color_list,
          xlab=NULL,ylab=NULL,
)
dev.off()


```


