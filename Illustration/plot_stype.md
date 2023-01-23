# this script summarise plot styles
# legend
# theme
# sample meta from google doc
# output: df
``` r
#### plot legend ####
## 1. population with shape and color
#### 3.1 read meta information from google docs ####
library(googlesheets4)
gs4_deauth()

meta_legend = read_sheet(
    "https://docs.google.com/spreadsheets/d/1s9wQQ85gLAFcWb1AI8uJ0eHqR_ScqwzCnorYMxcHAa8/edit#gid=1963352145",
    sheet="plot_legend")
rownames(meta_legend) = meta_legend$Pop
```
``` r
#### plot theme ####
# 1. pca
p_theme_pca = theme(
    panel.background = element_rect(fill = 'white', colour = "black", linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

# 2. 
```
``` r
#### sample meta ####
# from google doc
library(googlesheets4)
gs4_deauth()

sample_meta = read_sheet(
    "https://docs.google.com/spreadsheets/d/1s9wQQ85gLAFcWb1AI8uJ0eHqR_ScqwzCnorYMxcHAa8/edit#gid=1963352145",
    sheet="meta")
rownames(sample_meta) = sample_meta$id
```

