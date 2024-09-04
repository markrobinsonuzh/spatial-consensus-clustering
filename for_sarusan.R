
library(readr)
library(dplyr)

rt <- read_tsv("./data/C57BL6J-638850_7200/spaGCN/config_default/combined_domains.tsv",
         col_names = c("id", paste0("n_clusters_",4:20)),
         show_col_types = FALSE, skip = 1)

rt %>% pull("n_clusters_20") %>% table
  
rt %>% pull("n_clusters_18") %>% table

rt %>% pull("n_clusters_15") %>% table
