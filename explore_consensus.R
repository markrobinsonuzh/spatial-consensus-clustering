suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(clue)
  library(khroma)
  library(scran)
  library(limma)
  library(tibble)
  #library(googlesheets4)
  library(readr)
  library(ggrepel)
  library(mclust)
  library(pheatmap)
})

datadir <- "data"
spe <- readRDS(file.path(datadir,
                         "obs_Br8100_151673_spe.RDS"))
df <- data.frame(barcode = colnames(spe), label = spe$label)
# spe <- readRDS(file.path(datadir,
#                          "obs_Br8100_151673_spe_with_base_clusterings.RDS"))
spe <- readRDS(file.path(datadir,
                         "obs_Br8100_151673_spe_with_all_base_clusterings.RDS"))



calc_entropy <- function(u) {
  p <- u[u>0]
  p <- p/sum(p)
  -sum(p*log(p))
}


# given a matrix of labels, calculate all pairwise ARIs
calc_aris <- function(m) {
  require(mclust)
  a <- diag(ncol(m))
  for(i in 1:(ncol(m)-1))
    for(j in 2:ncol(m))
      a[i,j] <- a[j,i] <- mclust::adjustedRandIndex(m[,i], m[,j])
  rownames(a) <- colnames(a) <- colnames(m)
  a
}


cd <- colData(spe)[,-c(1:5)]
stopifnot(all(df$barcode==colnames(spe)))
cd$label <- df$label

cd <- cd[,-grep("_MajorityVote", colnames(cd))]
cd <- cd[,-grep("_KModes", colnames(cd))]
cd <- cd[,-grep("_LCA", colnames(cd))]

arism <- calc_aris(cd)

ph <- pheatmap(arism)


num_methods <- function(u) {
  s <- sapply(strsplit(u, "[_\\.]"), .subset, 1)
  g <- grep("label",s)
  if(length(g)) 
    s <- s[-g]
  s %>% table %>% length
}


calc_summary <- function(tree, arism, k = 3, gtcol = "label") {
  ct <- cutree(tree, k=k)
  tct <- table(ct)
  data.frame(group = names(tct), 
             number = as.integer(tct),
             ari_to_truth = sapply(names(tct),
                                   function(u) median(arism[ct==u & names(ct)!=gtcol,gtcol])),
             ari_block = sapply(names(tct),
                                   function(u) {z <- arism[ct==u,ct==u]; median(z[upper.tri(z)]) }),
             k = k,
             gt_in_group = (names(tct) == ct[gtcol]),
             number_distinct_methods = sapply(names(tct),
                                   function(u) num_methods(names(ct)[ct==u])))
}

stats <- lapply(3:10, function(u) calc_summary(ph$tree_col, 
                                      arism, k=u)) %>% bind_rows


# ggplot(stats, aes(x=number_distinct_methods, y=ari_to_truth, 
#                   colour = gt_in_group)) + 
#   geom_jitter(width=.1, height=.01, aes(size = number))

ggplot(stats, aes(x=number_distinct_methods, y=ari_block, 
                  shape = gt_in_group, colour = ari_to_truth)) + 
  geom_jitter(width=.1, height=.01, aes(size = number))
