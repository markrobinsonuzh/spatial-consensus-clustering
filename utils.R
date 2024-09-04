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


align_classes <- function(df, ref) {
  
  bcs <- df %>% as.data.frame
  levs <- lapply(bcs, levels)
  n_levs <- sapply(levs, length)
  
  cols_to_change <- setdiff(colnames(bcs), ref)
  
  for(i in cols_to_change) {
    if(n_levs[i] != n_levs[ref]) {
      message(paste0(i," has ",n_levs[i],
                     " levels; reference (", ref, 
                     ") has ",n_levs[ref], " (not modifying)"))
      next
    }
    hung <- clue::solve_LSAP(table(bcs[,ref], bcs[,i]), 
                             maximum = TRUE)
    lookup <- cbind(seq_along(hung), levs[[i]][hung])
    levels(bcs[,i]) <- lookup[order(as.integer(lookup[,2])),1]
    bcs[,i] <- as.factor(as.character(bcs[,i]))
  }
  
  bcs
  
}

spot_entropy <- function(spatial_coords, label, k=6) {
  require(dbscan)
  knns <- kNN(spatial_coords, k=k)
  label <- as.factor(label)
  neighb_labels <- apply(knns$id, 2, function(u) label[u])
  apply(neighb_labels, 
        1, function(u) calc_entropy(table(factor(u,
                                                 levels=levels(label)))))
}

