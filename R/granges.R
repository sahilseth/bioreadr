gr_intersect <- function(df_a, df_b){
  p_load(tictoc)
  tic()
  message("get genomic ranges")
  gr_a <- makeGRangesFromDataFrame(df_a, keep.extra.columns = T)
  gr_b <- makeGRangesFromDataFrame(df_b, keep.extra.columns = T)
  
  message("find overlaps")
  # using genomic ranges:
  #pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
  pairs <- findOverlaps(gr_a, gr_b, ignore.strand = TRUE)
  #ans <- pintersect(pairs, ignore.strand = TRUE)
  
  message("merge")
  df_merge = bind_cols(as.data.frame(gr_a[queryHits(pairs)]),
                       as.data.frame(gr_b[subjectHits(pairs)]))
  
  toc()
  df_merge
  
}
