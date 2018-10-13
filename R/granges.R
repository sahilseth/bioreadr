gr_intersect <- function(df_a, df_b){
  # get genomic ranges:
  gr_a <- makeGRangesFromDataFrame(df_a, keep.extra.columns = T)
  gr_b <- makeGRangesFromDataFrame(df_b, keep.extra.columns = T)
  
  # using genomic ranges:
  #pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
  pairs <- findOverlaps(gr_a, gr_b, ignore.strand = TRUE)
  #ans <- pintersect(pairs, ignore.strand = TRUE)
  
  df_merge = bind_cols(as.data.frame(gr_a[queryHits(pairs)]),
                       as.data.frame(gr_b[subjectHits(pairs)]))
  
  df_merge
  
}