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


# x = "bed/gatkcnv/AmpliSeqExome.20131001.designed.interval_list"
to_gr.interval_list <- function(x){
  # p_load(plyranges)
  df = data.table::fread(cmd = glue("grep -v '@' {x}"), 
                         data.table = F, col.names = c("seqnames", "start", "end", "strand", "name")) %>% 
    as_tibble()
  plyranges::as_granges(df)
}


to_gr.gatkcnv.cr <- function(x){
  p_load(plyranges)
  df = data.table::fread(cmd = glue("grep -v '@' {x}"), data.table = F) %>% 
    as_tibble() %>% clean_names() %>% 
    dplyr::rename(seqnames = contig)
  plyranges::as_granges(df)
}


to_gr.gatkcnv.baf <- function(x){
  
  df = data.table::fread(cmd = glue("grep -v '@' {x}"), data.table = F) %>% 
    as_tibble() %>% clean_names() %>% 
    dplyr::rename(seqnames = contig, start = position) %>% 
    dplyr::mutate(end = start)
  plyranges::as_granges(df)
  
}





# END

