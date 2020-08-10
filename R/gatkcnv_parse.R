to_gr_seg.gatkcna <- function(segfl, 
                              sample_id,
                              gencode_fl = "~/.rsrch3/home/iacs/sseth/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene_hgnc.bed",
                              lst_gen = NULL,
                              verbose = F){
  p_load(futile.logger)
  flog.debug("reading segfile")
  seg = data.table::fread(cmd = glue("grep -v '@' {segfl}"), data.table = F) %>% as_tibble() %>% 
    clean_names() %>% 
    mutate(sample_id = sample_id)
  seg
  seg$contig %>% table()
  
  # # http://seqanswers.com/forums/archive/index.php/t-32100.html
  # message("switch chr names")
  #   mutate(chr = gsub("23", "X", chr),
  #          chr = gsub("24", "Y", chr))
  
  # convert to GR
  lst = seg %>% to_gr_seg.seg(col_sample = "sample_id", 
                              col_chr = "contig", 
                              col_start = "start", col_end = "end", 
                              col_num_mark = "num_points_allele_fraction", 
                              col_seg_mean = "log2_copy_ratio_posterior_50")
  
  # convert to GR
  if(verbose)
    message("annotate")
  if(is.null(lst_gen))
    lst_gen <- to_grseg.gencode(gencode_fl = "~/.rsrch3/home/iacs/sseth/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene_hgnc.bed")

  source('~/Dropbox/projects/packs_dnaseq/R/to_mae.seg.R')
  length(lst$gr_seg)
  df_ann = annotate.gr_seg(lst$gr_seg, lst_gen = lst_gen)
  dim(df_ann)
  # df_ann = df_ann %>% group_by(gene_name) %>% add_count()
  # df_ann %>% filter(n>1) %>% View()
  
  list(df_ann = df_ann, df_igv_seg = lst$df_igv_seg, gr_seg = lst$gr_seg)
}
