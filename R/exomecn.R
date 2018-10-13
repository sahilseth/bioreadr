

read_gencode_gene <- function(
  gencode_gene_fl = "~/.rsrch2/iacs/iacs_dep/sseth/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene.bed"){
  df_gen = read_tsv(gencode_gene_fl) %>% clean_names() %>% 
    dplyr::rename(chr = number_seqnames) %>% 
    dplyr::mutate(chr = gsub("chr", "", chr))
  df_gen
}

exomecn_ann <- function(x, y, df_gen){
  
  # assuming x can be a vector of file names
  if(length(x) != length(y))
    stop("x and y, need to be of same length")
  
  df_cnv = lapply(seq_along(x), function(i){
    message(".", appendLF = F)
    df_cnv = read_tsv(x[i]) %>% clean_names() %>% 
      mutate(samplename = x[i])
  }) %>% bind_rows()
  
  df_cnv_ann = gr_intersect(df_gen, df_cnv)
  df_cnv_ann

}
  
  
  
  
  
  
  
  
