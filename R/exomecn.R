



exomecn_ann <- function(x, y, df_gen){
  pacman::p_load(tidyranges)
  
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


read_exomecn.trk.seg <- function(df_trk, 
                                col_fl = "exomcn_fl_full", 
                                col_samplename = "t_path_sample_id"){
  
  df_seg = lapply(1:nrow(df_trk), function(i) {
    samplename = df_trk[i, col_samplename] %>% unlist()
    seg_fl = df_trk[i, col_fl] %>% unlist()
    # print(class(seg_fl))
    message(samplename, appendLF = F)
    df_seg = read_tsv(seg_fl, 
                      col_types = cols(
                        ID = col_character(),
                        chrom = col_character(),
                        loc.start = col_double(),
                        loc.end = col_double(), 
                        num.mark = col_double(), 
                        seg.mean = col_double())) %>% clean_names()
    df_seg %<>% mutate(samplename = samplename) %>% 
      # this file is perfect for viewing in IGV
      dplyr::select(samplename, chrom, loc_start, loc_end,
                    num_mark, seg_mean, id, everything())
    df_seg
  }) %>% bind_rows()
  df_seg
}







