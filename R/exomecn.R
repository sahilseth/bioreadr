
# https://bioconductor.org/packages/release/bioc/vignettes/CNVRanger/inst/doc/CNVRanger.html
# 0: homozygous deletion (2-copy loss)
# 1: heterozygous deletion (1-copy loss)
# 2: normal diploid state
# 3: 1-copy gain
# 4: amplification (>= 2-copy gain)

get_state.exomecn <- function(x){
  case_when(
    x <= -0.6 ~ 0,
    x > -0.6 &  x <= -0.3 ~ 1,
    
    x > -0.3 & x < 0.3 ~ 2,

    x > 0.3 & x < 0.9 ~ 3,
    x >= 0.9 ~ 4)
}



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
                                col_samplename = "t_path_sample_id", 
                                cores = 10){
  
  wranglr::expect_columns(df_trk, c(col_fl, col_samplename))

  df_seg = mclapply(1:nrow(df_trk), function(i) {
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
  }, mc.cores = cores) %>% bind_rows()
  df_seg
}

# read all the files
# to GR
# annotate
# handle extemes
to_mae.exomecn.seg <- function(df_exomecn_seg,
                               gencode_fl = "/rsrch3/home/iacs/sseth/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene.bed"){

  source("~/projects/packs_dnaseq/R/to_mae.seg.R")
  gr_seg = df_exomecn_seg %>% to_gr_seg.seg(
                       col_sample = "samplename",
                       col_chr = "chrom",
                       col_start = "loc_start",
                       col_end = "loc_end",
                       col_num_mark = "num_mark",
                       col_seg_mean = "seg_mean",
                       save_fls = FALSE)

  lst_gen = to_grseg.gencode(gencode_fl)
  names(lst_gen)

  class(gr_seg)
  df_cnv_ann_lng = annotate.gr_seg(gr_seg$gr_seg, lst_gen = lst_gen)

  # create SE
  se = to_se.ann_seg.exomecn(df_cnv_ann_lng, lst_gen = lst_gen)
  
}






