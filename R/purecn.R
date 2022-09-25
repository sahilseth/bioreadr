
to_purecn_seg.gatkcnv <- function(seg_fl, out_segfl){
  flog.info(paste0("write out: ", out_segfl))
  p_load(dplyr)
  df_seg = read_counts.gatk(seg_fl) %>% 
    mutate(ID = sampleid, 
           contig = switch_chr_names.to_int(contig),
           contig = as.integer(contig)) %>% 
    dplyr::select(ID, chrom = contig, 
                  loc.start = start, loc.end = end, 
                  num.mark = num_points_copy_ratio, 
                  seg.mean = log2_copy_ratio_posterior_50)
  df_seg$chrom
  tail(df_seg)
  # need to write it out, since runabs expects a file :(
  # tmpfl = tempfile(fileext = ".txt");tmpfl
  write_tsv(df_seg, out_segfl)
  out_segfl
}


# read purecn -------
to_gr_seg.purecn <- function(segfl,
                            variantfl = NULL,
                            calc_minor_major,
                             gencode_fl = "~/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene_hgnc.bed",
                             lst_gen = NULL,
                             purity = NULL, ploidy = NULL,
                             verbose = FALSE){
  # segfl = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/runs/185_059/tmp/purecn2/IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_dnacopy.seg"  
  # library(my.ultraseq)
  seg = data.table::fread(segfl, data.table = F) %>% as_tibble() %>% 
    clean_names() %>% 
    mutate(chrom = chrom %>% my.ultraseq:::recode_chr_names.to_chr())
  seg$chrom %>% table()
  if(!is.null(purity) & !is.null(ploidy)){
    if(verbose)
      flog.info("adding adjust copy ratios")
    seg %<>% 
      mutate(seg_mean_adjusted = log2(PureCN:::.calcExpectedRatio(c, purity, ploidy)))
  }
  # calc minor major:
  if(!is.null(variantfl)){
    message("calculating major and minor alleles")
    df_var = read_csv(variantfl)
    # combine the two:
    # https://github.com/lima1/PureCN/issues/85
    # If it asks for major and minor copy number (maybe mentions b-allele or SNP array), 
    # it is almost certainly ML.C-ML.M.SEGMENT and ML.M.SEGMENT, respectively.
    # other option is genome_intersect
    seg2 = seg %>% mutate(chr = chrom, start = loc_start, end = loc_end)
    seg3 =  seg2 %>% 
      fuzzyjoin::genome_join(select(df_var, chr, start, end, ML.M, ML.C, ML.M.SEGMENT), 
        by = c("chr", "start", "end"), mode = "left")
    # now we have a lot more rows than expected:
    # 30,172, instead of 951
    # get columns seg
    seg3 = group_by(seg3, chrom, loc_start, loc_end) %>% 
      mutate(ML.M = mean(ML.M, na.rm = TRUE), 
             ML.C = mean(ML.C, na.rm = TRUE), 
             ML.M.SEGMENT = mean(ML.M.SEGMENT)) %>% 
      select(id:seg_mean_adjusted, ML.C, ML.M, ML.M.SEGMENT) %>% 
      unique() %>% 
      mutate(minor_cn = ML.M.SEGMENT, 
             major_cn = ML.C-ML.M.SEGMENT) %>% 
      filter(!is.na(major_cn)) %>% ungroup()
    # seg3
  }
  
  # convert to GR
  lst = seg3 %>% to_gr_seg.seg(col_sample = "id", 
                              col_chr = "chrom", 
                              col_start = "loc_start", col_end = "loc_end", 
                              col_num_mark = "num_mark", col_seg_mean = "seg_mean")
  
  # read igvseg
  
  # read subclonal seg for igv
  # segfl_sub = gsub(".segs.txt$", ".titan.txt", segfl)
  # df_sub_seg = data.table::fread(segfl_sub, data.table = F)
  # head(df_sub_seg)
  
  # convert to GR
  if(verbose)
    message("annotate")
  # use the supplied file
  if(is.null(lst_gen))
    lst_gen <- to_grseg.gencode(gencode_fl)
  
  source('~/.drop/projects/packs_dnaseq/R/to_mae.seg.R')
  df_ann = annotate.gr_seg(lst$gr_seg, lst_gen = lst_gen)
  dim(df_ann)
  # df_ann = df_ann %>% group_by(gene_name) %>% add_count()
  # df_ann %>% filter(n>1) %>% View()
  
  list(df_ann = df_ann, df_igv_seg = lst$df_igv_seg, gr_seg = lst$gr_seg)
  
}


to_gr_seg.purecnloh <- function(segfl,
                             gencode_fl = "~/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene_hgnc.bed",
                             lst_gen = NULL,
                             purity = NULL, ploidy = NULL,
                             verbose = FALSE){
  
  seg = data.table::fread(segfl, data.table = F) %>% as_tibble() %>% 
    clean_names() %>% 
    mutate(chr = chr %>% recode_chr_names.to_chr())
  seg$chr %>% table()
  if(!is.null(purity) & !is.null(ploidy)){
    if(verbose)
      flog.info("adding adjust copy ratios")
    seg %<>% 
      mutate(seg_mean_adjusted = log2(PureCN:::.calcExpectedRatio(c, purity, ploidy)))
  }
  
  # convert to GR
  lst = seg %>% to_gr_seg.seg(col_sample = "sampleid", 
                              col_chr = "chr", 
                              col_start = "start", col_end = "end", 
                              col_num_mark = "num_mark", col_seg_mean = "seg_mean_adjusted")
  
  # read igvseg
  
  # read subclonal seg for igv
  # segfl_sub = gsub(".segs.txt$", ".titan.txt", segfl)
  # df_sub_seg = data.table::fread(segfl_sub, data.table = F)
  # head(df_sub_seg)
  
  # convert to GR
  if(verbose)
    message("annotate")
  # use the supplied file
  if(is.null(lst_gen))
    lst_gen <- to_grseg.gencode(gencode_fl)
  
  # source('~/.drop/projects/packs_dnaseq/R/to_mae.seg.R')
  df_ann = annotate.gr_seg(lst$gr_seg, lst_gen = lst_gen)
  dim(df_ann)
  # df_ann = df_ann %>% group_by(gene_name) %>% add_count()
  # df_ann %>% filter(n>1) %>% View()
  
  list(df_ann = df_ann, df_igv_seg = lst$df_igv_seg, gr_seg = lst$gr_seg)
  
}



to_mut.purecn <- function(mutfl,
                          purity = NULL, ploidy = NULL,
                          verbose = FALSE){
  p_load(tidyverse)
  df = readr::read_csv(mutfl) %>% janitor::clean_names()
  df
}

purecn_read_rds <- function(rdsfl){
    p_load(PureCN)
    ret = read_rds(rdsfl)
    # The purity/ploidy combinations are sorted by likelihood and stored in ret$results.
    names(ret)
    # [1] "candidates" "results"    "input"     
    muts = predictSomatic(ret)
    
    dim(muts)
    ret$input$vcf@rowRanges
    
    filter(muts, gene.symbol == "TP53")
    filter(muts, gene.symbol == "PIK3CA")
    
    # ML.M: no. of chromosomes harboring the variant or mutation
    # 
    
}




           
