# code: https://github.com/mskcc/facets
# 


#fitrds = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/facets/x185_006_T0_fit.rds"
facets_purity <- function(fitrds){
  fit = read_rds(fitrds)
  round(fit$purity, 2)
}



facets_read_to_granges <- function(...){
  facets_read.segfl(...)
}


# out: data frame of segment summaries pre and post clustering of
#       segments. The columns are: 
# ‘chrom’ the chromosome to which the segment belongs; ‘seg’ the segment number; 
# ‘num.mark’ the number of SNPs in the segment; 
# ‘nhet’ the number of SNPs that are deemed heterozygous; 
# ‘cnlr.median’ the median log-ratio of the segment; 
# ‘mafR’ the log-odds-ratio summary for the segment; 
# ‘segclust’ the segment cluster to which segment belongs; 
# ‘cnlr.median.clust’ the median log-ratio of the segment cluster; 
# ‘mafR.clust’ the log-odds-ratio summary for the segment cluster;
# ‘cf’ the cellular fraction of the segment; 
# ‘tcn’ the total copy number of the segment;
# ‘lcn’ the minor copy number of the segment.
# https://user-images.githubusercontent.com/1105264/53222479-4ff3cd80-36b1-11e9-8137-6839027c4364.png
# cncf: dataframe consisting of the columns of segmentation output as
#       well as cellular fraction (cf), total (tcn) and lesser (lcn)
#       copy number of each segment and their em counterpart (with
#       .em suffix)
read_facets_seg_igv.trk <- function(df_trk, 
                                col_fl, 
                                col_samplename){
  
  df_trk = data.frame(df_trk, stringsAsFactors = F)
  df_seg = lapply(1:nrow(df_trk), function(i){
    
    samplename = df_trk[i, col_samplename]
    seg_fl = df_trk[i, col_fl]
    message(samplename, " ", appendLF = F)
    
    df_seg = read_tsv(seg_fl, 
                      col_types = cols(chrom = col_double(),
                                       seg = col_double(),
                                       num.mark = col_double(),
                                       nhet = col_double(),
                                       cnlr.median = col_double(),
                                       mafR = col_double(),
                                       segclust = col_double(),
                                       cnlr.median.clust = col_double(),
                                       mafR.clust = col_double(),
                                       start = col_double(),
                                       end = col_double(),
                                       cf.em = col_double(),
                                       tcn.em = col_double(),
                                       lcn.em = col_double())) %>% clean_names()
    df_seg %<>% mutate(samplename = samplename) %>% 
      mutate(major_cn = tcn_em - lcn_em) %>% 
      # this file is perfect for viewing in IGV
      dplyr::select(samplename, chrom, start, end,
                    num_mark = num_mark, 
                    seg_num = seg, 
                    
                    everything(), 
                    
                    log2_rd_r = cnlr_median, # log-ratio      (read depth ratio)      log2
                    log_baf_r = maf_r,       # log-odds-ratio (ratio of allele freqs) log_e  
                    
                    cluster = segclust,
                    minor_cn = lcn_em,
                    major_cn,
                    total_cn = tcn_em)
    df_seg
  }) %>% bind_rows()
  df_seg
  
}

# A SEG file (segmented data; .seg or .cbs) is a tab-delimited text file that lists loci and associated numeric values. 
# The first row contains column headings and each subsequent row contains a locus and an associated numeric value. 
# IGV ignores the column headings. 
# It reads the first four columns as track name, chromosome, start location, and end location. 
# It reads the last column as the numeric value for that locus (if the value is non-numeric, IGV ignores the row). 
# IGV ignores all other columns.
# 
# The segmented data file format is the output of the Circular Binary Segmentation algorithm (Olshen et al., 2004).


read_facets.trk.fits <- function(df_trk){
  
  df_facet_summ <- lapply(1:nrow(df_trk), function(i){
    message(i, appendLF = F)
    seg_fl = df_trk$facets_fl_full[i]
    fit_fl = gsub("_cncf.tsv", ".rds", seg_fl)
    
    fit = read_rds(fit_fl)
    data.frame(
      
      t_path_sample_id = df_trk$t_path_sample_id[i],
      n_path_sample_id = df_trk$n_path_sample_id[i],
      t_sample_type = df_trk$t_sample_type[i],
      path_patient_id = df_trk$path_patient_id[i],
      
      purity = fit$purity, 
      ploidy = fit$ploidy, 
      # if null, or absent, return NA
      loglik = ifelse(is.null(fit$loglik), NA, fit$loglik),
      dipLogR = fit$dipLogR)
  }) %>% bind_rows()
  df_facet_summ
}


