# devtools::install_github("https://github.com/Crick-CancerGenomics/ascat", subdir = "ASCAT")

# https://github.com/Crick-CancerGenomics/ascat/issues/19
# Hi Tommy,
# Unfortunately, we don’t have a polished pipeline to use sequencing data as ASCAT input just yet, but we are working on an integration for the ASCAT package.
# I can describe the workflow that we are planning to integrate which has worked quite well so far.
# Assuming that you start from the bam file, this is how we usually go about it:
#   We use alleleCount (https://github.com/cancerit/alleleCount) to run on bam files, with the SNP loci from the 1000 genomes project.
# To obtain the logR, we then use the read counts in the tumour sample (=countTumour) and in the normal (=countNormal).
# In a first step we just divide the tumour by the normal counts: tumour_logR=countTumour/countNormal.
# Then we transform the value into logarithm while at the same time normalising by the mean, as that helps centering it around 0:
#   log2(tumour_logR/mean(tumour_logR)). The logR of the matched normal is afterwards simply set to 0.
# As alleleCount returns all counts, the BAF is calculated as the ratio between the maximum count and the complete coverage of the locus.
# Afterwards you should have tumour and normal LogR and BAF values that you can use in an ASCAT analysis.
# 
# After running the steps mentioned above, you have logR and BAF files which you can use in ASCAT analogous to array data.
# So we just add some pre-processing to create the correct input from the sequencing data, the ASCAT functions stay the same. Just make sure
# to set the gamma parameter to 1 when working with sequencing data.
# 
# Best wishes,
# Kerstin

# this is good for multisample cases
to_ascat_logr.gatkcnv.cr <- function(fls, samps, outfile){
  
  i=2
  tmp = lapply(seq_along(fls), function(i){
    fl = fls[i];samp = samps[i]
    message(samp)
    df = data.table::fread(cmd = glue("grep -v '@' {fl}"), data.table = F) %>% 
      as_tibble() %>% clean_names() %>% 
      mutate(sampleid = samp)
  }) %>% bind_rows()
  tmp %<>% mutate(width = end - start,
                  pos = start + width/2,
                  region = glue("{contig}-{start}-{end}"),
                  log2_copy_ratio = round(log2_copy_ratio, 4))
  mat = tmp %>%
    pivot_wider(names_from = "sampleid", values_from = "log2_copy_ratio", id_cols = c("region", "contig", "pos")) %>% 
    select(region, chrs = contig, pos, everything()) %>% to_mat()
  head(mat)
  # tail(mat)
  write.table(mat, outfile, sep = "\t", quote = F)
  mat
}

to_ascat_baf.gatkcnv.hets <- function(fls, samps, outfile){
  source('~/Dropbox/public/github_wranglr/R/to_mat_df.R')
  
  i=2
  tmp = lapply(seq_along(fls), function(i){
    fl = fls[i];samp = samps[i]
    message(samp)
    df = data.table::fread(cmd = glue("grep -v '@' {fl}"), data.table = F) %>% 
      as_tibble() %>% clean_names() %>% 
      mutate(sampleid = samp)
  }) %>% bind_rows()
  tmp %<>% mutate(dp = ref_count + alt_count,
                  af = alt_count/dp) %>% 
    mutate(snp = glue("{contig}-{position}"))
  
  mat = tmp %>%
    pivot_wider(names_from = "sampleid", values_from = "af", id_cols = c("snp", "contig", "position")) %>% 
    select(snp, chrs = contig, position, everything()) %>% to_mat()
  head(mat)
  # tail(mat)
  write.table(mat, outfile, sep = "\t", quote = F)
  mat
}

# source("~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/granges.R")
ascat_prep_inputs <- function(trk){
  i=1
  tmp = mclapply(1:nrow(trk), function(i){
    
    # setwd(glue("{flowr_path}/{ind}/tmp"))
    # list.files()

    samp = trk$name[i]
    message(samp)
    gr_logr = to_gr.gatkcnv.cr(trk$gatkcnv_dn_cr[i])
    gr_baf = to_gr.gatkcnv.baf(trk$gatkcnv_acounts[i])
    
    # get allelic counts in our intervals
    gr_baf_logr = join_overlap_intersect(gr_baf, gr_logr)
    length(gr_baf_logr)
    # cont filter on counts, since we need SAME numbers
    # gr_baf_wes$cnts = gr_baf_wes$alt_count+gr_baf_wes$ref_count
    # keep = gr_baf_wes$cnts > 10;table(keep)
    # plot(density(gr_baf_wes$cnts))
    # gr_baf_wes = gr_baf_wes[keep, ]
    # length(gr_baf_wes)
    
    # baf
    df = gr_baf_logr %>% data.frame() %>% as_tibble() %>% 
      mutate(dp = ref_count+alt_count,
            af = alt_count/dp,
            baf = ifelse(af < 0.5,  1-af, af)) %>% 
      mutate(region = glue("{seqnames}-{start}")) %>% 
      dplyr::select(region, chrs = seqnames, pos = start, af)
    # use AF and not BAF, plots look weird
    mat_baf = df %>% wranglr::to_mat();
    colnames(mat_baf)[3] = samp
    head(mat_baf);dim(mat_baf)
    write.table(mat_baf, trk$ascat_baf[i], sep = "\t", quote = F)
    
    
    # logr
    # gr_logr_wes = join_overlap_left(gr_baf_wes, gr_logr)
    # length(gr_logr_wes)
    df = gr_baf_logr %>% data.frame() %>% as_tibble() %>% 
      mutate(region = glue("{seqnames}-{start}")) %>% 
      dplyr::select(region, chrs = seqnames, pos = start, log2_copy_ratio)
    mat_logr = df %>% wranglr::to_mat()
    colnames(mat_logr)[3] = samp
    head(mat_logr);dim(mat_logr)
    write.table(mat_logr, trk$ascat_logr[i], sep = "\t", quote = F)
  }, mc.cores = 6)

}

run_ascat <- function(){
  # add ascat files (need to move it out)
  trk %<>% mutate(
    ascat_logr = glue("ascat/{name}_logr.tsv"),
    ascat_baf = glue("ascat/{name}_baf.tsv"),
    ascat_in = glue("ascat/{name}_input.rds"),
    ascat_out = glue("ascat/{name}_output.rds")
  )

}

ascat_pipe_r <- function(){
  trk = read_tsv(glue("{flowr_path}/trk.tsv")) %>% 
    filter(individual == "185_003")
  p_load(ASCAT)
  p_load(glue, tidyverse, janitor, magrittr)
  
  # setwd("~/projects_git/ascat/ExampleData")
  # ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt","Germline_LogR.txt","Germline_BAF.txt")
  # dim(ascat.bc$Tumor_BAF)
  
  # read all logR data
  # read gatk rds
  # rds = read_rds("gatkcnv_v4/seg_results.rds")
  # trk_gtk_fl = glue("{wexpath}/gatkcnv_v4/df_trk_gtkcnv.tsv")
  # trk_gtk = read_tsv(trk_gtk_fl)
    
  # logR, use denoised
  trk_tum = filter(trk, normal == "NO")
  trk_norm = filter(trk, normal == "YES")
  
  ascat_prep_inputs(trk)

  # setwd("/rsrch3/scratch/iacs/sseth/flows/SS/sarco/mda/wex/ascat")
  as.c = as.character
  lapply(1:nrow(trk_tum), function(i){
    
    message(i)
    out = ascat.loadData(as.c(trk_tum[i, "ascat_logr"]),
                         as.c(trk_tum[i, "ascat_baf"]),
                         as.c(trk_norm[1, "ascat_logr"]),
                         as.c(trk_norm[1, "ascat_baf"]),
                         chrs = c(1:22,"X"), 
                         gender = "XX")
    # trk_tum$ascat_gender[i])
    
    # plotting raw data
    ascat.plotRawData(out, img.dir = "ascat")
    # default penalty 25
    out = ascat.aspcf(out, penalty = 100, out.dir = "ascat")
    ascat.plotSegmentedData(out, img.dir = "ascat")
    out2 = ascat.runAscat(out, gamma = 1, circos = "ascat/circos", img.dir = "ascat") #
    names(out);length(out)
    names(out2);length(out2)

    write_rds(out, as.c(trk_tum[i, "ascat_in"]))
    write_rds(out2, as.c(trk_tum[i, "ascat_out"]))
    
  })
  
  
}







