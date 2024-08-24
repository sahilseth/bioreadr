

sclust_read_seg <- function(x){
  df = read_tsv(x, 
                col_names = c("sample_id", "contig", "start", "end", "num_points", "log2_copy_ratio"), 
                col_types = cols(
                  sample_id = col_character(),
                  contig = col_character(),
                  start = col_double(),
                  end = col_double(),
                  num_points = col_double(),
                  log2_copy_ratio = col_double()))
}

# SCLUST has several types of file
# we will read all of them, and create a single object with all the information
# we will attempt to use a similar format as pureCN

# logical (e.g., TRUE, FALSE)
# integer (e.g,, 2L, as.integer(3))
# numeric (real or decimal) (e.g, 2, 2.0, pi)
# complex (e.g, 1 + 0i, 1 + 4i)
# character (e.g, "a", "swc")

# https://www.nature.com/articles/nprot.2018.033#Sec19
# this protocol includes all the information necessary
# to read sclust properly
slust_read <- function(sclust_cn_summary, sclust_allelic_states, sclust_subclonal_cn,
                       sclust_icn, sclust_clusters, sclust_mclusters, 
                       # sclust_corr_seg, 
                       sclust_uncorr_seg,
                       bms_sample_id){
  
  p_load(futile.logger)
  # this is summary of the run
  # sample_name      purity ploidy 
  # fraction_subclonal_cn: fraction of subclonal CN
  # sex_estimated: sex of the patient from CNV
  # status: optimum, invariant (no sufficient CNV changes; purity estimated from mutations), forced, OR failed
  # fraction_inconsistent_segs: if no unique CNV solution was found
  flog.info("read summary")
  df_summ = read_tsv(sclust_cn_summary, col_types = "cnnnccn") %>% 
    mutate(bms_sample_id = bms_sample_id)
  df_summ
  
  # we are assuming this is PER sample
  # fair assumption
  sample_name = df_summ$sample_name[1]
  futile.logger::flog.info(glue("working with {sample_name}"))
  
  
  # Sample Chromosome Start End 
  # Copy_Nr_Raw: uncorrected CNV
  # CopyNr: corrected integer CNV
  # A B: Major and minor allele CNV
  # LOH: 1/0 in case of LOH
  # Theta Theta_Exp: observed and model predicted allelic imbalance
  # n_SNPs: total number of SNPs
  # Is_Subclonal_CN: 1/0 subclonal? Subclonal_P_value: confidnece
  # Is_Inconsistent_State: 1 - no unique solution was found
  # this has HALF the segments as iCN and uncorrected
  flog.info("read allelic states")
  df_allelic_states = read_tsv(sclust_allelic_states, col_types = "cciiniiiinniiniini") %>% 
    mutate(bms_sample_id = bms_sample_id)
  
  # list of all subclonal CNV events
  # Sample Chromosome Start End 
  # Subclonal_CN:  total CNV
  # Clone1_A Clone1_B Clone1_Fraction; major minor in clone 1; estimated fraction of clone 1
  # Clone2_A Clone2_B Clone2_Fraction; major minor in clone 1; estimated fraction of clone 1
  # if we have more clones, this needs to be dynamic
  flog.info("read subclonal CN")
  # number of clusters would be different each time, cant specify the columns here
  df_subclonal_cn = read_tsv(sclust_subclonal_cn) %>% 
    mutate(bms_sample_id = bms_sample_id)
  
  # format good for IGV; with corrected CNV
  flog.info("read ICN")
  df_icn = read_tsv(sclust_icn, 
                    col_names = c("sample_name", "chr", "start", "stop", "probes", "cnv_corrected"), 
                    col_types = cols(
                      sample_name = col_character(),
                      chr = col_character(),
                      start = col_integer(),
                      stop = col_integer(),
                      probes = col_integer(),
                      cnv_corrected = col_double()
                    )) %>% 
    mutate(bms_sample_id = bms_sample_id)
  table(df_icn$cnv_corrected)
  
  # <sample>_ muts_expAF.txt
  # file missing
  # this is input mutations for clustering
  
  # <sample>_cn_profile.pdf
  # pg1: optimization process
  # For each ploidy (x axis); y: estimate of purity and ploidy estimates
  # lower panel: object func based in read ratio (L-CN); and biallelic freq (L-bi)
  #   red bar: optimal result by global min of L-cn within scanning range (b/w minp and maxp)
  # 
  # pg 2: 
  # observed allelic imbalance (theta)
  # in dependence of total CN
  # allelic states: depicted by different colors
  # 
  # pg3
  
  # cluster assignments
  # columns 1:3; same as muts_expAF:
  # CCF: raw cancer cell fraction and its coverage
  # cluster ID; and CCF of the cluster, and probability
  # P0-PXX: probability of belonging to all other clusters (these sum to 1)
  flog.info("read clusters")
  df_clusters = data.table::fread(sclust_clusters, fill = TRUE, data.table = FALSE) %>% janitor::clean_names() %>% as_tibble() %>% 
    mutate(bms_sample_id = bms_sample_id)
  #   Detected 11 column names but the data has 10 columns. Filling rows automatically. Set fill=TRUE explicitly to avoid this warning.
  # df_clusters = read_tsv(sclust_clusters) %>% clean_names()
  # problems(df_clusters)
  
  # summry of the clusters
  # cluster ID, its CCF, number of mutations in the cluster
  # 
  flog.info("read mcluster")
  df_mclusters = read_tsv(sclust_mclusters, 
                          col_types = cols(
                            Cluster_ID = col_double(),
                            CCF_Cluster = col_double(),
                            Cluster_Peak_Height = col_double(),
                            Mutations_In_Cluster = col_double())) %>% 
    mutate(bms_sample_id = bms_sample_id)
  
  # add sample ids
  # df_mclusters %>% mutate()
  
  # flog.info("read corrected seg")
  # df_corr = sclust_read_seg(sclust_corr_seg)
  flog.info("read uncorrected seg")
  if(is.na(sclust_uncorr_seg)){
    df_uncorr = data.frame()
  }else{
    df_uncorr = sclust_read_seg(sclust_uncorr_seg) %>% 
      mutate(bms_sample_id = bms_sample_id)
  }

  # df_uncorr = sclust_read_seg(sclust_uncorr_seg)
  
  list(df_summ = df_summ, 
       df_allelic_states = df_allelic_states,
       df_subclonal_cn = df_subclonal_cn, 
       df_icn = df_icn, 
       df_clusters = df_clusters, 
       df_mclusters = df_mclusters, 
       # df_corr = df_corr, 
       df_uncorr = df_uncorr)
}
















# END