# example -----------------

if(FALSE){
  library(pacman)
  p_load(tidyverse)
  p_load(glue, janitor)
  # detach("package:flowr")
  # detach("package:params")
  remotes::install_local("~/.drop/public/github_params", force = TRUE)
  remotes::install_local("~/.drop/public/github_flowr", force = TRUE)
  p_load(params, flowr)
  p_load(futile.logger)
  p_load(magrittr)
  source("~/.drop/public/flowr/my.ultraseq/my.ultraseq/R/gatkcnv.R")
  opts_flow = new_opts()
  
  bampath = "/rsrch3/home/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/bam"
  opts_flow$load_toml("~/.drop/public/flowr/my.ultraseq/projects/ss_tnbc/ms51/wex/flowr_conf.toml")
  # remotes::install_local("~/.drop/public/github_params")
  purecn_exe = opts_flow$get("purecn_exe")
  # purecndir = "purecn"
  normaldb_fl = opts_flow$get("purecn_normaldb_fl")
  mappingbias_fl = opts_flow$get("purecn_mappingbias_fl")
  snpblacklist_fl = opts_flow$get("purecn_snpblacklist_fl")
  cosmic_fl = opts_flow$get("purecn_cosmic_fl");
  intervals = opts_flow$get("purecn_intervals")
  # bed = "/rsrch3/home/iacs/sseth/ref/az_ref_beds/ss_downloads/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets_nochr_purecnv_optimized.bed"
  purecn_opts = opts_flow$get("purecn_opts")
  cpu = opts_flow$get("purecn_cpu")
  
  # p_load(rtracklayer)
  # df = import(snpblacklist_fl)
  
  # purecn
  runpath = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_noref/runs"
  trk_full = read_tsv(glue("{runpath}/trk_full.tsv"))
  covpath = "/rsrch3/home/iacs/sseth/flows/SS/tnbc/ms51_wex/purecn/coverage"
  trk_full %<>%  
    mutate(
      bamprefix = gsub(".bam", "", bam),
      purecn_cov_fl = glue("{covpath}/{bamprefix}_coverage_loess.txt.gz"),
      vcf_fl = glue("ssm/m1_m2/pcgr/{outprefix_paired}_combvcf.pcgr_acmg.grch37.vcf.gz"))
  # trk_full$purecn_cov_fl
  # trk_full$vcf_fl
  
  source('~/.drop/public/flowr/my.ultraseq/my.ultraseq/R/metadata_dnaseq.R')
  ind = "185_006"
  setwd(glue("{runpath}/{ind}/tmp"))
  wranglr::mkdir("purecn")
  trk = trk_full %>% filter(individual == ind)
  
  # create flowmat:
  p_load(RcppTOML)
  pipe_str <- parseToml("~/.drop/public/flowr/my.ultraseq/projects/ss_tnbc/ms51/wex/ss_ct_het/pipeline_str.toml")
  source("~/.drop/public/flowr/my.ultraseq/my.ultraseq/R/purecn.R")
  lst = purecn(trk, ind)
  system(lst$cmd_purecn[[3]])
  
  # extras -------
  sampleid = df$gatkcnv_prefix_paired %>% basename()
  bam = file.path(bampath, df$bam)
  # tumor_cov = 
  counts_hf5 = gsub("tsv", "hdf5", df$gatkcnv_cnts)
  m2_vcf = glue("{df$mutect2_outprefix_paired}_f2.vcf.gz")
  vcf_fl = 
    seg_fl = df$gatkcnv_final_seg
  logratio_fl = df$gatkcnv_dn_cr
  
  out_prefix = glue("purecn/{sampleid}")
  purecn_vcf = glue("purecn/{sampleid}_purecn.vcf.gz")
  purecn_cov = glue("purecn/{sampleid}_purecn.cov.loess.txt")
  
}



# funcs -------
purecn_prep_filenms <- function(trk, purecndir){
  trk %<>%
    mutate(
      # skip taking basename of bam!
      # bam = basename(bam),
      # oprefix = glue("{outpath}{name}"),
      purecn_oprefix = glue("{purecndir}/{outprefix_paired}"),
      purecn_segl_fl = glue("{purecn_oprefix}_gatkcnv.seg"),
      purecn_vcf_fl = glue("{purecn_oprefix}_purecn.vcf"),
      purecn_rds = glue("{purecn_oprefix}.rds")
    )
  trk
}

purecn_summary <- function(x){
  results = x$results
  
  length(results)
  i=1
  df_summ = lapply(seq_along(results), function(i){
    res = results[[i]]
    contamination <- res$SNV.posterior$posterior.contamination
    contamination <- if (is.null(contamination)) {
      0
    } else {
      contamination
    }
    log.ratio.offset <- mean(res$log.ratio.offset)
    log.ratio.sdev <- x$input$log.ratio.sdev
    num.segments <- nrow(res$seg)
    d.f.curation <- data.frame(
      solution = i,
      Sampleid = res$seg$ID[1], 
      Purity = res$purity, 
      Ploidy = res$ploidy, Sex = PureCN:::.getSexFromRds(x), 
      Contamination = contamination, 
      Flagged = res$flag, Failed = FALSE, Curated = FALSE, 
      Comment = res$flag_comment,
      
      log.ratio.offset = log.ratio.offset,
      log.ratio.sdev = log.ratio.sdev, 
      num.segments = num.segments
    )
    d.f.curation
  }) %>% bind_rows()
  df_summ
}

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

# create the pureCN cli
purecn <- function(trk, 
                   samplename,
                   # funr_exe = opts_flow$envir$R_funr_exe,
                   purecn_exe = opts_flow$envir$purecn_exe,
                   purecndir = pipe_str$purecn$dir,
                   
                   intervals = opts_flow$envir$purecn_intervals,
                   
                   normaldb_fl = opts_flow$envir$purecn_normaldb_fl,
                   mappingbias_fl = opts_flow$envir$purecn_mappingbias_fl,
                   snpblacklist_fl = opts_flow$envir$purecn_snpblacklist_fl,
                   cosmic_fl = opts_flow$envir$purecn_cosmic_fl,
                   
                   purecn_opts = opts_flow$envir$purecn_opts,
                   cpu = opts_flow$envir$purecn_cpu
){
  
  # print function name
  # message(match.call()[1], " ", appendLF = F)
  check_args()
  
  flog.debug(paste0("Generating a purecn_matched flowmat for sample: ", samplename))
  # check generic columns
  trk = metadata_for_dnaseq_tools(trk)
  # we would use
  # WEX-sarco14-T.denoisedCR.tsv instead of WEX-sarco14-T___matched.cr.seg
  # coz titan segments it
  expect_columns(trk, c(
    "outprefix", 
    "outprefix_paired",
    "gatkcnv_final_seg", 
    "gatkcnv_dn_cr",
    "purecn_cov_fl",
    "mutect2_vcf",
    "vcf_fl"
  ))
  
  source("~/.drop/public/flowr/my.ultraseq/my.ultraseq/R/purecn.R")
  # ceate a new trk with ALL reqd files
  flog.debug(paste0("Generating a filenames flowmat for sample"))
  trk = purecn_prep_filenms(trk, purecndir)
  
  trk_tum = filter(trk, normal == "NO")
  trk_norm = filter(trk, normal == "YES")
  
  # ** prep inputs ----------
  i <- 1
  # activate multiple sample segmentation where possible
  # deactivate segmentation
  tmp = lapply(1:nrow(trk_tum), function(i){
    
    # vcf_fl = trk_tum$vcf_fl[i]
    # vcf_fl = "ssm/m1_m2/IPCT-S4012-MOON0051-Cap2515-4-HTID295___pon_combvcf.vcf.gz"
    # vcf_fl = "ssm/m1_m2/IPCT-S4012-MOON0051-Cap2515-4-HTID295___pon_combvcf.vcf.gz"
    vcf_fl = trk_tum$mutect2_vcf[i]
    out_prefix = trk_tum$purecn_oprefix[i]
    sampleid = trk_tum$outprefix[i]
    logratio_fl = trk_tum$gatkcnv_dn_cr[i]
    seg_fl = trk_tum$gatkcnv_final_seg[i]
    cov_fl = trk_tum$purecn_cov_fl[i]
    cov_fls = trk_tum$purecn_cov_fl %>% setdiff(cov_fl)

    # get the cmds
    cmd = glue("Rscript {purecn_exe}/PureCN_v2.R ",
               "--out {out_prefix} ",
               "--sampleid {sampleid} ",
               "--tumor {cov_fl} ",
               "--normaldb {normaldb_fl} ",
               
               # only supported if counts were also from gatk!
               # "--logratiofile {logratio_fl} ",
               
               # "--segfile {seg_fl} ",
               # "--segfile_type gatkcnv ",
               
               "--mappingbiasfile {mappingbias_fl} ",
               "--vcf {vcf_fl} ",
               "--outvcf ", # this is a true/false
               
               "--cosmic.vcf.file {cosmic_fl} ",
               
               "--intervals {intervals} ",
               "--snpblacklist {snpblacklist_fl} ",
               "--parallel --cores {cpu} ",
               "{purecn_opts} ")
    if(length(cov_fls) > 0){
      cov_fls %<>% paste(collapse = ",")
      cmd = glue("{cmd} --additionaltumors {cov_fls}")
    }
    cmd
    # system(cmd)
    cmds <- list(cmd_purecn = cmd,
                 outfiles = list(trk_tum$purecn_rds))
  })
  lst = purrr::transpose(tmp)
  cmds = list(#purecn.prep = cmd_pileup, 
    purecn.mrg = c(unlist(lst$cmd_purecn)))
  # system(tmp[[2]]$cmd_purecn)
  # write_tsv(df_trk, path = "gatkcnv/df_trk.tsv")
  
  # cmds %>% unlist() %>% paste0(collapse = "\n") %>% cat()
  lst$flowmat = to_flowmat(cmds, samplename = samplename) %>%
    mutate(cmd = as.character(cmd))
  lst$trk = trk
  return(lst)
}



# example -------------------
if(FALSE){
  
  opts_flow$load_toml("~/.drop/public/flowr/my.ultraseq/projects/ss_tnbc/ms51/wex/flowr_conf.toml")
  setwd("/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gc/runs/185_023/tmp")
  setwd("/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/runs/185_003/tmp")
  as.c = as.character
  purecndir = pipe_str$purecn$dir
  df_trk_full = read_tsv("/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/runs/trk_full.tsv")
  trk = df_trk_full %>% filter(individual == "185_006") %>% 
    
    trk = purecn_prep_filenms(trk, purecndir) %>% 
    metadata_for_dnaseq_tools()
  trk_tum = filter(trk, normal == "NO")
  trk_norm = filter(trk, normal == "YES")
  i=1
  opt = list(
    # "purecn/IPCT-S4037-MOON0051-Cap2620-4-NT05_200720-A00422-0252-BHG3VMDSXY-2-ACCACTGT.bwa_recalibed_coverage_loess.txt.gz",
    tumor = as.c(trk_tum$purecn_cov_fl[i]),
    vcf  = as.c(trk_tum$mutect2_vcf[i]),
    out = as.c(trk_tum$purecn_oprefix[i]),
    # "purecn2/IPCT-S4013-MOON0051-Cap2529-4-HTID206___matched",
    # "IPCT-S4013-MOON0051-Cap2529-4-HTID206",
    sampleid = trk_tum$outprefix[i],
    
    # segfile = trk_tum$gatkcnv_final_seg[i],
    segfile = NULL,
    segfile_type = NULL,
    logratiofile = NULL,
    additionaltumors = NULL,
    
    normaldb  = opts_flow$envir$purecn_normaldb_fl,
    mappingbiasfile = opts_flow$envir$purecn_mappingbias_fl,
    snpblacklist = opts_flow$envir$purecn_snpblacklist_fl,
    intervals = opts_flow$envir$purecn_intervals,
    
    # vcf
    minaf = 0.03,
    dbinfoflag = formals(PureCN::runAbsoluteCN)$DB.info.flag,
    padding = formals(PureCN::filterVcfBasic)$interval.padding,
    
    # cosmic
    # popafinfofield = formals(PureCN::runAbsoluteCN)$POPAF.info.field,
    # cosmic.vcf.file = opts_flow$envir$purecn_cosmic_fl,
    cosmic.vcf.file = NULL,
    popafinfofield = "POPAF",
    mincosmiccnt = formals(PureCN::setPriorVcf)$min.cosmic.cnt,
    
    # purity
    minpurity = formals(PureCN::runAbsoluteCN)$test.purity[[2]],
    maxpurity = formals(PureCN::runAbsoluteCN)$test.purity[[3]],
    minploidy = formals(PureCN::runAbsoluteCN)$min.ploidy,
    maxploidy = formals(PureCN::runAbsoluteCN)$max.ploidy,
    
    # formals(PureCN::runAbsoluteCN)$test.num.copy
    maxcopynumber = max(eval(formals(PureCN::runAbsoluteCN)$test.num.copy)),
    maxnonclonal = formals(PureCN::runAbsoluteCN)$max.non.clonal,
    
    # loss
    modelhomozygousloss = FALSE,
    maxhomozygousloss = paste(formals(PureCN::runAbsoluteCN)$max.homozygous.loss[2:3], collapse=","),
    
    # model
    model = formals(PureCN::runAbsoluteCN)$model[[2]],
    
    # params
    alpha = formals(PureCN::segmentationCBS)$alpha,
    undosd = formals(PureCN::segmentationCBS)$undo.SD,
    maxsegments = formals(PureCN::runAbsoluteCN)$max.segments,
    intervalweightfile = NULL,
    
    outvcf = TRUE, 
    parallel = TRUE, 
    maxcandidatesolutions = 20,
    cores = 24,
    postoptimize = TRUE,
    seed = 123,
    genome = "hg19",
    
    # calc power:
    error = formals(PureCN::runAbsoluteCN)$error,
    
    # speedup
    speedupheuristics = max(eval(formals(PureCN::runAbsoluteCN)$speedup.heuristics)),
    logratiocalibration = formals(PureCN::runAbsoluteCN)$log.ratio.calibration,
    
    funsegmentation = "PSCBS",
    centromere_seq_style = "NCBI")
  opt
}



# read purecn -------
to_gr_seg.purecn <- function(segfl,
                            variantfl = NULL,
                            calc_minor_major,
                             gencode_fl = "~/.rsrch3/home/iacs/sseth/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene_hgnc.bed",
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
                             gencode_fl = "~/.rsrch3/home/iacs/sseth/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene_hgnc.bed",
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
  
  source('~/.drop/projects/packs_dnaseq/R/to_mae.seg.R')
  df_ann = annotate.gr_seg(lst$gr_seg, lst_gen = lst_gen)
  dim(df_ann)
  # df_ann = df_ann %>% group_by(gene_name) %>% add_count()
  # df_ann %>% filter(n>1) %>% View()
  
  list(df_ann = df_ann, df_igv_seg = lst$df_igv_seg, gr_seg = lst$gr_seg)
  
}



# mutfl = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/purecn_v2/185_011/IPCT-S4012-MOON0051-Cap2516-4-HTID343___matched_variants.csv"
to_mut.purecn <- function(mutfl,
                          purity = NULL, ploidy = NULL,
                          verbose = FALSE){
  p_load(tidyverse)
  df = readr::read_csv(mutfl) %>% janitor::clean_names()
  df
}

# 011 has a frameshift deletion which is NOT included!
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




# purecn outputs -------
# [1] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_amplification_pvalues.csv"
# [2] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_chromosomes.pdf"          
# [3] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_debug.csv"                
# [4] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_dnacopy.seg"              
# [5] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_full.csv"                 
# [6] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_genes.csv"                
# [7] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_local_optima.pdf"         
# [8] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_loh.csv"                  
# [9] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_multisample.seg"          
# [10] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_segmentation.pdf"         
# [11] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched_variants.csv"             
# [12] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched.csv"                      
# [13] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched.log"                      
# [14] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched.pdf"                      
# [15] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched.rds"                      
# [16] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched.vcf"                      
# [17] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched.vcf.gz"                   
# [18] "IPCT-S4013-MOON0051-Cap2531-4-HTID278___matched.vcf.gz.tbi"               
