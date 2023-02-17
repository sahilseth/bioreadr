# job::job({p_install(SRAdb)})

# generic --------
as_n <- function(...) {
  as.numeric(gsub("\\.", NA, unlist(...)))
}
as_i <- function(...) {
  as.integer(gsub("\\.", NA, unlist(...)))
}
as_c <- function(...) {
  as.character(gsub("\\.", NA, unlist(...)))
}

# transpose, if header info is NOT a matrix
t2 <- function(x){
  if(class(x)[1] == "matrix")
    x
  else
    t(x)
}

# split vals --------
splt_vcf_format <- function(x, format, prefix){
  x = as.character(unlist(x))
  format = as.character(unlist(format))
  
  lst <- lapply(1:length(x), function(i){
    flog.debug(i)
    xi = x[i];formati = format[i]
    splt = strsplit(xi, ":")[[1]]
    nms = tolower(strsplit(formati, ":")[[1]])
    len_splt = length(splt)
    len_nms = length(nms)
    
    # for muse, NORMAL has less fields, assume last one is missing
    # last var in TUMOR is variant type!
    flog.debug(glue("len_splt: {len_splt}, len_nms: {len_nms}"))
    if(len_splt < len_nms)
      names(splt) <- paste(prefix, nms[1:len_splt], sep = "")
    else
      names(splt) = paste(prefix, nms, sep = "")

    ret = as.data.frame(t(splt), stringsAsFactors = FALSE)
    return(ret)
  })
  mat = bind_rows(lst)
}


#' splt_vcf_info
#'
#' @param x a character vector
#' @param cores something
splt_vcf_info <- function(x, cores = 1){
  x = as.character(unlist(x))
  # for every record
  # split variables using ;
  # extract variable and value (separated by =)
  # 
  #profvis({
  lst <- mclapply(1:length(x), function(i){
    if(i %% 10^5 == 0) message(i)
    xi = x[i];
    splt = strsplit(xi, ";")[[1]]
    splt = gsub("(.?)\\=(.*)", "\\1\n\n\\2", splt)
    splt = strsplit(splt, "\n\n")
    vals = sapply(splt, tail, 1)
    #vals = sapply(splt, "[[", 2) # improves performance
    names(vals) = sapply(splt, "[[", 1)
    ret = as.data.frame(t(vals), stringsAsFactors = FALSE)
    return(ret)
  }, mc.cores = cores)
  #})
  mat = bind_rows(lst)
  colnames(mat) = tolower(colnames(mat))
  mat
}


split_vep_consequence <- function(df, 
                                  col,
                                  pattern = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID"){
  col_names = strsplit(pattern, "\\|")
  # select the first one (based on what??)
  
  
}
# funcs for specific fields -------
.get_value_type <- function(type){
  
  if(type == "Float")
    return(as_n)
  if(type == "Integer")
    return(as_i)
  if(type == "String")
    return(as_c)
  
}




.get_sample_names_mutect2 <- function(header){
  
  
}

# main funcs -------
#' Parse a VCF file into a TSV to be processed and annotated
#'
#' @param x a vcf file to be parsed
#' @param cores something
#'
#' @export
read_vcf <- function(x, cores = 1) {
  require(dplyr)
  require(parallel)
  require(futile.logger)
  require(bedr)
  require(futile.logger)
  require(bedr)

  # check if its a bedr object
  flog.info("Reading file...")
  if (is.list(x)) {
    if (length(x$header) > 1) {
      flog.info("using provided  bedr::vcf object")
      vcf <- x
    } else {
      stop("you supplied a list, however it is not a bedr object")
    }
  } else {
    # use bedr (trusting them :)
    vcf <- bedr::read.vcf(x)
  }
  # if(tools::file_ext(x) == "gz"){
  #   fl_con = gzfile(x)
  # }else{
  #   fl_con = file(x)
  # }
  #
  # # assuming that header is not longer than 100
  # hd = scan(fl_con, what = "character", sep = "\n", n = 1000, quiet = TRUE)
  # hd_row = grep("#CHROM", hd)
  # hd = tolower(strsplit(hd[hd_row], "\t")[[1]])
  # col_classes = paste(rep("c", length(hd)), collapse = "")
  #
  # #tab = data.table:::fread(x, data.table = FALSE, header = FALSE, skip = hd_row, sep = "\t", colClasses = "character")
  # tab = readr::read_tsv(x, col_names = FALSE, skip = hd_row, col_types = col_classes)
  # colnames(tab) = gsub("#chrom", "chrom", hd)
  # #tab = as_tibble(tab)

  flog.info("Parsing info columns...")
  # debug(splt_vcf_info)
  info_cols <- splt_vcf_info(vcf$vcf$INFO, cores = cores)
  colnames(info_cols) <- paste0("info_", colnames(info_cols))

  # total cols
  flog.info("Parsing format columns...")
  cols <- colnames(vcf$vcf)
  cols_samp <- seq(grep("FORMAT", cols) + 1, ncol(vcf$vcf))
  cols_samp <- cols[cols_samp]
  if (length(cols_samp) > 0) {
    message(" found ", length(cols_samp), " samples")
    # samp="334187-N"
    lst_tmp <- lapply(cols_samp, function(samp) {
      message(samp)
      format_cols <- NA
      format_cols <- splt_vcf_format(x = vcf$vcf[, ..samp], format = vcf$vcf$FORMAT, prefix = "")
      cols_format <- colnames(format_cols)
      # add fmt: to be sure where this column is coming from
      colnames(format_cols) <- glue::glue("{samp}_fmt_{cols_format}")
      format_cols
    })
    format_cols <- do.call(cbind, lst_tmp)
    # head(format_cols)
  } else {
    flog.info("no format columns")
  }
  # select vcf fixed
  colfixed <- !colnames(vcf$vcf) %in% c("FORMAT", cols_samp, "INFO")

  flog.info("Assembling data")
  tab2 <- as_tibble(cbind(vcf$vcf[, ..colfixed], info_cols, format_cols))
  head(tab2)
  colnames(tab2) <- tolower(colnames(tab2))
  return(list(vcf = vcf$vcf, tab = tab2, header = vcf$header))
}
parse_vcf <- read_vcf


#' Parse a somatic VCF, with two samples.
#'
#' @param x a vcf file
#' @param samp name of the 'tumor' sample
#' @param ref name of the 'reference' sample
#'
#' @export
read_vcf_somatic <- function(x, .vcf = NULL, samp = NULL, ref = NULL){
  
  
  # flog.info("Reading file...")
  # add a read_vcf function
  if(is.null(.vcf))
    vcf = read_vcf(x)
  else
    vcf = .vcf
  mat = vcf$tab
  head(mat)
  
  flog.info("get sample names...")
  # check if this is mutect
  if(length(vcf$header$tumor_sample) > 0 & is.null(samp)){
    flog.info("picking tumor_sample from vcf")
    samp = vcf$header$tumor_sample
    pair_type = "TO"
    flog.info(samp)
  }
  if(length(vcf$header$normal_sample) > 0 & is.null(ref)){
    flog.info("picking normal_sample from vcf")
    ref = vcf$header$normal_sample
    pair_type = "TN"
    flog.info(ref)
  }
  # check_args()
  
  colnms = colnames(mat)
  
  # ISSUE: 
  # replace all occurance of {tumor_name} with t
  # in case tumor name is a subset of the normal name, we will encounter a problem!
  # lets be more specific of what columns need to be changed
  colnms_new = gsub(glue("{tolower(samp)}_fmt"), "t_fmt", colnms)
  if(!is.null(ref))
    colnms_new = gsub(glue("{tolower(ref)}_fmt"), "n_fmt", colnms_new)
    # colnms_new %<>% gsub(tolower(ref), "n", .)
  colnames(mat) = colnms_new
  head(mat)
  
  flog.info("INFO: auto force column info type ...")
  # work on integer or float columns, which have a single value
  # no separators
  # - there are a few cases where we have "." are supposed to be translated as NA
  vcf_header_info = vcf$header$INFO %>% t2() %>% 
    data.frame(stringsAsFactors = F) %>% 
    dplyr::mutate(colnm = paste0("info_", tolower(ID)))
  vcf_header_info_f1 = tidylog::filter(vcf_header_info, 
                                       Number %in% c(1, "A"), 
                                       Type %in% c("Integer", "Float"), 
                                       colnm %in% colnames(mat))
  if(length(vcf_header_info_f1$colnm) > 0){
    flog.info(paste0("can auto-force following:\n", 
            paste0("\t", vcf_header_info_f1$colnm, collapse = "\n")))
    for(i in 1:nrow(vcf_header_info_f1)){
      df = vcf_header_info_f1[i, ]
      colfunc = .get_value_type(df$Type)
      colnm = df$colnm
      # mat[, colnm] = colfunc(mat[, colnm])
      mat %<>% dplyr::dplyr::mutate_at(.vars = .vars = colnm, .funs = .funs = colfunc)

    }
  }

  flog.info("auto split column type R ...")
  vcf_header_type_info_r <- vcf_header_info %>% filter(Number == "R")
  if (nrow(vcf_header_type_info_r) > 0) {
    flog.info(paste0(
      "can split following ALLELE related vars:\n",
      paste0("\t", vcf_header_type_info_r$colnm, collapse = "\n")
    ))
    for (i in 1:nrow(vcf_header_type_info_r)) {
      df <- vcf_header_type_info_r[i, ]
      colfunc <- .get_value_type(df$Type)

      # support upto three records
      colnms_new <- paste0(df$colnm, c("_ref", "_alt", "_alt2"))
      # drop multi-allelic records
      mat <- tidyr::separate(mat, col = df$colnm, colnms_new, sep = ",", extra = "drop", fill = "right", remove = F)
      # mat[, colnms_new] <- colfunc(mat[, colnms_new])
      mat %<>% dplyr::dplyr::mutate_at(.vars = .vars = colnms_new, .funs = .funs = colfunc)

    }
  }

  flog.info("FMT: auto force column fmt type ...")
  vcf_header_fmt = vcf$header$FORMAT %>% data.frame(stringsAsFactors = F) %>% 
    mutate(colnm_t = paste0("t_fmt_", "", tolower(ID)),
           colnm_n = paste0("n_fmt_", "", tolower(ID)))
  vcf_header_fmt_f1 = tidylog::filter(vcf_header_fmt, 
                                      Number %in% c(1, "A"), 
                                      Type %in% c("Integer", "Float"), 
                                      colnm_t %in% colnames(mat),
                                      colnm_n %in% colnames(mat))
  flog.info(paste0("can auto-force following:\n", nrow(vcf_header_fmt_f1),
          paste0("\t", vcf_header_fmt_f1$colnm_t, " ", vcf_header_fmt_f1$colnm_n, collapse = "\n")))
  if(nrow(vcf_header_fmt_f1) > 0){
    for(i in 1:nrow(vcf_header_fmt_f1)){
      df = vcf_header_fmt_f1[i, ]
      colfunc = .get_value_type(df$Type)
      # mat[, df$colnm_t] = colfunc(mat[, df$colnm_t])
      # mat[, df$colnm_n] = colfunc(mat[, df$colnm_n])
      mat %<>% mutate_at(.vars = .vars = c(df$colnm_t, df$colnm_n), .funs = .funs = colfunc)
    }
  }
  
  vcf_header_fmt_r = vcf_header_fmt %>% filter(Number == "R")
  if( nrow(vcf_header_fmt_r) > 0 ){
    flog.info(paste0("can split following ALLELE related vars:\n", 
            paste0("\t", vcf_header_fmt_r$colnm_t, collapse = "\n")))
    for(i in 1:nrow(vcf_header_fmt_r)){
      df = vcf_header_fmt_r[i, ]
      colfunc = .get_value_type(df$Type)
      
      # support upto three records
      colnms_t_new = paste0(df$colnm_t, c("_ref", "_alt", "_alt2"))
      mat = tidyr::separate(mat, col = df$colnm_t, colnms_t_new, sep = ",", extra = "drop", fill = "right", remove = F)
      # change col type:
      mat %<>% mutate_at(.vars = .vars = colnms_t_new, .funs = .funs = colfunc)            

      # repeat for normal sample
      if(!is.null(ref)){
        colnms_n_new = paste0(df$colnm_n, c("_ref", "_alt", "_alt2"))
        mat = tidyr::separate(mat, col = df$colnm_n, into = colnms_n_new, sep = ",", extra = "drop", fill = "right", remove = F)
        mat %<>% mutate_at(.vars = .vars = colnms_n_new, .funs = .funs = colfunc)
      }
      
      # mat[, colnms_t_new] = colfunc(mat[, colnms_t_new])
      # mat[, colnms_n_new] = colfunc(mat[, colnms_n_new])
    }
  }
  
  mat %<>% mutate(key = glue("{chrom}:{pos}_{ref}/{alt}"))
  
  head(mat)
  
  # cols_info_int =  %>% 
  #   filter(Number == 1, Type == "Integer")
  
  # mat = as_tibble(as.data.frame(mat,  stringsAsFactors = FALSE))
  # colnames(mat) = tolower(colnames(mat))
  # samplenames_vcf = colnames(mat)[grep("format", colnames(mat)) + 1:2]
  # samplenames_input = c(samp, ref)
  # message("found ", length(samplenames_vcf), " sample(s):\n", paste0(samplenames_vcf, collapse = "\n"))
  
  
  # assume first column
  # samp = colnames(mat)[grep("format", colnames(mat)) + 1]
  # ref = colnames(mat)[grep("format", colnames(mat)) + 2]
  
  # flog.info("Parsing format columns...")
  # tcols = splt_vcf_format(x = mat[,samp], format = mat$format, prefix = "t_")
  # ncols = splt_vcf_format(x = mat[,ref], format = mat$format, prefix = "n_")
  
  # flog.info("Parsing info columns...")
  # infocols = splt_vcf_info(mat$info)
  
  # not sure what this is
  # flog.info("Parsing func...")
  #funccols = splt_vcf_func(mat$func)
  #infocols[1:100,] %>% View()
  # mat$t_sample = samp
  # mat$n_sample = ref
  # colsel = !colnames(mat) %in% c("format", samp, ref, "info")
  
  # flog.info("Assembling data")
  #mat = as_tibble(cbind(mat[, colsel], tcols, ncols, infocols, funccols))
  # mat = as_tibble(cbind(mat[, colsel], tcols, ncols, infocols))
  return(list(vcf = vcf$vcf, tab = mat, header = vcf$header, samp = samp, ref = ref, pair_type = pair_type))
  
}

calc_taf_fwd_rev <- function(x){
  x %>% mutate(
    taf_fwd = t_fmt_f1r2_alt / (t_fmt_f1r2_alt + t_fmt_f1r2_ref),
    taf_rev = t_fmt_f2r1_alt / (t_fmt_f2r1_alt + t_fmt_f2r1_ref),
    taf_diff = abs(taf_fwd - taf_rev))

}

create_maf_key_ins_del <- function(x){
  x %>% 
    mutate(
      alt_maf = sub('.', '', alt),
      pos_maf = pos+1,
      ref_maf = sub('.', '', ref),
      variant_type = case_when(
        nchar(alt) > nchar(ref) ~ "INS",
        nchar(alt) < nchar(ref) ~ "DEL"),
      key_maf = case_when(
        variant_type == "INS" ~ glue("{chrom}:{pos}_-/{alt_maf}"),
        variant_type == "DEL" ~ glue("{chrom}:{pos_maf}_{ref_maf}/-"),
        TRUE ~ glue("{chrom}:{pos}_{ref}/{alt}")))
}
calc_taf_fwd_rev <- function(x){
  x %>% mutate(
    taf_fwd = t_fmt_f1r2_alt / (t_fmt_f1r2_alt + t_fmt_f1r2_ref),
    taf_rev = t_fmt_f2r1_alt / (t_fmt_f2r1_alt + t_fmt_f2r1_ref),
    taf_diff = abs(taf_fwd - taf_rev))

}

#' Parse a somatic VCF, with two samples.
#'
#' @param x a vcf file
#' @param samp name of the 'tumor' sample
#' @param ref name of the 'reference' sample
#'
#' @export
read_vcf_germline <- function(x, samp = NULL){
  
  # flog.info("Reading file...")
  # add a read_vcf function
  vcf = read_vcf(x)
  mat = vcf$tab
  
  if(is.null(samp)){
    samp = vcf$header$tumor_sample
    flog.info(glue("samp missing, extracting from vcf: {samp}"))
  }
  check_args()
  
  colnms = colnames(mat)
  colnms_new = gsub(paste0(tolower(samp), "_"), "", colnms)
  colnames(mat) = colnms_new
  head(mat)
  
  flog.info(">> VCF INFO column: auto force ...")
  vcf_header_info = vcf$header$INFO %>% data.frame(stringsAsFactors = F) %>% 
    mutate(colnm = paste0("info_", tolower(ID)))
  vcf_header_info_f1 = tidylog::filter(vcf_header_info, 
                                       Number %in% c(1, "A"), 
                                       Type %in% c("Integer", "Float"), 
                                       colnm %in% colnames(mat))
  flog.info(paste0("can auto-force following columns:\n", 
          paste0("\t", vcf_header_info_f1$colnm, collapse = "\n")))
  for(i in 1:nrow(vcf_header_info_f1)){
    df = vcf_header_info_f1[i, ]
    colfunc = .get_value_type(df$Type)
    colnm = df$colnm
    mat[, colnm] = colfunc(mat[, colnm])
  }
  
  flog.info(">> VCF FMT: auto force column fmt type ...")
  vcf_header_fmt = vcf$header$FORMAT %>% data.frame(stringsAsFactors = F) %>% 
    mutate(colnm = paste0("fmt_", "", tolower(ID)))
  vcf_header_fmt_f1 = tidylog::filter(vcf_header_fmt, 
                                      Number %in% c(1, "A"), 
                                      Type %in% c("Integer", "Float"), 
                                      colnm %in% colnames(mat))
  flog.info(paste0("can auto-force following:\n", 
            paste0("\t", vcf_header_fmt_f1$colnm, collapse = "\n")))
  for(i in 1:nrow(vcf_header_fmt_f1)){
    df = vcf_header_fmt_f1[i, ]
    colfunc = .get_value_type(df$Type)
    mat[, df$colnm] = colfunc(mat[, df$colnm])
  }
  
  flog.info(">> VCF split alleles: auto split column type R ...")
  vcf_header_type_r = bind_rows(vcf_header_info, 
                                vcf_header_fmt) %>% 
    tidylog::filter(Number == "R")
  flog.info(paste0("can split following ALLELE related vars:\n", 
            paste0("\t", vcf_header_type_r$colnm, collapse = "\n")))
  for(i in 1:nrow(vcf_header_type_r)){
    df = vcf_header_type_r[i, ]
    flog.info(df$colnm)
    colfunc = .get_value_type(df$Type)
    # support upto three records
    colnms_new = paste0(df$colnm, c("_ref", "_alt", "_alt2"))
    # drop multi-allelic records
    mat = tidyr::separate(mat, col = df$colnm, colnms_new, sep = ",", extra = "drop", fill = "right", remove = F)
    # convert value type them
    mat[, colnms_new] = apply(mat[, colnms_new], 2, colfunc)
  }
  mat %<>% mutate(key = glue("{chrom}:{pos}_{ref}/{alt}"))
  
  head(mat)
  return(list(vcf = vcf$vcf, tab = mat, header = vcf$header))
}

# example -------
if(FALSE){
  x ="/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/ponm/vcfs/WEX-1004-N.vcf.gz"
  
  library(pacman)
  p_load(dplyr, janitor, glue, magrittr)
  
  source('~/Dropbox/public/flow-r/ultraseq/ultraseq/R/parse_vcfs.R')
  x='/rsrch2/iacs/iacs_dep/sseth/flowr/runs/flowr_test/fastq_haplotyper-MS132-20150824-16-37-58-XScJT0OZ/tmp/GLizee-Pancreatic-MS132-MP013normalDNA.recalibed_1.haplotyper.vcf'
  x2 = parse_vcf(x)
  
  x='~/rsrch1_data/runs/sarcomatoid/jco/mutect/sarco10-T_1208XX_ST1374_073_H09WGADXX___sarco10-N_1208XX_ST1374_073_H09WGADXX_filt2_onco.vcf'
  x2 = parse_somatic_vcf(x)
  
  x = '/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/ssm/m1_m2_ir/WEX-334187-T___matched_combvcf.vcf.gz'
  debug(parse_somatic_vcf);debug(parse_vcf)
  df = parse_somatic_vcf(x, '334187-T', '334187-N')
  
  # ** m2: sarco
  df_ion = read_rds("/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/ionreporter/Exome_Sarcomatoid_334187-T.rds")
  fl_m2 = "/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/mutect2/WEX-334187-T___matched_f1_cleannms_snp.vcf"
  rds_m2 = "/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/mutect2/WEX-334187-T___matched_f1_cleannms_snp.rds"
  # debug(read_vcf_somatic);#debug(parse_vcf)
  vcf = read_vcf_somatic(fl_m2, 
                         samp = '334187-T', 
                         ref = '334187-N')
  write_rds(vcf, rds_m2)
  
  
  # get mutations with strand bias and high AF:
  vcf$tab %<>% mutate(
    taf_fwd = t_fmt_f1r2_alt / (t_fmt_f1r2_alt + t_fmt_f1r2_ref),
    taf_rev = t_fmt_f2r1_alt / (t_fmt_f2r1_alt + t_fmt_f2r1_ref),
    taf_diff = abs(taf_fwd - taf_rev),
    key = glue("chr{chrom}:{pos}:{ref}-{alt}"),
  )
  # test and install
  # get IGV snapshots of a few of these
  df_igv = arrange(vcf$tab, -taf_diff) %>% 
    mutate(
      taf_diff = round(taf_diff, 2),
      name = glue("{chrom}_{pos}_{ref}_{alt}_{taf_diff}_{filter}"), 
      end = pos) %>% 
    dplyr::rename(chr = chrom, start = pos) %>% 
    select(chr:alt, filter, taf_fwd, taf_rev, taf_diff, everything()) %>% 
    # apply usual filters
    tidylog::filter(t_fmt_dp > 20, n_fmt_dp > 20, 
                    t_fmt_af > 0.05, filter == "PASS") %>% 
    left_join(dplyr::select(df_ion, key, ion_taf_fwd=taf_fwd, ion_taf_rev = taf_rev), by = "key")
  nrow(df_igv)
  arrange(df_igv, -(ion_taf_fwd-taf_fwd)) %>% 
    select(chr:taf_diff, 
           t_fmt_sb,n_fmt_sb,
           ion_taf_fwd:ion_taf_rev, 
           n_fmt_f1r2:n_fmt_f2r1_alt2, 
           t_fmt_f1r2:t_fmt_f2r1_alt2) %>% 
    head() %>% data.frame()
  
  library(ggpubr)
  head(df_igv$key);head(df_ion$key)
  select(df_igv, taf_fwd, taf_rev, ion_taf_fwd, ion_taf_rev)
  ggscatter(df_igv, "taf_fwd", "ion_taf_fwd")
  ggscatter(df_igv, "taf_rev", "ion_taf_rev")
  
  # ** ggscatter 
  p_dp_full = tidylog::filter(df_ion, type == "snp") %>% 
    gghistogram("t_dp") + scale_x_log10()
  # peak cov in this one sample seems to be about 50
  p_dp_filt = tidylog::filter(df_ion, in_filtered == TRUE, type == "snp") %>% 
    gghistogram("t_dp") + scale_x_log10()
  cowplot::plot_grid(p_dp_full, p_dp_filt)
  
  
  sock <- SRAdb::IGVsocket()
  source('~/Dropbox/projects/packs_dnaseq/R/igv_snapshot.R')
  source('~/Dropbox/public/github_wranglr/R/expect_columns.R')
  bams = c("/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/bams/WEX-334187-N_nochr.bam",
           "/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/bams/WEX-334187-T_nochr.bam")
  # debug(get_igv_snap)
  get_igv_snap(df_igv, 
               input_type = "bed",
               bams = bams, 
               method = "socket",
               sock = sock,
               outdir = "~/Downloads/igv/sarco_334187_m2_strand_bias")
  
  
}


# EXTRA ---------------

splt_vcf_func <- function(x){
  x = as.character(unlist(x))
  
  splt_func <- function(i){
    str_replace_all(x[i], pattern = "\\]", "")
    xi = gsub("[{", "", gsub("}]", "", x[i], fixed = TRUE), fixed = TRUE)
    splt = strsplit(xi, ",")[[1]]
    splt = gsub("(.?)\\:(.*)", "\\1\n\n\\2", splt)
    splt = strsplit(splt, "\n\n")
    vals = sapply(splt, tail, 1)
    names(vals) = tolower(sapply(splt, head, 1))
    ret = as.data.frame(t(vals), stringsAsFactors = FALSE)
    return(ret)
  }
  
  lst <- lapply(1:length(x), function(i){
    y = try(splt_func(i))
    ifelse(class(y)=="try-error", "", y)
  })
  mat = bind_rows(lst)
}



format_vcf_info_ion <- function(x){
  x$freqt = gsub("Freq", "", x$freqt)
  x$freqn = gsub("Freq", "", x$freqn)
  
}



#' .depreciated
#'
#' @param key something
#' @param fl vcf fl 
#' 
#' @export
.read_vcf <- function (fl, key = TRUE){
  vcf <- readLines(fl)
  header <- unlist(strsplit(gsub("#", "", vcf[max(grep("#",
                                                       vcf))]), "\t"))
  vcfMat <- do.call("rbind", args = lapply(vcf[-grep("#", vcf)],
                                           function(x) t(unlist(strsplit(x, "\t")))))
  if (is.null(vcfMat)) {
    vcfMat <- matrix(rep(NA, length(header)), nrow = 1)
  }
  colnames(vcfMat) <- header
  if (key) {
    vcfMat <- cbind(vcfMat, key = paste(vcfMat[, "CHROM"],
                                        ":", vcfMat[, "POS"], ":", vcfMat[, "REF"], "-",
                                        vcfMat[, "ALT"], sep = ""))
  }
  return(vcfMat)
}




#' Parse a somatic VCF, with two samples.
#'
#' @param x a vcf file
#' @param samp name of the 'tumor' sample
#' @param ref name of the 'reference' sample
#'
#' @export
parse_somatic_vcf.oncotator <- function(x, samp, ref){
  message("Reading file...")
  
  # add a read_vcf function
  #mat = read_vcf(x)
  require(dplyr)
  require(parallel)
  
  message("Reading file...")
  if(tools::file_ext(x) == "gz"){
    fl_con = gzfile(x)
  }else{
    fl_con = file(x)
  }
  
  # assuming that header is not longer than 100
  hd = scan(fl_con, what = "character", sep = "\n", n = 1000, quiet = TRUE)
  hd_row = grep("#CHROM", hd)
  hd = tolower(strsplit(hd[hd_row], "\t")[[1]])
  col_classes = paste(rep("c", length(hd)), collapse = "")
  
  #tab = data.table:::fread(x, data.table = FALSE, header = FALSE, skip = hd_row, sep = "\t", colClasses = "character")
  tab = readr::read_tsv(x, col_names = FALSE, skip = hd_row, col_types = col_classes)
  colnames(tab) = gsub("#chrom", "chrom", hd)
  #tab = as_tibble(tab)
  
  #mat = as_tibble(as.data.frame(mat,  stringsAsFactors = FALSE))
  #colnames(mat) = tolower(colnames(mat))
  
  # assume first column is tumor
  samp1 = colnames(tab)[grep("format", colnames(tab)) + 1]
  # assume 2nd column in normal
  samp2 = colnames(tab)[grep("format", colnames(tab)) + 2]
  message("samp1: ", samp1, " samp2: ", samp2)
  
  message("Parsing format columns...")
  s1cols = splt_vcf_format(x = tab[,samp1], format = tab$format, prefix = "s1_")
  s2cols = splt_vcf_format(x = tab[,samp2], format = tab$format, prefix = "s2_")
  
  message("Parsing info columns...")
  infocols = splt_vcf_info(tab$info)
  
  #message("Parsing func...")
  #funccols = splt_vcf_func(mat$func)
  #infocols[1:100,] %>% View()
  tab$samp1 = samp1
  tab$samp2 = samp2
  
  colsel = !colnames(tab) %in% c("format", samp1, samp2, "info")
  
  message("Assembling data")
  #mat = as_tibble(cbind(mat[, colsel], tcols, ncols, infocols, funccols))
  tab = as_tibble(cbind(tab[, colsel], s1cols, s2cols, infocols))
  
  tab = clean_names(tab)
  return(tab)
}

