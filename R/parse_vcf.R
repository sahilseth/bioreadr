
#' Parse a VCF file into a TSV to be processed and annotated
#'
#' @param x a vcf file to be parsed
#' @param cores something 
#'
#' @export
parse_vcf <- function(x, cores = 1){
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
  #tab = tbl_df(tab)
  
  samp = colnames(tab)[grep("format", colnames(tab)) + 1]
  format_cols = NA
  if(length(samp) > 0){
    message("Parsing format columns...")
    format_cols = splt_vcf_format(x = tab[, samp], format = tab$format, prefix = "")
    colnames(format_cols) = paste0("format_", colnames(format_cols))
  }else{
    message("no format columns")
  }
  
  message("Parsing info columns...")
  info_cols = splt_vcf_info(tab$info, cores = cores)
  colnames(info_cols) = paste0("info_", colnames(info_cols))
  
  colsel = !colnames(tab) %in% c("format", samp, "info")
  message("Assembling data")
  tab2 = tbl_df(cbind(tab[, colsel], format_cols, info_cols))
  return(tab2)
  
}

read_vcf = parse_vcf


splt_vcf_format <- function(x, format, prefix){
  x = as.character(unlist(x))
  format = as.character(unlist(format))
  
  lst <- lapply(1:length(x), function(i){
    #message(".")
    xi = x[i];formati = format[i]
    splt = strsplit(xi, ":")[[1]]
    nms = tolower(strsplit(formati, ":")[[1]])
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



#' adapted from john
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
parse_somatic_vcf <- function(x, samp, ref){
  message("Reading file...")
  
  # add a read_vcf function
  mat = read_vcf(x)
  mat = tbl_df(as.data.frame(mat,  stringsAsFactors = FALSE))
  colnames(mat) = tolower(colnames(mat))
  ## assume first column
  samp = colnames(mat)[grep("format", colnames(mat)) + 1]
  ref = colnames(mat)[grep("format", colnames(mat)) + 2]
  message("Parsing format columns...")
  tcols = splt_vcf_format(x = mat[,samp], format = mat$format, prefix = "t_")
  ncols = splt_vcf_format(x = mat[,ref], format = mat$format, prefix = "n_")
  message("Parsing info columns...")
  infocols = splt_vcf_info(mat$info)
  message("Parsing func...")
  #funccols = splt_vcf_func(mat$func)
  #infocols[1:100,] %>% View()
  mat$t_sample = samp
  mat$n_sample = ref
  colsel = !colnames(mat) %in% c("format", samp, ref, "info")
  message("Assembling data")
  #mat = tbl_df(cbind(mat[, colsel], tcols, ncols, infocols, funccols))
  mat = tbl_df(cbind(mat[, colsel], tcols, ncols, infocols))
  return(mat)
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
  #tab = tbl_df(tab)

  #mat = tbl_df(as.data.frame(mat,  stringsAsFactors = FALSE))
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
  #mat = tbl_df(cbind(mat[, colsel], tcols, ncols, infocols, funccols))
  tab = tbl_df(cbind(tab[, colsel], s1cols, s2cols, infocols))
  
  tab = clean_names(tab)
  return(tab)
}

if(FALSE){
  
  source('~/Dropbox/public/flow-r/ultraseq/ultraseq/R/parse_vcfs.R')
  x='/rsrch2/iacs/iacs_dep/sseth/flowr/runs/flowr_test/fastq_haplotyper-MS132-20150824-16-37-58-XScJT0OZ/tmp/GLizee-Pancreatic-MS132-MP013normalDNA.recalibed_1.haplotyper.vcf'
  x2 = parse_vcf(x)
  
  x='~/rsrch1_data/runs/sarcomatoid/jco/mutect/sarco10-T_1208XX_ST1374_073_H09WGADXX___sarco10-N_1208XX_ST1374_073_H09WGADXX_filt2_onco.vcf'
  x2 = parse_somatic_vcf(x)
  
  
  
}
