

#' bam_set
#'
#' @param bam 
#' @param outprefix 
#' @param ref.fasta 
#' @param interval_split 
#' @param split_by 
#'
#' @return
#' @export
#'
#' @examples
bam_set <- function(bam, 
                    outprefix, 
                    ref.fasta, 
                    interval_split,
                    split_by = c("none", "chr_len", "interval_split")){
  
  split_by = match.arg(split_by)
  obj = list()
  
  # only done by MAIN pipeline function
  #bam = tools::file_path_as_absolute(bam)
  assert_that(has_extension(bam, "bam"))
  obj$bam = bam
  
  # determine output file name
  if(missing(outprefix)){
    message("using bamname as outprefix")
    obj$out_prefix <- gsub(".bam", "", basename(bam))
  }else{
    obj$out_prefix <- outprefix
  }
  
  # if file is available determine whether to split for faster processing
  tmp <- get_fasta_chrs(ref.fasta)
  obj$chr_names = tmp$chrs
  obj$chr_lengths = tmp$lens
  
  if(split_by == "chr_len"){
    # bam names
    obj$out_prefix_chr <- paste(obj$out_prefix , obj$chr_names, sep = "_") 
    
    # get intervals, using length
    intervals = get_gatk_intervals(obj)
    obj$out_prefix_interval = 1:length(intervals) # interval files
    intervals = lapply(intervals, function(x){
      paste(paste0(" -L ", x),collapse = "")
      })
    obj$intervals = unlist(intervals)
    
  }else if(split_by == "interval_split"){
    intervals = list.files(interval_split, pattern = ".interval_list", full.names = T)
    obj$intervals = paste0(" -L ", intervals)
    obj$out_prefix_interval = sprintf("%s_%03d", obj$out_prefix, 1:length(intervals))

  }else if(split_by == "none"){
    
    obj$out_prefix_chr = obj$out_prefix 
    obj$gatk_intervals = ""
    
  }
  
  
  obj
}

# ref_fasta_path = "~/Dropbox/public/flow-r/my.ultraseq/pipelines/dnaseq/Homo_sapiens_assembly19.dict"
# its BEST to keep the order same as in dict, so that sorting and indexing works as desired
get_gatk_intervals <- function(obj){
  
  if(any(is.na(obj$chr_lengths))){
    intervals = obj$chr_names
  }else{
    # get lengths
    lens = as.numeric(obj$chr_lengths);names(lens) = obj$chr_names
    #lens = sort(lens, decreasing = T);
    longest_len = max(lens)
    
    intervals = list(paste0(names(lens[1]), ":1+"))
    #temp_size = lens["12"]
    temp_size = lens[12]
    chr = 13
    for(chr in names(lens)[-1]){
      new_size = temp_size + lens[chr]
      if(new_size <= longest_len){
        #message(chr, appendLF = F)
        int = paste0(chr, ":1+")
        intervals[[length(intervals)]] = c(intervals[[length(intervals)]], int)
        temp_size = new_size
      }else{
        #message(chr)
        int = paste0(chr, ":1+")
        intervals = c(intervals, int)
        temp_size = lens[chr]
      }
    }
    # add it to the end
    intervals[[length(intervals)]] = c(intervals[[length(intervals)]], "unmapped")
  }
  intervals
}




if(FALSE){
  debug(get_gatk_intervals)
  get_gatk_intervals(obj)
}


# single: tumor/normal
#   |__: can create a out-prefix
#
# multiple: tumor/normal
#   |__: assumed one for each chr.
#   |__: length   -  of both files should be same
#   |__: sequence -  of files and those of chrs, is assumed to be the same
#   |__: length of chr and tumor normal files should be the same
#   | 
#   |__: DICT is available. 
#   |__: length and sequence of dict, is assumed (lexicographically sorted)
#   |__: 
#   |__: 
#   |__: 


# create a definitive out prefix
# samplename, may not always be unique, using bams is better
# determine output file name, if out_prefix is not provided


# this is shared, by ALL callers
#' Create a single object, to be used by multiple modules
#'
#' @param tumor_bam tumor bam file
#' @param normal_bam normal bam file
#' @param out_prefix output prefix for paired sample analysis, example tumorname___normalname
#' @param is_merged are bam files merged - single bam file for all chromosomes [TRUE/FALSE]
#' @param split_by_chr should downstream analysis be split by chr.
#' @param ref_fasta_path path to reference fasta file
#' 
#' @import assertthat
#' @import flowr
#'
paired_bam_set <- function(tumor_bam, 
                           normal_bam,
                           out_prefix,
                           is_merged,
                           split_by_chr,
                           ref_fasta_path = opts_flow$get('ref_fasta_path')
                           
){
  
  check_args(ignore = "out_prefix")
  
  obj = list()
  
  
  obj$split_by_chr = split_by_chr
  
  obj$out_prefix_chr = ""
  
  # check number of tumor and normal: should be same
  # if number of tumor/normal more than 1, is_merged should be FALSE
  
  sapply(sapply(tumor_bam, has_extension, "bam"), assert_that)
  sapply(sapply(normal_bam, has_extension, "bam"), assert_that)
  obj$bam$tumor = tumor_bam
  obj$bam$normal = normal_bam
  
  assert_that(length(tumor_bam) == length(tumor_bam))
  obj$bam$len = length(tumor_bam)
  
  
  if(split_by_chr & is_merged & length(tumor_bam) > 1){
    stop("multiple bams supplied, expected 1; perhaps is_merged should be FALSE?")
    
  }else if(split_by_chr & !is_merged & length(tumor_bam) == 1){
    stop("single bam supplied, expected multiple (one for each chromosome); perhaps is_merged should be TRUE")
  }
  
  
  if(obj$bam$len > 1 & missing(out_prefix))
    stop(">> multiple bams supplied, and out_prefix is missing")
  
  if(obj$bam$len == 1 & missing(out_prefix))
    obj$out_prefix <- paste0(gsub(".bam", "", basename(tumor_bam)), "__", 
                             gsub(".bam", "", basename(normal_bam)))
  else
    obj$out_prefix = out_prefix
  
  ## if file is available determine whether to split for faster processing
  tmp <- get_fasta_chrs(ref_fasta_path)
  obj$chr_names = tmp$chrs
  obj$chr_lengths = tmp$lens
  
  if(split_by_chr){
    obj$out_prefix_chr <- paste0(obj$out_prefix, "_", obj$chr_names) # bam names
    obj$gatk_intervals = paste0(" -L ", obj$chr_names)             ## interval files
    
  }else{
    obj$out_prefix_chr = obj$out_prefix 
    obj$gatk_intervals = ""
  }
  
  obj
  
}


#' Read the associated dictionary file and return a list of chromosome names
#'
#' @param fa a reference genome fasta file
#' @param length provide length of each chromosome as well [FALSE]
#'
#' @export
#'
#' @importFrom params read_sheet
#' @importFrom tools file_path_sans_ext
#'
get_fasta_chrs <- function(fa = opts_flow$get("ref_fasta_path"),
                           length = FALSE){
  check_args()
  
  #dict = gsub("fasta$", "dict", x)
  dict = paste0(tools::file_path_sans_ext(fa), ".dict")
  
  # use default chrs if they are already set
  if(!is.null(opts_flow$get("ref_fasta_chrs")))
    return(list(chrs = opts_flow$get("ref_fasta_chrs"), lens = NA))
  
  if(!file.exists(dict)){
    message(c("We need a dictionary (extension: .dict) for the reference fasta file to proceed. ",
              "Follow this link to learn more: http://lmgtfy.com/?q=create+dict+fasta"))
    warning("By default using hg19 chrs.")
    chrs = c(1:22, "X", "Y", "MT")
    lens = rep(NA, length(chrs))
    
  }else{
    seqs = read_sheet(dict, ext = "tsv", skip = 1, header = FALSE)
    chrs = gsub("SN:", "", seqs[, 2], fixed = TRUE)
    lens = gsub("LN:", "", seqs[, 3], fixed = TRUE)
    # 
    # if(length)
    #   return(list(chrs = chrs, lens = lens))
  }
  
  return(list(chrs = chrs, lens = lens))
}

# using samtools
.get_bam_chrs <- function(bam, samtools_exe = "samtools") {
  cmd <- sprintf("%s view -H %s | grep  '\\@SQ' | cut -f 2,3",
                 samtools_exe, bam)
  chrs_info <- system(cmd, intern = TRUE)
  chrs_info <- do.call(rbind, lapply(chrs_info, function(x){
    y = strsplit(x, "\t|:")[[1]]
    y[c(2,4)]
  }))
  return (chrs_info)
}

# using Rsamtools
#' @import Rsamtools
get_bam_chrs <- function(x){
  if(file.exists(x)){
    out = Rsamtools::scanBamHeader(x)
    chrs = names(out[[1]]$targets)
  }else{
    message("bam does not exists, returning hg19 chrs")
    chrs = c(1:22,"X","Y","MT")
  }
  return(chrs)
}

