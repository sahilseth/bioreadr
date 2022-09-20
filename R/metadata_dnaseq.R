


# STRUCTURE:
# REQUIRED (with values)
# individual  name  sampletype  normal  bam

# REQUIRED CALLS (missing acceptable)
# ssm_vcf  cnv_seg  sv_vcf
# this contains the final required columns for any downstream analysis:
# SSM/VCF: SSM simple somativ mutations
# SEG: copy number variantions
# SV: VCF of more complicated mutations
# 
# OPTIONAL:
# TERMS:
# 
# sampletype: 
#   t0, t1, ts: three time points
#   primary met
# 
# platform (check values):
#    wex, rnaseq, microarray, snparray, wgs
#    targetted, t200, ccp
# 
# batch:
#   ms51_wex_b1, ms51_wex_b2, 
#   (genohub) ms51_wex_g1, 
#             ms51_rna_g1

# ~/.rsrch1/iacs/iacs_dep/sseth/projects2/ss_tnbc/data/artemis/wex/ms51_dna_b1/bams/185_145_T0-D.bwa_recalibed.bam        NA      185_145 185_145_T0-D    T0      NO      ms51_dna_b1     2
# ~/.rsrch1/iacs/iacs_dep/sseth/projects2/ss_tnbc/data/artemis/wex/ms51_dna_b1/bams/185_145_T1-D.bwa_recalibed.bam        NA      185_145 185_145_T1-D    T1      NO      ms51_dna_b1     2
# ~/.rsrch1/iacs/iacs_dep/sseth/projects2/ss_tnbc/data/artemis/wex/ms51_dna_b1/bams/185_145_GB-D.bwa_recalibed.bam        NA      185_145 185_145_GB-D    GB      YES     ms51_dna_b1     2


# source('~/Dropbox/public/github_wranglr/R/expect_columns.R')

check_dnaseq_metadata <- function(...){
  .Deprecated("metadata_for_dnaseq")
  metadata_for_dnaseq(...)
}

#' metadata_for_dnaseq
#'
#' @param trk 
#' @param normalize_file_names 
#' @param basename_of_files 
#' @param stage 
#' @param check_files_exist 
#'
#' @export
#'
metadata_for_dnaseq <- function(trk, 
                                normalize_file_names = F,
                                basename_of_files = normalize_file_names,
                                stage = "none",
                                check_files_exist = F){
  
  pacman::p_load(testit)
  
  if(!is.data.frame(trk))
    trk = readr::read_tsv(trk) %>% data.frame(stringsAsFactors = F)
  
  # expect clean names (lower case, seperated by _)
  trk %<>% clean_names()
  
  # flexible to allow pairs, and multiple normals
  wranglr::expect_columns(trk, c("individual", "name", "normal", "bam"))
  # extra for longitudinal data
  # expect_columns(trk, "sampletype")
  wranglr::expect_columns(trk, "timepoint")
  
  if(stage == "variants_called"){
    .check_variant_columns(trk, normalize_file_names)
  }
  trk %<>% dplyr::mutate(gender = .resolve_gender.trk(gender), 
                         normal = .resolve_normal.trk(normal))
  
  attr(trk, "metadata_type") = "dnaseq"
  trk
}

.resolve_normal.trk <- function(x){
  vals = c("NO", "YES")
  if( !all(x %in% vals) ){
    stop("normal needs to have values: ", paste0(vals, collapse = " "))
  }
  x
}

.resolve_gender.trk <- function(x){
  sapply(x, function(xi){
    if (xi %in% c("male", "Male","MALE","m", "M")){
      xi = "m"
    }else{
      xi = "f"
    }
    xi
  })
}

# titan require full gender
.resolve_gender.titancna <- function(x){
  sapply(x, function(xi){
    if (xi %in% c("male", "Male","MALE","m", "M")){
    xi = "male"
  }else{
    xi = "female"
  }
  xi
  })
}


metadata_for_dnaseq_tools = metadata_for_dnaseq
as_dnaseq_metadata = metadata_for_dnaseq_tools
as_dnaseq_metadata = check_dnaseq_metadata
.check_variant_columns <- function(trk, basename_of_files = F){
  
  # bam file should already have been transferred
  if(basename_of_files){
    trk %<>% 
      mutate(bam = basename(bam),
             ssm_vcf = basename(ssm_vcf),
             cnv_seg = basename(cnv_seg))
  }
  
  # make sure files exist
  bams = trk$bam
  bais = gsub(".bam$", ".bam.bai", bams)
  muts = trk$ssm_vcf
  segs = trk$cnv_seg
  # # this will change based on input!
  # fits = gsub("_cncf.tsv", ".rds", segs)
  
  # ONLY REMOVE TEMPORARILY!!
  if(check_files_exist)
    testit::assert("check if all files exists", {
      # file.exists(bams) &
      #   file.exists(bais) &
      file.exists(segs) &
        # file.exists(fits) &
        file.exists(muts)
    })
  trk
}

.check_single_individual <- function(trk){
  # make sure it is a single individial
  inds = unique(trk$individual)
  if(length(inds) > 1)
    rlogging::warning("There are multiple individuals listed in the tracking sheet.\n", 
                      "This tool assumes all the tumor bams are from the same individual")
  
  trk
}


metadata_for_mutect1 <- function(trk, ...){
  metadata_for_dnaseq_tools(trk, ...)
  
  .check_single_individual(trk)
  
  
}

metadata_for_pyclone <- function(){
  
  
}


metadata_for_sfq <- function(){
  
}


metadata_for_phylowgs <- function(){
  
}





# END
