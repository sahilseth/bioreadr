


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


source('~/Dropbox/public/github_wranglr/R/expect_columns.R')

check_dnaseq_metadata <- function(trk, 
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
  expect_columns(trk, c("individual", "name", "sampletype", "normal", "bam"))
  
  if(stage == "variants_called"){
    .check_variant_columns(trk, normalize_file_names)
  }
  
  attr(trk, "metadata_type") = "dnaseq"
  trk
}

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



metadata_for_pyclone <- function(){
  
  
}


metadata_for_sfq <- function(){
  
}


metadata_for_phylowgs <- function(){
  
}





# END
