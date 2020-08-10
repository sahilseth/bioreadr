

# this calls the wrapper by brentp

#' bam_readcount_helper
#'
#' @param samplename 
#' @param vcf 
#' @param bam 
#' @param tumor_seq_sample_id 
#' @param ref_fasta 
#' @param ensembl_vep_setup 
#' @param vt_exe 
#' @param bam_readcount_helper_py 
#'
#' @return
#' @export
#'
#' @examples
bam_readcount_helper <- function(
  
  samplename,
  vcf,
  bam,
  tumor_seq_sample_id,
  
  ref_fasta = opts_flow$get("ref_fasta"),
  ensembl_vep_setup = opts_flow$get("ensembl_vep_setup"),
  
  vt_exe = opts_flow$get("vt_exe"),
  bam_readcount_helper_py = opts_flow$get("bam_readcount_helper_py")
){
  
  vep_setup_cmd = ensembl_vep_setup
  
  # Running vt decompose
  input_vcf = vcf
  output_vcf = gsub(".vcf", "_decomp.vcf", vcf)
  vt_decompose_cmd <- glue("vt decompose -s vcf -o {output_vcf}")
  input_vcf = output_vcf
  
  
  # runnning bam-readcount
  # bam-readcount -f <reference_fasta> -l <site_list> <bam_file> [-i] [-b 20]
  # The -i option must be used when running the indels site list in order to process indels in insertion-centric mode.
  # A minimum base quality of 20 is recommended which can be enabled using the -b 20 option.
  
  # The mgibio/bam_readcount_helper-cwl Docker container contains a bam_readcount_helper.py script that will create the snv and indel site list files from a VCF and run bam-readcount. Information on that Docker container can be found here: dockerhub mgibio/bam_readcount_helper-cwl.
  
  # bam_readcount_helper.py command
  # python /usr/bin/bam_readcount_helper.py <decomposed_vcf> <sample_name> <reference_fasta> <bam_file> <output_dir>
  bam_readcount_helper_cmd = glue("python {bam_readcount_helper_py} {input_vcf} {tumor_seq_sample_id} {ref_fasta} <bam_file> .")
  
  # Running the vcf-readcount-annotator
  # vcf-readcount-annotator <decomposed_vcf> <snv_bam_readcount_file> <DNA|RNA> -s <sample_name> -t snv -o <snv_annotated_vcf>
  vcf_ann_snv_cmd = glue("vcf-readcount-annotator <decomposed_vcf> <snv_bam_readcount_file> <DNA|RNA> -s <sample_name> -t snv -o <snv_annotated_vcf>")
  
  # vcf-readcount-annotator <snv_annotated_vcf> <indel_bam_readcount_file> <DNA|RNA> -s <sample_name> -t indel -o <annotated_vcf>
  vcf_ann_indel_cmd = glue("vcf-readcount-annotator <snv_annotated_vcf> <indel_bam_readcount_file> <DNA|RNA> -s <sample_name> -t indel -o <annotated_vcf>")
  
  
  
  
}










# END