



  
vcf_readcount_annotator <- function(){
  
  # decompose multiple allelic files  
  # vt decompose -s <input_vcf> -o <decomposed_vcf>
  
  # SNV file
  
  # indel file
  
  # bam-readcount -f <reference_fasta> -l <site_list> <bam_file> [-i] [-b 20]
  # The -i option must be used when running the indels site list in order to process indels in insertion-centric mode.
  
  # A minimum base quality of 20 is recommended which can be enabled using the -b 20 option.
  
  # The mgibio/bam_readcount_helper-cwl Docker container contains a bam_readcount_helper.py script that will create the snv and indel site list files from a VCF and run bam-readcount. Information on that Docker container can be found here: dockerhub mgibio/bam_readcount_helper-cwl.

  # Example bam_readcount_helper.py command
  
  # /usr/bin/python /usr/bin/bam_readcount_helper.py \
  # <decomposed_vcf> <sample_name> <reference_fasta> <bam_file> <output_dir>
  #   This will write two bam-readcount files to the <output_dir>: <sample_name>_bam_readcount_snv.tsv and <sample_name>_bam_readcount_indel.tsv, containing readcounts for the snv and indel positions, respectively.
  # 
  # Running the vcf-readcount-annotator
  # The readcounts for snvs and indels are then added to your VCF separately, by running the vcf-readcount-annotator twice.
  # 
  # Example vcf-readcount-annotator commands
  # 
  # vcf-readcount-annotator <decomposed_vcf> <snv_bam_readcount_file> <DNA|RNA> \
  # -s <sample_name> -t snv -o <snv_annotated_vcf>
  #   
  #   vcf-readcount-annotator <snv_annotated_vcf> <indel_bam_readcount_file> <DNA|RNA> \
  # -s <sample_name> -t indel -o <annotated_vcf>
  #   The data type DNA or RNA identifies whether you are annotating DNA or RNA readcount. DNA readcount annotations will be written to the AD/DP/AF format fields while RNA readcount annotations will be written to the RAD/RDP/RAF format fields. Please see the VAtools documentation for more information.  

}