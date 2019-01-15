

# Annotating your VCF with VEP
# Adding coverage data to your VCF
# Adding expression data to your VCF
# Creating a phased VCF of proximal variants

# maf>VCF-VEP -------



# adding read counts -------
# https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/readcounts.html
# Running the vcf-readcount-annotator
# If you have multiple files for SNVs and InDels you will first need to concatenate the two files together:
#   
#   cat snvs_bam_readcount_file indels_bam_readcount_file > bam_readcount_file
# You can now use the combined bam-readcount output file to add readcount information to your VCF.
# 
# vcf-readcount-annotator input_vcf bam_readcount_file DNA|RNA -s sample_name
# vcf-expression-annotator input_vcf expression_file kallisto|stringtie|cufflinks|custom gene|transcript


# phased variants -------
# https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/proximal_vcf.html
# Update sample names
# The sample names in the tumor.bam, the somatic.vcf, and the germline.vcf need to match. If they don’t you need to edit the sample names in the VCF files to match the tumor BAM file.
# 
# Combine somatic and germline variants using GATK’s CombineVariants
# /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar \
# -T CombineVariants \
# -R reference.fa \
# --variant germline.vcf \
# --variant somatic.vcf \
# -o combined_somatic_plus_germline.vcf \
# --assumeIdenticalSamples
# Sort combined VCF using Picard
# /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar SortVcf \
# I=combined_somatic_plus_germline.vcf \
# O=combined_somatic_plus_germline.sorted.vcf \
# SEQUENCE_DICTIONARY=reference.dict
# Phase variants using GATK’s ReadBackedPhasing
# /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar \
# -T ReadBackedPhasing \
# -R reference.fa \
# -I tumor.bam \
# --variant combined_somatic_plus_germline.sorted.vcf \
# -L combined_somatic_plus_germline.sorted.vcf \
# -o phased.vcf
# bgzip and index the phased VCF
# bgzip -c phased.vcf > phased.vcf.gz
# tabix -p vcf phased.vcf.gz
# The resulting phased.vcf.gz file can be used as the input to the --phased-proximal-variants-vcf option.
# 
# bgzip and index the input VCF
# In order to use the --phased-proximal-variants-vcf option you will also need to bgzip and index your main input VCF.
# 
# bgzip -c input.vcf > input.vcf.gz
# tabix -p vcf input.vcf.gz




pvacseq <- function(
  input_vcf,
  samplename,
  hla_alleles = "HLA-A*02:01,HLA-B*35:01,DRB1*11:01",
  mhc_tools = "MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign",
  pvacseq_opts = "-e 8,9,10",
  odir = "pvacseq",
  
  pvacseq_setup = "module load conda_/3.6;source activate pvactools;export KERAS_BACKEND=tensorflow"
  
  
){
  
  #hla_alleles = 
  #mhc_tools = "MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign"
  cmd_pvacseq1 = glue("{pvactools_setup}")
  cmd_pvacseq2 = glue("{pvacseq_py} run {input_vcf} {samplename} {hla_alleles} {mhc_tools} odir {pvacseq_opts}")
  
  
}

read_pvacseq <- function(x){
  df_all = read_tsv("pvacseq/MHC_Class_I/1004-t.all_epitopes.tsv")
  df_rnk = read_tsv("pvacseq/MHC_Class_I/1004-t.filtered.condensed.ranked.tsv")
  df_filt = read_tsv("pvacseq/MHC_Class_I/1004-t.filtered.tsv")
  
  
}














# END