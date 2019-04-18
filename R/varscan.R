

varscan_germline <- function(bam, 
                    samplename, 
                    ref_fasta = opts_flow$get("ref_fasta"),
                    ref_dict = opts_flow$get("ref_dict"),
                    samtools_exe = opts_flow$get("samtools_exe"),
                    varscan_jar = opts_flow$get("varscan_jar"),

                    java_exe = opts_flow$get("java_exe_jdk7"), #module load jdk/1.7.0_79;
                    gatk4_exe = opts_flow$get("gatk4_exe"),
                    java_mem = opts_flow$get("java_mem")
                    
){
  
  
  # cd /rsrch2/iacs/iacs_dep/sseth/flows/SS/tnbc/ms51_wex_b1/varscan_germ
  # module load samtools_/1.9 jdk/1.7.0_79
  library(pacman)
  p_load(glue, flowr, janitor, readr, dplyr, parallel)
  
  #bam = "~/.rsrch1/genomic_med/omics/Womens/MS51/bam/IPCT-S4002-MOON0051-Cap1716-4-HTID295_181220-A00422-0024-BH722GDSXX-1-AGCGTTGCACGG.bwa_recalibed.bam"
  # samplename = "185_006_T0-D"
  # bam="185_006_T0-D.bwa_recalibed.bam"
  # samtools_exe = "/rsrch2/iacs/iacs_dep/sseth/apps/samtools/1.9/bin/samtools"
  # ref_fasta = "~/.rsrch2/iacs/iacs_dep/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta"
  # ref_dict = "~/.rsrch2/iacs/iacs_dep/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.dict"
  # samtools_exe = "/rsrch2/iacs/iacs_dep/sseth/apps/samtools/1.9/bin/samtools"
  # varscan_jar = "/risapps/noarch/varscan/2.3.9/VarScan.v2.3.9.jar"
  # java_exe = "java" #module load jdk/1.7.0_79;
  # gatk4_exe = "module purge;module load jdk/1.8.0_45;/rsrch2/iacs/iacs_dep/sseth/apps/gatk/4.0.10.1/gatk"
  # java_mem = "-Xmx16g"
  # samtools mpileup -d 10000 -q 1 -Q 15 -A -f my/reference/genome.fasta sample.bam > sample.mpileup
  # java -jar VarScan.v2.3.6.jar mpileup2cns sample.mpileup --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.01
  bamset  = ultraseq:::bam_set(bam = bam, ref_fasta_path = ref_fasta, split_by_chr = T)
  
  # -l FILE 	BED or position list file containing a list of regions or sites where pileup or BCF should be generated
  # - r STR 	Only generate pileup in region STR
  chrs = bamset$chr_names[1:25]
  oprefix = bamset$out_prefix
  oprefix_chr = bamset$out_prefix_chr[1:25]
  out_pileup_chr = glue("{oprefix_chr}.pileup.gz")
  varscan_vcfs = glue("{oprefix_chr}.varscan.vcf")
  varscan_vcf = glue("{oprefix}.varscan.vcf")
  
  # A Note on Indel Positions
  # VarScan calls consensus bases, SNPs, and indels at the position reported by SAMtools in the pileup file. 
  # For SNPs and consensus bases, this is the 1-based position of the site or variant. 
  # Indels, however, are reported at the base immediately upstream of where they occur. 
  # Thus, the first inserted base occurs at (position + 1) and the first deleted base occurs at (position + 1). 
  # This is why, for deletions, the reference base (e.g. A) may not match the deleted sequence (e.g. */C).
  cmd_pileup_varscan = glue("{samtools_exe} mpileup -d 10000 -q 1 -Q 15 -A -r {chrs} -f {ref_fasta} {bam} | ", 
                            "{java_exe} -jar {varscan_jar} mpileup2cns ",
                            "--variants --strand-filter 0 --p-value 0.01 ",
                            "--min-var-freq 0.01 --output-vcf 1 > {varscan_vcfs}") %>% 
    as.character()
  # cmd_pileup_varscan
  
  # cmd_mpileup = glue("{samtools_exe} mpileup -d 10000 -q 1 -Q 15 -A -r {chrs} -f {ref_fasta} {bam} | gzip > {out_pileup_chr}");
  # out = mclapply(cmd_mpileup, system, mc.cores = 20)
  # module load varscan/2.3.9
  # cmd_varscan = glue("zcat {out_pileup_chr}|{java_exe} -jar {varscan_jar} mpileup2cns ",
                     # "--variants --strand-filter 0 --p-value 0.01 ",
                     # "--min-var-freq 0.01 --output-vcf 1 > {varscan_vcfs}");
  # out = mclapply(cmd_varscan, system, mc.cores = 20)
  
  # cmd_mergevcf = glue("{java_exe} -jar {varscan_jar} mpileup2cns {out_pileup_chr} --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.01 > {varscan_vcfs}")
  # could use gather vcfs instead of this
  # https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.2.0/picard_vcf_MergeVcfs.php
  varscan_vcfs_i = paste0(varscan_vcfs, collapse = " -I ")
  varscan_vcf = paste0(samplename, ".varscan.vcf") # final VCF
  cmd_mergevcf = glue("{gatk4_exe} --java-options {java_mem} MergeVcfs -I {varscan_vcfs_i} -O {varscan_vcf} -D {ref_dict}") %>% 
    as.character()
  # system(cmd_mergevcf)
  
  flowmat = list(
    varscan_pileup = cmd_pileup_varscan, 
    # varscan_varscan = cmd_varscan, 
    varscan_mergevcf = cmd_mergevcf) %>% to_flowmat(samplename)
  
  list(flowmat = flowmat, 
       outfiles = list(varscan_vcf = varscan_vcf))
}










# module load jdk/1.8.0_45;/risapps/rhel6/gatk/4.1.0.0/gatk --java-options -Xmx16g MergeVcfs -I 185_006_T0-D.bwa_recalibed_1.varscan.vcf -I 185_006_T0-D.bwa_recalibed_2.varscan.vcf -I 185_006_T0-D.bwa_recalibed_3.varscan.vcf -I 185_006_T0-D.bwa_recalibed_4.varscan.vcf -I 185_006_T0-D.bwa_recalibed_5.varscan.vcf -I 185_006_T0-D.bwa_recalibed_6.varscan.vcf -I 185_006_T0-D.bwa_recalibed_7.varscan.vcf -I 185_006_T0-D.bwa_recalibed_8.varscan.vcf -I 185_006_T0-D.bwa_recalibed_9.varscan.vcf -I 185_006_T0-D.bwa_recalibed_10.varscan.vcf -I 185_006_T0-D.bwa_recalibed_11.varscan.vcf -I 185_006_T0-D.bwa_recalibed_12.varscan.vcf -I 185_006_T0-D.bwa_recalibed_13.varscan.vcf -I 185_006_T0-D.bwa_recalibed_14.varscan.vcf -I 185_006_T0-D.bwa_recalibed_15.varscan.vcf -I 185_006_T0-D.bwa_recalibed_16.varscan.vcf -I 185_006_T0-D.bwa_recalibed_17.varscan.vcf -I 185_006_T0-D.bwa_recalibed_18.varscan.vcf -I 185_006_T0-D.bwa_recalibed_19.varscan.vcf -I 185_006_T0-D.bwa_recalibed_20.varscan.vcf -I 185_006_T0-D.bwa_recalibed_21.varscan.vcf -I 185_006_T0-D.bwa_recalibed_22.varscan.vcf -I 185_006_T0-D.bwa_recalibed_X.varscan.vcf -I 185_006_T0-D.bwa_recalibed_Y.varscan.vcf -I 185_006_T0-D.bwa_recalibed_MT.varscan.vcf -O 185_006_T0-D.bwa_recalibed.varscan.vcf --SEQUENCE_DICTIONARY ~/ref/human/b37/fastas/Homo_sapiens_assembly19.dict






# END