
#/risapps/src6/speedseq/bin/freebayes -f /scratch/genomic_med/mtang1/scratch/LOWPASS_WGS_LGG1/ref_genome/human_g1k_v37.fasta \
#--pooled-discrete \
#--min-repeat-entropy 1 \
#--genotype-qualities \
#--min-alternate-fraction 0.05 \
#--min-alternate-count 2 \
#--region $chrom:$start..$end \
#| somatic_filter 1e-5 18 0 \
#> somatic_temp/TCGA-02-2485_somatic_primaryVsnormal.$chrom:$start..$end.vcf


opts_flow$set(
  ref_fasta = "/scratch/genomic_med/mtang1/scratch/LOWPASS_WGS_LGG1/ref_genome/human_g1k_v37.fasta",
  freebayes_exe = "/risapps/src6/speedseq/bin/freebayes",
  freebayes_opts = "--pooled-discrete --min-repeat-entropy 1 --genotype-qualities --min-alternate-fraction 0.05 --min-alternate-count 2",
  tabix_exe = "/risapps/src6/speedseq/bin/tabix")


if(FALSE){

  x = "/scratch/genomic_med/mtang1/scratch/results/TCGA-02-2485-10A-01D-1494-08/TCGA-02-2485-10A-01D-1494-08.bam"
  y = "/scratch/genomic_med/mtang1/scratch/results/TCGA-02-2485-01A-01D-1494-08/TCGA-02-2485-01A-01D-1494-08.bam"

  library(ultraseq)
  out = ultraseq:::freebayes(x, y, samplename = "TCGA-02-2485-01A-01D-1494-08")
  out$cmd[24]
  #somatic_filter_exe = "~/Dropbox/public/github_ngsflows/inst/scripts/somatic_filter.sh")

}

freebayes <- function(tumor_bam,
                      normal_bam,
                      samplename,
                      outfile,

                      split_by_chr = TRUE,
                      ref_fasta = opts_flow$get("ref_fasta"),

                      freebayes_exe = opts_flow$get("freebayes_exe"),
                      freebayes_opts = opts_flow$get("freebayes_opts"),

                      somatic_filter_exe = system.file("scripts/somatic_filter.sh", package = "ultraseq")
){

  if(missing(outfile))
    vcf_prefix <- gsub(".bam", "", basename(x))
  else
    vcf_prefix <- gsub(".bam", "", basename(outfile))

  if(split_by_chr){
    chrs_info = get_fasta_chrs(ref_fasta)
    interval_opts = paste0("--region ", chrs_info)
    chrs_prefix <- paste(vcf_prefix, chrs_info, sep = "_") ## bam names
  }else{
    chrs_prefix = vcf_prefix
    intervals_opts = ""
  }

  ## confirm that all none of the arguments are null
  check_args(ignore = 'outfile')

  cmd_free = sprintf("source %s; %s -f %s %s %s %s %s | somatic_filter 1e-5 18 0 > %s.vcf",
                     somatic_filter_exe,
                     freebayes_exe, ref_fasta,
                     freebayes_opts, interval_opts, x, y, chrs_prefix)

  merged_vcf = vcf_prefix
  #cmd_tabix = sprintf("%s -f -p vcf %s",merged_vcf)

  cmds = list(freebayes = cmd_free)

  flowmat = to_flowmat(cmds, samplename)
  #echo -e "1\tTCGA-02-2485-10A-01D-1494-08\tNone\tNone\t0\t1\n1\tTCGA-02-2485-01A-01D-1494-08\tNone\tNone\t0\t2" >
  #  TCGA-02-2485_somatic_primaryVsnormal.ped

  return(list(flowmat = flowmat))

}
