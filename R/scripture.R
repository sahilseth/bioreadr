

opts_flow$set(java_exe = "module load jdk/1.8.0_45; java",
         peak_fdr = 0.05)

scripture <- function(x,
                      samplename,

                      chrs = c(1:22,"X","Y","M"),
                      peak_fdr = opts_flow$get("peak_fdr"),

                      samtools_exe = opts_flow$get("samtools_exe"),
                      java_exe = opts_flow$get("java_exe"),
                      scripture_exe = opts_flow$get("scripture_exe")

                      ){


  cmd_hg19table <- sprintf("%s idxstats %s | grep -E \"^chr[0-9XYM]{1,2}\" | awk '{OFS=\"\t\" ; print $1,$2}'  > %s.hg19.table",
                           samtools_exe, bam, sample_name)

  ## input: hg19 table
  cmd_scripture = sprintf("%s -Xmx2g -jar %s -task chip -alignment  -chr chr%s -out %s_chr%s.bed -windows 750,1500,5000 -trim -sizeFile %s.hg19.table -alpha %s",
                          java_exe, scripture_exe, samplename, chrs, samplename, chrs, samplename, peak_fdr)

  bedfile = sprintf("%s_scripture_out_%s_%s.bed", samplename, peak_fdr, time_tag)
  cmd_scripture_bed <- sprintf("cat %s_chr*.bed > %s",
                               samplename, bedfile)

  cmds = list(hg19table = cmd_hg19table,
              scripture = cmd_scripture,
              scripture_bed = cmd_scripture_bed)

  flowmat = to_flowmat(cmds, samplename = samplename)

  return(list(flowmat = flowmat, outfile = bedfile))
}
