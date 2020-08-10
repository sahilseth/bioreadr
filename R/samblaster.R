


samblaster <- function(x,
                       samblaster_exe = opts_flow$get("samblaster_exe"),
                       samblaster_opts = opts_flow$get("samblaster_opts"),
                       execute = FALSE
                       ){

  check_args()

  dedupbam <- paste0(bam_prefix, ".marked.bam")

  cmd = sprintf("samtools view -h %s | samblaster -M %s",
                samtools_exe, x, samblaster_exe, dedupbam)

  if(execute)
    system(cmd)

  flowmat = to_flowmat(list(samblaster = cmd))

  return(list(flowmat = flowmat, outfiles = dedupbam))
}
