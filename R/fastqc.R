

#' Fastqc
#'
#' @param fqs
#' @param fqpath
#' @param odir
#' @param fastqc_exe
#' @param cpu_fastqc
#' @param casava
#' @param fastqc_opts
#'
#' @export
#'
fastqc <- function(fqs,
									 samplename = opts_flow$get("samplename"),
         fqpath,
         odir,
         fastqc_exe = opts_flow$get("fastqc_exe"),
         cpu_fastqc = opts_flow$get("cpu_fastqc"),
         casava = TRUE,
         fastqc_opts = opts_flow$get("fastqc_opts")
         ){

  if(!missing(fqpath))
    fqs = list.files(fqpath, pattern = "fastq.gz", full.names = TRUE, recursive = TRUE)

  check_args(ignore = "fqpath")
  if(!mean(file.exists(fqs)))
    stop("Some files do not exist, please check")

  ## if input is casava: want to get a summary on fastqs
  if(casava){
    fastqc_opts = c(fastqc_opts, "--casava")
    fqs = paste(fqs, collapse = " ")
  }

  message("creating output directory")
  try(dir.create(odir, recursive = TRUE))
  
  fastqc_opts = paste(fastqc_opts, collapse = " ")
  cmds <- sprintf("%s -f fastq -o %s -t %s %s %s",
                  fastqc_exe, odir, cpu_fastqc, fastqc_opts, fqs)

  flowmat = to_flowmat(list(fastqc = cmds), samplename = samplename)

  return(list(flowmat = flowmat))
}
attr(fastqc, "type", "module" )
