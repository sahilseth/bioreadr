#' Creates a sample sheet file names in the provided folder
#'
#' This function would check for files ending in fastq.gz, fq.gz, fastq, fq.
#'
#' @param path path to the fastq files
#' @param project name of the project.  \emph{optional}
#' @param subproject name of the subproject \emph{optional}
#' @param runid name of the flowcell this data is from \emph{optional}
#' @param outfile name of the output csv files \emph{optional}
#' @param format the format for names of fastq files, we have defaults for CASAVA and miSeq
#' @param pattern extensions this function will look for \emph{optional}
#' @param fix.names change the sample names such that they are acceptable as column names for R
#' 
#' @keywords samplesheet fastq casava
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' create_sample_mat(levelipath)
#' 
#' }
create_sample_sheet <- function(path, project, subproject, runid, format,
                                fix.names = FALSE,  fix.names.char = "-",
                                out_sep = c("\t", ","),
                                include.undetermined = FALSE,
                                pattern = "fastq.gz|fq.gz|fastq|fq", outfile){
  
  message("Fetching path(s) for fastq files...\n")
  .Deprecated("create_fq_sheet")
  
  fqs <- unlist(lapply(path, list.files, pattern = pattern,full.names=TRUE,recursive=TRUE))
  
  if(!include.undetermined) fqs <- grep("Undetermined", fqs, value = TRUE, invert = TRUE)
  if(length(fqs) == 0) stop("No fastq files detected in this folder\n")
  if(missing(project)) {project = basename(path); cat("\nDetecting project name:", project)}
  if(missing(subproject)){subproject = substr(project, 1, 2); cat("\nDetecting subproject:", subproject)}
  if(missing(runid)){runid = basename(dirname(dirname(dirname(fqs[1])))) ; cat("\nDetecting runid:", runid)}## runid
  
  out_sep = match.arg(out_sep)
  
  if(missing(outfile)){
    outfile = sprintf("%s_%s_%s_sample_mat.%s", project, subproject, runid,
                      switch(out_sep,
                             "," = "csv",
                             "\t" = "tsv"))
    cat("\nDetecting outfile:", outfile)
  }## folder for samplemat
  
  if(missing(format)){
    format <- detect_fq_format(fqs[1])
  }
  
  fq_mat <- split_names_fastq(files = fqs, format = format)
  
  if(fix.names){
    fq_mat$sample_id_orig = fq_mat$sample_id
    fq_mat$sample_id = fix_names(fq_mat$sample_id, char = fix.names.char)
    if(fix.names.char == ".") ## . opt style 2, good for data.frames
      fq_mat$sample_id = make.names(fq_mat$sample_id)
    ## ------- cleanup things: use _ to seperate out other things
  }
  
  cat("\nThere are", length(unique(fq_mat$samplename)), "samples in this folder")
  
  out_basename <- sprintf("%s-%s-%s_%s", project, subproject, fq_mat$samplename, runid)
  
  ## sorted_bam <- sprintf("%s_rg.sorted.bam",out_basename)
  ## recal_bam <- sprintf("%s_rg.sorted.recalibed.bam", out_basename)
  
  fq_mat <- cbind(fq_mat, out_basename, runid, project, subproject)
  fq_mat = fq_mat[!grepl("Undetermined", fq_mat$samplename), ] ## remove undermined
  outpath = dirname(outfile)
  
  if(!file.exists(outpath) & outpath!='.') dir.create(outpath) ## is X exists and not 'blank'
  
  write.table(fq_mat, file = outfile, row.names=FALSE, sep = out_sep, quote = FALSE)
  
  return(fq_mat = fq_mat)
}
