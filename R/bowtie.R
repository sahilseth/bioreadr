#' @title parse_bowtie_metrics
#' @description parse_bowtie_metrics
#' @param files out files
#' @param wd working dir
#' 
#' @export
#' 
#' @examples
#' parse_bowtie_report(files = files, wd = wd)
parse_bowtie_metrics <- function(wd, files, aln.dir, cores = 1){
  
  if(missing(aln.dir))
    aln.dir <- grep("aln", list.dirs(wd, recursive = FALSE, full.names = FALSE), value = TRUE)
  aln.dir = file.path(wd, aln.dir)
  if(!file.exists(aln.dir)){
    message("aln.dir does not exist: ", basename(aln.dir))
    return()
  }
  if(length(list.files(aln.dir, pattern = "out"))==0){
    message("aln dir does not have *.out files yet, still running ?", basename(aln.dir))
    return()
  }
  
  message("Using aln.dir: ", basename(aln.dir))
  if(missing(files))
    files = list.files(aln.dir, pattern = ".out", full.names = TRUE)
  else
    message("Files supplied, will use those instead")
  
  fls = do.call(rbind, parallel::mclapply(files, function(x){
    y = scan(x, what = "character", sep = "\n", quiet = TRUE)
    y = y[grep("reads; of these:", y):grep("overall alignment rate", y)]
    #y = mypack:::trim_white(y) ## trim whitespace
    y = stringr:::str_trim(y) ## trim whitespace
    y = grep("^[0-9]", y, value = TRUE) ## get rows with numbers
    y = do.call(rbind, strsplit(gsub("(^[0-9]*)(.*)", "\\1,\\2", y), ","))
    ret = as.n(y[,1]);names(ret) = y[,2]
    names(ret) = stringr::str_trim(gsub(".*%[)]?(.*)", "\\1", names(ret)));ret
    return(ret)
  }, mc.cores = cores))
  aln_summary = colSums(fls)
  
  ## This is aln rate should be a percentage
  aln_summary["overall alignment rate"] = round(aln_summary["overall alignment rate"]/length(files), 2)
  
  ## These two may fail for PE
  aln_summary["bad_aln"] = round(100*aln_summary["aligned >1 times"] / aln_summary["reads; of these:"], 2)
  
  # detect is alignment is pe or se
  pe = sum(grepl("concordantly", names(aln_summary)))
  
  if(pe){
    uniq_map_reads = sum(c(aln_summary['aligned concordantly exactly 1 time'], 
                           aln_summary['aligned discordantly 1 time']))
  }else{
    uniq_map_reads = aln_summary["aligned exactly 1 time"]
  }
  
  aln_summary["good_aln"] = round(100*uniq_map_reads / aln_summary["reads; of these:"], 2)
  return(aln_summary)
}
