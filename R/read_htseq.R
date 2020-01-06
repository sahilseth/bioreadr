
#' Reads multiple HTSeq/QORTS files and creates a simple data.frame of all the counts. 
#'
#' @param x a vector of files
#' @param y a vector of sample names
#' 
#' @export
#' 
#' @importFrom parallel mclapply
#' 
read_htseq <- function(x, y, num_cores = 10){
  if (!length(y) == length(unique(y))) {
    stop("The supplied sample names (y) should be unique. Please check and run again.")
  }
  
  tmp <- mclapply(1:length(x), function(i) {
    fl = x[i]
    nm = make.names(y[i])
    dat = read.table(fl, sep = "\t", 
                     header = FALSE, 
                     col.names = c("feature", nm),
                     stringsAsFactors = FALSE)
    return(dat)
  }, mc.cores = num_cores)
  
  rows = sapply(tmp, nrow)
  
  # number of rows should be the same
  if (length(unique(rows)) > 1) {
    warning("Some of these files have different number of rows. Please check the library used for alignment:\n", 
            paste(y, rows, sep = "\t", collapse = "\n"))
    use_these = which(rows == min(rows))
    tmp = tmp[use_these]
    warning("Skipping these samples:\n", paste(y[-use_these], 
                                               collapse = "\n"))
  }
  
  dat = Reduce(function(...) merge.data.frame(..., by = "feature", all = FALSE), tmp)
  
  dat
  
}
