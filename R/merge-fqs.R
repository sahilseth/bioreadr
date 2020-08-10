# in-house


merge_fqs <- function (x, ...) {
  UseMethod("merge_fqs", x)
}

# https://www.biostars.org/p/81924/
# cat file1.gz file2.gz file3.gz > allfiles.gz
# The resulting hash/message digests should be identical.
# My experience is that the zcat method is around 40x slower, 
# but the cat method's resulting file is a few percent bigger depending on your the gzip parameters used in the methods.

merge_fqs.character <- function(x, samplename, outfile){
  
  cmd = sprintf("cat %s > %s", 
          paste(x, collapse = " "), 
          outfile)
  
  flowmat = to_flowmat(list(merge_fq = cmd), samplename)
  
  list(flowmat = flowmat, outfile = outfile)
}

merge_fqs.data.frame <- function(x, samplename, outfile){

  # required columns
  # samplename
  # merged_file
  # file
  message("> using merge_fqs.data.frame")
  reqd_cols = c("samplename", "file", "merged_file")
  if(!mean(reqd_cols %in% colnames(x)) == 1)
    stop("some required columns of fq sheet are missing for merging, we need:\n", 
         paste(reqd_cols, collapse = "\n"))
  
  tmp <- lapply(unique(x$samplename), function(samp){
    x.samp = filter(x, samplename == samp)
    
    # invoke vector function
    out = merge_fqs(x = x.samp$file, 
              samplename = x.samp$samplename[1],
              outfile = x.samp$merged_file)
    out$flowmat
  })
  
  flowmat = bind_rows(tmp)

  list(flowmat = flowmat)
}
