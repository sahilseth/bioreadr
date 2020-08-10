# transfer out  ----

#' copy_data_base
#'
#' @param fls 
#' @param out_dir 
#' @param transfer 
#'
#' @export
copy_data_base <- function(fls, 
                           out_dir, 
                           transfer = FALSE){
  if(!missing(path.file))
    paths = scan(path.file, what = "character", sep = "\t")
  samps = list.files(paths, full.names = TRUE)
  cat("We have",length(samps), "samples in this project\n", 
      paste(basename(samps), collapse = "\n"))
  
  if(transfer){
    cat("Transferring...\n")
    dir.create(file.path(out_dir, project.id, "orig.fqs"), 
               recursive = TRUE, showWarnings = FALSE)
    cmds = sprintf("rsync -a %s %s/%s/orig.fqs/", 
                   samps, out_dir, project.id)
    
    tmp <- sapply(cmds, system)
    cat("All done\n")
  }
  
}




if(FALSE){
  if(copy_results){
    # add run metrics to this cmd
    #metrics_exe = "/rsrch2/iacs/iacs_dep/sseth/flows/TH/run_metrics.R"
    metrics <- gsub(".bam$",".metrics.tsv", bam_name)
    
    cmd_mkdir <- glue("mkdir -p {out_path}")
    cmd_metrics <- "flowr my.hass:::get_metrics wd=../ force=TRUE"
    
    cmd_copyresults <- sprintf("rsync -av %s* %s/", 
                               gsub(".bam$","", bam_name), out_path, out_path)
    cmd_copymetrics <- sprintf("cp ../df_metrics.tsv %s/%s", out_path, metrics)
    cmd_perms <- glue("chmod -R ug+rw {out_path}/;chgrp -R ccct-group {out_path}/")
    
    cmd_transfer = paste(cmd_mkdir, cmd_metrics, cmd_copymetrics, 
                         cmd_copyresults, cmd_perms, sep = ";")
    
    cmds$copyresults = cmd_transfer
  }
}
