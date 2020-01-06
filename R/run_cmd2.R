#' run_cmd
#'
#' @param cmd cmd to run 
#' @param cmdname a descriptive name of the cmd
#' @param target target file, cmd is skipped if present
#' @param stderr redirect output (stderr, 2>) to a file
#' @param redo force rede
run_cmd <- function(cmd, cmdname, target, stderr = "/dev/stderr", redo = F){
  
  if(file.exists(target) & !redo){
    # if file exists, and no redo
    message("> ", cmdname, " is already complete, skipping. target: ", target, " exists")
  }else{
    message("> running ", cmdname, " creating target: ", target)
    cmd = paste0(cmd, " 2> ", stderr)
    # print(cmd)
    system(cmd)
  }
  
}
