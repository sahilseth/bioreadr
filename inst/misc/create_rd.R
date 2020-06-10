# get parsed source into roxygen-friendly format
env <- new.env(parent = globalenv())
#rfiles <- sapply(myfiles, function(f) file.path(mydir,f))

rfiles = list.files("inst/pipelines", ".R$", full.names = TRUE)[1:5]
blocks <- unlist(lapply(rfiles, roxygen2:::parse_file, env=env), recursive=FALSE)
parsed <- list(env=env, blocks=blocks)

# parse roxygen comments into rd files and output then into the "./man" directory
roc <- roxygen2:::rd_roclet()
results <- roxygen2:::roc_process(roc, parsed, "inst")
#roxygen2:::roc_output(roc, results, mydir, options=list(wrap=FALSE), check = FALSE)


run = "150806_SN1120_0347_AC7DL6ACXX"
