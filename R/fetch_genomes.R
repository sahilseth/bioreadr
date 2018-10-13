
#' @title fetch_genomes
#' @description fetch_genomes
#'
#'
#' @param genome_path output path to be used for downloading
#' @param verbose logical specifying verboseness
#' @param ... passed onto \link{fetch_genomes}
#'
#' @export
#' 
#' @importFrom curl curl_download
#'
#'
#' @examples \dontrun{
#' genomes_fetch(species = Homo_sapiens, src = 'NCBI', build = 'build37.2')
#' }
fetch_genomes <- function(genome_path = "~/flowr/genomes", verbose = FALSE,
                          ...){
  if(!file.exists(genome_path)) dir.create(genome_path, recursive = TRUE)
  setwd(genome_path)
  ## ---------- if any of these are missing call fetch to check
  tmp = avail_genomes(...)
  if(tmp$type == "tar"){
    url = tmp$url
    tarfl = tmp$lst
    message("Downloading tar...", tarfl, " from ", url)
    #tmp <- getURL(url, verbose=verbose)
    curl_download(file.path(url, tarfl), tarfl, quiet = !verbose)
    message("Extracting tar...", tarfl)
    message("All one with", tarfl)
    invisible()
  }
}




#' Use Illumina's iGenomes and get a list of available genomes
#'
#' additional details
#'
#' @param species species
#' @param src source
#' @param build build
#' @param from from
#' @param base_url base url to be used
#' @param verbose logical, specifying verboseness
#' @param ... not used
#'
#' @export
#' 
#' @importFrom RCurl getURL
#'
#' @examples
#' gen = avail_genomes(species = 'Homo_Sapiens', from = 'igenomes')
avail_genomes <- function(species, src, build,
                          from = "igenomes",
                          verbose = FALSE,
                          base_url = "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com", ...){
  url = paste(base_url)
  head = "\n################################################\n"
  msg = c(head, 'Available Species:',head)
  type = "list"
  if(!missing(species)){
    msg = c(head, 'Available Sources:',head)
    url = paste(base_url, species, sep = "/")
    if(!missing(src)){
      msg = c(head, 'Available builds:',head)
      url = paste(base_url, species, src, sep = "/")
      if(!missing(build)){
        msg = c(head, 'Available files:',head)
        url = paste(base_url, species, src, build, sep = "/")
        type = "tar"
      }
    }
  }
  message(msg)
  if(verbose) message(url)
  lst = getURL(sprintf("%s/", url), dirlistonly = TRUE)
  message(lst)
  message("Example:\nfetch_genomes species=Homo_sapiens src=NCBI build=build37.2")
  invisible(list(lst = lst, url = url, type = type))
}
