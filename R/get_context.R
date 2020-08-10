

# either use module load OR
# /scratch/rists/hpcapps/reference/human/broad_hg19/fastas/
#' Title
#'
#' @param x data.frame with chr start end
#' @param bedtools_exe 
#' @param tmpfile1 
#' @param tmpfile2 
#' @param reflib 
get_context <- function(x,
												bedtools_exe = "bedtools",
                        tmpfile1 = tempfile(fileext = ".bed"), #"filtered_wex_ccp_rna.bed",
                        tmpfile2 = tempfile(fileext = ".bed"), #"filtered_wex_ccp_rna2.bed",
                        reflib = "~/ref/human/fastas/Homo_sapiens_assembly19.fasta"){
  
  library(bit64)
  
	# creating a input bed file
	message("creating input bed file...")
	mybed = select(x, chr, start, end) %>%
    mutate(start = start-4, end = end+3)
	
	print(str(mybed))
  
  # fix chrnames
  mybed$chr = gsub("chr", "", mybed$chr)
  mybed$chr = gsub("23", "X", mybed$chr)
  mybed$chr = gsub("24", "Y", mybed$chr)
  
  message("writing input bed file... ", tmpfile1)
  mutate(mybed, 
         start = as.integer64(start), 
         end = as.integer64(end)) %>% 
    readr::write_tsv(file.path(tmpfile1), col_names = FALSE)
  
  # use bedtools to get the fasta
  message("using bedtools to get fasta...")
  cmd = sprintf("%s getfasta -tab -fi %s -bed %s -fo %s", 
  							bedtools_exe, reflib, tmpfile1, tmpfile2)
  tmp = system(cmd)
  
  message("reading/parsing bedtools output...")
  tb = params::read_sheet(tmpfile2, ext = "tsv", header = FALSE)
  #                first three   middle  last three
  x$context = gsub("([ATGC]{3})([ATGC]?)([ATGC]{3})", "\\1N\\3", tb$V2)
  x$fasta_ref_allele = gsub("([ATGC]{3})([ATGC]?)([ATGC]{3})", "\\2", tb$V2)
  
  return(x)
}
