
# First and foremost, they should not be followed by a PAM with either 5′-NGG or 5′-NAG sequences. 
# Second, their global sequence similarity to the target sequence should be minimized, and guide sequences with genomic off-target loci that have fewer than three mismatches should be avoided. 
# Third, at least two mismatches should lie within the PAM-proximal region of the off-target site. 
# Fourth, a maximal number of mismatches should be consecutive or spaced less than four bases apart. 
# Finally, the amount of SpCas9 and sgRNA can be titrated to optimize on- to off-target cleavage ratio

# get position
#grep("[a-z]", strsplit(s, "")[[1]])

run_blast_fasta_pam <- function(fa, blast_db, cores = 20, out_prefix){
  
  message("> checking blast db")
  bl <- blast(db = blast_db)
  print(class(bl))
  print(bl)
  
  # read in the fasta file
  message("> reading fasta")
  seqs = readDNAStringSet(fa)
  
  nts = c("A", "T", "G", "C")
  pam1 = paste0(nts, "GG")
  pam2 = paste0(nts, "AG")
  pam = c(pam1, pam2)
  
  #seq = seqs[[1]]
  message("> adding the PAM sequences")
  
  pb <- txtProgressBar(min = 1, max = length(seqs), style = 3)
  tmp = lapply(seq_along(seqs), function(i){
    pb$up(i)
    seq = seqs[[i]]
    tmp = paste0(as.character(seq), pam)
    names(tmp) = rep(names(seqs[i]), length(pam))
    tmp
  })
  close(pb)
  
  message("> merging PAM sequences")
  seqs_pam = DNAStringSet(unlist(tmp))
  
  seqs_pam2 = seqs_pam[which(names(seqs_pam) == "sgUPP2_GB_2_65760_a")]
  
  # query a sequence using BLAST
  message("> mapping using blast")
  df_blast <- predict(bl, seqs_pam2, BLAST_args="-task blastn-short")
  colnames(df_blast) = tolower(colnames(df_blast))
  
  # add the strand
  df_blast = mutate(df_blast, strand = ifelse(s.start < s.end, "+", "-"))
  
  message("> writing output ..")
  write_sheet(df_blast, paste0(out_prefix, ".tsv"))
  
  df_blast
}


# -------- shRNA annotation functions -------- #
# this function is for shRNA
#' Title
#'
#' @param mp mapping file with columns: antisense_seq sense_seq
#' @param blast_db 
#' @param out_prefix 
#' @param do_revcomp 
#'
run_blast_map_seq <- function(mp, 
                              blast_db = "/rsrch2/iacs/iacs_dep/sseth/ref/human/annotations/gencode/v19/gencode.v19.pc_transcripts.fa", 
                              out_prefix, 
                              do_revcomp = TRUE){
  
  message("> checking blast db")
  bl <- blast(db = blast_db)
  #print(class(bl))
  #print(bl)
  
  df = read_sheet(mp)# %>% head()
  if(do_revcomp){
    # do rev comp
    message("> reading antisense")
    seqs = DNAStringSet(df$antisense_seq) %>% Biostrings:::reverseComplement() %>% as.character()
  }else{
    # read in the fasta file
    message("> reading sense map")
    seqs = df$sense_seq
  }
  
  names(seqs) = df$contig_id
  seqs = DNAStringSet(seqs)

  # query a sequence using BLAST
  message("> mapping using blast")
  df_blast <- predict(bl, seqs, BLAST_args="-task blastn-short")
  colnames(df_blast) = tolower(colnames(df_blast))
  
  # add the strand
  df_blast = mutate(df_blast, strand = ifelse(s.start < s.end, "+", "-"))
  
  message("> writing output ..")
  write_sheet(df_blast, paste0(out_prefix, ".tsv"))
  
  out_prefix
}


#grep '>' gencode.v19.pc_transcripts.fa > gencode.v19.pc_transcripts.txt
# ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/_README.TXT
# The pc_transcripts.fa files provide the cDNA sequences of coding transcripts in FASTA format with 
# the header listing:
#   transcript-id|
#   gene-id|
#   Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
#   Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
#   transcript-name|
#   gene-name|
#   sequence-length|
#   5'-UTR (3'-UTR if reverse strand) location in the transcript|
#   CDS location in the transcript|
#   3'-UTR (5'-UTR if reverse strand) location in the transcript
read_fa_hdr <- function(){
  fa_fl = "/rsrch2/iacs/iacs_dep/sseth/ref/human/annotations/gencode/v19/gencode.v19.pc_transcripts.txt"
  df_fa = read_sheet(fa_fl)
  df_fa = data.frame(fa_hdr = df_fa, stringsAsFactors = FALSE)
  df_fa
}


#' Split DNA sequences into smaller files.
#'
#' Given a vector of sequences or DNAStringSet or a FASTA filename, the function 
#' splits it into smaller pieces as denoted by totalFiles parameter.
#'
#' @param x a DNAStringSet object, or a FASTA filename.
#' @param totalFiles an integer indicating how many files to create. Default is 4.
#' @param suffix a word to add to each file created. Default is "tempy".
#' @param filename name of the file if x is a DNAStringSet object. 
#' Default is "queryFile.fa".
#' @param outDir directory to write the output file. Default is current directory.
#'
#' @return a vector of filename names created.
#'
#' @seealso \code{\link{blatSeqs}}
#'
#' @export
#'
#' @examples
#' seqs <- DNAStringSet(sapply(sample(c(100:1000), 500),
#' function(size) paste(sample(DNA_BASES, size, replace = TRUE), collapse = ""))) 
#' splitSeqsToFiles(seqs, 5, "tempyQ", "myDNAseqs.fa", tempdir())
splitSeqsToFiles <- function(x, totalFiles = 4, suffix = "tempy",
                             filename = "queryFile.fa", outDir = getwd()) {
  if(is.atomic(x)) {
    message("Splitting file ",x)
    totalSeqs <- length(fasta.info(x, use.names = FALSE))
    chunks <- round(totalSeqs/totalFiles)
    ## incase totalSeqs is lower than number of files to be created!
    chunks <- ifelse(chunks > 0, chunks, totalSeqs) 
    
    starts <- seq(0, totalSeqs, by = chunks) ## create chunks of starts    
    for(skippy in starts[starts != totalSeqs]) {
      filename.out <- 
        file.path(outDir, paste(x, skippy, runif(1), suffix, sep = "."))
      ## no need to read the entire file...save memory by reading in N lines
      query.tmp <- readBStringSet(x, nrec = chunks, skip = skippy) 
      writeXStringSet(query.tmp, filepath = filename.out, format = "fasta")            
    }
    return(list.files(path = outDir,
                      pattern = paste0(basename(x), ".*", suffix, "$"),
                      full.names = TRUE))
  } else if (class(x) == "DNAStringSet") {
    message("Splitting Reads.")
    totalSeqs <- length(x)
    chunks <- round(totalSeqs / totalFiles)
    starts <- seq(1, totalSeqs, by = chunks)
    stops <- unique(c(seq(chunks, totalSeqs, by = chunks), totalSeqs))
    stopifnot(length(starts) == length(stops))        
    for(skippy in 1:length(starts)) {
      filename.out <- 
        file.path(outDir, paste(filename, skippy, runif(1), suffix, sep = "."))
      writeXStringSet(x[starts[skippy]:stops[skippy]], filepath = filename.out,
                      format = "fasta")            
    }            
    return(list.files(path = outDir,
                      pattern = paste(filename, ".*", suffix, "$", sep = ""),
                      full.names = TRUE))
  } else {
    stop("Dont know what is supplied in parameter x.")
  }
}

parse_gencode_fa_header <- function(df_fa, 
                                    fa_hdr_col = "fa_hdr"){
  
  library(tidyr)
  fa_hdr_cols = c("transcript_id", "gene_id", 
                  "havana_gene_id", "havana_transcript_id", 
                  "transcript_name", "gene_name", 
                  "sequence_length",
                  "5_utr_fwd", "cds", "3_utr_fwd", "extra")
  
  # blast headers
  #mutate(df_fa, subjectid = gsub("\\s|\\d+", "", subjectid)) %>%
  # fa_hdr_col
  #fa_hdr_col = "subjectid"
  df_fa = separate(df_fa, col = "subjectid", into = fa_hdr_cols, 
                   sep = "\\|", remove = FALSE)
  # df_fa2 = separate(df_fa, col = fa_hdr_col, into = fa_hdr_cols, 
  #                   sep = "\\|", remove = FALSE)
  
  
  # further cleaning
  df_fa = mutate(df_fa,
         transcript_id = gsub(">", "", transcript_id))
  df_fa
}


# message("> writing output ..")
# write_sheet(df_blast, paste0(out_prefix, ".tsv"))
# 
# # will save all results
# tmp = mclapply(seq_along(seqs), function(i){
#   message(".", appendLF = FALSE)
#   #print(seqs[i]);print(class(seqs[i]))
#   
#   message("_", appendLF = FALSE)
#   write_sheet(cl, sprintf("%s_%s.tsv", out_prefix, i))
# 
# }, mc.cores = cores)
# 
# df_blast = bind_rows(tmp)

