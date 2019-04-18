

# build index -------
#' cd ~/rsrch2_home/ref/human/b37/indexes/kallisto/0.43.0/
#' kallisto=/rsrch2/iacs/iacs_dep/sseth/apps/kallisto/0.43.0/kallisto
#' $kallisto index -i gencode.v19.pc_transcripts ../../../annotations/gencode/v19/gencode.v19.pc_transcripts.fa
#' 
# kallisto index -i gencode.v19.vM11.pc_transcripts gencode.v19.vM11.pc_transcripts.fa

# kallisto quant -i index -o output pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq
# for ALL fq based tools, needs output of create_fq_sheet

# https://github.com/pachterlab/sleuth/issues/55
#"call": "/projects/ps-yeolab/software/kallisto/kallisto quant -i /projects/ps-yeolab/genomes/hg19/kallisto/gencode.v19.pc_transcripts.only_protein_coding_transcripts.fa.k31 --threads=4 -o /home/obotvinnik/projects/autism_brain_rnaseq/analysis/kallisto_concatenated/Sample_10A_Full_Set_All_Runs_Concatenated -l 320 --single /projects/ps-yeolab/seqdata/20150422_heather_all_runs_autism_brain_postmortem_data/final_sym_links_to_files/Sample_10A_Full_Set_All_Runs_Concatenated.fastq.gz"


#fq1=/rsrch2/iacs/iacs_dep/sseth/projects2/av_pdac_clonal_tracking/201606_rnaseq_e012/wrangl_01_merge_fq/flowname-RNA-20160627-18-20-13-pV51YJZ7/tmp/GM_AV_0134_R1.fastq.gz
#fq2=/rsrch2/iacs/iacs_dep/sseth/projects2/av_pdac_clonal_tracking/201606_rnaseq_e012/wrangl_01_merge_fq/flowname-RNA-20160627-18-20-13-pV51YJZ7/tmp/GM_AV_0134_R2.fastq.gz
#index=/rsrch2/iacs/iacs_dep/sseth/ref/human/indexes/kallisto/0.42.4/gencode.v19.pc_transcripts

#/rsrch2/iacs/apps/conda/3.18.9/bin/kallisto quant -i $index -o . $fq1 $fq2 > GM_AV_0134_kallisto.txt




#' kallisto PE fqsheet
#' 
#' @param fqsheet needs columns: samplename, lane, num, read
#' @param kallisto_exe kallisto_exe
#' @param kallisto_index_path kallisto_index_path
#' 
#' 
kallisto_pe.fqsheet <- function(fqsheet, 
                                kallisto_exe = "/rsrch2/iacs/apps/conda/3.18.9/bin/kallisto", 
                                kallisto_index_path = "/rsrch2/iacs/iacs_dep/sseth/ref/human/indexes/kallisto/0.42.4/gencode.v19.pc_transcripts"){
  
  samples = unique(fqsheet$samplename)
  
  #s = "AV-Heterogeneity-RNA205-PATC53P178"
  
  lst = lapply(samples, function(s){
    fqsheet.s = filter(fqsheet, samplename == s)
    
    outfile = sprintf("%s_kallisto.txt", s)
    
    # fq input string
    fq_str = arrange(fqsheet.s, lane, num, read) %>% select(file) %>% unlist() %>% paste(collapse = " ")
    
    # create kallisto cmd
    cmd_kal = sprintf("%s quant -i %s -o . %s", 
                      kallisto_exe, kallisto_index_path, fq_str, outfile)
    
    flowmat = to_flowmat(list(kallisto = cmd_kal), s)
    return(list(flowmat = flowmat, outfile = outfile))
  })
  
  flowmat  = bind_rows(lapply(lst,  "[[", 'flowmat'))
  outfiles = sapply(lst,  "[[", 'outfile')
  
  return(list(flowmat = flowmat, outfiles = outfiles))
  
}


#' Run a paired end kallisto pipeline
#'
#' @param fq1 something
#' @param fq2 something 
#' @param samplename something 
#' @param out_prefix something 
#' @param kallisto_exe something 
#' @param kallisto_index_path something 
#' @param kallisto_opts something 
#'
#' @export
kallisto_pe <- function(fq1, fq2, 
                        samplename,
                        
                        out_prefix,
                        
                        kallisto_exe = opts_flow$get("kallisto_exe"),
                        kallisto_index_path = opts_flow$get("kallisto_index_path"),
                        #"/rsrch2/iacs/iacs_dep/sseth/ref/human/indexes/kallisto/0.42.4/gencode.v19.pc_transcripts", 
                        kallisto_opts = opts_flow$get("kallisto_opts")
                        #"-b 30"
                        
){
  
  if(missing(out_prefix))
    out_prefix = paste0(samplename, "_kallisto")
  
  check_args()
  
  # create kallisto cmd
  cmd_kal = sprintf("%s quant -i %s %s -o %s %s %s", 
                    kallisto_exe, kallisto_index_path, kallisto_opts, out_prefix, fq1, fq2)
  
  flowmat = to_flowmat(list(kallisto = cmd_kal), samplename)
  
  return(list(flowmat = flowmat, outfiles = out_prefix))
  
}

#' Run a paired end kallisto pipeline
#'
#' @param fq1 something
#' @param fq2 something 
#' @param samplename something 
#' @param out_prefix something 
#' @param kallisto_exe something 
#' @param kallisto_index_path something 
#' @param kallisto_opts something 
#'
#' @export
kallisto_se <- function(fq1,
                        samplename,
                        
                        out_prefix,
                        
                        kallisto_exe = opts_flow$get("kallisto_exe"),
                        kallisto_index_path = opts_flow$get("kallisto_index_path"),
                        #"/rsrch2/iacs/iacs_dep/sseth/ref/human/indexes/kallisto/0.42.4/gencode.v19.pc_transcripts", 
                        kallisto_opts = opts_flow$get("kallisto_opts")
                        #"-b 30"
                        
){
  
  if(missing(out_prefix))
    out_prefix = paste0(samplename, "_kallisto")
  
  check_args()
  
  # create kallisto cmd
  cmd_kal = sprintf("%s quant -i %s %s -o %s %s", 
                    kallisto_exe, kallisto_index_path, kallisto_opts, out_prefix, fq1)
  
  flowmat = to_flowmat(list(kallisto = cmd_kal), samplename)
  
  return(list(flowmat = flowmat, outfiles = out_prefix))
  
}

read_kallisto <- function(x, reader = read_sheet, ...){
  df = reader(x, ...) %>% tbl_df()
  
  # summarize gene level counts
  df = separate(df, col = target_id, 
                into = c("transcript_id", "gene_id"), sep = "\\|", remove = FALSE, extra = "drop")
  # head(df) %>% as.data.frame()
  
  # sum est_counts and transcripts
  df = group_by(df, gene_id) %>%
    mutate(gene_est_counts = sum(est_counts), 
           gene_tpm = sum(tpm))
  df
}


#' parse_kallisto_metrics
#'
#' @param pattern something 
#' @param wd path for flowr dir
parse_kallisto_metrics <- function(wd, pattern = "kallisto"){
  
  if(missing(pattern))
    kal.dir <- grep("kallisto", list.dirs(wd, recursive = FALSE, full.names = FALSE), value = TRUE)
  
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
}
# wd = "~/.rsrch1/iacs/iacs_dep/sseth/projects2/ss_tnbc/run_rnaseq/kallisto_pipe/flows/genohub1827102/flowname-20170719-16-26-19-ODdFEWlt"


#' Title
#'
#' @param x one metric file
.parse_kallisto_metrics <- function(x){
  
  if(!file.exists(x))
    return(NA)
  
  y = scan(x, what = "character", sep = "\n", quiet = TRUE)
  
  # extract relevant columns
  y1 = grep("[quant]", y, value = T, fixed = T)
  y2 = grep("[index]", y, value = T, fixed = T)
  y3 = grep("[   em]", y, value = T, fixed = T)
  y = c(y1, y2, y3)
  
  y = stringr:::str_trim(y) ## trim whitespace
  
  # extract specific values
  vars = c(kmer_length = "k-mer length: ", 
           est_avg_frag_length = "estimated average fragment length: ",
           num_kmers = "number of k-mers: ")
  
  tmp = lapply(seq_along(vars), function(i){
    var = vars[i]
    df = grep(var, y, value = T) %>% strsplit(":") %>% 
      do.call(rbind, .)
    df[, 2] = stringr::str_trim(df[, 2])
    df[, 1] = names(var)
    data.frame(df, stringsAsFactors = F)
  }) %>% bind_rows()
  names(tmp) = c("var", "value")
  
  tmp2 = gsub(".* processed (.*) reads, (.*) reads pseudoaligned", "\\1;\\2", y[5]) %>% 
    strsplit(";")
  tmp2 = data.frame(var = c("total_reads", "pseudoaligned"), value = tmp2[[1]], stringsAsFactors = F)
  
  df_metrics = bind_rows(tmp2, tmp)
  
  df_metrics = mutate(df_metrics, 
                      value = gsub(",", "", value), 
                      value = as.numeric(value))
  df_metrics
}

# x="~/.rsrch1/iacs/iacs_dep/sseth/projects2/ss_tnbc/run_rnaseq/kallisto_pipe/flows/genohub1827102/flowname-20170719-16-26-19-ODdFEWlt/flowname-1732R-02-01-20170719-16-26-19-Ua0ci9ZH/002.kallisto/kallisto_cmd_1.out"


cp_kallisto <- function(flow_base_path = "~/flows/AV/20180313_CT_RNA/flowname-20180716-17-04-03-5WlTgRNo", 
                       out_path = "~/.rsrch1/iacs/iacs_dep/sseth/data/runs/20180313_CT_RNA/kallisto", 
                       runid = basename(flow_base_path)
                       ){
  #library(glue)
  #kallisto_dirs = list.files();kallisto_dirs
  cmd = glue("rsync -avP {flow_base_path}/*/tmp/*_kallisto {out_path}/{runid}/")
  system(cmd)
  glue("{out_path}/{runid}")
}

#' create a tracking sheet for kallisto
#'
#' @param trk something 
#' @param base_out something 
#' @param project_id something 
#' @param sequencing_id something 
#'
#' @export
create_trk.kallisto <- function(trk, 
                                col_sequencing_id = "sample_id",
                                kallisto_path = "~/.rsrch1/iacs/iacs_dep/sseth/data/runs/20180313_CT_RNA/kallisto/flowname-20180716-17-04-03-5WlTgRNo"){
  
  seq_id = trk[, col_sequencing_id] %>% unlist() %>% as.character();seq_id
  trk$kallisto_fl = glue("{kallisto_path}/{seq_id}_kallisto/abundance.tsv")
  trk$kallisto_fl %>% file.exists() %>% table() %>% print()
  trk
}


if(FALSE){
  # get a set of kallisto out files
  wd = "~/.rsrch1/iacs/iacs_dep/sseth/projects2/ss_tnbc/run_rnaseq/kallisto_pipe/flows/genohub1827102/flowname-20170719-16-26-19-ODdFEWlt"
  wd = "~/.rsrch1/iacs/iacs_dep/sseth/projects2/av_pdac_ct/run_rnaseq/kallisto/wrangl_01_kallisto/RNA180/ iacs/iacs_dep/sseth/projects2/ss_tnbc/run_rnaseq/kallisto_pipe/flows/genohub1827102/flowname-20170719-16-26-19-ODdFEWlt"
  
  all_dirs = list.dirs(wd, recursive = T, full.names = FALSE)
  kal_dirs <- grep("002.kallisto", all_dirs, value = TRUE, fixed = T)
  kal_dirs = file.path(wd, kal_dirs)
  
  # get kallisto out files
  message("Using kal_dirs: ", basename(kal_dirs))
  if(missing(files))
    files = list.files(kal_dirs, pattern = ".out", full.names = TRUE)
  else
    message("Files supplied, will use those instead")
  
  
  df_metrics = lapply(files, function(x){
    tmp = .parse_kallisto_metrics(x)
    data.frame(wd = basename(dirname(dirname(x))), tmp, stringsAsFactors = F)
  }) %>% bind_rows()
  
  dcast(df_metrics, wd ~ var, value = "value") %>% View()
  
  
  
  
  
}













# END

