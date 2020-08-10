


phlat <- function(samplename, 
                  fqs,
                  paired_end = T,
                  seq_sample_id,
                  phlat_sif = opts_flow$get("phlat.sif"),
                  phlat_dir = "/opt/phlat-release",
                  phlat_index_dir = opts_flow$get("phlat.index_dir"),
                  phat_opts = opts_flow$get("phlat.opts"),
                  odir = "phlat",

                  funr_exe = opts_flow$get("R_funr_exe")
                     
                     ){
  
  fq1=fqs$fq1;fq2=fqs$fq2;fq3=fqs$fq3
  check_args()
  
  # singularity exec /rsrch3/home/iacs/sseth/singularity-images/phlat_1.1.sif python -O /opt/phlat-release/dist/PHLAT.py
  # Usage: PHLAT.py -1 fastq1 [-2 fastq2] -index indexdir -b2url b2url[-orientation pairedReadOrientation] [-tag samplename] [-p nthreads] [-e phlatdir] [-o outdir] [-pe 1]
  # -1: fastq file of the reads if single-end, or the first reads if paired-end
  # -2: fastq file of the second reads if paired-end;ignored if single-end
  # -index: url to the index files for Bowtie2 [default: b2folder subfolder in phlat-release packge]
  # -b2url: url to Bowtie2 executable
  # -orientation: --fr, --rf, --ff etc as defined in bowtie2 [default --fr], this parameter is not used for single-end input
  # -tag: name label for the sample associated with the fastq files
  # -p: number of threads for running Bowtie2 [default 8]
  # -e: url to the home folder of phlat-release package
  # -o: url to a directory where results shall be stored
  # -pe: flag indicating whether the data shall be treated as paired-end(1) or single-end(0) [default 1]
  # -tmp: flag for whether temporary folder should be kept[default 0]
  if(paired_end)
    cmd_phlat = glue("singularity exec {phlat_sif} python -O {phlat_dir}/dist/PHLAT.py -1 {fq1} -2 {fq2} -index {phlat_index_dir} -b2url bowtie2 -orientation '--fr' -tag {seq_sample_id} -e {phlat_dir} -o phlat -tmp 0 {phat_opts}")
  else
    cmd_phlat = glue("singularity exec {phlat_sif} python -O {phlat_dir}/dist/PHLAT.py -1 {fq3} -pe 0 -index {phlat_index_dir} -b2url bowtie2 -tag {seq_sample_id} -e {phlat_dir} -o phlat -tmp 0 {phat_opts}")
  
  cmd_phlat = glue("mkdir {odir};{cmd_phlat}")
  
  # seq_sample_id = seq_sample_id_wex_n
  phlat_fl = glue("{odir}/{seq_sample_id}_HLA.sum")
  phlat_pvac_mhc1 = glue("{odir}/{seq_sample_id}_phlat_pvac_mhc1.txt")
  phlat_pvac_mhc2 = glue("{odir}/{seq_sample_id}_phlat_pvac_mhc2.txt")
  cmd_hla_pvac = glue("{funr_exe} my.ultraseq::to_pvacseq.phlat phlat_fl={phlat_fl} outfile1={phlat_pvac_mhc1} outfile2={phlat_pvac_mhc2}")

  cmds = list(phlat = c(cmd_phlat, cmd_hla_pvac))
  flowmat = to_flowmat(cmds, samplename)  
  
  outfiles = list(phlat_pvac_mhc1 = phlat_pvac_mhc1,
                  phlat_pvac_mhc2 = phlat_pvac_mhc2)
  
  list(flowmat = flowmat, outfiles = outfiles)
}



#' to_pvacseq.phlat
#' 
#' convert optitype for pvacseq
#'
#' @param path 
#' @param outfile1 outfile
#' @param outfile2 outfile
#'
#' @export
to_pvacseq.phlat <- function(phlat_fl, outfile1, outfile2){
  
  # setwd("~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/neoantigen/runs/185_003/tmp")
  # list.files(path)
  # fl = list.files(path, pattern = "result.tsv", recursive = T) %>% tail(1)
  df_phlat = data.table::fread(phlat_fl, data.table = F) %>% 
    mutate(gene = gsub("HLA_", "", Locus)) %>% 
    clean_names() %>% 
    mutate(allele1 = hla_trim_d4(allele1),
           allele2 = hla_trim_d4(allele2))
  
  hla_pvac_mhc1 = df_phlat %>% 
    dplyr::filter(gene %in% LETTERS) %>%
    mutate(allele1 = paste0("HLA-", allele1), 
           allele2 = paste0("HLA-", allele2)) %>% 
    select(allele1, allele2) %>% unlist() 
  
  # DRB1*11:01 OR DRB1*11:01-DRB1*11:01
  hla_pvac_mhc2 = df_phlat %>% 
    dplyr::filter(!gene %in% LETTERS) %>%
    mutate(allele = paste0(allele1, "-", allele2)) %>% 
    select(allele1, allele2, allele) %>% unlist() 
  
  # hla_pvac = c(hla_pvac_mhc1, hla_pvac_mhc2)
  
  paste0(hla_pvac_mhc1, collapse = ",") %>% write(outfile1)
  paste0(hla_pvac_mhc2, collapse = ",") %>% write(outfile2)
  
  c(outfile1, outfile2)

}