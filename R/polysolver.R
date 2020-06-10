


polysolver <- function(samplename, 

                       seq_sample_id_wex_n = NULL,
                       seq_sample_id_wex_t = NULL,
                       bam_wex_t = NULL,
                       bam_wex_n = NULL, 
                       
                       singularity_exec_opts = "--bind /rsrch3",
                       polysolver_sif = opts_flow$get("polysolver_sif"),
                       polysolver_scripts = "~/apps/polysolver/4.1/scripts",
                       
                       funr_exe = opts_flow$get("R_funr_exe"),
                       paired_end
                     
){
  
  # why polysolver
  # https://www.biostars.org/p/365614/
  # In order to understand what has changed in HLA in the tumour cells, one has to know the HLA genotype of the normal, so, they go 'hand in hand' - you need both. 
  # Whatever program you use should have options to perform somatic HLA genotyping.
  # I know that people at UCL Cancer Institute and Francis crick Institute (both London) are using Polysolver and xHLA. Check all of the options for these programs.
  
  check_args()

  # /home/polysolver/scripts/shell_call_hla_type_test bam race includeFreq build format insertCalc outDi
  cmd_call_hla = glue("singularity exec {singularity_exec_opts} {polysolver_sif} bash {polysolver_scripts}/shell_call_hla_type {bam_wex_n} Unknown 1 hg19 STDFQ 0 polysolver")
  # convert to mhc type alleles
  cmd_hla2netmhc = glue("singularity exec {singularity_exec_opts} {polysolver_sif} /home/polysolver/scripts/convertToNetmhc.pl polysolver/winners.hla.txt polysolver/winners.hla_netmhc.txt /home/polysolver")
  
  
  # 1.2 POLYSOLVER-based mutation detection
  # another installation: 
  # https://github.com/jason-weirather/hla-polysolver
  # This tool detects mutations in HLA genes using POLYSOLVER-inferred HLA alleles as input.
  # Input parameters:
  #     -normal_bam_hla: path to the normal BAM file
  #     -tumor_bam_hla: path to the tumor BAM file
  #     -hla: inferred HLA allele file from POLYSOLVER (winners.hla.txt or winners.hla.nofreq.txt)
  #     -build: reference genome used in the BAM file (hg18, hg19 or hg38)
  #     -format: fastq format (STDFQ, ILMFQ, ILM1.8 or SLXFQ; see Novoalign documentation)
  #     -outDir: output directory
  #     -indiv: individual ID, used as prefix for output files
  # normal_bam_hla tumor_bam_hla hla build format outDir indiv
  cmd_hla_mut = glue("singularity exec {singularity_exec_opts} {polysolver_sif} bash {polysolver_scripts}/shell_call_hla_mutations_from_type {bam_wex_n} {bam_wex_t} polysolver/winners.hla.txt hg19 STDFQ polysolver {seq_sample_id_wex_t}")
  
  # 1.3 Annotation of mutations
  # This tool annotates POLYSOLVER detected HLA mutations with gene compartment and amino acid change information.
  # Script: shell_annotate_hla_mutations
  # Input parameters:
    # -indiv: individual ID, used as prefix for output files
    # -tarZipDir: tar zipped directory containing the raw call files (Mutect: call_stats*, Strelka: *all.somatic.indels.vcf)
    # -outDir: the output directory  
  cmd_mut_ann = glue("singularity exec {singularity_exec_opts} {polysolver_sif} bash /home/polysolver/scripts/shell_annotate_hla_mutations {seq_sample_id_wex_t} hla_mut.tar.gz polysolver")
  
  # seq_sample_id = seq_sample_id_wex_n
  polysolver_pvac_mhc1 = glue("polysolver/{seq_sample_id_wex_n}_polysolver_pvac_mhc1.txt")
  cmd_polysolver_pvac = glue("{funr_exe} my.ultraseq::to_pvacseq.polysolver x=polysolver/winners.hla_netmhc.txt outfile={polysolver_pvac_mhc1}")

  cmds = list(polysolver = c(cmd_call_hla, cmd_hla2netmhc,
                             cmd_hla_mut, cmd_mut_ann))
  flowmat = to_flowmat(cmds, samplename)  
  
  outfiles = list(polysolver_pvac_mhc1 = polysolver_pvac_mhc1,
                  hla_mut_tar = "hla_mut.tar.gz",
                  polysolver_dir = "polysolver")
  
  list(flowmat = flowmat, outfiles = outfiles)
}


#' to_pvacseq.polysolver
#' 
#' convert optitype for pvacseq
#'
#' @param path 
#' @param outfile outfile
#'
#' @export
to_pvacseq.polysolver <- function(x, outfile){
  
  # x = "~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/neoantigen/runs/185_003/tmp/polysolver/winners.hla_netmhc.txt"
  df_poly = data.table::fread(x, data.table = F, header = F)
  hla_alleles = df_poly$V1 %>% unlist() %>% paste0(collapse = ",")
  
  write(hla_alleles, outfile)
  outfile
  
}