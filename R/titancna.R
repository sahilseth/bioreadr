
# singularity pull docker://quay.io/bcbio/bcbio-vc


# https://github.com/gavinha/TitanCNA/tree/master/scripts/R_scripts
# numClusters=3
# numCores=4
# for ploidy in [2, 3, 4]:
# for num_clusters in [1, 2, 3]:
# ## run TITAN for each ploidy (2,3,4) and clusters (1 to numClusters)
# echo "Maximum number of clusters: $numClusters";
# for ploidy in $(seq 2 4)
# do
# echo "Running TITAN for $i clusters.";
# outDir=run_ploidy$ploidy
# mkdir $outDir
# for numClust in $(seq 1 $numClusters)
# do
# echo "Running for ploidy=$ploidy";
# Rscript titanCNA_v1.10.1.R --id test --hetFile test.het.txt --cnFile test.cn.txt \
# --numClusters $numClust --numCores $numCores --normal_0 0.5 --ploidy_0 $ploidy \
# --chrs "c(1:22, \"X\")" --estimatePloidy TRUE --outDir $outDir
# done
# echo "Completed job for $numClust clusters."
# done
# 
# ## select optimal solution
# Rscript selectSolution.R --ploidyRun2=run_ploidy2 --ploidyRun3=run_ploidy3 --ploidyRun4=run_ploidy4 --threshold=0.05 --outFile optimalClusters.txt

.run_titan_example <- function(){
  
  # test trk file:
  trk = "~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex_b1/superfreq/df_sfq_metadata_ms51_dna_b1.tsv"
  
  
  
  
}

#' Title
#'
#' @param trk 
#' @param samplename 
#' @param num_cores 
#' @param test 
#'
titancna_somatic_matched <- function(trk,
                                    samplename,
                                    # outpath = "gatkcnv/",
                                    num_cores = nrow(df_trk), 
                                    test = T
){
  
  # check generic columns
  trk = metadata_for_dnaseq_tools(trk)
  expect_columns(trk, c("outprefix", "outprefix_paired"))
  
  # ceate a new trk with ALL reqd files
  trk %<>% 
    mutate(
      # skip taking basename of bam!
      # bam = basename(bam),
      # oprefix = glue("{outpath}{name}"),
      gatkcnv_counts = glue("{outprefix}.counts.tsv"), 
      gatkcnv_std_cr = glue("{outprefix}.standardizedCR.tsv"),
      gatkcnv_dn_cr = glue("{outprefix}.denoisedCR.tsv"),
      gatkcnv_acounts = glue("{outprefix}.allelicCounts.tsv"),
      
      ## ModelSegments outputs
      # oprefix_final = glue("{outpath}{name}_clean"),
      gatkcnv_hets = glue("{outprefix_paired}.hets.tsv"),
      gatkcnv_cr_seg = glue("{outprefix_paired}.cr.seg"),
      gatkcnv_final_seg = glue("{outprefix_paired}.modelFinal.seg"),
      ## CallCopyRatioSegments outputs
      gatkcnv_called_seg = glue("{outprefix_paired}.called.seg")
    )
  trk
  
  trk_tum = filter(trk, normal == "NO")
  trk_norm = filter(trk, normal == "YES")
  
  # wranglr::mkdir("gatkcnv")
  
  
  # ** .collectreadcounts -------
  # need for ALL samples  i=1
  out_readcnts = mclapply(1:nrow(trk), function(i){
    .collect_read_counts(trk$bam[i], 
                         trk$gatkcnv_counts[i],
                         test = test, redo = F)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** denoise -------
  out_dr = mclapply(1:nrow(trk), function(i){
    .denoise_counts(counts = trk$gatkcnv_counts[i],
                    std_cr = trk$gatkcnv_std_cr[i],
                    dn_cr = trk$gatkcnv_dn_cr[i],
                    test = test)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** plot denoise -------
  out_plot_dr = mclapply(1:nrow(trk), function(i){
    .plot_denoise_cr(std_cr = trk$gatkcnv_std_cr[i],
                     dn_cr = trk$gatkcnv_dn_cr[i],
                     outprefix = trk$outprefix[i], 
                     test = test, redo = F)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** .collectallelicounts -------
  i=1
  out_acounts = mclapply(1:nrow(trk), function(i){
    .collect_allelic_counts(bam = trk$bam[i], 
                            acounts = trk$gatkcnv_acounts[i],
                            test = test)
  }, mc.cores = 1)
  # tmp
  warnings()
  
  # ** modelSegments -------
  # womn
  i=3
  out_seg = mclapply(1:nrow(trk_tum), function(i){
    .model_segments(bam = trk_tum$bam[i], 
                    dn_cr = trk_tum$gatkcnv_dn_cr[i],
                    outprefix = trk_tum$outprefix_paired[i],
                    acounts = trk_tum$gatkcnv_acounts[i],
                    normal_acounts = trk_norm$gatkcnv_acounts[1],
                    
                    hets = trk_tum$gatkcnv_hets[i],
                    cr_seg = trk_tum$gatkcnv_cr_seg[i],
                    final_seg = trk_tum$gatkcnv_final_seg[i],
                    gatkcnv_modelseg_params = "--minimum-total-allele-count-case 30 --maximum-number-of-smoothing-iterations 25 --number-of-smoothing-iterations-per-fit 1 --kernel-variance-allele-fraction 0.8 --kernel-variance-copy-ratio 0.8",
                    called_seg = trk_tum$gatkcnv_called_seg[i],
                    test = test, redo = F)
  }, mc.cores = num_cores)
  # tmp
  
  # Number of segmentation-smoothing iterations per MCMC model refit. 
  # (Increasing this will decrease runtime, but the final number of segments may be higher. 
  # Setting this to 0 will completely disable model refitting between iterations.)
  
  # --minimum-total-allele-count-case:Integer
  # --minimum-total-allele-count-case 30: min depth of seq in case
  # (but then how would you call COMPLETE LOSS?)
  # Minimum total count for filtering allelic counts in the case sample, if available.  The
  # default value of zero is appropriate for matched-normal mode; increase to an appropriate
  # value for case-only mode.  Default value: 0. 
  
  # write_tsv(df_trk, path = "gatkcnv/df_trk.tsv")
  cmds <- list(gtkcnv.splt = unlist(out_readcnts), 
               gtkcnv.splt = unlist(out_acounts),
               gtkcnv.mrg = unlist(out_dr),
               gtkcnv.mrg = unlist(out_plot_dr),
               gtkcnv.mrg = unlist(out_seg))
  # mutect_gather_bams = cmd_mutect_gather_bams
  
  flowmat = to_flowmat(cmds, samplename = samplename) %>% 
    mutate(cmd = as.character(cmd))
  
  list(flowmat = flowmat, outfiles = list(), trk = trk)
}

