

# funcs --------
.collect_read_counts.split <- function(bam, 
                                       counts.auto,
                                       counts.allo,
                                       gatkcnv_intervals.auto = opts_flow$get("gatkcnv_intervals_auto"),
                                       gatkcnv_intervals.allo = opts_flow$get("gatkcnv_intervals_allo"),
                                       gatk4_exe = opts_flow$get("gatk4_exe")){
  
  check_args()
  cmd_auto = glue("mkdir -p `dirname {counts.auto}`;{gatk4_exe} --java-options '-Xmx32g' CollectReadCounts -I {bam} -L {gatkcnv_intervals.auto} -O {counts.auto} --format TSV --interval-merging-rule OVERLAPPING_ONLY")
  cmd_allo = glue("mkdir -p `dirname {counts.allo}`;{gatk4_exe} --java-options '-Xmx32g' CollectReadCounts -I {bam} -L {gatkcnv_intervals.allo} -O {counts.allo} --format TSV --interval-merging-rule OVERLAPPING_ONLY")
  
  cmds = c(cmd_auto, cmd_allo)
  list(cmds = cmds, outfiles = c(counts.auto, counts.allo))
}


.denoise_counts.split <- function(counts.auto, counts.allo,
                                  dn_cr.auto, dn_cr.allo,
                                  std_cr.auto, std_cr.allo,
                                  
                                  # merged
                                  dn_cr,
                                  std_cr,
                                  
                                  gatkcnv_pon_auto_hdfs = opts_flow$get("gatkcnv_pon_auto_hdfs"),
                                  gatkcnv_pon_allo_hdfs,
                                  
                                  gatkcnv_cnts_concat_exe = opts_flow$get("gatkcnv_cnts_concat_exe"),
                                  gatk4_exe = opts_flow$get("gatk4_exe")){
  
  check_args()
  
  # get cmds
  cmd_auto = glue("{gatk4_exe} --java-options '-Xmx12g' DenoiseReadCounts ", 
                  "-I {counts.auto} --count-panel-of-normals {gatkcnv_pon_auto_hdfs} ",
                  "--standardized-copy-ratios {std_cr.auto} ",
                  "--denoised-copy-ratios {dn_cr.auto}")
  cmd_allo = glue("{gatk4_exe} --java-options '-Xmx12g' DenoiseReadCounts ", 
                  "-I {counts.allo} --count-panel-of-normals {gatkcnv_pon_allo_hdfs} ",
                  "--standardized-copy-ratios {std_cr.allo} ",
                  "--denoised-copy-ratios {dn_cr.allo}")
  # merge both
  # this from glass throws errors in flowr due to curly braces etc...
  # cmd_dn_std = glue("awk 'FNR==1 && NR!=1 {{ while (/^@|^CONTIG/) getline; }} 1 {{ print }}' ",
  #                   "{std_cr.auto} {std_cr.allo} 1> {std_cr}")
  # cmd_dn_cr = "awk 'FNR==1 && NR!=1 {{ while (/^@|^CONTIG/) getline; }} 1 {{ print }}' \
  #         {dn_cr.auto} {dn_cr.allo} 1> {dn_cr}"
  cmd_dn_std = glue("{gatkcnv_cnts_concat_exe} {std_cr.auto} {std_cr.allo} {std_cr}")
  cmd_dn_cr = glue("{gatkcnv_cnts_concat_exe} {dn_cr.auto} {dn_cr.allo} {dn_cr}")
  
  cmd = c(cmd_auto, cmd_allo, cmd_dn_cr, cmd_dn_std)
  
  list(cmds = cmd, outfiles = dn_cr)
  
}




gatkcnv_somatic_matched.split <- function(trk,
                                          samplename,
                                          # outpath = "gatkcnv/",
                                          num_cores = 1, 
                                          gatkdir = pipe_str$gatkcnv$dir,
                                          
                                          gatkcnv_interval_mode = opts_flow$get("gatkcnv_pon_type"),
                                          gatkcnv_pon_m_hdfs = opts_flow$get("gatkcnv_pon_type"),
                                          gatkcnv_pon_f_hdfs = opts_flow$get("gatkcnv_pon_type"),
                                          
                                          run_cmds = F){
  # print function name
  # message(match.call()[1], " ", appendLF = F)
  
  check_args(ignore = c("gatkcnv_pon_allo_hdfs"))
  warnings()
  
  trk = gatkcnv_prep_filenms(trk, gatkdir, gatkcnv_interval_mode)
  warnings()
  
  trk_tum = dplyr::filter(trk, normal == "NO")
  trk_norm = dplyr::filter(trk, normal == "YES")
  gender = trk$gender[1]
  # wranglr::mkdir("gatkcnv")
  
  
  # ** .collectreadcounts -------
  # need for ALL samples  i=1
  out_readcnts = lapply(1:nrow(trk), function(i){
    .collect_read_counts.split(trk$bam[i], 
                               counts.auto = trk$gatkcnv_cnts_auto[i],
                               counts.allo = trk$gatkcnv_cnts_allo[i])
  })
  # tmp
  warnings()
  
  # ** denoise -------
  gatkcnv_pon_allo_hdfs = opts_flow$get(glue("gatkcnv_pon_{gender}_hdfs"))
  
  out_dr = mclapply(1:nrow(trk), function(i){
    .denoise_counts.split(counts.auto = trk$gatkcnv_cnts_auto[i],  
                          counts.allo = trk$gatkcnv_cnts_allo[i],  
                          dn_cr.auto = trk$gatkcnv_dn_cr_auto[i],
                          dn_cr.allo = trk$gatkcnv_dn_cr_allo[i],
                          std_cr.auto = trk$gatkcnv_std_cr_auto[i], 
                          std_cr.allo = trk$gatkcnv_std_cr_allo[i],
                          std_cr = trk$gatkcnv_std_cr[i],
                          dn_cr = trk$gatkcnv_dn_cr[i],
                          gatkcnv_pon_allo_hdfs = gatkcnv_pon_allo_hdfs)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** plot denoise -------
  out_plot_dr = mclapply(1:nrow(trk), function(i){
    .plot_denoise_cr(std_cr = trk$gatkcnv_std_cr[i],
                     dn_cr = trk$gatkcnv_dn_cr[i],
                     outprefix = trk$gatkcnv_prefix[i], 
                     run_cmds=run_cmds, redo = F)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** .collectallelicounts -------
  i=1
  out_acounts = mclapply(1:nrow(trk), function(i){
    .collect_allelic_counts(bam = trk$bam[i], 
                            acounts = trk$gatkcnv_acounts[i],
                            run_cmds=run_cmds)
  }, mc.cores = 1)
  # tmp
  warnings()
  
  # ** modelSegments -------
  # womn
  i=3
  out_seg = mclapply(1:nrow(trk_tum), function(i){
    .model_segments(bam = trk_tum$bam[i], 
                    dn_cr = trk_tum$gatkcnv_dn_cr[i],
                    outprefix = trk_tum$gatkcnv_prefix_paired[i],
                    acounts = trk_tum$gatkcnv_acounts[i],
                    normal_acounts = trk_norm$gatkcnv_acounts[1],
                    
                    hets = trk_tum$gatkcnv_hets[i],
                    cr_seg = trk_tum$gatkcnv_cr_seg[i],
                    final_seg = trk_tum$gatkcnv_final_seg[i],
                    gatkcnv_modelseg_params = "--minimum-total-allele-count-case 30 --maximum-number-of-smoothing-iterations 25 --number-of-smoothing-iterations-per-fit 1 --kernel-variance-allele-fraction 0.8 --kernel-variance-copy-ratio 0.8",
                    called_seg = trk_tum$gatkcnv_called_seg[i],
                    run_cmds=run_cmds, redo = F)
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
  
  # since this is tumor normal, we expect multiple cmds in each case
  out_readcnts = purrr::transpose(out_readcnts)
  out_acounts = purrr::transpose(out_acounts)
  out_dr = purrr::transpose(out_dr)
  out_plot_dr = purrr::transpose(out_plot_dr)
  out_seg = purrr::transpose(out_seg)
  
  cmds <- list(gatkcnv.splt = unlist(out_readcnts$cmds), 
               gatkcnv.splt = unlist(out_acounts$cmds),
               gatkcnv.mrg = unlist(out_dr$cmds),
               gatkcnv.mrg = unlist(out_plot_dr$cmds),
               gatkcnv.mrg = unlist(out_seg$cmds))
  # mutect_gather_bams = cmd_mutect_gather_bams
  
  flowmat = to_flowmat(cmds, samplename = samplename) %>% 
    dplyr::mutate(cmd = as.character(cmd))
  
  list(flowmat = flowmat, outfiles = list(), trk = trk)
}













# END
