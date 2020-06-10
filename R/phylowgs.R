# phylowgs
# 
# 
# install:
# # instalation:
# https://bitbucket.org/aroth85/pyclone/wiki/Installation
# module load conda_/2.7
# conda create --name pyclone #python=2
# conda create --name phylowgs phylowgs


# # https://github.com/morrislab/phylowgs
# module purge
# module load python/2.7.15-bio
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda create -n phylowgs_py27 phylowgs python=2.7

# source activate phylowgs_py27



# For cnv_data.txt, the parser requires CNV calls. Battenberg and TITAN are supported. As the process of going from CNV calls to PhyloWGS input is complex -- the a and d values used by PhyloWGS for each CNV depend on the size of the CNV, the read depth, the total copy-number change, and other factors -- we perform CNV parsing as a two-step process:
#   
#   Run parse_cnvs.py to create the intermediate cnvs.txt file, which will contain for each CNV its chromosome, start and end coordinates, major and minor copy numbers, and cellular prevalence (i.e., fraction of cells in sample containing the CNV, not just the fraction of tumor cells containing the CNV).
# 
# Run create_phylowgs_inputs.py with the --cnvs <sample_name>=cnvs.txt parameter.
# <sample_name> is a string that uniquely identifies the given sample. This isn't terribly useful when running on only a single sample; in such cases, you can just name your sample S1, for example. Supporting other CNV callers is simply a matter of converting their results to our intermediate cnvs.txt format. For an example of how this file should be formatted, please see cnvs.txt.example.
# 
# Please note that your listing of CNVs must list regions of normal copy number as well. Any SSMs given in your VCF that fall outside regions listed in cnvs.txt will be ignored by the parser, as it assumes your CNV caller could not make a proper determination of copy number (whether normal or abnormal) for unlisted regions, meaning that PhyloWGS will be unable to correct for copy-number changes for the SSMs in question.




merge_all_calls_mutect_pindel <- function(trk, 
                                          somatic_mode = c("matched", "pon")){
  
  # https://github.com/morrislab/phylowgs/blob/262325b219e6d31f672791a05c6f927a18963ded/parser/create_phylowgs_inputs.py#L253
  # trk = data.frame(trk, check.names = F, stringsAsFactors = F)
  
  message("reading mutect ...")
  df_mutect = gm_mutect.read(trk, 
                             col_fl = "mutect_fl",
                             col_samp = "NAME")
  message("reading pindel ...")
  df_pindel = gm_pindel.read(trk, 
                             col_fl = "pindel_fl",
                             col_samp = "NAME")
  
  df_variants = bind_rows(df_mutect, df_pindel)
  
  # create well annotated bed, judgement is always KEEP
  # message("\ncreating uniq bed ...")
  # df_variants_bed = dplyr::select(df_variants, chr, start, end, ref_allele, alt_allele, 
  #                               context, key:entrez_gene_id, 
  #                               # should add aaannotation
  #                               aaannotation) %>% unique()
  # write_tsv(df_variants_bed, "ssm/df_mutect_pindel_uniq.bed")
  
  # convert to VCF -----
  source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/mutect_ann_to_maf.R')
  # debug(mutect_ann_to_maf)
  combined_maf = "ssm/merged_variants.maf"
  combined_vcf = "ssm/merged_variants.vcf"
  maf = mutect_ann_to_maf(df_variants, 
                          gene = "gene", 
                          func = "exonicfunc_knowngene", 
                          aa_change = "aaannotation", 
                          sample_name = "samplename",
                          ref_name = "ref_name",)
  write_tsv(maf, combined_maf)
  
  # convert to VCF
  maf2vcf_exe = "/rsrch3/home/iacs/sseth/.conda/envs/pcgr_py36/bin/maf2vcf.pl"
  cmd_maf2vcf = glue("{maf2vcf_exe} --input-maf {combined_maf} --output-dir ssm --per-tn-vcfs --ref_fasta {ref_fasta}")
  cmd_maf2vcf;
  run_cmd(cmd_maf2vcf, target = combined_vcf, cmdname = "maf2vcf")
  
  # somatic_calling_mode = detect_somatic_calling_mode(trk)
  
  # call mutect2 on all samples
  # get final calls for ALL samples!
  
  
  
  
  
}


if(FALSE){
  
  p_load(tidylog)
  outpath="/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex_b1/ss_cl_het/185_145"
  setwd(outpath)
  
  # mutectpath = "/rsrch1/iacs/iacs_dep/sseth/projects2/ss_tnbc/data/artemis/wex/2019_b1/mutect"
  # pindelpath = "/rsrch1/iacs/iacs_dep/sseth/projects2/ss_tnbc/data/artemis/wex/2019_b1/pindel"
  
  trk = read_tsv("gatkcnv/df_trk.tsv")
  trk %<>% mutate(
    # add mutect/pindel calls
    mutect_fl = glue("ssm/{NAME}___185_145_GB-D_mutect.annotated.tsv"), 
    pindel_fl = glue("ssm/{NAME}___185_145_GB-D_pindel.annotated.tsv")
  ) # %>% filter(file.exists(mutect_fl), file.exists(pindel_fl))
  
  
}


phylowgs_parse_cnvs <- function(){
  
  # ./parse_cnvs.py -f titan -c 0.81 --cnv-output cnvs1.txt samp1_segs.txt
  # ./parse_cnvs.py -f titan -c 0.74 --cnv-output cnvs2.txt samp2_segs.txt
  
  # ./parse_cnvs.py -f titan -c 0.81 cnv_calls_segs.txt
  # ./create_phylowgs_inputs.py --cnvs sample1=cnvs.txt --vcf-type sample1=vardict sample1=sample.vcf
  # Valid VCF types are strelka,mutec
  # t_pcawg,dkfz,muse,vardict,mutect_smchet,mutect_tcga,sa
  # nger,pcawg_consensus.
  
  
}


phylowgs.titancna <- function(trk, samplename,
                     vcf_type = opts_flow$envir$phylowgs_vcf_type,
                     singularity_bind_opts = opts_flow$envir$singularity_bind_opts,
                     pwgsdir = pipe_str$phylowgs$dir,
                     cores = resources$pwgs$cpu,
                     # normal_and_abnormal_cn
                     input_regions = opts_flow$envir$phylowgs_input_regions
                     # pwgsdir = pipe_
                     ){
  
  check_args()
  phylowgs_sif = "$HOME/singularity-images/phylowgs.sif"
  
  # read in titancna trk
  # ./parse_cnvs.py -f titan -c 0.81 cnv_calls_segs.txt
  # ./create_phylowgs_inputs.py --cnvs sample1=cnvs.txt --vcf-type sample1=vardict sample1=sample.vcf
  # convert to phylowgs inputs:
  
  source('~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/metadata_dnaseq.R')
  
  if(nrow(trk) == 0)
    stop("no rows in trk")

  trk = metadata_for_dnaseq_tools(trk)
  expect_columns(trk, c("outprefix", "outprefix_paired",
                        "titan_optimal_cluster_file",
                        "titan_phylowgs_input",
                        "filtered_vcf"))
  trk = filter(trk, normal == "NO")

  titan_optimal_cluster_file = trk$titan_optimal_cluster_file
  titan_phylowgs_input = trk$titan_phylowgs_input
  filtered_vcf = trk$filtered_vcf
  cmd_titan_vars = glue("eval $(tsv-to-env {titan_optimal_cluster_file} titancna)")
  parse_cnvs_py="$HOME/Dropbox/public/github_phylowgs/parser/parse_cnvs_ss.py"
  cmd_cnv = glue("mkdir -p {pwgsdir};singularity run {singularity_bind_opts} {phylowgs_sif} python {parse_cnvs_py} -f titan -c $titancna_purity --cnv-output {titan_phylowgs_input} $titancna_path.segs.txt")
  cmd_cnv
  
  # python ${create_phylowgs_inputs_py}
  # vcf="../../../ssm/m1_m2_ir/filtered/WEX-1004-T___matched_mrgd_pcgr_edit.vcf"
  # vcf2="../../../ssm/m1_m2_ir/filtered/not_annotated/WEX-1004-T___matched_combvcf.vcf.gz"
  # zcat ../../../ssm/m1_m2_ir/filtered/not_annotated/WEX-1004-T___matched_combvcf.vcf.gz | grep -v "##" | wc -l
  # [--regions {normal_cn,normal_and_abnormal_cn,all}]
  # this may remove variants which were not in any of the segments, so we are left with 403 variants
  create_phylowgs_inputs_py="$HOME/Dropbox/public/github_phylowgs/parser/create_phylowgs_inputs_ss.py"
  cmd_inp = glue("singularity run {singularity_bind_opts} {phylowgs_sif} python {create_phylowgs_inputs_py} ",
                 "--cnvs S1={titan_phylowgs_input} ",
                 "--vcf-type S1={vcf_type} S1={filtered_vcf} ",
                 "--output-cnvs {pwgsdir}/cnv_data.txt --output-variants {pwgsdir}/ssm_data.txt ",
                 "--regions {input_regions} --verbose")
  cmd_inp
  # wc -l phylowgs/ssm_data.txt
  # wc -l phylowgs/cnv_data.txt
  
  # for some weird reason, the opt file works well
  multievolve_py="/opt/phylowgs/multievolve.py"
  cmd_phylo = glue("singularity run {singularity_bind_opts} {phylowgs_sif} python {multievolve_py} --num-chains {cores} --ssms {pwgsdir}/ssm_data.txt --cnvs {pwgsdir}/cnv_data.txt -O {pwgsdir}/chains")
  cmd_phylo
  
  # results
  write_results_py="/opt/phylowgs/write_results.py"
  cmd_results = glue("mkdir -p {pwgsdir}/results;singularity run {singularity_bind_opts} {phylowgs_sif} python {write_results_py} {samplename} {pwgsdir}/chains/trees.zip ",
                     "{pwgsdir}/results/{samplename}.summ.json.gz {pwgsdir}/results/{samplename}.muts.json.gz {pwgsdir}/results/{samplename}.mutass.zip")
  cmd_results
  
  # Multi-modal solutions won't be properly characterized by the single solution reported by this code. This code assumes that there is a dominant mode amongst the sampled solutions.
  cmd_report = glue("mkdir -p {pwgsdir}/report;PYTHONPATH=/opt/phylowgs singularity run {singularity_bind_opts} {phylowgs_sif} python /opt/smchet-challenge/create-smchet-report/write_report.py ",
                    "{pwgsdir}/results/{samplename}.summ.json.gz ",
                    "{pwgsdir}/results/{samplename}.muts.json.gz ",
                    "{pwgsdir}/results/{samplename}.mutass.zip {pwgsdir}/report")
  cmd_report
  
  # a brown script
  create_phylowgs_output_py="$HOME/Dropbox/public/github_phylowgs/parse_pwgs_outputs.py"
  # usage: parse_pwgs_outputs.py [-h] --cnv_input CNV_INPUT --ssm_input SSM_INPUT
  # --summary_file SUMMARY_FILE --mutasgn_path
  # MUTASGN_PATH --output_folder OUTPUT_FOLDER
  # [--k K]
  cmd_output = glue("mkdir -p outputs;singularity run {singularity_bind_opts} {phylowgs_sif} python {create_phylowgs_output_py} ",
                    "--ssm_input {pwgsdir}/ssm_data.txt --cnv_input {pwgsdir}/results/cnv_data.txt ",
                    "--summary_file {pwgsdir}/results/{samplename}.summ.json.gz --mutasgn_path {pwgsdir}/results/{samplename}.mutass.zip ",
                    " --output_folder outputs")
  # cmd_output
  
  cmds = c(cmd_titan_vars, cmd_cnv, cmd_inp, 
           cmd_phylo,
           cmd_results, cmd_report)
  cmds = list(phylowgs = cmds)
  flowmat = to_flowmat(cmds, samplename)  
  
  outfiles = list(trees = "{pwgsdir}/chains/trees.zip",
                  results = "{pwgsdir}/results",
                  report = "{pwgsdir}/report")
  
  list(flowmat = flowmat, outfiles = outfiles)
  
}  

calc_u.ccf <- function(ccf, m){
  ccf*m
}

calc_ccf.u <- function(u, m){
  return (u/m)
}

calc_m.u <- function(u){
  if(u >=1)
    m = round(abs(u), 0)
  if(u < 1)
    m = 1
  m
}


calc_u  <- function(f, p, cn){
  (f * 1/p) * (p*cn + (1-p)*2)
  
}

phylowgs_parse_report.single_samp <- function(path = "/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/runs/1004/tmp/phylowgs/tesla", rds){

  # list.files(path)
  # setwd(path)
  
  #  Inferred cellularity of the sample
  purity = read.table(glue("{path}/report/1A.txt", header = F)) %>% as.numeric();purity
  # number of subclones
  subclones = read.table(glue("{path}/report/1B.txt"), header = F) %>% as.numeric();subclones
  
  # https://github.com/morrislab/phylowgs/issues/42
  # https://cloud.githubusercontent.com/assets/1753326/16813882/84cf7838-4901-11e6-861f-e5627f5c5adf.png
  # CCF: = cellular prevalence of node + cellular prevalence of children
  df_size = read.table(glue("{path}/report/1C.txt"), col.names = c("clone", "snv", "phi"));
  
  # ssm map
  df_ssm_map = read.table(glue("{path}/report/2A.txt"));df_ssm_map
  mat_co_cluster = read.table(glue("{path}/report/2B.txt.gz"));dim(mat_co_cluster)
  
  df_tree1 = read.table(glue("{path}/report/3A.txt"), col.names = c("clone", "parent")) %>% 
    left_join(df_size, by = "clone") %>% as_tibble()
  # 3B.txt: NxN Ancestor-decedent matrix. Entry i,j = The probability that i is in node that is an ancestor of node containing j.
  mat_ances = read.table(glue("{path}/report/3B.txt.gz"));dim(mat_ances)
  
  # read ssm-data
  ssm_data = read_tsv(glue("{path}/ssm_data.txt"));dim(ssm_data)
  ssm_data$clone = df_ssm_map$V1
  ssm_data$purity = purity
  ssm_data$subclones = subclones
  ssm_data %<>% left_join(df_tree1, by = "clone")
  
  # ** read annotated tsv --------------
  # p_load(VariantAnnotation, glue)
  df_mut1 = read_rds(rds)
  df_mut1$df_mut_ann_filt %<>% mutate(key2 = glue("{chrom}_{pos}"))
  setdiff(ssm_data$gene, df_mut1$df_mut_ann_filt$key2)
  # extra in ssm
  # "7_12376499"
  setdiff(df_mut1$df_mut_ann_filt$key2, ssm_data$gene)
  # missing in SSM
  # [1] "1_5925156"
  
  # add mutation info
  ssm_data = ssm_data %>%
    left_join(dplyr::select(df_mut1$df_mut_ann_filt, -id), by = c("gene" = "key2"))

  list(ssm_data = ssm_data, df_tree1 = df_tree1)
}


# ** summarize loop -------
if(FALSE){
  wexpath = "/rsrch3/scratch/iacs/sseth/flows/SS/sarco/mda/wex"
  runpath = glue("{wexpath}/runs")
  # pwgsdir = "phylowgs2"
  pwgsdir = "phylowgs"
  
  # read manually curated file
  df_mut_curated = read_sheet("/rsrch3/scratch/iacs/sseth/flows/SS/sarco/mda/wex/ssm/m1_m2_ir/filtered/df_mut_ssm_f2_edit.xlsx")
  df_trk = read_tsv(glue("{runpath}/trk_ann.tsv"))
  df_trk %<>% 
    mutate(pwgs_path = glue("{runpath}/{individual}/tmp/{pwgsdir}"),
           pwgs_res_fl = glue("{pwgs_path}/report/1A.txt"),
           filtered_vcf_rds = glue("{wexpath}/ssm/m1_m2_ir/filtered/{outprefix_paired}_mrgd_pcgr.rds")) %>% 
    filter(normal == "NO")
  # pwgs complete
  df_trk = df_trk %>% filter(file.exists(pwgs_res_fl))# %>% 
    #filter(!individual %in% c(2016, 768));dim(df_trk)
  df_trk$individual
  # file.exists(df_trk$filtered_vcf_rds)
  
  # add tree.zip
  source('~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/phylowgs.R')
  # i=2
  for(i in 1:nrow(df_trk)){
    
    ind = df_trk$individual[i]
    path = df_trk$pwgs_path[i]
    rds = df_trk$filtered_vcf_rds[i]
    summpath = file.path(df_trk$pwgs_path[i], "summary")
    wranglr::mkdir(summpath)
    
    message("working on:, ", df_trk$individual[i])
    
    lst = phylowgs_parse_report.single_samp(path, rds = rds)
    lst$ssm_data %>% dim()
    # lst$ssm_data %>% filter(f2)
    write_rds(lst, glue("{summpath}/{ind}_lst_ssm.rds"))
    try(phylowgs_plots(lst, summpath, ind = ind, df_mut_curated = df_mut_curated))
    
  }
  system(glue("rsync -avP /rsrch3/scratch/iacs/sseth/flows/SS/sarco/mda/wex/runs/*/tmp/{pwgsdir}/summary/ /rsrch3/scratch/iacs/sseth/flows/SS/sarco/mda/wex/{pwgsdir}/summary/"))
  
  
  system("rsync -avP /rsrch3/scratch/iacs/sseth/flows/SS/sarco/mda/wex/runs/*/tmp/phylowgs/results/ /rsrch3/home/iacs/sseth/apps/phylowgs_witness/data/sarco-mda-v1.1/")
  system("rsync -avP /rsrch3/scratch/iacs/sseth/flows/SS/sarco/mda/wex/runs/*/tmp/phylowgs2/results/ /rsrch3/home/iacs/sseth/apps/phylowgs_witness/data/sarco-mda-v1.2/")
  
  
}

phylowgs_plots <- function(lst, summpath, ind, df_mut_curated){
  
  wranglr::mkdir(summpath)
  p_load(clonevol)
  
  # ind = lst$ssm_data$ind[1]
  # cnv is missing in this report
  message("getting F2 mutations")
  df_f2 = filter(lst$ssm_data, f2) %>% 
    dplyr::select(ind, name,
                  in_ir, in_m1, in_m2,
                  gene, symbol, protein_change, 
                  taf, naf, tdp, ndp, 
                  clone, parent, phi, f1, f2, cds_change, 
                  consequence, existing_variation, fathmm_dbnsfp, fathmm_mkl_dbnsfp, 
                  impact, intogen_driver, intogen_driver_mut, 
                  mutationassessor_dbnsfp,	mutationtaster_dbnsfp,	mutation_hotspot,	mutation_hotspot_cancertype,	
                  tier, tier_description, mutation_hotspot_transcript,
                  everything()) %>% 
    dplyr::filter(consequence != "synonymous_variant") %>% 
    arrange(-phi)
    write_sheet(df_f2, glue("{summpath}/{ind}_df_mut_f2.xlsx"))
    
    i = ind
    dim(df_f2)
    df_mut_curated_ind = filter(df_mut_curated, ind == i, in_tree == 1);dim(df_mut_curated_ind)
    df_f2.2 = df_f2 %>% filter(key %in% df_mut_curated_ind$key);dim(df_f2.2)
    
    message("clone plot")
    p_load(cowplot, ggplot2, ggsci, forcats)
    p_clone_v1 = df_f2.2 %>% 
      mutate(gene_change = glue("{symbol} {protein_change}"),
             gene_change = fct_reorder(gene_change, -phi)) %>% 
      ggplot(aes(gene_change, phi, fill = consequence)) +
      geom_col() + 
      theme_cowplot() +
      scale_fill_npg() + ylab("cellularity") + xlab("") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      ggtitle(ind)
    save_plot(glue("{summpath}/{ind}_p_clone_v1.pdf"), p_clone_v1)
    #
    
    treefl = tempfile()
    variantfl = tempfile()
    df_tree1 = lst$df_tree1 %>% 
      mutate(sample.with.nonzero.cell.frac.ci = round(phi, 2),
             parent = ifelse(parent == 0, -1, parent))
    write_tsv(df_tree1, treefl)
    df_vars = lst$ssm_data %>% 
      mutate(is.driver = key %in% df_f2.2$key) %>% 
      dplyr::select(cluster = clone, is.driver, gene = symbol)
    # add cnv for missing snv
    missing_clusters = setdiff(df_tree1$clone, df_vars$cluster)
    if(length(missing_clusters) > 0){
      df_vars = data.frame(cluster = missing_clusters, is.driver = T, gene = "cnv event") %>% 
        bind_rows(df_vars)
    }
    write_tsv(df_vars, variantfl)

    # df_f2.2
    
    message("tree plot")
    y = import.tree(treefl, variantfl)
    if(sum(y$variants$is.driver) > 0)
      y <- transfer.events.to.consensus.trees(y, 
                                              events = y$variants[y$variants$is.driver,], 
                                              cluster.col.name = 'cluster', event.col.name = 'gene')
    y = convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')
    # clonevol::plot.tree.clone.as.branch(y)
    # clonevol:::plot.clonal.models(y, out.dir = "~/Downloads/clone")
    # plot.all.trees.clone.as.branch(y)
    
    # ** plot tree  -------
    pdf(glue("{summpath}/{ind}_p_tree_clonevol_v1.pdf"))
    node_sz = lst$df_tree1$phi*8 + 1;node_sz
    plot.all.trees.clone.as.branch(y, 
                                   branch.width = 1,
                                   
                                   # ccf label
                                   node.text.size = 0,
                                   
                                   # events
                                   branch.text.size = 1,
                                   
                                   node.size = node_sz,
                                   node.label.size = 1,
                                   
                                   tree.label = ind, main = ind
                                   )
    
    dev.off()
}






# END













# END
