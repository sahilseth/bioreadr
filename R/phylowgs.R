# !diagnostics off
# https://stackoverflow.com/questions/39041115/fixing-a-multiple-warning-unknown-column

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






phylowgs_parse_cnvs <- function(){
  
  # ./parse_cnvs.py -f titan -c 0.81 --cnv-output cnvs1.txt samp1_segs.txt
  # ./parse_cnvs.py -f titan -c 0.74 --cnv-output cnvs2.txt samp2_segs.txt
  
  # ./parse_cnvs.py -f titan -c 0.81 cnv_calls_segs.txt
  # ./create_phylowgs_inputs.py --cnvs sample1=cnvs.txt --vcf-type sample1=vardict sample1=sample.vcf
  # Valid VCF types are strelka,mutec
  # t_pcawg,dkfz,muse,vardict,mutect_smchet,mutect_tcga,sa
  # nger,pcawg_consensus.
  
  
}



phylowgs_parse_report.single_samp <- function(path = "~/flows/1004/", rds){

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
