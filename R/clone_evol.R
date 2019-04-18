

cloneevol_r <- function(pyclone_run_path = "~/.rsrch2/iacs/iacs_dep/sseth/flows/SS/tnbc/ms51_wex_b1/pyclone/185_198"){
  pyclone_run_path = "~/.rsrch2/iacs/iacs_dep/sseth/flows/SS/tnbc/ms51_wex_b1/pyclone/185_057"
  setwd(pyclone_run_path)
  
  #devtools::install_github("hdng/clonevol")
  library(pacman)
  p_load(tidyverse, reshape2)
  p_load(clonevol, timescape)
  
  message("read loc ann")
  df_loci_ann_f = read_rds("tables/loci_ann_f.rds") %>% 
    tbl_df()
  
  # create wide format; cloneevol ------
  df_ccf = dplyr::rename(df_loci_ann_f, 
                         cluster = cluster_id, 
                         ccf = cellular_prevalence) %>% 
    mutate(cluster = as.numeric(cluster), 
           ccf = ccf*100) %>% 
    #dplyr::count(mutation_id, cluster, gene, aachange, sample_group) %>% arrange(-n)
    dcast(mutation_id + cluster + gene + aachange ~ sample_id, value.var = "ccf")
  df_ccf = arrange(df_ccf, cluster)
  
  df_vaf = dplyr::rename(df_loci_ann_f, 
                         cluster = cluster_id, 
                         vaf = variant_allele_frequency) %>% 
    mutate(cluster = as.numeric(cluster), 
           vaf = vaf*100) %>% 
    #dplyr::count(mutation_id, cluster, gene, aachange, sample_group) %>% arrange(-n)
    dcast(mutation_id + cluster + gene + aachange ~ sample_id, value.var = "vaf")
  head(df_vaf)
  df_vaf = arrange(df_vaf, cluster)
  #df_vaf$cluster
  
  sample.groups = df_loci_ann_f$sample_id %>% unique() %>% sort()
  names(sample.groups) = sample.groups
  vaf.col.names = sample.groups
  
  
  # https://github.com/hdng/clonevol/issues/7
  # Sorry for the late reply. This is a known bug when your cluster ID is not consecutive and 
  # founding clone is not 1. I'll release a fix soon. 
  # The workaround is to add the following code before infer.clonal.models, to make sure founding clone is 1:
  # infer.clonal.models ------
  message("infer.clonal.models")
  clone.colors <- c('#999793', '#8d4891', '#f8e356', '#fe9536', '#d7352e')
  y = infer.clonal.models(variants = df_ccf,
                          cluster.col.name = 'cluster',
                          ccf.col.names = vaf.col.names,
                          #vaf.col.names = vaf.col.names,
                          sample.groups = sample.groups,
                          cancer.initiation.model = 'monoclonal',
                          subclonal.test = 'bootstrap',
                          subclonal.test.model = 'non-parametric',
                          num.boots = 1000,
                          founding.cluster = 1,
                          cluster.center = 'mean',
                          ignore.clusters = NULL,
                          #clone.colors = clone.colors,
                          min.cluster.vaf = 0.01,
                          # min probability that CCF(clone) is non-negative
                          sum.p = 0.05,
                          # alpha level in confidence interval estimate for CCF(clone)
                          alpha = 0.05)
  
  
  # get drivers
  # y <- transfer.events.to.consensus.trees(y,
  #                                         x[x$is.driver,],
  #                                         cluster.col.name = 'cluster',
  #                                         event.col.name = 'gene')
  
  y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')
  
  
  tmp = plot.clonal.models(y,
                           # box plot parameters
                           box.plot = TRUE,
                           fancy.boxplot = TRUE,
                           fancy.variant.boxplot.highlight = 'is.driver',
                           fancy.variant.boxplot.highlight.shape = 21,
                           fancy.variant.boxplot.highlight.fill.color = 'red',
                           fancy.variant.boxplot.highlight.color = 'black',
                           fancy.variant.boxplot.highlight.note.col.name = 'gene',
                           fancy.variant.boxplot.highlight.note.color = 'blue',
                           fancy.variant.boxplot.highlight.note.size = 2,
                           fancy.variant.boxplot.jitter.alpha = 1,
                           fancy.variant.boxplot.jitter.center.color = 'grey50',
                           fancy.variant.boxplot.base_size = 12,
                           fancy.variant.boxplot.plot.margin = 1,
                           fancy.variant.boxplot.vaf.suffix = '.VAF',
                           # bell plot parameters
                           clone.shape = 'bell',
                           bell.event = TRUE,
                           bell.event.label.color = 'blue',
                           bell.event.label.angle = 60,
                           clone.time.step.scale = 1,
                           bell.curve.step = 2,
                           # node-based consensus tree parameters
                           merged.tree.plot = TRUE,
                           tree.node.label.split.character = NULL,
                           tree.node.shape = 'circle',
                           tree.node.size = 30,
                           tree.node.text.size = 0.5,
                           merged.tree.node.size.scale = 1.25,
                           merged.tree.node.text.size.scale = 2.5,
                           merged.tree.cell.frac.ci = FALSE,
                           # branch-based consensus tree parameters
                           merged.tree.clone.as.branch = TRUE,
                           mtcab.event.sep.char = ',',
                           mtcab.branch.text.size = 1,
                           mtcab.branch.width = 0.75,
                           mtcab.node.size = 3,
                           mtcab.node.label.size = 1,
                           mtcab.node.text.size = 1.5,
                           # cellular population parameters
                           cell.plot = TRUE,
                           num.cells = 100,
                           cell.border.size = 0.25,
                           cell.border.color = 'black',
                           clone.grouping = 'horizontal',
                           #meta-parameters
                           scale.monoclonal.cell.frac = TRUE,
                           show.score = FALSE,
                           cell.frac.ci = TRUE,
                           disable.cell.frac = FALSE,
                           # output figure parameters
                           out.dir = 'output',
                           out.format = 'pdf',
                           overwrite.output = TRUE,
                           width = 14,
                           height = 8,
                           # vector of width scales for each panel from left to right
                           panel.widths = c(3,4,2,4,2))
  
  # plots  ------
  # ** fishplot -------
  library(fishplot)
  f = generateFishplotInputs(results=y)
  
  fishes = createFishPlotObjects(f)
  #plot with fishplot
  pdf('fish.pdf', width=8, height=5)
  #i=1
  for (i in 1:length(fishes)){
    fish = layoutClones(fishes[[i]])
    fish = setCol(fish,f$clonevol.clone.colors)
    fishPlot(fish, shape="spline", title.btm=partid, cex.title=0.5,
             vlines=seq(1, length(sample.groups)), vlab=sample.groups, pad.left=0.5)
  }
  dev.off()
  
  
  # ** timescape -------
  df_tree_edges = y$matched$merged.trees[[1]] %>% 
    dplyr::select(source = parent, target = lab) %>% dplyr::filter(source != -1)
  df_tree_edges
  
  df_clonal_prev = select(df_loci_ann_f, timepoint = sample_group, clone_id = cluster_id, clonal_prev = cellular_prevalence) %>% 
    unique() %>% group_by(timepoint, clone_id) %>% summarise(clonal_prev = mean(clonal_prev)) %>% 
    arrange(timepoint)
  df_clonal_prev
  
  timescape(df_clonal_prev, df_tree_edges)
  
  group_by(df_clonal_prev, timepoint) %>% 
    summarise(clonal_prev = sum(clonal_prev))
  
  # totals to 1
  group_by(clonal_prev, timepoint) %>% 
    summarise(clonal_prev = sum(clonal_prev))
  
  # extract clone phylogeny
  
  # try timescape:
  p_load(timescape)
  
  # clonal_prev  is a data frame consisting of clonal prevalences for each clone at each time point. The columns in this data frame are:
  # character() timepoint - time point
  # character() clone_id - clone id
  # numeric() clonal_prev - clonal prevalence.
  
  # tree_edges is a data frame describing the edges of a rooted clonal phylogeny. The columns in this data frame are:
  # character() source - source node id
  # character() target - target node id.
  
  
}


if(FALSE){
  library(clonevol)
  
  
  # pyclone_run_path = "~/rsrch2_home/tmp/test_clone_conv/pyclone/run_v1/185_057"
  # 
  # 
  # min_clust_size = 7
  
  # df_clust = read_tsv(glue("{pyclone_run_path}/tables/cluster.tsv")) %>% 
  #   tbl_df()
  # df_loci = read_tsv(glue("{pyclone_run_path}/tables/loci.tsv")) %>% 
  #   tbl_df()
  # df_loci_ann = read_tsv(glue("{pyclone_run_path}/tables/loci_ann.tsv")) %>% 
  #   tbl_df()
  # df_loci_ann_f = dplyr::filter(df_loci_ann, size > min_clust_size) %>% 
  #   separate(sample_id, c("trial", "partid", "sample_group", "suffix"), "_", extra = "merge", remove = F)
  # df_clust_f = dplyr::filter(df_clust, size > min_clust_size)
  
  
  # prepare data -------
  # variants:
  # data frame of the variants.
  # At least cluster column and VAF or CCF columns are required.
  # Cluster column should contain cluster identities as continuous integer values starting from 1.
  # Either c or variants parameter is required.
  # 
  # founding.cluster: not sure?
  # 
  # ignore.clusters: can ask it to ignore ?
  # 
  # cancer.initiation.model
  
  # cluster 1, 2, 3....
  # df_var = dplyr::select(df_clust_f, 
  #                        cluster = cluster_id, 
  #                        ccf = mean, 
  #                        sample_group = sample_id) %>% 
  #   mutate(cluster = factor(cluster), 
  #          cluster = as.numeric(cluster)) %>% 
  #   dcast(cluster ~ sample_group, value.var = "ccf")
  # head(df_var)
  
  # order according to max prevelance, at T0
  # df_ccf_wd = dcast(df_clust_f, cluster_id ~ sample_id, value.var = "mean")
  # df_ccf_wd = df_ccf_wd[order(df_ccf_wd[, 2], decreasing = T), ]
  # df_ccf_wd
  
  # https://github.com/BradnerLab/TNBC
  
  
  
  

  
}



# END