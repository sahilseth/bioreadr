

# https://bitbucket.org/dranew/citup/issues/4/install-issue-exception-no-submit-queue
# /risapps/rhel6/citup/anaconda2-4.2/bin/run_citup_qip.py freq.txt clusters.txt results.h5

# install --------

# https://shahlab.ca/projects/citup/

# ** new conda env 2.7 --------
# module load conda_/2.7
# conda create --name citup



# ** cplex (final) --------
# student website
# https://ibm.onthehub.com/WebStore/OfferingDetails.aspx?o=733c3d21-0ce1-e711-80fa-000d3af41938&pmv=00000000-0000-0000-0000-000000000000
# installing cplex using the recipe:
# cd ~/apps/dranew_conda_recipes/cplex
# mkdir src;cd src
# cp -r ~/apps/cplex/12.7.1/ .
# nano build.sh
# nano meta.yaml

# https://github.com/spacetelescope/cubeviz/wiki/Testing-Conda-Recipes-Locally
# conda install conda-build
# conda build cplex
# /rsrch2/iacs/iacs_dep/sseth/apps/conda/2.7/envs/citup/conda-bld/cplex_1551541229318/work/tmp
# # conda install --use-local cplex, not sure if this reqd
# try this: (based on: https://github.com/conda/conda/issues/7024)
# conda install -c ~/apps/conda/2.7/envs/citup/conda-bld cplex


# ** install citup --------
# conda install -c dranew citup




# conda config --add channels http://conda.anaconda.org/dranew
# conda install citup
# for some reason, this does not seem to work.
# conda install -c dranew citup 



# try using python 2.7, instead:
# conda install -c dranew citup 
# conda install -c dranew cplex

# try using a clean conda 3.6 session
# module load conda_/3.6
# cplex installs, and it seems it automatically install py3.6
# which may be messing up things!!
# conda install -c ibmdecisionoptimization cplex
# this fails in 3.6, i gues it needs 2.7
# conda install -c dranew citup


# FINAL try ------!
# conda install -c ibmdecisionoptimization cplex
# conda install -c dranew citup
# citup-0.1.1-py27_1


# TEST -----
# run_citup_qip.py freq.txt clusters.txt results.h5 --submit local
# STILL FAILS:
# which run_citup_qip.py
# ~/apps/conda/2.7/envs/citup/bin/run_citup_qip.py

# ATTEMPT #4 -------
# conda create --name citup6 python=2.7
# conda activate citup6
# conda install -c ibmdecisionoptimization cplex
# conda install citup
# conda list citup
# packages in environment at /home/pulintz/miniconda2/envs/citup6:
# #
# citup                     0.1.1                    py27_1    dranew
# run_citup_iter.py --submit local --max_nodes 5 freq.txt clusters.txt results.h5
# min_nodes: 1, max_nodes: 5
# log file: ./tmp/log/20180803-124215/pipeline.log

# ** test -----
# (citup6) [sseth@sharklogin2 ~]$ run_citup_qip.py freq.txt clusters.txt results.h5 --submit local
# it works!!!!


# Citup Iter
# To run the iterative version of citup, provide a table of frequencies. The input format is tab/whitespace separated, with each row
# representing the frequency of a mutation and each column is a tumour sample. No header is required.
# 
# For example, the following would be input for 3 mutations in 2 samples.
# 
# 0.2 0.1
# 0.4 0.3
# 0.5 0.1


# Citup QIP
# To run the QIP version of citup, provide a table of mutation frequencies, and a table of mutation clusters.
# The input format for mutation frequencies is described above. The mutation clusters is a single line per mutation,
# containing a 0 based integer cluster index for that mutation. For example, the following specifies 3 mutations, the
# first two in the same cluster the last in a different cluster.

# we would use these;
# now the question is use, freq of mutation, or CCF of the mutation??

# could try with both, actually.


# OUTPUT:
# Output Format
# The output of citup is pandas hdf5 format.
# 
# The results store contains the following pandas series, with an entry per tree solution:
# 
# /results/bic
# /results/error_rate
# /results/likelihood
# /results/num_mutations
# /results/num_nodes
# /results/num_samples
# /results/objective_value
# /results/optimal
# /results/tree_id
# /results/tree_index
# /results/tree_string
# The store also contains the following pandas series, with one series per tree solution:
# 
# /trees/{tree_solution}/cluster_assignment (QIP only)
# /trees/{tree_solution}/variant_assignment (Iter only)
# /trees/{tree_solution}/objective_value
# Finally, the store contains the following pandas dataframes, with one frame per tree solution:
# 
# /trees/{tree_solution}/adjacency_list
# /trees/{tree_solution}/clade_freq
# /trees/{tree_solution}/objective_value
# Finally, the store contains the following pandas dataframes, with one frame per tree solution:
# 
# /trees/{tree_solution}/adjacency_list
# /trees/{tree_solution}/clade_freq
# /trees/{tree_solution}/clone_freq
# /trees/{tree_solution}/gamma_matrix
# Recent activity 
# Repo activity is currently unavailable.



citup_prep_input <- function(pyclone_path){
  
  setwd(pyclone_path)
  
  # read in combined df
  df_loci_ann = read_tsv("tables/loci_ann.tsv", col_types = cols(.default = col_character())) %>% 
    mutate(cellular_prevalence = as.numeric(cellular_prevalence), 
           size = as.integer(size), 
           cellular_prevalence = round(cellular_prevalence, 2),
           variant_allele_frequency = as.numeric(variant_allele_frequency),
           sample_id = gsub("_pyclone_inp", "", sample_id))
  
  # ** simple filtering -----
  df_loci_ann_f = filter(df_loci_ann, size > 7)
  
  # conv to wide
  df_vaf = reshape2::dcast(df_loci_ann_f, key + cluster_id ~ sample_id, value.var = "variant_allele_frequency")
  df_ccf = reshape2::dcast(df_loci_ann_f, key + cluster_id ~ sample_id, value.var = "cellular_prevalence")
  
  # 0.2 0.1
  # 0.4 0.3
  # 0.5 0.1
  
  mat_vaf = dplyr::select(df_vaf, -cluster_id, -key)
  mat_ccf = dplyr::select(df_ccf, -cluster_id, -key)
  mat_clus = dplyr::select(df_ccf, cluster_id)
  
  write_tsv(mat_vaf, "mat_vaf.tsv", col_names = F)
  write_tsv(mat_ccf, "mat_ccf.tsv", col_names = F)
  write_tsv(mat_clus, "mat_clus.tsv", col_names = F)
  
}

#' pyclone_citup
#' 
#' this step would use 16 cores
#'
#' @param pyclone_path 
#' @param citup_exe 
#' @param citup_params 
#'
#' @export
pyclone_citup <- function(pyclone_path, 
                          citup_exe = "module load conda_/2.7;source activate citup;run_citup_qip.py",
                          citup_params = "--submit local --loglevel DEBUG --maxjobs 8"
                          ){
  
  citup_prep_input(pyclone_path)
  
  odir_vaf = "citup_vaf";wranglr::mkdir(odir_vaf)
  odir_ccf = "citup_ccf";wranglr::mkdir(odir_ccf)
  
  # cmd = "/risapps/rhel6/citup/anaconda2-4.2/bin/run_citup_qip.py mat_ccf.tsv mat_clus.tsv results.h5 --submit local"
  # cmd = "module load conda_/2.7;source activate citup;run_citup_qip.py mat_ccf.tsv mat_clus.tsv results.h5 --submit local"
  # system(cmd)
  
  cmd_vaf = glue("cd {odir_vaf};{citup_exe} ../mat_vaf.tsv ../mat_clus.tsv results.h5 {citup_params} > citup.log 2>&1")
  cmd_ccf = glue("cd {odir_ccf};{citup_exe} ../mat_ccf.tsv ../mat_clus.tsv results.h5 {citup_params} > citup.log 2>&1")

  cmds = c(cmd_ccf, cmd_vaf)
  tmp = mclapply(c(cmds), system, mc.cores = 2)
  

}

# if we load the conf, prior to running the cmd, all is set!


parse_citup <- function(fl, verbose = F){
  p_load(reticulate)
  # load some libs
  pd = import("pandas")
  np = import("numpy")
  h5 = import("h5py")
  os = import("os")
  pathlib = import("pathlib")
  
  # source_python()
  
  fl_tree_adj = gsub("results.h5", "tree_adj.tsv", fl)
  fl_clone_freq = gsub("results.h5", "clone_freq.tsv", fl)
  fl_tree = gsub("results.h5", "timescape.html", fl)
  
  pd_arr_opt = pd$read_hdf(fl, key = 'results/optimal')
  if(verbose)
    print(pd_arr_opt)
  opt_tree = names(pd_arr_opt)[1] %>% stringr::str_trim()
  message("> optimal tree: ", opt_tree, ".")
  # pd_arr_opt = pd.read_hdf(fl, key = 'results/optimal')
  # pd_arr_opt.head()
  # get the opt solution:
  # pd_arr_opt.name
  # opt_tree = pd_arr_opt.index[0]
  
  print("> read various pd.dfs... ")
  tree_adj = pd$read_hdf(fl, key = glue("trees/{opt_tree}/adjacency_list"))
  clade_freq = pd$read_hdf(fl, key = glue("trees/{opt_tree}/clade_freq"))
  clone_freq = pd$read_hdf(fl, key = glue("trees/{opt_tree}/clone_freq"))
  gamma_matrix = pd$read_hdf(fl, key = glue("trees/{opt_tree}/gamma_matrix"))
  
  colnames(tree_adj) = c("source", "target")
  
  
  print("> print out ...")
  # The adjacency list can be written as a TSV with the column names  source, target to be input into E-scape
  # source, target
  # write_tsv
  write_tsv(tree_adj, fl_tree_adj)
  
  # clone frequencies should be reshaped such that each row represents a clonal frequency in a 
  # specific sample for a specific clone, with the columns representing the time or 
  # space ID, the clone ID, and the clonal prevalence.
  clone_freq %>% t() %>% wranglr::to_df() %>% 
    write_tsv(fl_clone_freq)
  
}


citup_timescape <- function(fl){
  fl_tree_adj = gsub("results.h5", "tree_adj.tsv", fl)
  fl_clone_freq = gsub("results.h5", "clone_freq.tsv", fl)
  fl_tree = gsub("results.h5", "timescape.html", fl)
  
  
  message("read data")
  df_tree_adj = read_tsv(fl_tree_adj, 
                         col_types = cols(
                           source = col_double(),
                           target = col_double()))
  
  
  df_clone_freq = read_tsv(fl_clone_freq, col_types = cols(.default = col_character()))
  df_clone_freq = df_clone_freq %>% 
    dplyr::rename(clone_id = id) %>% 
    gather(key = "timepoint", "clonal_prev", -clone_id) %>% 
    mutate(clonal_prev = as.numeric(clonal_prev))
  
  
  # gather
  # df_tree_adj
  # df_clone_freq
  
  message("timescape")
  tmp = timescape(df_clone_freq, df_tree_adj)
  htmlwidgets::saveWidget(widget = tmp, file = fl_tree)
}


# END
