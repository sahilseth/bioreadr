



# https://bioconductor.org/packages/release/bioc/vignettes/timescape/inst/doc/timescape_vignette.html

# reqd params:
# The required parameters for TimeScape are as follows:
#   
# clonal_prev is a data frame consisting of clonal prevalences for each clone at each time point. The columns in this data frame are:
#   
#   character() timepoint - time point
#   character() clone_id - clone id
#   numeric() clonal_prev - clonal prevalence.
# 
# tree_edges is a data frame describing the edges of a rooted clonal phylogeny. The columns in this data frame are:
#   
#   character() source - source node id
#   character() target - target node id.

# mutations:
# mutations is a data frame consisting of the mutations originating in each clone. The columns in this data frame are:
#       character() chrom - chromosome number
#       numeric() coord - coordinate of mutation on chromosome
#       character() clone_id - clone id
#       character() timepoint - time point
#       numeric() VAF - variant allele frequency of the mutation in the corresponding timepoint.
#       If this parameter is provided, a mutation table will appear at the bottom of the view.


# E-scape takes as input a clonal phylogeny and clonal prevalences per clone per sample.
# At the time of submission many methods have been proposed for obtaining these values,
# and accurate estimation of these quantities is the focus of ongoing research.
# 
# We describe a method for estimating clonal phylogenies and clonal prevalence using
# PyClone (Roth et al., 2014; source code available at https://bitbucket.org/aroth85/pyclone/wiki/Home) and
# citup (Malikic et al., 2016; source code available at https://github.com/sfu-compbio/citup).
# In brief, PyClone inputs are prepared by processing fastq files resulting from a targeted deep sequencing experiment.
# 
# Using samtools mpileup (http://samtools.sourceforge.net/mpileup.shtml), the number of nucleotides
# matching the reference and non-reference are counted for each targeted SNV.
# Copy number is also required for each SNV.
# 
# We recommend inferring copy number from whole genome or whole exome sequencing of samples taken
# from the same anatomic location / timepoint as the samples to which targeted deep sequencing was applied.
# 
# Copy number can be inferred using Titan (Ha et al., 2014; source code available at https://github.com/gavinha/TitanCNA).
# 
# Sample specific SNV information is compiled into a set of TSV files, one per sample.
# The tables includes mutation id, reference and variant read counts, normal copy number,
# and major and minor tumour copy number (see PyClone readme). PyClone is run on these files using the
# PyClone run_analysis_pipeline subcommand, and produces the tables/cluster.tsv in the working directory.
# 
# Citup can be used to infer a clonal phylogeny and clone prevalences from the cellular prevalences produced by PyClone.
# # so not VAF, rather CCF
# 
# The tables/cluster.tsv file contains per sample, per SNV cluster estimates of cellular prevalence.
# 
# The table is reshaped into a TSV file of cellular prevalences with rows as clusters and columns as samples,
# and the mean of each cluster taken from tables/cluster.tsv for the values of the table.
# The iterative version of citup is run on the table of cellular frequencies, producing an hdf5 output results file.
# (OK we need to run another version of citup, with only cellular freq, and clusters)
#
# cluster   T0      T1     TS
# 
# 
# Within the hdf5 results, the /results/optimal can be used to identify the id of the optimal tree solution.
# 
# The clonal phylogeny as an adjacency list is then the /trees/{tree_solution}/adjacency_list entry and the
# clone frequencies are the /trees/{tree_solution}/clone_freq entry in the hdf5 file.
# 
# The adjacency list can be written as a TSV with the column names source, target to be input into E-scape,
# and the clone frequencies should be reshaped such that each row represents a clonal frequency in a
# specific sample for a specific clone, with the columns representing the time or space ID, the clone ID, and the clonal prevalence.

timescape <- function(){
  library(timescape)
  
  # EXAMPLE 1 - Acute myeloid leukemia patient, Ding et al., 2012
  # genotype tree edges
  tree_edges <- read.csv(system.file("extdata", "AML_tree_edges.csv", 
                                     package = "timescape"))
  
  # clonal prevalences
  clonal_prev <- read.csv(system.file("extdata", "AML_clonal_prev.csv",
                                      package = "timescape"))
  
  # targeted mutations
  mutations <- read.csv(system.file("extdata", "AML_mutations.csv", 
                                    package = "timescape"))
  
  # perturbations
  perturbations <- data.frame( pert_name = c("Chemotherapy"), 
                               prev_tp = c("Diagnosis"))
  
  # run timescape
  timescape(clonal_prev = clonal_prev, tree_edges = tree_edges, 
            perturbations = perturbations, mutations = mutations)
  
}







# END
