# download reference ----
# To run MutSigCV in this way, please first download the following four reference files:
#   
#   genome reference sequence:   chr_files_hg18.zip    or    chr_files_hg19.zip
# wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/chr_files_hg19.zip
# unzip this file to yield a directory (chr_files_hg18/ or chr_files_hg19/) of chr*.txt files

# mutation_type_dictionary_file.txt
# wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/mutation_type_dictionary_file.txt

# exome_full192.coverage.txt.zip
# wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/exome_full192.coverage.zip

# unzip this file to yield exome_full192.coverage.txt
# gene.covariates.txt
# wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/gene.covariates.txt


# inputs ---------
# MUTATION TABLE
# 
# This is a list of all the mutations to be analyzed, with all patients concatenated together.  It should be a tab-delimited file with a header row.  The positions of the columns are not important, but the header names are.  The program will look up the information it needs based on the names of the headers.  
# 
# The file can be a MAF file, and some of the columns come directly from the standard columns of the MAF format.  However, there are a few additional columns required.
# 
# The four columns required by the program are:
#   
#   "gene" = name of the gene that the mutation was in.  (can also be called "Hugo_Symbol")
# "patient" = name of the patient that the mutation was in.  (can also be called "Tumor_Sample_Barcode")
# "effect" = what broad class of effect does the mutation exert on the gene?  This should be one of "nonsilent" (it changes the protein sequence or splice-sites), "silent" (it is a synonymous change), or "noncoding" (it is intronic or otherwise in a flanking noncoding region.).  (Note, in version1.0.0, this information was split into two columns, "is_coding" and "is_silent", and this method is still supported.)
# "categ" = mutation category.  MutSigCV splits mutations into different categories depending on their DNA context (e.g. was the mutation in a CpG dinucleotide? in another C:G basepair?  in an A:T basepair?).  For each category, there is a different set of "bases at risk" for that category of mutation.  There is also a special category outside this system called "null/indel", into which all truncating mutations and indels are placed.  By convention, the "bases at risk" for a null/indel mutation is taken to be the entire territory of the gene.
# The standard set of categories we have used for many sequencing projects is as follows:
#   
#   1. CpG transitions
# 2. CpG transversions
# 3. C:G transitions
# 4. C:G transversions
# 5. A:T transitions
# 6. A:T transversions
# 7. null+indel mutations
# To compute the "categ" column, information from the reference genome is required: specifically, the identity of the nucleotides directly on each side of the mutation site.  Together with that information and the Variant_Classification, Reference_Allele, and Tumor_Seq_Allele1+2 columns, "categ" can be computed.
# 
# Starting with MutSigCV version 1.3,  an integrated preprocessing module assists with the calculation of these categ numbers, and also enables the automated determination of the optimal set of categories to be used for a given dataset.
# 
# 
# 
# COVERAGE TABLE
# 
# The coverage table tells how many nucleotides were sequenced to adequate depth for mutation calling.  Coverage is tabulated for each patient, and for each gene.  It is also broken down by "categ" and "effect" (as listed above).  Again, "effect" can be either "noncoding" (this refers to the flanking territory outside exons), "nonsilent" (this refers to bases which, when mutated, yield a change in the protein sequence--including splice-sites), or "silent" (this refers to bases which give a synonymous change when mutated).  Note, some coding positions can contribute fractionally to the "silent" and "nonsilent" zones, in a ratio of 1/3 to 2/3 (or vice versa), depending on the consequences of mutating to each of the three possible alternate bases.
# 
# The columns required in the coverage table are: 
#   
#   "gene": the gene name, corresponding to the "gene" column in the mutation_table.
# "effect": whether this row tabulates the "silent", "nonsilent", or "noncoding" territory for this gene.  (Note, in version 1.0.0 column was called "zone", and this name is still allowed.)
# "categ": which category this row tabulates-- should be same as in mutation_table
# <patient_name_1>: number of sequenced bases for patient#1 in this gene and effect/categ bin
# <patient_name_2>: number of sequenced bases for patient#2 in this gene and effect/categ bin
# <patient_name_3>
#   (etc.)
# We recognize that detailed coverage information is not always available.  In such cases, a reasonable approach is to carry out the computation assuming full coverage.  We have prepared a file that can be used for this purpose: it is a "full coverage" file, or more accurately a "territory" file: the only information it contributes is a tabulation of how the reference sequence of the human exome breaks down by gene, categ, and effect.  To download this file, see the section below about MutSigCV v1.3.
# 
# 
# 
# COVARIATES TABLE
# 
# This table lists genomic parameters for each gene being analyzed.  They are called covariates because they co-vary with mutation rate.  They will be used to calculate distances between pairs of genes in a "covariate space" in order to determine the nearest neighbors of each gene, in order to pool information among them about the local background mutation rate (BMR).
# 
# The columns of this file are:
#   
#   "gene": the gene name, should match those used in the first two tables.
# <covariate_name_1>:  the value of the first covariate for each gene
# <covariate_name_2>:  the value of the second covariate for each gene
# <covariate_name_3>:  the value of the third covariate for each gene
# etc.
# The covariates table provided in the Example Data has proven useful for analyzing many cancer types.  The table contains one value per gene for:  (1) global expression, derived from RNA-Seq data and summed across the 91 cell lines in the CCLE (Barretina et al.).  (2) DNA replication time (from Chen et al.).  (3) the HiC statistic, a measure of open vs. closed chromatin state (from Lieberman-Aiden et al.).
# 


# running -----
# If you have a license for Matlab, you can run MutSigCV from its source code file: MutSigCV.m
# Open Matlab and type the following command at the ">>" Matlab prompt.
# MutSigCV('mutations.maf','coverage.txt','covariates.txt','output.txt')
# 
# 
# If you do not have a license for Matlab, you can run the compiled version of MutSigCV using the free Matlab MCR:
#   run_MutSigCV.sh <path_to_MCR> mutations.maf coverage.txt covariates.txt output.txt
# The algorithm will load the three input files, process them using the MutSigCV algorithm, and then finally write the output table to the file 'output.txt'.
# Be sure to replace mutations.maf, etc with the actual paths to the input files.

# mutsig2CV: forum instructions: -------
# https://groups.google.com/a/broadinstitute.org/g/gdac-users/c/CGrKr8ncq6s?pli=1
# Krutika Gaonkar
# unread,
# Jul 27, 2018, 9:56:22 AM
# to beth.bou...@gmail.com, Gdac-users
# Hi Beth,
# I was pointed to download the latest version of MutSig2CV using the link below ( this includes the scripts and the reference files needed including the mutation_type_dictionary_file.txt) by Julian Hess <jh...@broadinstitute.org>. 
# https://archive.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/MutSig2CV.tar.gz
# https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/MutSig2CV.tar.gz
# I download the tar file and was able to install and run the tool successfully using the README file. 
# Hope that helps!
# Thanks,
# Krutika
# The coverage files are no longer required, and the covariates table is baked into the files you downloaded 
# (reference/covariates_transformed.v5a.txt).


# MutSig 1.0, 1.5, 2.0, and CV ------------
# What are the differences between MutSig 1.0, 1.5, 2.0, and CV?
# A: MutSig relies on several sources of evidence in the data in order to estimate the amount of positive selection a gene underwent during tumorigenesis.
# 
# The three main sources are:
#   1) Abundance of mutations relative to the background mutation rate (BMR)
# 2) Clustering of mutations in hotspots within the gene
# 3) Conservation of the mutated positions (i.e. did the mutation happen at a position that is conserved across vertebrates?)
# The first line of evidence, Abundance, goes into the core significance calculation performed in all versions of MutSig. In MutSig1.0, this is simply called "p". MutSig1.0 assumes a constant BMR across all genes in the genome and all patients in the patient cohort. In MutSig1.5, this is also called "p", but MutSig1.5 uses information from synonymous mutations to roughly estimate gene-specific BMRs. Later versions of MutSig (MutSigS2N and MutSigCV) have increasingly sophisticated procedures for treating the heterogeneity in per-gene, per-patient, and per-context BMRs, but they are all answering essentially the same question about Abundance of mutations above the background level.
# The other lines of evidence, Conservation and Clustering, are examined by a separate part of MutSig (historically called "MutSig2.0") that carries out many permutations, comparing the distributions of mutations observed to the null distribution from these permutations. The output of this permutation procedure is a set of additional p-values: p_clust is the significance of the amount of clustering in hotspots within the gene. p_cons is the significance of the enrichment of mutations in evolutionarily conserved positions of the gene. Finally, p_joint is the joint significance of these two signals (Conservation and Clustering), calculated according to their joint distribution. The reason for calculating p_joint is to ensure there is no double-counting of the significance due, for example, to clustering in a conserved hotspot.
# Combining all three lines of evidence: In order to make a full accounting of the signals of positive selection in a given gene, we combine all three lines of evidence. This is done by using the Fisher method of combining p-values. The two p-values combined are the "p" (or "p_classic") from the analysis of mutation Abundance (performed by MutSig 1.0/1.5/S2N/CV), and the p_joint from the analysis of Conservation and Clustering (performed by MutSig2.0).

# summary: 
# In order to make a full accounting of the signals of positive selection in a given gene, 
# we combine all three lines of evidence. 
# 1. This is done by using the Fisher method of combining p-values. 
# 2.  The two p-values combined are the "p" (or "p_classic") from the analysis of mutation Abundance 
#   (performed by MutSig 1.0/1.5/S2N/CV), 
# 3. and the p_joint from the analysis of Conservation and Clustering (performed by MutSig2.0).

# example:
# http://gdac.broadinstitute.org/runs/analyses__2016_01_28/reports/cancer/BRCA-TP/MutSigNozzleReport2CV/nozzle.html
# http://gdac.broadinstitute.org/runs/analyses__2016_01_28/reports/cancer/BRCA-TP/MutSigNozzleReport2CV/BRCA-TP_coMut.pdf
# PIK3CA
# TP53
#....
# VS 
# http://gdac.broadinstitute.org/runs/analyses__2016_01_28/reports/cancer/BRCA-TP/MutSigNozzleReportCV/nozzle.html
# MAP3K1
# PIK3CA
# PTEN


# papers -----
# Discovery and saturation analysis of cancer genes across 21 tumour types
# https://www.nature.com/articles/nature12912
# 2014
# launch of tumorportal
# Three significance metrics were calculated for each gene, using the previously described methods MutSigCV, MutSigCL, 
# and MutSigFN. 
# MutSigCL and MutSigFN measure the significance of the positional clustering of the mutations observed, 
# as well as the significance of the tendency for mutations to occur at positions that are highly evolutionarily conserved (using conservation as a proxy for probably functional impact).
# The three MutSig tests described above (MutSigCV, MutSigCL and MutSigFN) were combined into a single final P value as 
# follows. First, a joint P value (CL + FN) for the observed clustering and conservation was calculated from the 
# joint probability distribution of the random permutations. Next, this was combined with the MutSigCV P value using 
# two methods: the Fisher method of combining P values from independent tests (http://en.wikipedia.org/wiki/Fisher's_method); 


# example:
# mcr_root=/risapps/rhel7/MATLAB/R2020a
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mcr_root/bin/glnxa64/
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mcr_root/sys/java/jre/glnxa64/jre/lib/amd64
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mcr_root/sys/java/jre/glnxa64/jre/lib/amd64/server
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mcr_root/sys/java/jre/glnxa64/jre/lib/amd64/native_threads
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mcr_root/sys/os/glnxa64
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mcr_root/bin/glnxa64
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mcr_root/runtime/glnxa64
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mcr_root/lib

# module load matlab_/v81
# ./MutSig2CV BRCA-TP.final_analysis_set.maf brca

# 1.4 was actually updated on: -rw-rw-r-- 1 sseth iacs 118K Oct 13  2016 MutSigCV.m




# functions ------
# Q: WHAT DO THE DIFFERENT FIELDS FOR SIGNIFICANTLY MUTATED GENES MEAN?
#   A: Many of these fields depend on what version of MutSig was used. The following table covers the majority of them:
#   Fields	
# Description
# MutSig_1.5	MutSig_2.0	MutSig_CV	MutSig_2CV
# gene	gene	gene	gene	HUGO Symbol
# description	description	 	longname	Full description/name of the gene
# N	N	 	 	number of sequenced bases in this gene across the individual set
# n	n	 	 	number of (nonsilent) mutations in this gene across the individual set
# nnon	nnon	number of nonsense mutations
# npat	npat	npat	npat	number of patients (individuals) with at least one nonsilent mutation
# nsite	nsite	nsite	nsite	number of unique sites having a nonsilent mutation
# nsil	nsil	nsil	nsil	number of silent mutations in this gene across the individual set
# n1	n1	 	 	number of nonsilent mutations of type "*CpG->T"
# n2	n2	 	 	number of nonsilent mutations of type "*Cp(A/C/T)->T*"
# n3	n3	 	 	number of nonsilent mutations of type "A->G"
# n4	n4	 	 	number of nonsilent mutations of type "transver"
# n5	n5	 	 	number of nonsilent mutations of type "indel+null"
# n6	n6	 	 	number of nonsilent mutations of type "double_null"
# p_ns_s	p_ns_s	 	 	p-value for the observed nonsilent/silent ratio being elevated in this gene
# p	p	p	p	p-value (overall)
# q	q	q	q	q-value, False Discovery Rate (Benjamini-Hochberg procedure)
# p_classic	 	 	p-value for the observed amount of nonsilent mutations being elevated in this gene
# p_clust	 	pCL	Clustering. Probability that recurrently mutated loci in this gene have more mutations than expected by chance. While pCV assesses the gene's overall mutation burden, pCL assesses the burden of specific sites within the gene. This allows MutSig to differentiate between genes with uniformly distributed mutations and genes with localized hotspots.
#  	p_cons	 	pFN	Conservation. Probability that mutations within this gene occur disproportionately at evolutionarily conserved sites. Sites highly conserved across vertebrates are assumed to have greater functional impact than weakly conserved sites.
#  	p_joint	 	 	p-value for joint model of clustering and conservation
#  	 	 	pCV	Abundance. Probability that the gene's overall nonsilent mutation rate exceeds its inferred background mutation rate (BMR), which is computed based on the gene's own silent mutation rate plus silent mutation rates of genes with similar covariates. BMR calculations are normalized with respect to patient-specific and sequence context-specific mutation rates.
#  	 	 	codelen	the gene's coding length
# nncd	number of noncoding mutations
# nmis	number of missense mutations
# nstp	number of readthrough mutations
# nspl	number of splice site mutations
# nind	number of indels

read_mutsigcv <- function(path = "~/projects2/ss_tnbc/analysis/art/integrative/paper_t0/ssm/mutsig/maf_t0_gb_gc"){
  patient_counts_and_rates = read_tsv(glue("{path}/patient_counts_and_rates.txt"))
  n_samples = nrow(patient_counts_and_rates)

  # https://broadinstitute.atlassian.net/wiki/spaces/GDAC/pages/844334036/FAQ
  sig_genes = read_tsv(glue("{path}/sig_genes.txt"))
  sig_genes %<>% 
    clean_names() %>% 
    mutate(frac_pat = npat/n_samples, 
           log10_q = -log10(q),
           ratio_non_sil = (nnon+1)/(nsil + 1),
           pct_pat = round(frac_pat*100, 0),
           lbl = glue("{gene} ({pct_pat}%)"))

  
  
  list(sig_genes = sig_genes, patient_counts_and_rates = patient_counts_and_rates, n_samples = n_samples)
}

run_mutsigcv <- function(x){
  
}






# END