

# https://github.com/mskcc/lohhla
# this is part of
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5720478/pdf/main.pdf

# Running LOHHLA with no arguments prints the usage information.
# 
# USAGE: Rscript /location/of/LOHHLA/script [OPTIONS]
# 
# OPTIONS:
#   
#   -id CHARACTER, --patientId=CHARACTER
# patient ID
# 
# -o CHARACTER, --outputDir=CHARACTER
# location of output directory
# 
# -nBAM CHARACTER, --normalBAMfile=CHARACTER
# normal BAM file
# can be FALSE to run without normal sample
# 
# -BAM CHARACTER, --BAMDir=CHARACTER
# location of all BAMs to test
# 
# -hla CHARACTER, --hlaPath=CHARACTER
# location to patient HLA calls
# 
# -hlaLoc CHARACTER, --HLAfastaLoc=CHARACTER
# location of HLA FASTA [default=~/lohhla/data/hla_all.fasta]
# 
# -cn CHARACTER, --CopyNumLoc=CHARACTER
# location to patient purity and ploidy output
# can be FALSE to only estimate allelic imbalance
# 
# -ov CHARACTER, --overrideDir=CHARACTER
# location of flagstat information if already run [default= FALSE]
# 
# -mc CHARACTER, --minCoverageFilter=CHARACTER
# minimum coverage at mismatch site [default= 30]
# 
# -kmer CHARACTER, --kmerSize=CHARACTER
# size of kmers to fish with [default= 50]
# 
# -mm CHARACTER, --numMisMatch=CHARACTER
# number of mismatches allowed in read to map to HLA allele [default= 1]
# 
# -ms CHARACTER, --mappingStep=CHARACTER
# does mapping to HLA alleles need to be done [default= TRUE]
# 
# -fs CHARACTER, --fishingStep=CHARACTER
# if mapping is performed, also look for fished reads matching kmers of size kmerSize [default= TRUE]
# 
# -ps CHARACTER, --plottingStep=CHARACTER
# are plots made [default= TRUE]
# 
# -cs CHARACTER, --coverageStep=CHARACTER
# are coverage differences analyzed [default= TRUE]
# 
# -cu CHARACTER, --cleanUp=CHARACTER
# remove temporary files [default= TRUE]
# 
# -no CHARACTER, --novoDir=CHARACTER
# path to novoalign executable [default= ]
# 
# -ga CHARACTER, --gatkDir=CHARACTER
# path to GATK executable [default= ]
# 
# -ex CHARACTER, --HLAexonLoc=CHARACTER
# HLA exon boundaries for plotting [default=~/lohhla/data/hla.dat]
# 
# -w CHARACTER, --ignoreWarnings=CHARACTER
# continue running with warnings [default= TRUE]
# 
# -h, --help
# Show this help message and exit           