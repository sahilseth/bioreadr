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

mutectann2maf <- function(){
  
}

# combine_mutect1_pindel <- function(mutect_ann, 
#                                    pindel_ann, 
#                                    ref_fasta = opts_flow$get("ref_fasta")
#                                    ){
#   # we would like this to create a file compatible with 
#   # maf2vcf.pl --input-maf maf.tsv --output-dir vcfs --per-tn-vcfs --ref_fasta ~/ref/human/vep/homo_sapiens/95_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
#   mutect_ann = "ssm/185_145_T0-D___185_145_GB-D_mutect.annotated.tsv"
#   mutect_maf = "ssm/185_145_T0-D___185_145_GB-D_mutect.annotated.maf"
#   mutect_vcf = "ssm/185_145_T0-D___185_145_GB-D_mutect.annotated.vcf"
#   # mutect_ann = "ssm/IPCT-S4002-MOON0051-Cap1725-4-HTID347_181220-A00422-0024-BH722GDSXX-2-CAGAAACTATTT--S4002-Cap1713-8-HTID328_190103-A00422-0027-AHGMW2DMXX-1-ATACATACCATT.011120191547202598.annotated.tsv"
#   # mutect_ann = "/rsrch3/home/iacs/sseth/projects2/ss_tnbc/data/artemis/wex/ms51_dna_b1/mutect/185_145_T0-D___185_145_GB-D_mutect.annotated.tsv"
#   
#   # mutect2maf
#   source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/mutect_ann_to_maf.R')
#   # debug(mutect_ann_to_maf)
#   maf = mutect_ann_to_maf(mutect_ann, 
#                     gene = "gene", 
#                     func = "exonicfunc.knowngene", 
#                     aa_change = "aaannotation")
#   write_tsv(maf, mutect_maf)
#   
#   # convert to VCF
#   maf2vcf_exe = "/rsrch3/home/iacs/sseth/.conda/envs/pcgr_py36/bin/maf2vcf.pl"
#   cmd_maf2vcf = glue("{maf2vcf_exe} --input-maf {mutect_maf} --output-dir ssm --per-tn-vcfs --ref_fasta {ref_fasta}")
#   run_cmd(cmd_maf2vcf, target = mutect_vcf, cmdname = "maf2vcf")
#   
#   # conv pindel2maf
#   pindel_ann = "ssm/185_145_T0-D___185_145_GB-D_pindel.annotated.tsv"
#   pindel_maf = "ssm/185_145_T0-D___185_145_GB-D_pindel.annotated.maf"
#   pindel_vcf = "ssm/185_145_T0-D___185_145_GB-D_pindel.annotated.vcf"
#   source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/mutect_ann_to_maf.R')
#   # debug(mutect_ann_to_maf)
#   debug(getTumorRef)
#   maf = mutect_ann_to_maf(pindel_ann, 
#                           gene = "gene", 
#                           mutation_score = "power",
#                           func = "exonicfunc.knowngene", 
#                           aa_change = "aaannotation")
#   write_tsv(maf, pindel_maf)
#   cmd_maf2vcf = glue("{maf2vcf_exe} --input-maf {pindel_maf} --output-dir ssm --per-tn-vcfs --ref_fasta {ref_fasta}")
#   run_cmd(cmd_maf2vcf, target = pindel_vcf, cmdname = "maf2vcf")
#   
# }




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

# ** example run --------
# samplename="
# create_phylowgs_inputs.py with the --cnvs <sample_name>=cnvs.txt

phylowgs <- function(df_trk){
  # read in titancna trk
  # ./parse_cnvs.py -f titan -c 0.81 cnv_calls_segs.txt
  # ./create_phylowgs_inputs.py --cnvs sample1=cnvs.txt --vcf-type sample1=vardict sample1=sample.vcf
  
  # run using non-matched data
  # .combine_mutect1_pindel
  
  
  
  
}
