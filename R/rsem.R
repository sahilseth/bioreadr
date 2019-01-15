
# Following tutorial from Bo Li
# https://github.com/bli25ucb/RSEM_tutorial


#  Build References

# gunzip ref/Mus_musculus.GRCm38.dna.toplevel.fa.gz
# gunzip ref/Mus_musculus.GRCm38.82.chr.gtf.gz
# software/RSEM-1.2.25/rsem-prepare-reference --gtf ref/Mus_musculus.GRCm38.82.chr.gtf \
#     --bowtie2 --bowtie2-path software/bowtie2-2.2.6 \
#     ref/Mus_musculus.GRCm38.dna.toplevel.fa ref/mouse_ref

# software/RSEM-1.2.25/rsem-calculate-expression -p 8 --paired-end \
#     --bowtie2 --bowtie2-path software/bowtie2-2.2.6 \
#     --estimate-rspd \
#     --append-names \
#     --output-genome-bam \
#     data/SRR937564_1.fastq data/SRR937564_2.fastq \
#     ref/mouse_ref exp/LPS_6h

# /Users/sahilseth/tech/02_bio/references/01_gencode.md

# mouse ref:
# ref="/scratch/rists/hpcapps/reference/mouse/mm10/indexes/rsem/1.2.12bowtie2/mm10.refgenes"

rsem_star_pe <- function(fq1, fq2, samplename,
                    
                    out_prefix,
                    #rsem_exp_exe = "/scratch/rists/hpcapps/x86_64/rsem/1.2.12/rsem-calculate-expression", 
                    rsem_exp_exe =  "/rsrch2/iacs/apps/conda/3.18.9/bin/rsem-calculate-expression",
                    rsem_opts = "-p 8", 
                    
                    #bowtie_exe = "/scratch/rists/hpcapps/x86_64/bowtie2/2.2.2/bowtie2", 
                    star_path = "/risapps/rhel6/STAR/2.4.2a/bin",
                    rsem_ref = "~/ref/human/indexes/rsem/1.2.28_star_gencode_v19/Genome"
                    
                    ){
  
  
  if(missing(out_prefix))
    out_prefix = paste0(samplename, "_rsem")
  
  check_args()
  
  # paired end:
  # cmd_rsem <- sprintf("%s %s --bowtie2 --bowtie2-path %s --paired-end %s %s %s %s", 
  #                     rsem_exp_exe, rsem_opts, bowtie_exe, 
  #                     fq1, fq2, rsem_ref, out_prefix)
  cmd_rsem <- sprintf("%s %s --star --star-path %s --paired-end %s %s %s %s", 
                      rsem_exp_exe, rsem_opts, star_path, 
                      fq1, fq2, rsem_ref, out_prefix)
  
  out_isoforms <- paste(out_prefix, '.isoforms.results',sep='')
  out_genes <- paste(out_prefix, '.genes.results',sep='')
  
  flowmat = list(rsem = cmd_rsem) %>% to_flowmat(samplename)
  
  return(list(flowmat = flowmat, 
              outfiles = list(isoform_counts = out_isoforms, gene_counts = out_genes)))
}

# x = "/Users/sahilseth/rsrch1_data/public_data/tnbc/art/rnaseq/rsem/IPCT-FC217-MOON0051-Cap1678-6-ID48_180910-STJ00106-0553-AHWNGNBBXX-4-TCGGCA.genes.results"
read_rsem <- function(x, reader = read_tsv, ...){
  df = reader(x, ...) %>% tbl_df() %>% clean_names()
  
  # summarize gene level counts
  # df = separate(df, col = target_id, 
  #               into = c("transcript_id", "gene_id"), sep = "\\|", remove = FALSE, extra = "drop")
  # head(df) %>% as.data.frame()
  
  # sum est_counts and transcripts
  df = group_by(df, gene_id) %>%
    mutate(gene_expected_count = sum(expected_count), 
           gene_tpm = sum(tpm), 
           gene_fpkm = sum(fpkm))
  df
}





