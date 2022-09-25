

# purecen can predict and assoc cell fraction of CN as well
# https://www.bioconductor.org/packages/release/bioc/vignettes/PureCN/inst/doc/PureCN.pdf


if(FALSE){
  
  # debug(to_purecn_vcf)
  to_purecn_vcf(tsv, gatkcnv_hets, tumor_name, normal_name, outfile)
  vcf <- PureCN:::.readAndCheckVcf(outfile)

  # , genome = genome, 
  #                         DB.info.flag = DB.info.flag, POPAF.info.field = POPAF.info.field,
  #                         min.pop.af = min.pop.af)
  
}

# this is initially geared towards sarco tsv
to_purecn_vcf <- function(tsv, gatkcnv_hets, 
                          tumor_name,
                          normal_name,
                          outfile, 
                          verbose = T){
  library(pacman)
  p_load(magrittr, tidyverse, glue, janitor)
  
  df_mut =  read_rds(tsv)$df_mut_ann_filt

  if(!any(class(tsv) %in% c("data.frame")))
    message("supply a df")
  
  df_hets = data.table::fread(data.table = F, cmd = glue("grep -v '@' {gatkcnv_hets}")) %>% as_tibble() %>% clean_names()
  df_hets_n = data.table::fread(data.table = F, cmd = glue("grep -v '@' {gatkcnv_hets}")) %>% as_tibble() %>% clean_names()
  head(df_hets)
  
  # create df VCF:
  # AD format field
  # and a DB info flag for membership in germline databases2
  df_mut_vcf = df_mut %>% 
    dplyr::rename(CHROM = chrom, 
                  POS = pos,
                  REF = ref,
                  ALT = alt) %>% 
    mutate(FILTER = "PASS", ID = ".", QUAL = ".") %>% 
    mutate(t_alt_count = round(taf*tdp, 0),
           t_ref_count = tdp - t_alt_count,
           n_alt_count = round(naf*ndp, 0),
           n_ref_count = ndp - n_alt_count) %>% 
    # no SB fields
    mutate(
      INFO = glue("TAF={taf};TDP={tdp};NAF={naf};NDP={ndp}"), 
      FORMAT = "GT:GQ:AD:AF:DP:SB",
      tumor = glue("0/1:.:{t_ref_count},{t_alt_count}:{taf}:{tdp}"),
      # should not matter as SB from normal is NEVER used....
      # it seems combine variants from GATK is skipping this value if its missing!
      normal = glue("0/0:.:{n_ref_count},{n_alt_count}:{naf}:{ndp}"))
  
  df_het_vcf = df_hets %>% 
    dplyr::rename(CHROM = contig, 
                  POS = position,
                  REF = ref_nucleotide,
                  ALT = alt_nucleotide) %>% 
    mutate(FILTER = "REJECT", ID = ".", QUAL = ".") %>% 
    mutate(t_alt_count = alt_count,
           t_ref_count = ref_count,
           tdp = t_alt_count + t_ref_count,
           taf = round(t_alt_count/tdp, 2),
           ndp = 10) %>% 
    mutate(
      INFO = glue("TAF={taf};TDP={tdp};NAF=.;NDP={ndp};DB"), 
      FORMAT = "GT:GQ:AD:AF:DP",
      tumor = glue("0/0:.:{t_ref_count},{t_alt_count}:{taf}:{tdp}"),
      # it seems combine variants from GATK is skipping this value if its missing!
      normal = glue("0/0:.:.,.:.:{ndp}"))
  head(df_het_vcf$FORMAT)
  head(df_het_vcf$tumor)
  head(df_het_vcf$normal)
  
  df_vcf = bind_rows(df_mut_vcf, df_het_vcf);dim(df_vcf)
  dplyr::count(df_vcf, is.na(tdp))
  
  # tumor is always assumed to be the last columns
  df_vcf_final = dplyr::select(df_vcf, CHROM, POS, ID, REF, ALT, QUAL, 
                               FILTER, INFO, FORMAT, 
                               normal, tumor)
  colnames(df_vcf_final) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                             "FILTER", "INFO", "FORMAT", 
                             normal_name, tumor_name)
  
  # dplyr::count(df_vcf_final, is.na(tdp))

  if(verbose)
    message("writing out a gz file: ", outfile)
  # write header
  header = .gatkcnv_header()
  # gz1 <- gzfile(outfile, "w")
  cat(header, file = outfile, sep = "\n")
  write_tsv(df_vcf_final, outfile, append = T, col_names = T)
  # close(gz1)
  # check it
  
  
  
}


.gatkcnv_header <- function(){
  c(
    '##fileformat=VCFv4.2',
    '##fileDate=20200601',
    '##source=tsv+gatkcnv hets??',
    '##reference=file:///seq/references/Homo_Sapiens.fasta',
    '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">',
    '##FILTER=<ID=REJECT,Description="Rejected as a confident somatic mutation">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    '##INFO=<ID=TAF,Number=1,Type=Float,Description="VAF in tumor">',
    '##INFO=<ID=NAF,Number=1,Type=Float,Description="VAF in normal">',
    '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Depth in tumor">',
    '##INFO=<ID=NDP,Number=1,Type=Integer,Description="Depth in normal">',
    '##INFO=<ID=DB,Number=0,Type=Flag,Description="Presence in GATK-hets file">',
    
    # picking up from mutect2 vcf
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fishers Exact Test to detect strand bias.">')
}  
  
  
  
  
