### HLAVBSeq - a recent algorithm supposed to be good. 

# Extract a list of read name that were aligned to HLA loci (HLA-A, B, C, DM, DO, DP, DQ, DR, E, F, G, H, J, K, L, P, V, MIC, and TAP)
# samtools view NA12878.sorted.bam chr6:29907037-29915661 chr6:31319649-31326989 chr6:31234526-31241863
# chr6:32914391-32922899 chr6:32900406-32910847 chr6:32969960-32979389 chr6:32778540-32786825
# chr6:33030346-33050555 chr6:33041703-33059473 chr6:32603183-32613429 chr6:32707163-32716664
# chr6:32625241-32636466 chr6:32721875-32733330 chr6:32405619-32414826 chr6:32544547-32559613
# chr6:32518778-32554154 chr6:32483154-32559613 chr6:30455183-30463982 chr6:29689117-29699106
# chr6:29792756-29800899 chr6:29793613-29978954 chr6:29855105-29979733 chr6:29892236-29899009
# chr6:30225339-30236728 chr6:31369356-31385092 chr6:31460658-31480901 chr6:29766192-29772202
# chr6:32810986-32823755 chr6:32779544-32808599 chr6:29756731-29767588
# | awk '{print $$1}' | sort | uniq > NA12878_partial_reads.txt
# Build read name index and search read pairs and their sequences on HLA loci
# java -jar -Xmx32g -Xms32g bamNameIndex.jar index NA12878.sorted.bam --indexFile NA12878.sorted.bam.idx
# java -jar bamNameIndex.jar search NA12878.sorted.bam --name NA12878_partial_reads.txt --output NA12878_partial.sam
# java -jar SamToFastq.jar I=NA12878_partial.sam F=NA12878_partial_1.fastq F2=NA12878_partial_2.fastq
# 
# If for some reason bamNameIndex.jar doesn't work, please use bedtools to extract reads from bam files:
#     http://bedtools.readthedocs.org/en/latest/content/tools/bamtofastq.html
# 
#     Or, alternatively, below is a python script that uses pysam to extract reads by read name from a bam file:
#     http://timoast.github.io/2015/10/12/ExtractReads/
# 
#     Or, please try samgrep:
#     https://github.com/lindenb/jvarkit/wiki/SamGrep
# Extract unmapped reads
#     samtools view -bh -f 12 NA12878.sorted.bam > NA12878.sorted_unmapped.bam
#     java -jar SamToFastq.jar I=NA12878.sorted_unmapped.bam F=NA12878_unmapped_1.fastq F2=NA12878_unmapped_2.fastq
# Combine reads in FASTQ format
#     cat NA12878_partial_1.fastq NA12878_unmapped_1.fastq > NA12878_part_1.fastq
#     cat NA12878_partial_2.fastq NA12878_unmapped_2.fastq > NA12878_part_2.fastq
# For analysis with v2 HLA database:
# 
# Alignment by BWA-MEM allowing multiple alignments for each read
#     bwa index hla_all_v2.fasta
#     bwa mem -t 8 -P -L 10000 -a hla_all_v2.fasta NA12878_part_1.fastq NA12878_part_2.fastq > NA12878_part.sam
# Estimation of HLA types by HLA-VBSeq
#     For paired-end read data:
#     java -jar HLAVBSeq.jar hla_all_v2.fasta NA12878_part.sam NA12878_result.txt --alpha_zero 0.01 --is_paired
# 
#     For single-end read data:
#     java -jar HLAVBSeq.jar hla_all_v2.fasta NA12878_part.sam NA12878_result.txt --alpha_zero 0.01
# 
#     Here, alpha_zero is a hyperparameter as described in the paper and we recommend to use 0.01.


#' runHLAVBSeq
#'
#' @param bamName 
#' @param odir 
#' @param oprefix 
#' @param region 
#' @param paired 
#' 
#' @details follow this link for more details:
#' 
#' http://nagasakilab.csml.org/hla/
#' 
#' @export
#'
hlavbseq <- function(bam, 
                     samplename,
                     # region = "6:29690552-33111102", 
                     paired_end = TRUE,
                     
                     odir = "hlavbseq",
                     fqs,
                     seq_sample_id,
                     
                     bwa_exe = opts_flow$get("bwa.exe"),
                     
                     hla_all_fa = opts_flow$get("hlavbseq.hla_v2_fa"),
                     hla_allele_list = opts_flow$get("hlavbseq.hla_allele_list"),
                     hla_allele_db = opts_flow$get("hlavbseq.allele_db"),
                     
                     java_exe = opts_flow$get("java.exe"),
                     hlavbseq_jar = opts_flow$get("hlavbseq.jar"),
                     hlavbseq_call_digits_py = opts_flow$get("hlavbseq.call_digits_py"),
                     
                     hla_call_digits_opts = opts_flow$get("hlavbseq.call_digits_opts"),
                     funr_exe = opts_flow$get("R.funr_exe")
                     
){
  
  check_args()
  # get FQs
  # out_hla_fq = hla_fqs(bam = bam, samplename = samplename, region = region)
  # fqs = out_hla_fq$outfiles;
  fq1=fqs$fq1;fq2=fqs$fq2;fq3=fqs$fq3
  
  # fl nms
  # bamnm = basename(bam)
  oprefix = file.path(odir, gsub("_1.fq", "", basename(fq1)))
  samnm = glue("{oprefix}.sam")
  hlavb_out_fl = glue("{oprefix}_hlavb_out.txt")
  hlavb_d4_fl = glue("{oprefix}_hlavb_d4.txt")
  hla_pvacseq_fl = glue("{oprefix}_hla_pvac.txt")
  
  # aln
  if(paired_end){
    cmd_hla_align = glue("mkdir {odir};{bwa_exe} mem -t 8 -P -L 10000 -a {hla_all_fa} {fq1} {fq2} > {samnm}")
  }else{
    cmd_hla_align = glue("mkdir {odir};{bwa_exe} mem -t 8 -P -L 10000 -a {hla_all_fa} {fq3} > {samnm}")
  };cmd_hla_align
  
  # hlavbseq
  if(paired_end){
    # For paired-end read data:
    cmd_hlavbseq = glue("{java_exe} -jar {hlavbseq_jar} {hla_all_fa} {samnm} {hlavb_out_fl} --alpha_zero 0.01 --is_paired")
  }else{
    # single end
    cmd_hlavbseq = glue("{java_exe} -jar {hlavbseq_jar} {hla_all_fa} {samnm} {hlavb_out_fl} --alpha_zero 0.01")
  };cmd_hlavbseq
  
  # call digits
  # v xxxxx_result.txt : Need to set the output file from the HLA-VBSeq
  # a Allelelist.txt : IMGT HLA Allelelist
  # r 90 : mean single read length (mean_rlen)
  # d 4 : HLA call resolution i4 or 6 or 8
  # ispaired : if set, twice the mean rlen for depth calculation (need to specify when the sequenced data is paired-end protocol)
  if(paired_end){
    cmd_hla_call_digits = glue("{hlavbseq_call_digits_py} -v {hlavb_out_fl} -a {hla_allele_list} {hla_call_digits_opts} --ispaired > {hlavb_d4_fl}")
  }else{
    cmd_hla_call_digits = glue("{hlavbseq_call_digits_py} -v {hlavb_out_fl} -a {hla_allele_list} {hla_call_digits_opts} > {hlavb_d4_fl}")
  };cmd_hla_call_digits
  
  #to_hla_pvacseq.hlavbseq(hlavb_d4_fl, hla_pvacseq_fl)
  hlavbseq_pvac_mhc1 = glue("{odir}/{seq_sample_id}_phlat_pvac_mhc1.txt")
  hlavbseq_pvac_mhc2 = glue("{odir}/{seq_sample_id}_phlat_pvac_mhc2.txt")
  cmd_hla_pvac = glue("{funr_exe} my.ultraseq::to_pvacseq.hlavbseq x={hlavb_d4_fl} outfile1={hlavbseq_pvac_mhc1} outfile2={hlavbseq_pvac_mhc2} allele_db_fl={hla_allele_db}")
  
  cmds = list(hlavbseq = c(cmd_hla_align, cmd_hlavbseq,
                           cmd_hla_call_digits, cmd_hla_pvac))
  
  # temporarily run the last cmd only
  # cmds = list(cmd_hla = cmd_hla_pvac)
  
  
  flowmat = to_flowmat(cmds, samplename)
  
  list(flowmat = flowmat, 
       # returning these as a list if perfect for things downstream. 
       # combining them becomes easy as well
       outfiles = list(hlavb_d4_fl = hlavb_d4_fl, hlavb_out_fl = hlavb_out_fl, 
                       hla_pvacseq_fl = hla_pvacseq_fl))
}

# call digits -----
# Prediction Results:
#   Output file format
# 1st column:
#   HLA allele ID (e.g. HLA00001)
# 2nd column:
#   genomic locus length (in bp) in fasta format (e.g. 3503 bp)
# 3rd column:
#   Z, the number of reads that the algorithm assigned to the HLA allele (e.g. 300 reads)
# 4th column:
#   Normalized number of reads (Fragments Per Kilobase per Million mapped fragments)
# 5th column:
#   Relative abundance, theta, which is Z divided by the number of total mapped reads.
# How HLA typing result looks like (e.g. HLA-A)
# For v2 HLA database:
#   ./parse_result.pl Allelelist_v2.txt NA12878_result.txt | grep "^A\*" | sort -k2 -n -r > HLA_A.txt
# less HLA_A.txt
# A*01:01:01:01   17.4022266628604
# A*11:01:01      12.0376819868684
# 
# 1st column:
#   HLA allele name
# 
# 2nd column:
#   Average depth of coverage
# 
# Here, the perl script, parse_result.pl, calculates the average depth of coverage for each HLA allele.
# Please modify "parse_result.pl" according to your data. This perl script assumes 100bp x 2 data as in the paper.





#' to_hla_pvacseq.hlavbseq
#'
#' @param x d4 file
#' @param outfile1 
#'
#' @details more details here:
#' http://nagasakilab.csml.org/hla/
#' 
#' @export
to_pvacseq.hlavbseq <- function(x, 
                                outfile1, outfile2,
                                #allele_db_fl = "~/.rsrch3/home/iacs/sseth/ref/human/b37/annotations/hlavbseq/all_alleles.txt"
                                allele_db_fl = opts_flow$get("hla_allele_db")
){
  pacman::p_load(janitor, readr, dplyr)
  # x="WEX-1004-N_nochr_hla_umap_hlavb_d4.txt"
  # x = "~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/neoantigen/runs/185_003/tmp/IPCT-S2006-MOON0051-Cap2023-8-ID01_190503-A00728-0034-BHJ3W3DMXX-1-ATCACG.bwa_recalibed_hla_umap_hlavb_d4.txt"
  # hla_nethmc_db = read_tsv("nethmc_alleles.txt", col_names = "gene")$gene
  
  # pvac/netmhc recognized alleles
  hla_db = read_tsv(allele_db_fl, col_names = "gene")$gene
  
  df_hla = read_tsv(x, col_types = cols(
    Gene = col_character(),
    Allele1 = col_character(),
    Allele2 = col_character()
  )) %>% clean_names()
  head(df_hla)
  
  # "HLA-A*02:01,HLA-B*35:01",
  hla_pvac_mhc1 = df_hla %>% 
    dplyr::filter(gene %in% LETTERS) %>%
    mutate(allele1 = paste0("HLA-", allele1), 
           allele2 = paste0("HLA-", allele2)) %>% 
    dplyr::select(allele1, allele2) %>% unlist() 
  
  # DRB1*11:01 OR DRB1*11:01-DRB1*11:01
  hla_pvac_mhc2 = df_hla %>% 
    dplyr::filter(!gene %in% LETTERS) %>%
    mutate(allele = paste0(allele1, "-", allele2)) %>% 
    select(allele1, allele2, allele) %>% unlist() 
  
  # filter mhc2:
  hla_pvac_mhc2 = hla_pvac_mhc2[hla_pvac_mhc2 %in% hla_db]
  # hla_pvac = c(hla_pvac_mhc1, hla_pvac_mhc2)
  # # keep those recognized by pvacseq
  # hla_pvac = hla_pvac[hla_pvac %in% hla_db]
  # hla_pvac = hla_pvac %>% paste(collapse = ",")
  
  paste0(hla_pvac_mhc1, collapse = ",") %>% write(outfile1)
  paste0(hla_pvac_mhc2, collapse = ",") %>% write(outfile2)
  
}

to_hla_pvacseq.hlavbseq = to_pvacseq.hlavbseq




# create ref -------
# wget http://nagasakilab.csml.org/hla/hla_all_v2.fasta
# wget http://nagasakilab.csml.org/hla/Allelelist_v2.txt



# END
