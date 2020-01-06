# URL: https://software.broadinstitute.org/cancer/cga/mutect_run 
# 
# ---------------------- example run -----------------------------------
# Require Java 6 Runtime
# 
# java -Xmx2g -jar muTect-XXXX-XX-XX.jar
# --analysis_type MuTect
# --reference_sequence <reference>
#   --cosmic <cosmic.vcf>
#   --dbsnp <dbsnp.vcf>
#   --intervals <intervals_to_process>
#   --input_file:normal <normal.bam>
#   --input_file:tumor <tumor.bam>
#   --out <call_stats.out>
#   --coverage_file <coverage.wig.txt> 
#   These parameters are based upon your genome build for your alignments:
#   
#   For HG18
# 
# <reference> - Homo_sapiens_assembly18.fasta
# <dbsnp.vcf> - dbsnp_132.hg18.vcf
# <cosmic.vcf> - hg18_cosmic_v54_120711.vcf
# For HG19/GRC37
# 
# <reference> - Homo_sapiens_assembly19.fasta
# <dbsnp.vcf> - dbsnp_132_b37.leftAligned.vcf
# <cosmic.vcf> - hg19_cosmic_v54_120711.vcf
# For Mouse MM9
# 
# <reference> - Mus_musculus_assembly9.fasta
# <dbsnp.vcf> - dbsnp_128_mm9.vcf
# <cosmic.vcf> - there is no cosmic VCF available for mouse, this entire parameter can be eliminated
# Whereas these parameters are related to the sample/BAM:
#   
#   <intervals_to_process> - either a literal list of "chrom:start-end" separated by semicolons  (e.g. chr1:1500-2500; chr2:2500-3500) or a file of such entries with one entry per line
# <normal.bam> - BAM file for the Normal (positional, this must be before the tumor BAM file)
# <tumor.bam> - BAM file for the Tumor
# <call_stats.out> - filename to write detailed caller output
# <coverage.wig.txt> - filename for coverage output
# 
# ---------------------- output params -----------------------------------
# contig - the contig location of this candidate
# position - the 1-based position of this candidate on the given contig
# ref_allele - the reference allele for this candidate
# alt_allele - the mutant (alternate) allele for this candidate
# tumor_name - name of the tumor as given on the command line, or extracted from the BAM
# normal_name - name of the normal as given on the command line, or extracted from the BAM
# score - for future development
# dbsnp_site - is this a dbsnp site as defined by the dbsnp bitmask supplied to the caller
# covered - was the site powered to detect a mutation (80% power for a 0.3 allelic fraction mutation)
# power - tumor_power * normal_power
# tumor_power - given the tumor sequencing depth, what is the power to detect a mutation at 0.3 allelic fraction
# normal_power - given the normal sequencing depth, what power did we have to detect (and reject) this as a germline variant
# total_pairs - total tumor and normal read depth which come from paired reads
# improper_pairs - number of reads which have abnormal pairing (orientation and distance)
# map_Q0_reads - total number of mapping quality zero reads in the tumor and normal at this locus
# init_t_lod - deprecated
# t_lod_fstar - CORE STATISTIC: Log of (likelihood tumor event is real / likelihood event is sequencing error )
# tumor_f - allelic fraction of this candidated based on read counts
# contaminant_fraction - estimate of contamination fraction used (supplied or defaulted)
# contaminant_lod - log likelihood of ( event is contamination / event is sequencing error )
# t_ref_count - count of reference alleles in tumor
# t_alt_count - count of alternate alleles in tumor
# t_ref_sum - sum of quality scores of reference alleles in tumor
# t_alt_sum - sum of quality scores of alternate alleles in tumor
# t_ins_count - count of insertion events at this locus in tumor
# t_del_count - count of deletion events at this locus in tumor
# normal_best_gt - most likely genotype in the normal
# init_n_lod - log likelihood of ( normal being reference / normal being altered )
# n_ref_count - count of reference alleles in normal
# n_alt_count - count of alternate alleles in normal
# n_ref_sum - sum of quality scores of reference alleles in normal
# n_alt_sum - sum of quality scores of alternate alleles in normal
# judgement - final judgement of site KEEP or REJECT (not enough evidence or artifact)

# to_vcf --------

if(FALSE){
  # mutect2vcf
  library(pacman)
  p_load(dplyr, readr, janitor, magrittr)
  x = "/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/mutect/WEX-334187-T__WEX-334187-N_merged.muTect_call_stats.txt"
  outfile = "/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/mutect/WEX-334187-T__WEX-334187-N_merged.muTect_call_stats2.vcf"
  
  chrom = "contig"
  start_pos = "position"
  end_pos = "position"
  
  ref_allele = "ref_allele"
  alt_allele = "alt_allele"
  
  filter = "judgement"
  
  # AF
  t_alt_count = "t_alt_count"
  t_ref_count = "t_ref_count"
  
  sample_name = "tumor_name"
  ref_name = "sample_name"
  
  # mut qual
  mutation_score = "t_lod_fstar"
  
}




#' to_vcf.mutect_call_stats
#'
#' Assumes fixed coloumns of VCF
#' 
#' @param x 
#' @param outfile ouput VCF file
#' 
#' @import pacman
#' 
#' @export
to_vcf.mutect_call_stats <- function(x, outfile){
  
  pacman::p_load(rlogging)
  rlogging::message("loading file")
  if(is.data.frame(x))
    df = as.data.frame(x, stringsAsFactors = FALSE)
  else
    df = read_tsv(x, col_types = cols(.default = col_character()))
  
  # MAIN VCF cols:
  # "\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  
  # need to encode NA, as .
  rlogging::message("parsing info/format fields")
  df_vcf = df %>% 
    select(CHROM = contig, 
           POS = position,
           REF = ref_allele, 
           ALT = alt_allele, 
           FILTER = judgement, 
           everything()) %>% 
    mutate(ID = ".", 
           QUAL = ".",
           FILTER = ifelse(FILTER == "KEEP", "PASS", FILTER))
  # df_vcf
  # FORMAT
  # INFO
  df_vcf %<>% mutate(
    t_alt_max_mapq,
    t_ref_count = as.integer(t_ref_count),
    t_alt_count = as.integer(t_alt_count),
    n_ref_count = as.integer(n_ref_count),
    n_alt_count = as.integer(n_alt_count),
    
    t_dp = t_ref_count + t_alt_count,
    n_dp = n_ref_count + n_alt_count,
    t_lod_fstar = round(as.integer(t_lod_fstar), 2),
    init_n_lod = round(as.integer(init_n_lod), 2),
    normal_f = n_alt_count/n_dp,
    
    INFO = glue("TLOD={t_lod_fstar}:NLOD={init_n_lod}:TAF={tumor_f}:TDP={t_dp}:NAF={normal_f};NDP={n_dp}"), 
    FORMAT = "GT:GQ:AD:AF:DP",
    tumor = glue("0/1:.:{t_ref_count},{t_alt_count}:{tumor_f}:{t_dp}"),
    normal = glue("0/0:.:{n_ref_count},{n_alt_count}:{normal_f}:{n_dp}"))

  # extract tumor normal name
  tumor_name = df_vcf$tumor_name[1]
  normal_name = df_vcf$normal_name[1] 
  df_vcf_final = dplyr::select(df_vcf, 
                               CHROM, POS, ID, REF, ALT, QUAL, 
                               FILTER, INFO, FORMAT, 
                               tumor, normal)
  colnames(df_vcf_final) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                             "FILTER", "INFO", "FORMAT", 
                             tumor_name, normal_name)
  
  rlogging::message("writing out")
  # write header
  header = .mutect1_header()
  cat(header, file = outfile, sep = "\n")
  write_tsv(df_vcf_final, outfile, append = T, col_names = T)
  
}


#' to_vcf.mutect_call_stats
#'
#' Assumes fixed coloumns of VCF
#' 
#' @param x 
#' @param outfile ouput VCF file
#' 
#' @import pacman
#' 
#' @export
to_vcf.mutect1.7_extended_tsv <- function(x, outfile){
  
  ext = tools::file_ext(outfile)
  if(ext == "gz")
    stop("writing compressed files is not supported")
  
  # x = "/rsrch3/home/iacs/sseth/flows/SS/sarco/jco/wex/mutect1/WEX-sarco10-T___matched_mutect.merged.tsv.gz"
  pacman::p_load(rlogging)
  rlogging::message("loading file")
  if(is.data.frame(x))
    df = as.data.frame(x, stringsAsFactors = FALSE)
  else
    df = read_tsv(x, col_types = cols(.default = col_character()))
  
  # MAIN VCF cols:
  # "\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  
  # need to encode NA, as .
  rlogging::message("parsing info/format fields")
  # table(df$failure_reasons, df$judgement, useNA = "always")
  # NA: PASS, otherwise we have a reason
  tmp = df$failure_reasons %>% strsplit(split = ",") %>% unlist()
  table(tmp) %>% names()
  df_vcf = df %>% 
    dplyr::select(CHROM = contig, 
           POS = position,
           REF = ref_allele, 
           ALT = alt_allele, 
           FILTER = failure_reasons, 
           everything()) %>% 
    dplyr::mutate(ID = ".", 
           QUAL = ".",
           FILTER = ifelse(is.na(FILTER), "PASS", FILTER))
  # df_vcf
  # FORMAT
  # INFO
  df_vcf %<>% mutate(
    t_alt_max_mapq,
    t_ref_count = as.integer(t_ref_count),
    t_alt_count = as.integer(t_alt_count),
    n_ref_count = as.integer(n_ref_count),
    n_alt_count = as.integer(n_alt_count),
    
    t_dp = t_ref_count + t_alt_count,
    n_dp = n_ref_count + n_alt_count,
    t_lod_fstar = round(as.integer(t_lod_fstar), 2),
    init_n_lod = round(as.integer(init_n_lod), 2),
    normal_f = n_alt_count/n_dp,
    sb = gsub("\\(", "", strand_bias_counts),
    sb = gsub("\\)", "", sb),
    
    INFO = glue("TLOD={t_lod_fstar}:NLOD={init_n_lod}:TAF={tumor_f}:TDP={t_dp}:NAF={normal_f};NDP={n_dp}"), 
    FORMAT = "GT:GQ:AD:AF:DP:SB",
    tumor = glue("0/1:.:{t_ref_count},{t_alt_count}:{tumor_f}:{t_dp}:{sb}"),
    # should not matter as SB from normal is NEVER used....
    # it seems combine variants from GATK is skipping this value if its missing!
    normal = glue("0/0:.:{n_ref_count},{n_alt_count}:{normal_f}:{n_dp}:-99,-99,-99,-99"))
  
  # extract tumor normal name
  tumor_name = df_vcf$tumor_name[1]
  normal_name = df_vcf$normal_name[1] 
  df_vcf_final = dplyr::select(df_vcf, 
                               CHROM, POS, ID, REF, ALT, QUAL, 
                               FILTER, INFO, FORMAT, 
                               tumor, normal)
  colnames(df_vcf_final) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                             "FILTER", "INFO", "FORMAT", 
                             tumor_name, normal_name)
  
  rlogging::message("writing out a gz file: ", outfile)
  # write header
  header = .mutect1_header()
  # gz1 <- gzfile(outfile, "w")
  cat(header, file = outfile, sep = "\n")
  write_tsv(df_vcf_final, outfile, append = T, col_names = T)
  # close(gz1)
}

.mutect1_header <- function(){
  c(
    '##fileformat=VCFv4.2',
    '##fileDate=20090805',
    '##source=mutect v1.7??',
    '##reference=file:///seq/references/Homo_Sapiens.fasta',
    #'##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>',
    #'##phasing=partial',
    '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">',
    '##FILTER=<ID=REJECT,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=alt_allele_in_normal,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=clustered_read_position,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=fstar_tumor_lod,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=germline_risk,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=nearby_gap_events,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=normal_lod,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=poor_mapping_region_alternate_allele_mapq,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=poor_mapping_region_mapq0,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=possible_contamination,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=strand_artifact,Description="Rejected as a confident somatic mutation">',
    '##FILTER=<ID=triallelic_site,Description="Rejected as a confident somatic mutation">',
    
    #'##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">',
    # '##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">',
    # '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">',
    # '##INFO=<ID=VT,Number=1,Type=String,Description="Variant type, can be SNP, INS or DEL">',
    #'##INFO=<ID=TLOD,Number=1,Type=String,Description="CORE STATISTIC: Log of (likelihood tumor event is real / likelihood event is sequencing error )">',
    # copy over values from mutect2 (gatk)
    '##INFO=<ID=TLOD,Number=A,Type=Float,Description="Log odds ratio score for variant">',
    '##INFO=<ID=NLOD,Number=A,Type=Float,Description="Normal LOD score">',
    '##INFO=<ID=TAF,Number=1,Type=Float,Description="VAF in tumor">',
    '##INFO=<ID=NAF,Number=1,Type=Float,Description="VAF in normal">',
    '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Depth in tumor">',
    '##INFO=<ID=NDP,Number=1,Type=Integer,Description="Depth in normal">',
    
    # picking up from mutect2 vcf
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fishers Exact Test to detect strand bias.">'
      
    
    # '##FORMAT=<ID=BQ,Number=A,Type=Float,Description="Average base quality for reads supporting alleles">',
    # '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
    # '##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">',
    # '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">',
    # '##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">'
    
  )
}

# https://raw.githubusercontent.com/ruping/seqare/master/mutect2vcf.pl
# perl $mutect2vcfBin $inMutect >$outVCF

# use strict;
# 
# my $mutect = shift;
# my $realSamples = shift;
# 
# my $vcfheader .= '##fileformat=VCFv4.1'."\n";
# $vcfheader .= '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">'."\n";
# $vcfheader .= '##FILTER=<ID=REJECT,Description="Rejected as a confident somatic mutation">'."\n";
# $vcfheader .= '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'."\n";
# $vcfheader .= '##FORMAT=<ID=BQ,Number=A,Type=Float,Description="Average base quality for reads supporting alleles">'."\n";
# $vcfheader .= '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'."\n";
# $vcfheader .= '##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">'."\n";
# $vcfheader .= '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'."\n";
# $vcfheader .= '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'."\n";
# $vcfheader .= '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">'."\n";
# $vcfheader .= '##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">'."\n";
# $vcfheader .= '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">'."\n";
# $vcfheader .= '##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">'."\n";
# $vcfheader .= '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">'."\n";
# $vcfheader .= '##INFO=<ID=VT,Number=1,Type=String,Description="Variant type, can be SNP, INS or DEL">'."\n";
# 
# my $samples = 'SRP';
# my %colindex;
# if ($mutect =~ /\.gz$/) {
#   open MU, "gzip -dc $mutect |";
# } else {
#   open MU, "$mutect";
# }
# while ( <MU> ) {
#   chomp;
#   next if /^#/;   #skip comments
#     my @cols = split /\t/;
#     if ( /^contig\t/ ){
#       for (my $i = 0; $i <= $#cols; $i++) {
#            $colindex{$cols[$i]} = $i;
#     }
#     next;
# }
# #my ($contig, $position, $context, $ref_allele, $alt_allele, $tumor_name, $normal_name, $score, $dbsnp_site, $covered, $power, $tumor_power, $normal_power, $normal_power_nsp, $normal_power_wsp, $total_pairs, $improper_pairs, $map_Q0_reads, $init_t_lod, $t_lod_fstar, $t_lod_lqs, $t_lod_fstar_forward, $t_lod_fstar_reverse, $tumor_f, $tumor_f_lb, $contaminant_fraction, $contaminant_lod, $t_q20_count, $t_ref_count, $t_alt_count, $t_ref_sum, $t_alt_sum, $t_ref_max_mapq, $t_alt_max_mapq, $t_ins_count, $t_del_count, $normal_best_gt, $init_n_lod, $n_lod_fstar, $normal_f, $normal_f_quals, $normal_artifact_lod_tf, $normal_artifact_lod_low_tf, $normal_artifact_lod_nf, $normal_artifact_lod_nfq, $n_q20_count, $n_ref_count, $n_alt_count, $n_ref_sum, $n_alt_sum, $power_to_detect_positive_strand_artifact, $power_to_detect_negative_strand_artifact, $strand_bias_counts, $tumor_alt_fpir_median, $tumor_alt_fpir_mad, $tumor_alt_rpir_median, $tumor_alt_rpir_mad, $alt_fpir, $alt_rpir, $powered_filters, $normal_artifact_power_tf, $normal_artifact_power_low_tf, $normal_artifact_power_nf, $normal_global_qll, $normal_local_qll, $normal_qmodel_lod, $observed_in_normals_count, $failure_reasons, $judgement) = split /\t/;
# 
# if ($samples eq 'SRP'){
#   $vcfheader .= "\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
#   my $tumor_name = $cols[$colindex{'tumor_name'}];
#   my $normal_name = $cols[$colindex{'normal_name'}];
#   if ($realSamples ne '') {
#     my @realSamples = split (/\,/, $realSamples);
#     $samples = join("\t", @realSamples);
#     $vcfheader .= "\t$samples\n";
#   } elsif ($tumor_name ne '' and $normal_name ne ''){
#     $vcfheader .= "\t$tumor_name\t$normal_name\n";
#     $samples = "$tumor_name\t$normal_name";
#   }
#   print "$vcfheader";
# }
# 
# my $chr = $cols[$colindex{'contig'}];
# my $pos = $cols[$colindex{'position'}];
# my $id = '.';
# my $ref = $cols[$colindex{'ref_allele'}];
# my $alt = $cols[$colindex{'alt_allele'}];
# my $qual = '.';
# my $judgement = $cols[$colindex{'judgement'}];
# my $filter = ($judgement eq 'KEEP')? 'PASS':$judgement;
# my $t_lod_fstar = $cols[$colindex{'t_lod_fstar'}];
# my $init_n_lod = $cols[$colindex{'init_n_lod'}];
# my $info .= 't_lod_fstar='.$t_lod_fstar.';';
# $info .= 'init_n_lod='.$init_n_lod.';';
# my $format = 'GT:AD:BQ:DP:FA';
# 
# my $t_ref_count = $cols[$colindex{'t_ref_count'}];
# my $t_alt_count = $cols[$colindex{'t_alt_count'}];
# my $n_ref_count = $cols[$colindex{'n_ref_count'}];
# my $n_alt_count = $cols[$colindex{'n_alt_count'}];
# my $t_alt_max_mapq = $cols[$colindex{'t_alt_max_mapq'}];
# my $t_dp = $t_ref_count + $t_alt_count;
# my $n_dp = $n_ref_count + $n_alt_count;
# my $tumor_f = $cols[$colindex{'tumor_f'}];
# my $normal_f = (exists($colindex{'normal_f'}))? $cols[$colindex{'normal_f'}]: (($n_dp > 0)? sprintf("%.5f", $n_alt_count/$n_dp):0);
# my $tumor = '0/1:'."$t_ref_count,$t_alt_count\:"."$t_alt_max_mapq\:"."$t_dp\:"."$tumor_f";
# my $normal = '0:'."$n_ref_count,$n_alt_count\:".'.:'."$n_dp\:"."$normal_f";
# printf("%s\n", join("\t", $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $tumor, $normal));
# 
# }
# close MU;
# 
# exit 0;






# END
