# install:
# conda create -n theta_py27 python=2.7 theta2
# https://bioconda.github.io/recipes/theta2/README.html

# get the repo as well:
# git clone https://github.com/raphael-group/THetA.git

# positional arguments:
#   QUERY_FILE            Interval file

# optional arguments:
#   -h, help            show this help message and exit
#   TUMOR_FILE TUMOR_FILE
#                         File containing allelic counts for tumor sample SNPs.
#   NORMAL_FILE NORMAL_FILE
#                         File containing allelic counts for normal samlpe SNPs.
#   -n N, N N           Number of subpopulations
#   -k K, MAX_K K       The maximum value allowed for entries in C
#   -t TAU, TAU TAU     Expected number of copies in normal genome
#   -d DIR, DIR DIR     Directory where result file is written to
#   -p PRE, OUTPUT_PREFIX PRE
#                         Prefix for output files created. By default, it will
#                         be the beginning of the input file name (i.e.if input
#                         filename were example.input, output filed would be
#                         example.output and example.withbounds)
#   -m MAX_NORMAL, MAX_NORMAL MAX_NORMAL
#                         The maximum fraction to consider for normal. Only
#                         enforced for n=2
#   NUM_PROCESSES NUM_PROCESSES
#                         The number of processes to be used
#   NUM_INTERVALS NUM_INTERVALS
#                         The maximum number of intervals used by automatic
#                         interval selection.
#   BOUND_HEURISTIC BH
#   NORMAL_BOUND_HEURISTIC NBH
#   HEURISTIC_LB LB
#   HEURISTIC_UB UB
#   BOUNDS_ONLY
#   NO_MULTI_EVENT
#   RESULTS filename
#   FORCE
#   GET_VALUES
#   NO_INTERVAL_SELECTION
#   READ_DEPTH_FILE FILENAME
#   GRAPH_FORMAT GRAPH_FORMAT
#                         Options are .pdf, .jpg, .png, .eps
#   BAF                 Option to run the BAF model.
#   RATIO_DEV RATIO_DEV
#                         The deviation away from 1.0 that a ratio must be to be
#                         considered a potential copy number aberration.
#   MIN_FRAC MIN_FRAC   The minimum fraction of the genome that must have a
#                         potential copy number aberration to be a valid sample
#                         for THetA analysis.
#   NO_CLUSTERING       Option to run THetA without clustering.

theta <- function(trk, samplename,
                  funr_exe = opts_flow$get("R_funr_exe"),
                  thetadir = pipe_str$theta$dir,
                  theta_opts = opts_flow$envir$theta_opts


){

    # print function name
    # message(match.call()[1], " ", appendLF = F)
    check_args()

    flog.debug(paste0("Generating a titancna_somatic_matched flowmat for sample: ", samplename))


    # check generic columns
    trk = metadata_for_dnaseq_tools(trk)
    # we would use
    # WEX-sarco14-T.denoisedCR.tsv instead of WEX-sarco14-T___matched.cr.seg
    # coz titan segments it
    expect_columns(trk, c(
        "outprefix", "outprefix_paired",
        "cnr_file", "hets_file",
        "gender"
    ))

    # scenario 1:
    #     gatk is already complete, we just need to run this
    # scenario 2:
    #   cnr file is coming from some other tool

    # ceate a new trk with ALL reqd files
    trk %<>%
        mutate(
            # skip taking basename of bam!
            # bam = basename(bam),
            # oprefix = glue("{outpath}{name}"),
            theta_oprefix = glue("{thetadir}/{outprefix}"),
            theta_cnr_file = glue("{thetadir}/{outprefix}.denoisedCR_theta.tsv"),
            theta_hets_file = glue("{thetadir}/{outprefix_paired}.hets_theta.tsv"),
            theta_hets_n_file = glue("{thetadir}/{outprefix_paired}.hets.normal_theta.tsv")
        )
    trk

    trk_tum = filter(trk, normal == "NO")
    trk_norm = filter(trk, normal == "YES")

    # ** prep inputs -
    # opts_flow$load_toml("../flowr_conf.toml")
    cnr_file = trk_tum$gatkcnv_dn_cr
    gatkcnv_cnts = trk_tum$gatkcnv_cnts
    theta_cnr_file = trk_tum$theta_cnr_file
    cmd_prep_cnr = glue("{funr_exe} my.ultraseq::to_theta_cn_file.gatkcnv cnr_file={cnr_files} theta_cnr_file={theta_cnr_file} gatkcnv_cnts={gatkcnv_cnts}")
    
    hets_files = c(trk_tum$gatkcnv_hets, trk_tum$gatkcnv_hets_n)
    theta_hets_files = c(trk_tum$theta_hets_file, trk_tum$theta_hets_n_file)
    cmd_prep_hets = glue("{funr_exe} my.ultraseq::to_theta_het_file.gatkcnv hets_file={hets_files} theta_hets_file={theta_hets_files}")


    # install.my.ultra
    cmds_prep = c(cmd_prep_cnr, cmd_prep_hets)
    # lapply(cmds_prep, system)

    # RunTHetA.py
    theta_exe = "conda activate theta_py27;RunTHetA.py"
    theta_opts = "--NUM_PROCESSES 5"
    oprefix = trk_tum$theta_oprefix
    cmd_theta = glue("{theta_exe} {theta_cnr_file} --TUMOR_FILE {theta_hets_files[1]} --NORMAL_FILE {theta_hets_files[2]} -p {oprefix} {theta_opts}")
    cmd_theta

    if(FALSE){
        to_theta_cn_file.gatkcnv(cnr_files, gatkcnv_cnts, theta_cnr_file)
        to_theta_het_file.gatkcnv(hets_files[1], theta_hets_files[1])
        to_theta_het_file.gatkcnv(hets_files[2], theta_hets_files[2])
        # RunTHetA.py theta/IPCT-S4013-MOON0051-Cap2531-4-HTID254.denoisedCR_theta.tsv --TUMOR_FILE theta/IPCT-S4013-MOON0051-Cap2531-4-HTID254___matched.hets_theta.tsv --NORMAL_FILE theta/IPCT-S4013-MOON0051-Cap2531-4-HTID254___matched.hets.normal_theta.tsv -p theta/IPCT-S4013-MOON0051-Cap2531-4-HTID254 --NUM_PROCESSES 5

    }

}

# 1. Create an interval_count_file:
# The interval_count_file must contain one line per genomic interval.  Each line
# is tab delimited containing the following fields:
# 	1. interval ID (String) - an identifier for the interval
# 	2. chromosome number (String) - chromosome on which the interval
# 	   occurs. Valid formats include (1, 2, ..., X, Y), (chr1, chr2, ..., X, Y), 
#        or (chrm1, chrm2, ..., X, Y)
# 	3. interval start coordinate (Integer) - starting position for the
# 	   interval
# 	4. interval end coordinate (Integer) - ending position for the interval
# 	5. tumorCount (Integer) - number of reads contained within interval for
# 	   tumor BAM
# 	6. normalCount (Integer) - number of reads contained within interval for
# 	   normal BAM
# 	7. * upperBound (Integer) [Optional] - the maximum copy number to consider for the
# 		interval.
# 	8. * lowerBound (Integer) [Optional]- the minimum copy number to consider for the
# 		interval.
# 	The interval_count_file may contain an optional header line, starting with
# 	a "#" symbol.
# 	ADDITIONAL NOTES for more than one tumor population:
# 	After running for a single tumor population, THetA now produces a .bash
# 	script [prefix].RunN3.bash that can be run to consider more than one tumor
# 	population. This is the recommended way to run THetA for multiple tumor
# 	populations.
# 	NOTE: See the section on WHOLE-GENOME PREPROCESSING SCRIPTS to see how to convert
# 		  BIC-Seq output into THetA input.
# 	NOTE: See the section on WHOLE-EXOME PREPROCESSING SCRIPTS to see how to convert
# 		 Whole-exome data (with specific support for ExomeCNV and EXCAVATOR)
# 		 into THetA input.

#' to_theta_cn_file.gatkcnv
#'
#' @param hets_file
#' @param titan_hets_file
#'
#' @export
to_theta_cn_file.gatkcnv <- function(cnr_file,
                                     gatkcnv_cnts,
                                     theta_cnr_file) {
    
    # these are centered around 0
    # summary(df_cnr$log2_copy_ratio)
    #  Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    # -30.00634  -0.21394   0.00006   0.02714   0.22351   4.79311 
    df_cnr = read_counts(cnr_file)
    df_cnts = read_counts(gatkcnv_cnts)
    df_cnr = left_join(df_cnr, df_cnts, by = c("contig", "start", "end"))
    df_cnr %<>% 
        mutate(count, normal_count = count*(1/(2^log2_copy_ratio)),
            # conv DP to actual count

               normal_count = as.integer(normal_count),
               # interval_id = glue("{contig}:{start}-{end}"),
               id = 1:n(), cr_calc = log2(count/normal_count)) 
    head(df_cnr, 100)
    # cor.test(df_cnr$log2_copy_ratio, df_cnr$cr_calc)
    df_cnr %<>% select(id, contig, start, end, tumor_count = count, normal_count)
    colnames(df_cnr)[1] = "#id"

    write_tsv(df_cnr, theta_cnr_file)
}

# Format 1: A tab (or comma) delimited file with the following columns
# (1) chromosome (Integer) - Chromosome on which a germline SNP may occur.
# (2) position (Integer) - Genomic coordinate of potential SNP position.
# (3) ref_allele - The number of reads aligning to the position with the reference allele.
# (4) mut_allele - The number of reads aligning to the position with the mutant allele.

#' to_theta_het_file.gatkcnv
#'
#' @param hets_file
#' @param titan_hets_file
#'
#' @export
to_theta_het_file.gatkcnv <- function(hets_file, theta_hets_file) {

    # cat SAMPLE.hets.tsv | awk '/^[^@]/ { print $1,$2,$5,$3,$6,$4 }' | tr ' ' '\t' > tmp.hets.tsv

    # POSITON BASED
    df_hets = read_counts(hets_file)
    df_hets %<>% select(
        chromosome = contig, 
        pos = position,
        ref_allele = ref_count,
        mut_allele = alt_count)
    colnames(df_hets)[1] = "#chromosome"

    write_tsv(df_hets, theta_hets_file)
}

#  * bash bin/CreateExomeInput -s segment_file -t tumor_bam -n normal_bam --FA fasta_file --EXON_FILE exon_bed_file [Options}
#         INPUT:
#         * segment_file: Tab delimited file with a set of intervals of the reference
#                         genome returned by running an exome segmenter/copy number program.
#                         This file has the the following columns:
#                         (1) chrm
#                         (2) start
#                         (3) end
#                 Note: If using ExomeCNV this file is: FILE_NAME.segment.copynumber.txt
#                       If using EXCAVATOR this file is: FastCallResults_FILE_NAME.txt
#         * tumor_bam: the tumor bam file (should be indexed and sorted)
#         * normal_bam: the matched normal bam file (should be indexed and sorted)
#         * fasta_file: reference genome file used in mapping the raw reads
#         * exon_bed_file: bed file of all exons (data/ contains bed files for hg18 and hg19)
#                 data/hg18.exons.bed
#                 data/hg19.exons.bed
#         * [Options] - see later section on all optional parameters.
#         OUTPUT:
#         * theta_input: a file that is properly formatted as input to the THetA algorithm
# this operates on already segmented input from GATKCNV, which is perfect
theta_create_exome_input <- function(){
    #  -s gatkcnv/IPCT-S4013-MOON0051-Cap2531-4-HTID254___matched.called.igv.seg -t IPCT-S4013-MOON0051-Cap2531-4-HTID254_ngs-pipe-2-GCTACTGA.bwa_recalibed.bam -n IPCT-S4012-MOON0051-Cap2509-8-HTID373_ngs-pipe-2-GGCAGTCTCTGG.bwa_recalibed.bam --FA /rsrch3/home/iacs/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta --EXON_FILE /rsrch3/home/iacs/sseth/ref/az_ref_beds/ss_downloads/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets_pad250.bed 
    theta_exome_input_exe = "~/apps/theta/THetA/bin/CreateExomeInput"
    segfile = 
    glue("{theta_exome_input_exe}")
}