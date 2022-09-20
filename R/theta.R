# install:
# conda create -n theta_py27 python=2.7 theta2
# https://bioconda.github.io/recipes/theta2/README.html

# get the repo as well:
# git clone https://github.com/raphael-group/THetA.git

# example ---------
# --NO_MULTI_EVENT: should be false, allow multi events at same intervals
# RunTHetA.py example/Example.intervals --TUMOR_FILE example/TUMOR_SNP.formatted.txt --NORMAL_FILE example/NORMAL_SNP.formatted.txt --BAF -p example2 -n 2 --NUM_PROCESSES 8
# bash example2.RunN3.bash
# python ~/apps/theta/THetA/python/RunBAFModel.py example/TUMOR_SNP.formatted.txt example/NORMAL_SNP.formatted.txt example/Example.intervals example2.BEST.results -O example/

theta_prep_filenms <- function(trk, thetadir){
    trk %<>%
        mutate(
            # skip taking basename of bam!
            # bam = basename(bam),
            # oprefix = glue("{outpath}{name}"),
            theta_oprefix = glue("{thetadir}/{outprefix_paired}"),
            theta_cnr_file = glue("{theta_oprefix}.input"),
            theta_hets_file = glue("{theta_oprefix}.hets_theta.tsv"),
            theta_hets_n_file = glue("{theta_oprefix}.hets.normal_theta.tsv"),
            theta_results2 = glue("{theta_oprefix}.n2.igv.seg"),
            theta_results3 = glue("{theta_oprefix}.n3.igv.seg"),
            theta_fitfl = glue("{theta_oprefix}.df_fit_summary.tsv"),
            theta_results = glue("{theta_oprefix}.BEST.results")
        )
    trk
}

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

theta_matched <- function(trk, samplename,
                  funr_exe = opts_flow$get("R_funr_exe"),
                  theta_exe = opts_flow$envir$theta_exe,
                  thetadir = pipe_str$theta$dir,
                  theta_opts = opts_flow$envir$theta_opts
){

    # print function name
    # message(match.call()[1], " ", appendLF = F)
    check_args()

    flog.debug(paste0("Generating a theta_matched flowmat for sample: ", samplename))
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
    flog.debug(paste0("Generating a filenames flowmat for sample"))
    trk = theta_prep_filenms(trk, thetadir)

    trk_tum = filter(trk, normal == "NO")
    trk_norm = filter(trk, normal == "YES")

    # ** pileup cmds -------
    cmd_pileup = theta_create_mpileup(trk$bam)
    cmd_rm_pileup = theta_rm_pileup(trk$bam)

    # ** prep inputs ----------
    i <- 1
    tmp = lapply(1:nrow(trk_tum), function(i){
        # opts_flow$load_toml("../flowr_conf.toml")
        theta_oprefix = trk_tum$theta_oprefix[i]
        cnr_file = trk_tum$gatkcnv_dn_cr[i]
        gatkcnv_cnts = trk_tum$gatkcnv_cnts[i]
        theta_cnr_file = trk_tum$theta_cnr_file[i]
        segfile = trk_tum$gatkcnv_cr_seg[i]
        # cmd_prep_cnr = glue("{funr_exe} my.ultraseq::to_theta_cn_file.gatkcnv cnr_file={cnr_files} theta_cnr_file={theta_cnr_file} gatkcnv_cnts={gatkcnv_cnts}")
        cmd_prep_cnr = theta_create_exome_input(
            segfile = segfile,
            tumor_bam = trk_tum$bam[i],
            normal_bam = trk_norm$bam[1], 
            theta_oprefix = theta_oprefix,
            thetadir = thetadir)
        cmd_prep_cnr = glue("{cmd_prep_cnr};cp {segfile} {thetadir}/")

        hets_files = c(trk_tum$gatkcnv_hets[i], trk_tum$gatkcnv_hets_n[i])
        theta_hets_files = c(trk_tum$theta_hets_file[i], trk_tum$theta_hets_n_file[i])
        cmd_prep_hets = glue("{funr_exe} my.ultraseq::to_theta_het_file.gatkcnv hets_file={hets_files} theta_hets_file={theta_hets_files}")
        cmd_prep_hets

        # install.my.ultra
        cmds_prep = c(cmd_prep_cnr, cmd_prep_hets)
        # lapply(cmds_prep, system)

        # RunTHetA.py
        # theta_exe = "conda activate theta_py27;RunTHetA.py"
        # theta_opts = "--NUM_PROCESSES 5"
        cmd_theta = glue("{theta_exe} {theta_cnr_file} --TUMOR_FILE {theta_hets_files[1]} --NORMAL_FILE {theta_hets_files[2]} -p {theta_oprefix} {theta_opts} --BAF")
        cmd_theta
        cmd_parse_theta = glue("{funr_exe} my.ultraseq::read_theta.results oprefix={theta_oprefix} path=.")

        # theta_baf_model_exe = "python ~/apps/theta/THetA/python/RunBAFModel.py"
        # cmd_baf = glue("{theta_baf_model_exe} {theta_hets_files[1]} {theta_hets_files[2]} {theta_cnr_file} {oprefix.BEST.results} -O {thetadir}")
        # to_theta_cn_file.gatkcnv(cnr_files, gatkcnv_cnts, theta_cnr_file)
        # to_theta_het_file.gatkcnv(hets_files[1], theta_hets_files[1])
        # to_theta_het_file.gatkcnv(hets_files[2], theta_hets_files[2])

        cmds <- list(cmds_prep = cmds_prep,
                     cmd_theta = cmd_theta, 
                     cmd_parse_theta = cmd_parse_theta,
                     outfiles = list(trk_tum$theta_results))
    })
    lst = purrr::transpose(tmp)
    cmds = list(theta.prep = cmd_pileup, 
                theta.mrg = c(unlist(lst$cmds_prep), 
                              unlist(lst$cmd_theta), 
                              unlist(lst$cmd_parse_theta), 
                              cmd_rm_pileup))
    # write_tsv(df_trk, path = "gatkcnv/df_trk.tsv")

    # cmds %>% unlist() %>% paste0(collapse = "\n") %>% cat()
    lst$flowmat = to_flowmat(cmds, samplename = samplename) %>%
        mutate(cmd = as.character(cmd))
    lst$trk = trk
    return(lst)
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

# etadata (current_repodata.json): done
# https://github.com/biod/sambamba/wiki/%5Bsambamba-mpileup%5D-documentation
theta_create_mpileup <- function(bam, odir = pipe_str$theta$dir, 
                                quality = 30,
                                sambamba_exe = "/rsrch3/home/iacs/sseth/apps/sambamba/sambamba-0.7.1-linux-static",
                                ref_fasta = opts_flow$envir$ref_fasta,
                                cores = opts_flow$envir$theta_pileup_cores,
                                bed = opts_flow$envir$gatkcnv.intervals_bed
                                ){
    # bam = trk$bam[2] %>% basename()
    check_args()
    prefix = basename(bam) %>% tools::file_path_sans_ext()
    pileup = glue("{odir}/{prefix}.pileup")
    log = glue("{pileup}.log")

    cmd = glue("module load samtools/1.10 bcftools/1.10.2;",
            "{sambamba_exe} mpileup -t {cores} -L {bed} {bam} --samtools '-f {ref_fasta} -q {quality}' > {pileup};echo $? > {log}")
    cmd
       
}

theta_rm_pileup <- function(bam, odir = pipe_str$theta$dir){

    prefix <- basename(bam) %>% tools::file_path_sans_ext()
    pileup <- glue("{odir}/{prefix}.pileup")
    cmd = glue("rm -rf {pileup}")
    cmd
}

theta_create_exome_input <- function(
    segfile = trk_tum$gatkcnv_cr_seg,
    tumor_bam = trk_tum$bam[1],
    normal_bam = trk_norm$bam[1],
    theta_oprefix, thetadir,
    theta_exome_input_exe = opts_flow$envir$theta_exome_input_exe,
    ref_fasta = opts_flow$envir$ref_fasta,
    gatkcnv_intervals_bed = opts_flow$envir$gatkcnv_intervals_bed
){
    check_args()
    #  --FA /rsrch3/home/iacs/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta --EXON_FILE /rsrch3/home/iacs/sseth/ref/az_ref_beds/ss_downloads/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets_pad250.bed 
    # only uses start stop positions
    theta_oprefix = basename(as.character(theta_oprefix))
    cmd = glue("{theta_exome_input_exe} -s {segfile} -t {tumor_bam} -n {normal_bam} --FA {ref_fasta} --EXON_FILE {gatkcnv_intervals_bed} --OUTPUT_PREFIX {theta_oprefix} --DIR {thetadir} --QUALITY 30")
    return(cmd)
}

# the results from this are NOT the same, unfortunately
theta_create_exome_input.gatk <- function(){
    setwd("~/flows/SS/tnbc/ms51_wex/ss_cl_het/runs/185_003/tmp")
    t_cnts = "gatkcnv/IPCT-S2006-MOON0051-Cap2023-8-ID01.counts.tsv" %>% 
        read_counts() %>% mutate(width = end-start)
    n_cnts = "gatkcnv/IPCT-S4008-MOON0051-Cap2043-4-ID01.counts.tsv" %>%
        read_counts() %>% mutate(width = end - start)
    seg = "gatkcnv/IPCT-S4008-MOON0051-Cap2043-4-ID01___matched.called.igv.seg" %>%
        read_counts()

    # for every seg, merge counts:
    i=2
    seg2 = lapply(1:nrow(seg), function(i){
        segi = seg[i,]
        tumor_count = filter(t_cnts, 
                            contig == segi$chromosome,
                            start >= segi$start,
                            end <= segi$end) %>% 
                    summarize(count = sum(count))
        normal_count = filter(n_cnts, 
                            contig == segi$chromosome,
                            start >= segi$start,
                            end <= segi$end) %>% 
                    summarize(count = sum(count))
        segi$tumorCount = tumor_count$count
        segi$normalCount = normal_count$count
        segi
    }) %>% bind_rows()
    seg2 %>% mutate(ID = glue("start_{chromosome}_{start}:end_{chromosome}_{end}")) %>% 
        select(chrm = chromosome, start, end, tumorCount, normalCount, segment_mean) %>% 
        data.frame() %>% head()
    

}

.read_theta.results <- function(x, y, df_seg, oprefix){
    res = data.table::fread(x, data.table = F) %>% 
        clean_names()
    # res2 = read_tsv(y) %>% clean_names

    # neg log lik of the fit
    nll = res$number_nll

    # purity
    mu = res$mu %>% stringr::str_split(pattern = ",") %>% unlist()
    mu_normal = mu[1] %>% as.numeric()
    mu_tumor = mu[-1] %>% as.numeric()

    # inferred interval CN
    cn = res$c %>% stringr::str_split(pattern = ":") 
    df_seg$cn = cn[[1]]
    df_seg %<>% tidyr::separate(cn, c("cn_cl1", "cn_cl2"), sep = ",", fill = "right") %>% 
        dplyr::mutate(samplename = oprefix, cn_cl1 = as.integer(cn_cl1), cn_cl2 = as.integer(cn_cl2)) %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(total_cn = sum(cn_cl1, cn_cl2, na.rm = T), num_marks = NA, num_clones = length(mu_tumor)) %>% 
        dplyr::select(samplename, chr = chrm, start, end, num_marks, total_cn, everything())

    # multinomial param:
    p = res$p %>% stringr::str_split(pattern = ",") %>% unlist() %>% as.numeric()
    df_seg$p = p

    readr::write_tsv(df_seg, y)

    list(df_seg = df_seg, pct_purity = mu_tumor, pct_normal = mu_normal, nll = nll)
}

# results files: the set of all maximum likelihood solutions found.
# 	  Each solution will be a single line in this file with fields:
# 		1. NLL (Double) - solution neg. log likelihood (w/o constants)
# 		2. mu  (Double,Double) - comma delimited list:  %normal, %tumor
# 		3. C_2 (Integer: ... :Integer) - colon delimited list of
# 		   inferred tumor interval copy numbers.  When multiple tumor
# 		   populations are inferred this list is in the form
# 		   (Integer,Integer: ... : Integer, Integer).
# 		4. p* (Double, ... ,Double) - comma delimted list of inferred
# 		   multinomial parameter \widehat{C\mu}
# 		The chosen solution will be in a file with postfix ".BEST.results".
# 		THetA will also output the individual solutions for n=2 and n=3 in
# 		in files ending in ".n3.results" and ".n3.results" respectively

#' to_theta_cn_file.gatkcnv
#'
#' @param oprefix expecting files named `{path}/{oprefix}.input"`,  `{path}/{oprefix}.n2.results"` etc.
#' @param path output dir, can be ".", or ""
#'
#' @export
read_theta.results <- function(oprefix, path){
    # path = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het/runs/185_057/tmp/theta"
    # oprefix = "IPCT-S4013-MOON0051-Cap2531-4-HTID254___matched"
    # list.files(path)

    # setwd(path)
    input = glue("{path}/{oprefix}.input")
    flog.debug(glue("working on: {input}"))

    # https://software.broadinstitute.org/software/igv/sites/cancerinformatics.org.igv/files/linked_files/example.seg
    clnms = c("n2", "n3")
    xs = glue("{path}/{oprefix}.{clnms}.results") %>% as.character()
    ys = glue("{path}/{oprefix}.{clnms}.igv.seg") %>% as.character()
    fitfl = glue("{path}/{oprefix}.df_fit_summary.tsv") %>% as.character()

    df_seg <- read_counts(input) %>% clean_names()
    lst = purrr::map2(xs, ys, .read_theta.results, df_seg, oprefix)
    lst = purrr::transpose(lst)
    lst$df_seg %>% bind_rows()

    df_mod = data.frame(n_clones = c(2, 3),
                        pct_normal = lst$pct_normal %>% unlist(),
                        nll = lst$nll %>% unlist(), 
                        n_segments = nrow(df_seg)) %>% as_tibble() %>% 
            mutate(pct_normal = round(pct_normal, 3), 
                   seg = basename(ys),
                   oprefix = oprefix, 
                   final_model = nll == max(nll)) %>% arrange(-nll)
    df_mod
    write_tsv(df_mod, fitfl)

    list(df_mod = df_mod, df_seg = df_seg)
}

if(FALSE){
    con = db_connect_art(type = "postgres")
    dbListTables(con)
    p_load(dbx)
    p_load(purrr)

    # 23 samples fail because of lack of enough CNV
    df_trk_tum = tidylog::filter(df_trk, normal == "NO") %>% 
        mutate(theta_results_full = glue("{individual}/tmp/{theta_results}")) %>% 
        tidylog::filter(file.exists(theta_results_full))
    getwd()
    flog.threshold(DEBUG)
    lst = future_map2(df_trk_tum$outprefix_paired, file.path(df_trk_tum$individual, "tmp/theta"), read_theta.results, .progress = TRUE)
    lst2 = transpose(lst)
    df_mod = lst2$df_mod %>% bind_rows()
    df_mod2 = df_mod %>% filter(final_model)

    # add these to the DB
    df_mod2 %>% mutate(tool = "theta2", version = "1.1", purity = 1-pct_normal) %>% 
        left_join(df_trk_tum, by = c("oprefix" = "outprefix_paired")) %>% 
        select(path_patient_id = individual, sample_type, seq_sample_id = outprefix, outprefix_paired = oprefix, 
                tool, version,
                n_clones, purity, n_segments, score = nll) %>% 
        dbxUpsert(con, "wex_cnv_insights", ., where_cols = c("outprefix_paired","tool", "version"))
    write_rds(lst, "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het/theta_v1.1/lst_theta_results.rds")
    # then add titan to DB as well
}



# EXTRA ------
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