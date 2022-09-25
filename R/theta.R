# install:
# conda create -n theta_py27 python=2.7 theta2
# https://bioconda.github.io/recipes/theta2/README.html

# get the repo as well:
# git clone https://github.com/raphael-group/THetA.git



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