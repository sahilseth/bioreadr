
# singularity pull docker://quay.io/bcbio/bcbio-vc


# https://github.com/gavinha/TitanCNA/tree/master/scripts/R_scripts
# numClusters=3
# numCores=4
# for ploidy in [2, 3, 4]:
# for num_clusters in [1, 2, 3]:
# ## run TITAN for each ploidy (2,3,4) and clusters (1 to numClusters)
# echo "Maximum number of clusters: $numClusters";
# for ploidy in $(seq 2 4)
# do
# echo "Running TITAN for $i clusters.";
# outDir=run_ploidy$ploidy
# mkdir $outDir
# for numClust in $(seq 1 $numClusters)
# do
# echo "Running for ploidy=$ploidy";
# Rscript titanCNA_v1.10.1.R --id test --hetFile test.het.txt --cnFile test.cn.txt \
# --numClusters $numClust --numCores $numCores --normal_0 0.5 --ploidy_0 $ploidy \
# --chrs "c(1:22, \"X\")" --estimatePloidy TRUE --outDir $outDir
# done
# echo "Completed job for $numClust clusters."
# done
# 
# ## select optimal solution
# Rscript selectSolution.R --ploidyRun2=run_ploidy2 --ploidyRun3=run_ploidy3 --ploidyRun4=run_ploidy4 --threshold=0.05 --outFile optimalClusters.txt

.run_titan_example <- function(){
  
  # test trk file:
  trk = "~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex_b1/superfreq/df_sfq_metadata_ms51_dna_b1.tsv"
}


# https://github.com/broadinstitute/gatk/blob/master/scripts/unsupported/combine_tracks_postprocessing_cnv/combine_tracks.wdl
# https://github.com/bcbio/bcbio-nextgen/blob/520fad7da7fa46635e6514b8bb6b12b988cc20ac/bcbio/structural/titancna.py#L162
# Convert GATK hets from ModelSegments to Titan input:
# cat SAMPLE.hets.tsv | awk '/^[^@]/ { print $1,$2,$5,$3,$6,$4 }' | tr ' ' '\t' > tmp.hets.tsv 

#' .titan_cn_file
#' 
#' Convert CNVkit or GATK4 normalized input into TitanCNA ready format.
#'
#' @param cnr_file 
#' 
#' @export
to_titan_cn_file <- function(cnr_file, 
                             titan_cnr_file,
                             file_type = c("gatkcnv", "cnvkit")){
  
  df_cnr = read_tsv(cnr_file, comment = "@", col_types = cols(.default = col_character()))
  
  file_type = match.arg(file_type)
  support_cols = list("cnvkit" = c("chromosome", "start", "end", "log2"),
                      "gatkcnv" = c("CONTIG", "START", "END", "LOG2_COPY_RATIO"))
  cnr_cols = support_cols[[file_type]]
  df_cnr = df_cnr[, support_cols[[file_type]]]
  # std col names
  # NAME BASED
  colnames(df_cnr) = c("chr", "start", "end", "log2_copy_ratio")
  # library(GenomicRanges)
  # library(TitanCNA)
  # df_cnr$chr <- setGenomeStyle(df_cnr$chr, genomeStyle = "NCBI")
  # gr_cnr <- data.frame(df_cnr, stringsAsFactors = F) %>% as("GRanges")		
  
  write_tsv(df_cnr, titan_cnr_file)
  
}


#' to_titan_hets_file
#'
#' @param hets_file 
#' @param titan_hets_file 
#'
#' @export
to_titan_hets_file <- function(hets_file, titan_hets_file){
  
  # cat SAMPLE.hets.tsv | awk '/^[^@]/ { print $1,$2,$5,$3,$6,$4 }' | tr ' ' '\t' > tmp.hets.tsv 
  
  # POSITON BASED
  df_hets = read_tsv(hets_file, comment = "@", col_types = cols(.default = col_character()))
  df_hets %<>% select(chr = CONTIG, pos = POSITION, 
                      ref_allele = REF_NUCLEOTIDE, 
                      ref_count = REF_COUNT,
                      alt_allele = ALT_NUCLEOTIDE, 
                      alt_count = ALT_COUNT)
  
  write_tsv(df_hets, titan_hets_file)
}


#' Title
#'
#' @param trk 
#' @param samplename 
#' @param num_cores 
#' @param test 
#'
titancna_somatic_matched <- function(trk,
                                     samplename,
                                     # num_cores = nrow(df_trk), 
                                     run_cmds = F,
                                     funr_exe = opts_flow$get("R_funr_exe"),
                                     num_cores = 4,
                                     rscript_exe = opts_flow$envir$R.sing_rscript,
                                     titandir = pipe_str$ttn$dir
){
  
  # check generic columns
  trk = metadata_for_dnaseq_tools(trk)
  # we would use 
  # WEX-sarco14-T.denoisedCR.tsv instead of WEX-sarco14-T___matched.cr.seg
  # coz titan segments it
  expect_columns(trk, c("outprefix", "outprefix_paired",
                        "cnr_file", "hets_file",
                        "gender"))
  
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
      titan_cnr_file = glue("{titandir}/{outprefix}.denoisedCR_titan.tsv"),
      titan_hets_file = glue("{titandir}/{outprefix_paired}.hets_titan.tsv")
    )
  trk
  
  trk_tum = filter(trk, normal == "NO")
  trk_norm = filter(trk, normal == "YES")
  
  # ** prep inputs -------
  # opts_flow$load_toml("../flowr_conf.toml")
  cnr_files = trk_tum$cnr_file
  titan_cnr_files = trk_tum$titan_cnr_file
  cmd_prep_cnr = glue("{funr_exe} my.ultraseq::to_titan_cn_file cnr_file={cnr_files} titan_cnr_file={titan_cnr_files} file_type=gatkcnv")
  
  hets_files = trk_tum$hets_file
  titan_hets_files = trk_tum$titan_hets_file
  cmd_prep_hets = glue("{funr_exe} my.ultraseq::to_titan_hets_file hets_file={hets_files} titan_hets_file={titan_hets_files} file_type=gatkcnv")
  # install.my.ultra
  cmds_prep = c(cmd_prep_cnr, cmd_prep_hets)
  
  # ttn.split
  titancna="$HOME/Dropbox/public/github_titancna/scripts/R_scripts/titanCNA.R"
  centromere_gap_fl = "$HOME/ref/human/b37/annotations/broad-somatic-b37/GRCh37.p13_centromere_UCSC-gapTable.txt"
  
  # --chrs 'c(1:22, \"X\")'
  ploidies = c(2:4)
  clusters = c(1:10)
  
  # tic("> starting titan cna")
  # cpu_titan=24 # assumming we have 24 cores
  
  i = 1
  ttn_cmds = lapply(1:nrow(trk_tum), function(i){
    
    name = trk_tum$name[i]
    target = glue("{titandir}/{name}_optimalClusters.txt")
    titan_hets_file = trk_tum$titan_hets_file[i]
    titan_cnr_file = trk_tum$titan_cnr_file[i]
    gender = trk_tum$gender[i]
    tmp = mclapply(ploidies, function(ploidy){
      mclapply(clusters, function(cluster){
        # cluster = 2; ploidy = 2
        # estimateNormal: method to estimate normal map 
        outdir = glue("{titandir}/run_ploidy{ploidy}")
        # wranglr::mkdir(outdir)
        cmd_titancna = glue("mkdir -p {outdir};{rscript_exe} {titancna} --id {name} --hetFile {titan_hets_file} --cnFile {titan_cnr_file} ", 
                            "--numClusters {cluster} --numCores {num_cores} --normal_0 0.5 --ploidy_0 {ploidy} ",
                            "--estimateNormal map --estimatePloidy TRUE --estimateClonality TRUE --outDir {outdir} ", 
                            # Hyperparameter on Gaussian variance; for WES, use 2500; for WGS, use 10000; float (Default: 10000
                            "--gender {gender} --alphaK 2500 --alphaKHigh 2500 ", 
                            "--centromere {centromere_gap_fl} >> {outdir}_{name}.out")
        
        if(run_cmds){
          cluster = sprintf("%02d", cluster)
          cmdname = glue("{outdir}/{name}_cluster{cluster}")
          tmp = run_cmd(cmd_titancna, 
                        stderr = glue("{cmdname}.log"),
                        cmdname = cmdname, 
                        target = glue("{cmdname}.params.txt"), redo = redo)
        }
        cmd_titancna
      }, mc.cores = 1)
    }, mc.cores = 1)
  })
  ttn_splt_cmds = unlist(ttn_cmds)
  
  
  # ttn.mrg
  titancna_selectsolution = "$HOME/Dropbox/public/github_titancna/scripts/R_scripts/selectSolution.R"
  cmd_titancna_select_solutions = glue("{rscript_exe} {titancna_selectsolution} ",
                                       "--ploidyRun2={titandir}/run_ploidy2 --ploidyRun3={titandir}/run_ploidy3 ",
                                       "--ploidyRun4={titandir}/run_ploidy4 --threshold=0.05 ",
                                       "--outFile {titandir}/optimalClusters.txt --outFileFull {titandir}/allClusters.txt")
  # system(cmd_titancna_select_solutions)
  
  
  # write_tsv(df_trk, path = "gatkcnv/df_trk.tsv")
  cmds <- list(titancna.prep = cmds_prep,
               titancna.splt = unlist(ttn_splt_cmds),
               titancna.mrg = unlist(cmd_titancna_select_solutions))
  
  # cmds %>% unlist() %>% paste0(collapse = "\n") %>% cat()
  flowmat = to_flowmat(cmds, samplename = samplename) %>%
    mutate(cmd = as.character(cmd))
  
  list(flowmat = flowmat, outfiles = list(), trk = trk)
}


titancna_merge_opt_cluster <- function(df_trk){
  i=1
  tmp = lapply(1:nrow(df_trk), function(i){
    opt_clus = df_trk$ttn_optimal_fl[i] %>% 
      read_tsv(col_types = cols(
        Phi = col_double(),
        id = col_character(),
        barcode = col_character(),
        numClust = col_double(),
        cellPrev = col_character(),
        purity = col_double(),
        norm = col_double(),
        ploidy = col_double(),
        loglik = col_double(),
        sdbw = col_double(),
        path = col_character()
      )) %>% 
      clean_names()
  }) %>% bind_rows()
  # tmp
  df_trk = left_join(df_trk, tmp, by = c("name" = "barcode"))
}

titan_cna_states <- function(){
  # based on discussions here
  # https://github.com/gavinha/TitanCNA/issues/8
  # read.delim(pipe("pbpaste")) %>% dput()
  
  df_titan_states = structure(list(titan_state = -1:24, 
                                   titan_genotype = c("NULL", 
                                                      "NULL", "A", "AA", "AB", "AAA", "AAB", "AAAA", "AAAB", "AABB", 
                                                      "AAAAA", "AAAAB", "AAABB", "AAAAAA", "AAAAAB", "AAAABB", "AAABBB", 
                                                      "AAAAAAA", "AAAAAAB", "AAAAABB", "AAAABBB", "AAAAAAAA", "AAAAAAAB", 
                                                      "AAAAAABB", "AAAAABBB", "AAAABBBB"), 
                                   titan_total_copy_number = c("NULL", 
                                                               "0", "1", "2", "2", "3", "3", "4", "4", "4", "5", "5", "5", "6", 
                                                               "6", "6", "6", "7", "7", "7", "7", "8", "8", "8", "8", "8"), 
                                   titan_call = c("OUT", "HOMD", "DLOH", "NLOH", "HET", "ALOH", 
                                                  "GAIN", "ALOH", "ASCNA", "BCNA", "ALOH", "ASCNA", "UBCNA", 
                                                  "ALOH", "ASCNA", "UBCNA", "BCNA", "ALOH", "ASCNA", "UBCNA", 
                                                  "UBCNA", "ALOH", "ASCNA", "UBCNA", "UBCNA", "BCNA"), 
                                   titan_corrected_call = c("OUT", 
                                                            "HOMD", "DLOH", "NLOH", "HET", "ALOH", "GAIN", "ALOH", "ASCNA", 
                                                            "BCNA", "ALOH", "ASCNA", "UBCNA", "ALOH", "ASCNA", "UBCNA", 
                                                            "BCNA", "ALOH", "ASCNA", "UBCNA", "UBCNA", "HLAMP", "HLAMP", 
                                                            "HLAMP", "HLAMP", "HLAMP"),
                                   titan_corrected_call2 = c("OUT", 
                                                             "HOMD", "HETD", "NLOH", "NEUT", "ALOH", "GAIN", "ALOH", "ASCNA", 
                                                             "BCNA", "ALOH", "ASCNA", "UBCNA", "ALOH", "ASCNA", "UBCNA", 
                                                             "BCNA", "ALOH", "AMP", "AMP", "AMP", "HLAMP", "HLAMP", 
                                                             "HLAMP", "HLAMP", "HLAMP"),
                                   titan_call_description = c("Outlier state", 
                                                              "Homozygous deletion", "Hemizygous deletion LOH", "Copy neutral LOH", 
                                                              "Diploid heterozygous", "Amplified LOH", "Gain/duplication of 1 allele", 
                                                              "Amplified LOH", "Allele-specific copy number amplification", 
                                                              "Balanced copy number amplification", "Amplified LOH", "Allele-specific copy number amplification", 
                                                              "Unbalanced copy number amplification", "Amplified LOH", 
                                                              "Allele-specific copy number amplification", "Unbalanced copy number amplification", 
                                                              "Balanced copy number amplification", "Amplified LOH", "Allele-specific copy number amplification", 
                                                              "Unbalanced copy number amplification", "Unbalanced copy number amplification", 
                                                              "Amplified LOH", "Allele-specific copy number amplification", 
                                                              "Unbalanced copy number amplification", "Unbalanced copy number amplification", 
                                                              "Balanced copy number amplification"), 
                                   titan_corrected_call_description = c("Outlier state", 
                                                                        "Homozygous deletion", "Hemizygous deletion LOH", "Copy neutral LOH", 
                                                                        "Diploid heterozygous", "Amplified LOH", "Gain/duplication of 1 allele", 
                                                                        "Amplified LOH", "Allele-specific copy number amplification", 
                                                                        "Balanced copy number amplification", "Amplified LOH", "Allele-specific copy number amplification", 
                                                                        "Unbalanced copy number amplification", "Amplified LOH", 
                                                                        "Allele-specific copy number amplification", "Unbalanced copy number amplification", 
                                                                        "Balanced copy number amplification", "Amplified LOH", "Allele-specific copy number amplification", 
                                                                        "Unbalanced copy number amplification", "Unbalanced copy number amplification", 
                                                                        "High-level amplification", "High-level amplification", "High-level amplification", 
                                                                        "High-level amplification", "High-level amplification")), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                      -26L))
  df_titan_states %>% dplyr::mutate(titan_total_copy_number = as.integer(titan_total_copy_number))
}

# ** read final titiancna-------
if(FALSE){
  
  wexpath = "/rsrch3/scratch/iacs/sseth/flows/SS/sarco/mda/wex"
  titanpath = "titancna"
  inds = df_trk$individual %>% unique()
  cmds = glue("rsync -a {wexpath}/runs/{inds}/tmp/{titanpath}/ {wexpath}/{titanpath}/{inds}/");cmds
  # system(cmds)
  
  df_trk_ttn = df_trk %>% filter(normal == "NO") %>% 
    mutate(ttn_optimal_fl = glue("{wexpath}/{titanpath}/{individual}/optimalClusters.txt")) %>% 
    titancna_merge_opt_cluster()
  write_tsv(df_trk_ttn, glue("{wexpath}/{titanpath}/df_trk_ttn.tsv"))
  
  # read seg files:
  df_trk_ttn = read_tsv(glue("{wexpath}/{titanpath}/df_trk_ttn.tsv"))
  
  i = 1
  
  lapply(1:nrow(df_trk_ttn), function(i){
    
    segfl = df_trk_ttn$path[i]
    segfl_ploidy = dirname(segfl) %>% basename()
    segfl_cluster = segfl %>% basename()
    ind = df_trk_ttn$individual[i]
    segfl = glue("{wexpath}/{titanpath}/{ind}/{segfl_ploidy}/{segfl_cluster}.segs.txt")
    file.exists(segfl)
    
    source('~/Dropbox/projects/packs_dnaseq/R/to_mae.seg.R')
    
    lst = to_seg.titan(segfl = segfl)
    invisible(lst)
  })
  # END
  
  
}


to_seg.titan <- function(segfl){
  
  seg = data.table::fread(segfl, data.table = F) %>% as_tibble() %>% 
    clean_names()
  seg
  seg$chromosome %>% table()
  
  # # http://seqanswers.com/forums/archive/index.php/t-32100.html
  # message("switch chr names")
  # seg <- seg %>% mutate(major_cn = tcn_em-lcn_em) %>% 
  #   dplyr::rename(chr = chrom, start = start, end = end, total_cn =  tcn_em, 
  #                 major_cn = major_cn, minor_cn = lcn_em, everything()) %>% 
  #   mutate(chr = gsub("23", "X", chr),
  #          chr = gsub("24", "Y", chr))
  # gr_seg <- makeGRangesFromDataFrame(seg, keep.extra.columns = T)
  
  lst = seg %>% to_gr_seg.seg(col_sample = "sample", 
                              col_chr = "chromosome", 
                              col_start = "start_position_bp", col_end = "end_position_bp", 
                              col_num_mark = "length_snp", col_seg_mean = "corrected_copy_number")
  message("annotate")
  source('~/Dropbox/projects/packs_dnaseq/R/to_mae.seg.R')
  length(lst$gr_seg)
  df_ann = annotate.gr_seg(lst$gr_seg, gencode_fl = "~/.rsrch3/home/iacs/sseth/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene_hgnc.bed") %>% 
    dplyr::mutate(corrected_copy_number = seg_mean)
  dim(df_ann)
  df_ann = df_ann %>% group_by(gene_name) %>% add_count()
  # df_ann %>% filter(n>1) %>% View()
  
  list(df_ann = df_ann, df_igv_seg = lst$df_igv_seg, gr_seg = lst$gr_seg)
  
}

