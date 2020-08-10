


# The depth.ratio column is the GC content normalized ratio. a depth ratio of 1 means it has copy number of 2 (the same as the normal blood control in my case).
# segment file:
# (1) Sample (sample name)
# (2) Chromosome (chromosome number)
# (3) Start Position (segment start position, in bases)
# (4) End Position (segment end position, in bases)
# (5) Num markers (number of markers in segment)
# (6) Seg.CN (log2() -1 of copy number)
to_gistic.exomecn <- function(df_exomecn_seg){
    # mean seems to be 0, so it seems we have already done
    # log2(CN)
    # https://www.genepattern.org/modules/docs/GISTIC_2.0
    # Seg.CN: log2(CN) -1)
    # Seg.CN: log2(TDP/NDP)
    # 0: no diff
    # -1: 1 copy loss
    # 1: 1 copy gain
    summary(df_exomecn_seg$seg_mean)
    fl = tempfile(fileext = ".tsv")
    dplyr::select(df_exomecn_seg, 
        samplename, chrom, loc_start, loc_end, num_mark, seg_mean) %>%
        write_tsv(fl)
    fl
}


# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
# gistic2 \
# -b <base_directory> \
# -seg <segmentation_file> \
# -mk <marker_file> \
# -refgene <reference_gene_file> \
# -ta 0.1 \
# -armpeel 1 \
# -brlen 0.7 \
# -cap 1.5 \
# -conf 0.99 \
# -td 0.1 \
# -genegistic 1 \
# -gcm extreme \
# -js 4 \
# -maxseg 2000 \
# -qvt 0.25 \
# -rx 0 \
# -savegene 1 \
# (-broad 1)
run_gistic <- function(input, output){
    output = tools::file_path_as_absolute(output)
    wranglr::mkdir(output)
    # #!/bin/sh
    # ## set MCR environment and launch GISTIC executable

    # ## NOTE: change the line below if you have installed the Matlab MCR in an alternative location
    # MCR_ROOT=/scratch/genomic_med/apps/Matlab_Complier_runTime
    # MCR_VER=v83

    # echo Setting Matlab MCR root to $MCR_ROOT

    # ## set up environment variables
    # LD_LIBRARY_PATH=$MCR_ROOT/$MCR_VER/runtime/glnxa64:$LD_LIBRARY_PATH
    # LD_LIBRARY_PATH=$MCR_ROOT/$MCR_VER/bin/glnxa64:$LD_LIBRARY_PATH
    # LD_LIBRARY_PATH=$MCR_ROOT/$MCR_VER/sys/os/glnxa64:$LD_LIBRARY_PATH
    # export LD_LIBRARY_PATH
    # XAPPLRESDIR=$MCR_ROOT/$MCR_VER/MATLAB_Component_Runtime/v83/X11/app-defaults
    # export XAPPLRESDIR

    # ## launch GISTIC executable
    # ./gp_gistic2_from_seg $@
    gistic_exe = "module load gistic/2.0;gistic2"
    # -mk markers_gistic.txt
    refgenemat = "/risapps/rhel7/gistic/2.0/refgenefiles/hg19.mat"
    gistic_opts = "-genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme"
    gisticdir = "/risapps/rhel7/gistic/2.0"
    cmd = glue("cd {gisticdir};{gistic_exe} -b {output} -seg {input} -refgene {refgenemat} {gistic_opts}")
    system(cmd)

    lst = read_gistic(output)
    write_rds(lst, glue("{output}/lst_maftools.rds"))

}

read_gistic <- function(output){

    p_load(maftools)
    lst = readGistic(gisticAllLesionsFile = glue("{output}/all_lesions.conf_90.txt"), 
                    gisticAmpGenesFile = glue("{output}/amp_genes.conf_90.txt"), 
                    gisticDelGenesFile = glue("{output}/del_genes.conf_90.txt"), 
                    gisticScoresFile = glue("{output}/scores.gistic"))

    pdf(glue("{output}/p_chromplot.pdf"), width = 20)
    gisticChromPlot(gistic = lst, markBands = "all")
    dev.off()

    pdf(glue("{output}/p_bubble.pdf"), width = 8, height = 8)
    gisticBubblePlot(gistic = lst)
    dev.off()

    pdf(glue("{output}/p_onco.pdf"), width = 20, height = 20)
    gisticOncoPlot(gistic = lst)
    dev.off()
    lst

}