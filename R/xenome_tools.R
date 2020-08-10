## sseth@mdanderson.org

if(FALSE){

    ## /scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java -Xmx4g -XX:-UseGCOverheadLimit -XX:-UseParallelGC -Djava.io.tmpdir=java_tmp -jar /scratch/rists/hpcapps/x86_64/picard/1.112/SamToFastq.jar INPUT=/scratch/iacs/gcc/levelii/140603_SN208_0517_BC4JL7ACXX/LYNDA-FORD-HF-HW-24-60-01R_140603_SN208_0517_BC4JL7ACXX_s_6_CATACCAA.rg.sorted.recalibed.bam FASTQ=read1.fq SECOND_END_FASTQ=read2.fq VALIDATION_STRINGENCY=LENIENT

    ##  /scratch/rists/hpcapps/x86_64/xenome/1.0.1-r/xenome classify --num-threads 24 -P /IACS1/NGS/hg19-mm9-xenome-index/hg19-mm9-xenome-index --pairs --graft-name hg19 --host-name mm9 --output-filename-prefix HF-HW24-60-01D -i read1.fq -i read2.fq -v

    ## ~/projects/tools_xenome/xenome.to.fastq.sh HF-HW24-60-01D_HF-HW24-60-01D_mm9_2.fastq > HF-HW24-60-01D_HF-HW24-60-01D_mm9_end2.fq

    ##  export PYTHONPATH=$PYTHONPATH:/IACS1/lib/python2.7/site-packages/; ~/projects/tools_xenome/filter.xeno.mmuReads.sseth.py ./ HF-HW24-60-01D /scratch/iacs/gcc/levelii/140603_SN208_0517_BC4JL7ACXX/LYNDA-FORD-HF-HW-24-60-01R_140603_SN208_0517_BC4JL7ACXX_s_6_CATACCAA.rg.sorted.recalibed.bam /scratch/iacs/gcc/levelii/140603_SN208_0517_BC4JL7ACXX/LYNDA-FORD-HF-HW-24-60-01R_140603_SN208_0517_BC4JL7ACXX_s_6_CATACCAA.rg.sorted.recalibed.xeno_bam

    ## The following samples in 0284,

    patterns = c('Cap454-101029-5', 'Cap454-101030-6',
        'Cap453-100458-2', 'Cap453-100466-1',
        'Cap453-101024-3','Cap454-101028-4')
    runpath = "/scratch/iacs/ngs_runs/141013_SN959_0284_AC5JR0ACXX"
    odir = file.path(runpath, "xenome")
    files = list.files(file.path(runpath, "bams"), pattern = "recalibed.bam$", full.names = TRUE)
    bams = sapply(patterns, grep, files, value = TRUE)
    # require(SaturnV)
    ##source("~/iacsSVN/RPacks/SaturnV/R/xenome_tools.R")
    ##debug(runXenome)
    require(tools)
    require(parallel)
                                        #cmds <- runXenome(bam = bams[i], fq1 = fq1, fq2 = fq1, odir = odir, execute = TRUE, verbose = TRUE)

    sys_cmds <- sapply(bams, function(b){
        fq1 = gsub("star_recalibed.bam", "R1.fastq", b);file.exists(fq1)
        fq2 = gsub("star_recalibed.bam", "R2.fastq", b);file.exists(fq1)
        ## RNA example:
        sprintf("run_module runXenome bam_C='%s' fq1_C=%s fq2_C=%s odir_C=%s execute_L=TRUE",
                b, fq1, fq2, odir)
    })

    require(flow)
    qobj <- queue(queue="long", type="torque")
    jobj <- job(cmds=sys_cmds[2], q_obj = qobj, cpu = 24, submission_type = "scatter", memory = "24g", walltime = "72:00:00")
    fobj <- flow(jobs=list(jobj))
    ## debug(flow:::.submit_job)
    fobj2 <- flow:::.submit_flow(fobj, execute = TRUE)


}


#' @title runXenome
#' @description runXenome runXenome
#' @param bam bam
#' @param fq1 fq1
#' @param fq2 something
#' @param outbam something
#' @param odir something
#' @param picard_dir something
#' @param xenome_exe something
#' @param xen_index something
#' @param java something
#' @param javaMem something
#' @param xen_cores something
#' @param xenome_fastq_exe something
#' @param python_exe something
#' @param filter_reads_exe something
#' @param execute something
#' @param verbose something
#' @export
#' @examples \dontrun{
#' runXenome(bam = bam, fq1 = fq1, fq2 = fq2, outbam = outbam, odir = odir, picard_dir = /scratch/rists/hpcapps/x86_64/picard/1.112, xenome_exe = /scratch/rists/hpcapps/rhel6/xenome/1.0.1-r/xenome, xen_index = /IACS1/NGS/hg19-mm9-xenome-index/hg19-mm9-xenome-index, java = /scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java, javaMem = -Xmx8g, xen_cores = 24, xenome_fastq_exe = /scratch/iacs/bin/xenome/xen_to_fastq.sh, python_exe = /scratch/rists/hpcapps/x86_64/Python/2.7.6-new/bin/python, filter_reads_exe = /scratch/iacs/bin/xenome/filter_xeno_mmReads.py, execute = TRUE, verbose = FALSE)
#' }
runXenome <- function(bam, fq1, fq2,
                      outbam, odir,
                      picard_dir = "/scratch/rists/hpcapps/x86_64/picard/1.112", #optsaturn$get("ngs.picard_dir"),
                      xenome_exe = "/scratch/rists/hpcapps/rhel6/xenome/1.0.1-r/xenome", #optsaturn$get("ngs.xenome"),
                      xen_index = "/IACS1/NGS/hg19-mm9-xenome-index/hg19-mm9-xenome-index", #optsaturn$get("ngs.xen.index"),
                      java = "/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java", #optsaturn$get("ngs.java"),
                      javaMem = "-Xmx8g", #optsaturn$get("ngs.xen.javaMem"),
                      xen_cores = 24, ## optsaturn$get("ngs.xen.cores")
                      xenome_fastq_exe = "/scratch/iacs/bin/xenome/xen_to_fastq.sh",
                      ## optsaturn$get("ngs.xen.to.fastq"),
                      python_exe = "/scratch/rists/hpcapps/x86_64/Python/2.7.6-new/bin/python", ## hpcc python
                      filter_reads_exe = "/scratch/iacs/bin/xenome/filter_xeno_mmReads.py",
                      samtools_dir = "/risapps/rhel5/samtools/1.1/bin",
                      ##optsaturn$get("ngs.xen.filter.bam")
                      execute = TRUE,
                      verbose = FALSE){
    ## get missing params
    if(missing(odir)) odir  = dirname(bam)
    if(missing(outbam)) outbam  = file.path(odir, gsub(".bam$", ".xeno.bam", basename(bam)))
    ## evaluating where to start from bam/merged fastqs
    fq_available = sum(file.exists(c(fq1, fq2)));
    bam_fastq_cmd = "" ## initializing command as empty
    ## get a tmp_dir to work in
    oprefix = file_path_sans_ext(basename(bam))
    tmpPath = file.path(odir, oprefix); dir.create(tmpPath, recursive = TRUE, showWarnings = FALSE)
    ## fq1=file.path(tmpPath,gsub(".bam","read1.fq",basename(bam)));file.path(tmpPath,fq2=gsub(".bam","_read2.fq",basename(bam)))

    ## Get FASTQ files to work with
    if(!fq_available){
        fq1 = sprintf("%s/%s_read1.fq", odir, oprefix)
        fq2 = sprintf("%s/%s_read2.fq", odir, oprefix)
    	## fq1 = file.path(tmpPath, "read1.fq");fq2 = file.path(tmpPath, "read2.fq")
    	## Parameters for BAM to fastq
    	## javaMem="-Xmx12g"
    	params = "-XX:-UseGCOverheadLimit -XX:-UseParallelGC"
    	tool = "SamToFastq.jar"
    	extra_opts = ""
    	bam_fastq_cmd = sprintf("time %s %s -Djava.io.tmpdir=%s -jar %s/%s INPUT=%s FASTQ=%s SECOND_END_FASTQ=%s VALIDATION_STRINGENCY=LENIENT %s ",
            java, javaMem, tmpPath, picard_dir, tool, bam, fq1, fq2, extra_opts)##2>>%s ,logFile
    	if(verbose) cat("Converting bam to fastq\n", bam_fastq_cmd, "\n")
    	if(execute) system(bam_fastq_cmd)
    }

    ## RUN xenome
    ## xen.index='/IACS1/NGS/hg19-mm9-xenome-index/hg19-mm9-xenome-index'
    xenome_prefix= sprintf("%s/%s", odir, oprefix)
    xenome_cmd <- sprintf("%s classify --num-threads %s -P %s --pairs --graft-name hg19 --host-name mm9 --output-filename-prefix %s -i %s -i %s -v ",
                          xenome_exe, xen_cores, xen_index, xenome_prefix, fq1, fq2) #2>> %s , logFile
    if(verbose) cat("\nXENOME CMD:\n", xenome_cmd, "\n")
    if(execute) system(xenome_cmd)

    ## convert the xenome fastqs
    xen_out_suffix <- c("mm9_1.fastq","mm9_2.fastq","neither_1.fastq","neither_2.fastq",
                        "ambiguous_1.fastq","ambiguous_2.fastq","both_1.fastq","both_2.fastq")
    fqs <- paste(xenome_prefix, xen_out_suffix, sep = "_")
    ## fqs=list.files(tmpPath,pattern="fastq$",full.name=TRUE)
    ## func <- sapply; if(execute) func <- mclapply
    func = mclapply ## or any other parallelization
    conv_fastq_cmd <- func(fqs, function(x){
                                        #cat(x, "\n")
        cmd = sprintf("%s %s > %s", xenome_fastq_exe, x, gsub("(.*)_([1-2]{1}).fastq","\\1_end\\2.fq",x))
        if(verbose) cat("\nCONVERT CMD:\n", cmd, "\n")
        if(execute) system(cmd)
        return(cmd)
    },mc.cores = xen_cores)

    ## Filtering of the bam file
                                        #env="export PYTHONPATH=$PYTHONPATH:/IACS1/lib/python2.7/site-packages/"
    env = "" # dep: pysam, Bio::SeqIO, [os, time, sys]
    filter_cmd <- sprintf("%s %s %s %s %s %s", python_exe, filter_reads_exe, odir, oprefix, bam, outbam)
    if(verbose) cat("\nFILTER CMD:\n", filter_cmd, "\n")
    if(execute) system(filter_cmd)

    if(execute) indexBAM(outbam, samtoolPath = samtools_dir)
    cmds = list(bam_fastq_cmd = bam_fastq_cmd,
        xenome_cmd = xenome_cmd,
        conv_fastq_cmd = unlist(conv_fastq_cmd),
        filter_cmd = filter_cmd)
    return(list(cmds = cmds, outfile = outbam))
}


format_xen_fastq <- function(fq){


}
