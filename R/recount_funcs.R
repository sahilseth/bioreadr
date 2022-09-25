

if(FALSE){
    x = tab
}


## get big intersect of all
## idea is to re-evaluate significnat mutations in the whole cohort
## This is implemented on a variable 'group'
## key
mk_mut_superset <- function(x, key, group, sample){

}


mk_intervals <- function(x, ## a table with 1st col for chr followed by start and stop position; then ref allele and ALT allele
                         reflib,
                         with_dict = FALSE,
                         outfile = "interval_list"){
    x <- as.matrix(x)
    if(with_dict){
        dict <- gsub("fasta", "dict", reflib)
        if(!file.exists(dict))
            stop("Dictionary file of reflib does not exist, please check")
        bed <- cbind(str_trim(x[, 1]),
                     str_trim(x[, 2]),
                     str_trim(x[, 3]),
                     strand = "+", target = 1:nrow(x))
        if(ncol(x)>3)
            out = cbind(out, x[,4:ncol(x)])
        tmp <- file.copy(dict, outfile, overwrite = TRUE)
        write.table(out, file = outfile, append = TRUE, col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
    }else{
        out = paste0(str_trim(x[,1]), ":", str_trim(x[,2]), "-", str_trim(x[,3]))
        cat(out, file = outfile, sep = "\n")
        write.table(x, sprintf("%s_map.txt", tools::file_path_sans_ext(outfile)), quote = FALSE, row.names = FALSE)
    }
    invisible(outfile)
}


## This used gatk depth of coverage to get coverage and base count using bam and bed file
# "module load jdk/1.7.0_79"
count_bases <- function(
    bam,
    intervals,
    outfile,
    cores = 12,
    java_exe = "java",
    gatk_jar = "/scratch/rists/hpcapps/x86_64/gatk/3.1-1/GenomeAnalysisTK.jar",
    reflib = "/scratch/iacs/gcc/Sarcomatoid/referenceLibrary/hg19/hg19.fasta"){
    params <- sprintf("--printBaseCounts --omitLocusTable --omitIntervalStatistics --omitPerSampleStats --outputFormat table -nt %s",
                      cores)
    cov_tmp <- flowr:::get_unique_id(tools::file_path_sans_ext(outfile))
    cmd <- sprintf("%s -jar %s -T DepthOfCoverage -R %s -L %s -o %s -I %s %s",
                   java_exe, gatk_jar, reflib, intervals, cov_tmp, bam, params)
    message("Running command:\n", cmd)
    tmp <- system(cmd)
                                        #parse_cov(cov_tmp, outfile)
    invisible(outfile)
}


parse_cov <- function(x, outfile,
                      n = 1000 ## read in batches of 1000 lines
                      ){ ## dop_table
    require(dplyr)
                                        #n = 1000; #lines to read
                                        #x <- "/rsrch2/iacs/iacs_dep/sseth/tmp/tmp.table"
                                        #outfile = "/rsrch2/iacs/iacs_dep/sseth/tmp/cov.txt"
    require(stringr)
    con  <- file(x, open = "r")
    line <- readLines(con, n = 1)
    cols = c("chrom", "pos", "A", "C", "G", "T", "N")
    outcols = c("chrom", "pos", "base", "value")
    cat(c(outcols, "\n"), sep = "\t", file = outfile)
    db_path = sprintf("%s.sqlite3", tools::file_path_sans_ext(outfile))
    if(db){
        require(RSQLite)
        my_db <- src_sqlite(db_path, create = T)
        mydb = dbConnect(SQLite(), dbname = db_path)
                                        #tmp <- try(db_create_table(mydb, 'counts', types = c(chrom = "string", pos = "numeric", base = "string", value = "numeric")))
                                        #counts <- data.frame(chrom = NA, pos = 0, base = NA, count = 0)
        db_drop_table(con = mydb, table = "counts")
                                        #copy_to(dest = my_db, counts)
        qry <- "create table counts (chrom string not null, pos INTEGER not null, base string not null, count integer not null,
				PRIMARY KEY ( chrom, pos, base));"
        tmp <- dbSendQuery(conn = mydb, qry)
    }
                                        #
    while (length(line <- readLines(con, n = n)) > 0) {
        cat(".")
        mat = do.call(rbind, str_split(line, "\t"))
        count = matrix(na.omit(as.numeric(unlist(str_split(mat[,5], "[ATGCN: ]")))), ncol = 5, byrow = TRUE)
        count = data.table::data.table(str_match(mat[, 1, drop = FALSE], "(.*):(.*)")[,-1, drop = FALSE], count)
        colnames(count) = cols
        counts <- data.table::melt(as.data.frame(count), measure.vars = c("A", "T", "G", "C", "N"),
                                   variable.name = "base", value.name = "count")
        tmp <- dbWriteTable(conn = mydb, name = "counts", value = counts, append = TRUE)
        write.table(counts, file = outfile, append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

    }
    close(con)
    dbDisconnect(mydb)
}


if(FALSE){
    ## accesing the DB
    tab <- tbl(my_db, "counts")


}
