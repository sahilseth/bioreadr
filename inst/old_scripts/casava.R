
####### Move the Data path
initialize_casava_paths <- function(runid,intPath=intPath){
    ## if intensities is not a symbolic link
    ass <- ngs.createSampleSheet(runid,intPath) ## outside nodes
    if(system(sprintf("test -L %s; echo $?",intPath),intern=TRUE)!=0){ #0 means true
        ## move intensity folder and place RunInfo.xml
        dir.create(file.path(dirname(intPath),"Data"),showWarnings=FALSE)
        file.rename(file.path(intPath),file.path(dirname(intPath),"Data","Intensities"))
        system(sprintf("ln -s %s %s",file.path(dirname(intPath),"Data","Intensities"),file.path(intPath)))
    }
    return(ass)
}

create_casava_sample_sheet <- function(runid,intPath=getOption("ngs.intPath")){
    ass <- getOption("ngs.NGSAssignment")
    if(is.null(ass)){options(ngs.NGSAssignment=getAssignment(runid));ass <- getOption("ngs.NGSAssignment")}
    mask <- getOption("ngs.mask");indLen <- gsub(".*I([0-9]*).*","\\1",mask)#index,length
    ass[,"Flowcell"] <- substr(ass[,"Flowcell"],2,10) #remove first char
    cols <- c("FCID","Lane","SampleID","SampleRef","Index","Description","Control","Recipe","Operator","SampleProject")
    ## ass=ass[ass[,c("Index")]!="NA",]
    ##if we are getting a shorter length to accomodate...
    ass[,c("Index")]=ifelse(ass[,c("Index")]=="NA","",ass[,c("Index")])
    ## reorder columns
    ass2 <- cbind(ass[,c("Flowcell","Lane","Sample.Name")],"NA",ass[,c("Index")],"NA","N","NA","NA",
                  apply(ass[, c("Project","Sub.Project")],1, function(x) paste(x,collapse = "-")))
    colnames(ass2) <- cols;
    ## trim assignment
    ass3 <- ass2;ass3[,c("Index")]=ifelse(ass[,c("Index")] %in% c("NA",""),"",substr(ass[,c("Index")],1,indLen))
    ## fill assignment with 16 dummy tags
    ass4 <- ngs.fill.samplesheet(ass2,newlength=indLen)
    write.csv(ass4,file.path(intPath,"BaseCalls","SampleSheet_fill.csv"),quote=FALSE,row.names=FALSE)
    write.csv(ass2,file.path(intPath,"BaseCalls","SampleSheet_all.csv"),quote=FALSE,row.names=FALSE)
    write.csv(ass3,file.path(intPath,"BaseCalls","SampleSheet.csv"),quote=FALSE,row.names=FALSE)
    return(ass2)
}


create_casava_plain_sheet <- function(outPath,lanes,flowcellid){
    cols <- c("FCID","Lane","SampleID","SampleRef","Index","Description","Control","Recipe","Operator","SampleProject")
    ass <- data.frame(flowcellid,lanes,paste("Sample",lanes,sep=""),"","","","","","","SampleProject")
    write.csv(ass,file.path(outPath,"SampleSheetPlain.csv"),quote=FALSE,row.names=FALSE)
}


## sets up seqLen also
get_casava_mask <- function(path=getOption("ngs.bclPath")){
    require(XML)
    ##ass <- getOption("ngs.NGSAssignment")
    xml <- file.path(path,"RunInfo.xml")
    reads <- xmlRoot(xmlTreeParse(xml))[["Run"]][["Reads"]]
    out <- sapply(1:length(reads),function(i){
        len <- xmlAttrs(reads[[i]])[["NumCycles"]]
        type <- xmlAttrs(reads[[i]])[["IsIndexedRead"]]
        t <- ifelse(type=="Y","I","Y")
        options(ngs.seqLen=len)
        sprintf("%s%s",t,len)
    })
    mask <- paste(out,collapse=",")
    options(ngs.mask=mask)
    return(mask)
}

