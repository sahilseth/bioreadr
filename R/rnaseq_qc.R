rnaseqcRNASeQC<-function(samplebamlist,outwkdir,
                         java="/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                         bwa = "/scratch/rists/hpcapps/x86_64/bwa-0.7.5a/bwa",
                         rnaSeqcjar="/scratch/rists/hpcapps/x86_64/RNA_SeQC/RNA-SeQC_v1.1.7.jar",
                         reffa="/scratch/rists/hpcapps/reference/human/hg19/fastas/Homo_sapiens_assembly19.fasta",
                         rRnaFa="/scratch/rists/hpcapps/reference/human/broad_hg19/annotations/gencode/human_all_rRNA.fasta",
                         gtf="/scratch/rists/hpcapps/reference/human/broad_hg19/annotations/gencode/gencode.v7.annotation_goodContig.gtf",
                         execute=FALSE,verbose=FALSE) {
  
  cmd <-paste(java, '-Xmx8g -jar',rnaSeqcjar, '-n 1000 -s',samplebamlist,'-bwa',bwa,'-t',gtf,
              '-r',reffa,'-BWArRNA', rRnaFa, '-o',outwkdir,sep=' ')
  cmd <- sprintf("WD=$(pwd);cd %s; %s; cd $WD", outwkdir, cmd) 
  if (execute) {
    system(cmd);
  }
  if (verbose) {
    cat(cmd,"\n");
  }
  return(cmd);
}

rnaseqc2samplelist<-function(runName,outfile,proj,indir,bampostfix="recalibed.bam",
                             bamlistpl="/scratch/iacs/iacs_dep/xsong3/check_new_info/runcodes/rnaseqc_run2samplebamlist.pl",
                             levelii= "/scratch/iacs/gcc/levelii",
                             execute=FALSE,verbose=FALSE) {
  if (missing(outfile)) {
    outfile=paste(runName,'.samplelist',sep='')
    if (!missing(proj)) {
      outfile=paste(runName,proj,'samplelist',sep='.')
    }
  }
  if (missing(proj)) {
    proj='';
  }
  if (missing(indir)) {
    indir=file.path(levelii,runName);
  }
  cmd <- paste('perl', bamlistpl,'-r',runName, '-indir', indir, '-o', outfile, sep=' ');
  if (proj !='') {
    cmd <- paste('perl', bamlistpl,'-r',runName,'-p',proj,'-indir', indir, '-o', outfile, sep=' ');
  } 
  if (execute) {
    system(cmd);
  }
  if (verbose) {
    cat(cmd,"\n");
  }
  return(cmd);
}

rnaseqcRunone<-function(runName,outdir,indir,proj='',bampostfix="recalibed.bam",
                        orootdir="/scratch/iacs/gcc/leveliii/rnaqc",
                        bamlistpl="/scratch/iacs/iacs_dep/xsong3/check_new_info/runcodes/rnaseqc_run2samplebamlist.pl",
                        java="/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                        bwa = "/scratch/rists/hpcapps/x86_64/bwa-0.7.5a/bwa",
                        rnaSeqcjar="/scratch/rists/hpcapps/x86_64/RNA_SeQC/RNA-SeQC_v1.1.7.jar",
                        reffa="/scratch/rists/hpcapps/reference/human/hg19/fastas/Homo_sapiens_assembly19.fasta",
                        rRnaFa="/scratch/rists/hpcapps/reference/human/broad_hg19/annotations/gencode/human_all_rRNA.fasta",
                        gtf="/scratch/rists/hpcapps/reference/human/broad_hg19/annotations/gencode/gencode.v7.annotation_goodContig.gtf",
                        levelii= "/scratch/iacs/gcc/levelii",
                        execute=FALSE,verbose=FALSE) {
  if (missing(outdir)) {
    outdir<-file.path(orootdir,runName)
  }
  if(missing(indir)) {
    indir <-file.path(levelii,runName)
  }
  samplelistfile <-file.path(outdir,paste(runName,'.samplelist',sep=''))
  if (proj =='') {
    outdir<-file.path(outdir,'qcdir')
  } else {
    outdir<-file.path(outdir,proj);
  }
  
  if (!file.exists(outdir)) {
    dir.create(outdir,recursive=TRUE)
  }
  samplelistfile<- file.path(outdir,paste(runName,proj,'samplelist',sep='.'))
  
  cmd1<-rnaseqc2samplelist(runName=runName,outfile=samplelistfile,proj=proj,indir=indir,
                           bampostfix=bampostfix,bamlistpl=bamlistpl,levelii=levelii,
                           execute=execute,verbose=verbose);
  cat("cmd1: ",cmd1,"\n");
  
  cmd2<-rnaseqcRNASeQC(samplebamlist=samplelistfile,outwkdir=outdir,
                       java=java,bwa=bwa,rnaSeqcjar=rnaSeqcjar,reffa=reffa,rRnaFa=rRnaFa,gtf=gtf,
                       execute=execute,verbose=verbose);
  cat("cmd2: ",cmd2,"\n");
  
  oreport<-file.path(outdir,'report.html')
  odindexhtmpl <- file.path(outdir,'index.html')
  cmds<-c(cmd1,cmd2)
  outfiles<-c(oreport,odindexhtmpl);
  rtrns <-c( cmds,outfiles)
  return(rtrns)
}
rnaseqcRunoneMouse<-function(runName,outdir,indir,proj='',bampostfix="recalibed.bam",
                             orootdir="/scratch/iacs/gcc/leveliii/rnaqc",
                             bamlistpl="/scratch/iacs/iacs_dep/xsong3/check_new_info/runcodes/rnaseqc_run2samplebamlist.pl",
                             java="/scratch/rists/hpcapps/x86_64/jdk/jdk1.7.0/bin/java",
                             bwa = "/scratch/rists/hpcapps/x86_64/bwa-0.7.5a/bwa",
                             rnaSeqcjar="/scratch/rists/hpcapps/x86_64/RNA_SeQC/RNA-SeQC_v1.1.7.jar",
                             reffa="/scratch/rists/hpcapps/reference/mouse/mm10/fastas/mm10.fa",
                             rRnaFa="/scratch/rists/hpcapps/reference/mouse/mm10/annotations/rnaseqc/mouse_all_rRNA.fasta",
                             gtf="/scratch/rists/hpcapps/reference/mouse/mm10/annotations/ensembl/mm10.GRCm38.75.tr.gtf",
                             levelii= "/scratch/iacs/gcc/levelii",
                             execute=FALSE,verbose=FALSE) {
  if (missing(outdir)) {
    outdir<-file.path(orootdir,runName)
  }
  if(missing(indir)) {
    indir <-file.path(levelii,runName)
  }
  samplelistfile <-file.path(outdir,paste(runName,'.samplelist',sep=''))
  if (proj=='') {
    outdir<-file.path(outdir,'qcdir')
  } else {
    samplelistfile<- file.path(outdir,paste(runName,proj,'samplelist',sep='.'))
    outdir<-file.path(outdir,proj);
  }
  if (!file.exists(outdir)) {
    dir.create(outdir,recursive=TRUE)
  }
  
  cmd1<-rnaseqc2samplelist(runName=runName,outfile=samplelistfile,proj=proj,indir=indir,
                           bampostfix=bampostfix,bamlistpl=bamlistpl,levelii=levelii,
                           execute=execute,verbose=verbose);
  cat("cmd1: ",cmd1,"\n");
  
  cmd2<-rnaseqcRNASeQC(samplebamlist=samplelistfile,outwkdir=outdir,
                       java=java,bwa=bwa,rnaSeqcjar=rnaSeqcjar,reffa=reffa,rRnaFa=rRnaFa,gtf=gtf,
                       execute=execute,verbose=verbose);
  cat("cmd2: ",cmd2,"\n");
  
  oreport<-file.path(outdir,'report.html')
  odindexhtmpl <- file.path(outdir,'index.html')
  cmds<-c(cmd1,cmd2)
  outfiles<-c(oreport,odindexhtmpl);
  rtrns <-c( cmds,outfiles)
  return(rtrns)
}

