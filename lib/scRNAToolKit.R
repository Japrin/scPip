loginfo <- function(msg) {
  timestamp <- sprintf("%s", Sys.time())
  msg <- paste0("[",timestamp, "] ", msg,"\n")
  cat(msg)
}

entrezToXXX<-function(x,type="SYMBOL",species="human")
{
  ret=c()
  x=as.character(x)
  if(species=="human")
  {
      suppressPackageStartupMessages(require("org.Hs.eg.db"))
  }else if(species=="mouse")
  {
      suppressPackageStartupMessages(require("org.Mm.eg.db"))
  }else
  {
      stop("species must be human or mouse\n")
  }
  if(species=="human")
  {
      if(type=="SYMBOL")
      {
        ret=unlist(as.list(org.Hs.egSYMBOL))[x]
      }else if(type=="ENSG")
      {
        ret=unlist(as.list(org.Hs.egENSEMBL))[x]
      }else if(type=="GENENAME")
      {
        ret=unlist(as.list(org.Hs.egGENENAME))[x]
      }
  }else if(species=="mouse")
  {
      if(type=="SYMBOL")
      {
        ret=unlist(as.list(org.Mm.egSYMBOL))[x]
      }else if(type=="ENSG")
      {
        ret=unlist(as.list(org.Mm.egENSEMBL))[x]
      }
  }
  return(as.character(ret))
}

XXXToEntrez<-function(x,type="SYMBOL",species="human")
{
  ret=c()
  x=as.character(x)
  if(species=="human")
  {
      suppressPackageStartupMessages(require("org.Hs.eg.db"))
  }else if(species=="mouse")
  {
      suppressPackageStartupMessages(require("org.Mm.eg.db"))
  }else
  {
      stop("species must be human or mouse\n")
  }
  if(species=="human")
  {
      if(type=="SYMBOL")
      {
        ret=unlist(as.list(org.Hs.egSYMBOL2EG))[x]  
      }else if(type=="ENSG")
      {
        ret=unlist(as.list(org.Hs.egENSEMBL2EG))[x]
      }
  }else if(species=="mouse")
  {
      if(type=="SYMBOL")
      {
        ret=unlist(as.list(org.Mm.egSYMBOL2EG))[x]  
      }else if(type=="ENSG")
      {
        ret=unlist(as.list(org.Mm.egENSEMBL2EG))[x]
      }
  }
  return(as.character(ret))
}

run.SIMLR <- function(inputData,n.clusters=4,myseed=NULL,my.title="",out.prefix,legend.txt,col.points,col.legend,
                      pch=16,width.pdf=10,height.pdf=8,margin.r=10,legend.inset=-0.25,legend.txt2=NULL,...)
{
    require(Matrix)
    require(parallel)
    require(igraph)
    # load the palettes for the plots
    require(grDevices)
    # load the SIMLR R package
    SIMLR.DIR <- "/Share/BP/zhenglt/01.bin/scRNASeq/SIMLR"
    source(sprintf("%s/R/SIMLR.R",SIMLR.DIR))
    source(sprintf("%s/R/compute.multiple.kernel.R",SIMLR.DIR))
    source(sprintf("%s/R/network.diffusion.R",SIMLR.DIR))
    source(sprintf("%s/R/utils.simlr.R",SIMLR.DIR))
    source(sprintf("%s/R/tsne.R",SIMLR.DIR))
    # load the C file
    dyn.load(sprintf("%s/src/projsplx_R.so",SIMLR.DIR))
    if(is.null(myseed)){ myseed <- as.integer(Sys.time()) }
    set.seed(myseed)
    loginfo(sprintf("set.seed(%d) for SIMLR\n",myseed))

    res.SIMLR <- SIMLR(X=inputData,c=n.clusters)

    pdf(file=sprintf("%s.SIMLR.default.pdf",out.prefix),width=width.pdf,height=height.pdf)
    par(mar=c(5,5,4,margin.r),cex.lab=1.5,cex.main=1.5)
    ### color by predefined sample type
    plot(res.SIMLR$ydata, t='n', main=my.title,xlab="SIMLR component 1",ylab="SIMLR component 2")
    points(res.SIMLR$ydata,col=col.points,pch=pch,...)
    legend("right",legend=legend.txt,fill = NULL,inset = legend.inset,xpd = NA,cex=1.5,pch=16,border =NA,col = col.legend)
    if(!is.null(legend.txt2)) { 
           legend("topright",legend=legend.txt2,fill = NULL,inset = c(legend.inset,0),xpd = NA,cex=1.5,pch=unique(pch),border =NA,col = col.legend[1]) 
    }

    ## color by clustering result
    clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(res.SIMLR$y$cluster))),
                              names=as.character(unique(sort(res.SIMLR$y$cluster))))
    ##clusterColor["0"] <- "gray"
    print(clusterColor)
    plot(res.SIMLR$ydata, col=clusterColor[as.character(res.SIMLR$y$cluster)], 
         pch=16,main=my.title,xlab="SIMLR component 1",ylab="SIMLR component 2")
    legend("right",legend=sprintf("cluster%s",names(clusterColor)),
           fill = NULL,inset = legend.inset,xpd = NA,cex=1.5,pch=16,border =NA, 
           col = clusterColor)
    dev.off()

    return(res.SIMLR)
}


run.limma.from.matrixFile <- function(in.dat,myDesign,infile=NULL,designFile=NULL,T.fdr=0.1,T.logFC=1,withPatient=T,verbose=F,gid.mapping=NULL,
                                      do.voom=F)
{
	suppressPackageStartupMessages(require("limma"))
	suppressPackageStartupMessages(require("dplyr"))
	suppressPackageStartupMessages(require("BiocParallel"))
    if(is.null(in.dat) && !is.null(infile) && file.exists(infile)){
	    inTable <- read.table(infile,row.names=1,check.names = F,header = T,stringsAsFactors = F)
        gid.mapping <- inTable[,1]
        names(gid.mapping) <- rownames(inTable)
        f.na <- which(is.na(gid.mapping))
        gid.mapping[f.na] <- names(gid.mapping)[f.na]
	    inMatrix <- as.matrix(inTable[,-1])
    }else{
        inMatrix <- in.dat
    }
	###myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
    if(is.null(myDesign) && !is.null(designFile) && file.exists(designFile)){
        myDesign<-read.table(designFile,header=T,check.names=F,stringsAsFactors = F)
        rownames(myDesign) <- myDesign$sample
    }
    f.sample <- intersect(colnames(inMatrix),rownames(myDesign))
    #print(str(inMatrix))
    #print(str(myDesign))
    #print(str(f.sample))
    inMatrix <- inMatrix[,f.sample]
    myDesign <- myDesign[f.sample,]
    myDesign$sampleType <- factor(myDesign$sampleType)
    if(verbose){
        #print(str(myDesign))
        #print(head(inMatrix))
    }
    print("XXXX")
    if(withPatient){
	    design <- model.matrix(~patient+sampleType,data=myDesign)
    }else{
	    design <- model.matrix(~sampleType,data=myDesign)
    }
	colnames(design)[length(colnames(design))]<-"II"
    print("XXXX")

	register(MulticoreParam(4))
    if(do.voom){
        pdf(sprintf("./voom.pdf"),width = 7,height = 7)
        v <- voom(inMatrix, design, plot=TRUE)
        dev.off()
	    fit <- lmFit(v, design)
    }else{
	    fit <- lmFit(inMatrix, design)
    }
	fit <- eBayes(fit)
	all.table  <- topTable(fit, coef = "II", n = Inf, sort = "p", p = 1)
	#out.resSig <- topTable(fit, coef = "II", n = Inf, sort = "p", p = fdr)
	
	all.table <- cbind(data.frame(geneID=rownames(all.table),
                                  geneSymbol=if(is.null(gid.mapping)) rownames(all.table) else gid.mapping[rownames(all.table)],
                                  stringsAsFactors = F),
                       all.table)

    .getMean <- function(inDat,str.note=NULL){
        .Grp.mean <- t(apply(inDat,1,function(x){
                                 .mexp <- aggregate(x ~ myDesign$sampleType, FUN = mean)
                                 structure(.mexp[,2],names=sprintf("mean.Grp%s",.mexp[,1])) }))
        if(!is.null(str.note)){
            colnames(.Grp.mean) <- sprintf("%s.%s",colnames(.Grp.mean),str.note)
        }
        #print(head(.Grp.mean))
        .Grp.mean.df <- data.frame(geneID=rownames(.Grp.mean),stringsAsFactors = F)
        .Grp.mean.df <- cbind(.Grp.mean.df,.Grp.mean)
    }
    .Grp.mean.df <- .getMean(inMatrix)
    all.table <- left_join(all.table,.Grp.mean.df)

    if(withPatient){
        inMatrix.rmBatchEff <- removeBatchEffect(inMatrix, myDesign$patient)
        .Grp.mean.df <- .getMean(inMatrix.rmBatchEff,str.note="rmBatchEff")
        all.table <- left_join(all.table,.Grp.mean.df)
    }


    print(head(all.table))
    out.resSig <- subset(all.table,adj.P.Val<T.fdr & abs(logFC)>T.logFC)
	#out.resSig <- cbind(data.frame(Gene=rownames(out.resSig),stringsAsFactors=F),
	#      data.frame(Entrez=XXXToEntrez(rownames(out.resSig)),stringsAsFactors=F),
	#      out.resSig)

	#write.table(all.table,paste0(output.prefix,".all"),row.names = F,quote = F,sep = "\t")
	#write.table(out.resSig,paste0(output.prefix,".sig"),row.names = F,quote = F,sep = "\t")
	ret.dat <- list(all=all.table,sig=out.resSig)
    if(verbose){
        ret.dat[["in.dat"]] <- inMatrix
        ret.dat[["myDesign"]] <- myDesign
    }
    return(ret.dat)
}

readGMT<-function(file) 
{
    f <- readLines(file)
    lst = sapply(f, function(x) unlist(strsplit(x, "\t", fixed = TRUE)))
    names(lst) = sapply(lst, function(x) x[1])
    gSet = lapply(lst, function(x) x[-(1:2)])
    ###gSet = toupper(gSet)
    gLink = unlist(lapply(lst, function(x) x[2]))
    list(gSet=gSet,gLink=gLink)
}

gageAnalysis<-function(exp.fc,outDir,fdr=0.1)
{
    suppressPackageStartupMessages(require("gage"))
    suppressPackageStartupMessages(require("pathview"))
    suppressPackageStartupMessages(require("ReportingTools"))
    suppressPackageStartupMessages(require("hwriter"))
    data(kegg.gs)
    #c2.cp.gs=readGMT("/lustre1/zeminz_pkuhpc/00.database/MSigDB/msigdb_v4.0_files_to_download_locally/msigdb_v4.0_GMTs/c2.cp.v4.0.entrez.gmt")
    c2.cp.gs=readGMT("/Share/BP/zhenglt/00.database/MSigDB/msigdb_v4.0_files_to_download_locally/msigdb_v4.0_GMTs/c2.cp.v4.0.entrez.gmt")
    #### gage using kegg.gs
    gagePathway<-function(exp.fc,outDir,gs,gs.name,fdr,gl=NULL)
    {
	# 1 direction
	fc.p <- gage(exp.fc, gsets = gs, ref = NULL, samp = NULL)
	fc.p$greater=as.data.frame(fc.p$greater)
	fc.p$less=as.data.frame(fc.p$less)
	sel.h <- fc.p$greater[, "q.val"] < fdr & !is.na(fc.p$greater[, "q.val"])
	path.ids.h <- rownames(fc.p$greater)[sel.h]
	sel.l <- fc.p$less[, "q.val"] < fdr & !is.na(fc.p$less[,"q.val"])
	path.ids.l <- rownames(fc.p$less)[sel.l]
	# 2 direction
	fc.p.2d <- gage(exp.fc, gsets = gs, ref = NULL, samp = NULL, same.dir = F)
	fc.p.2d$greater=as.data.frame(fc.p.2d$greater)
	fc.p.2d$less=as.data.frame(fc.p.2d$less)
	sel.2d <- fc.p.2d$greater[, "q.val"] < fdr & !is.na(fc.p.2d$greater[, "q.val"])
	path.ids.2d <- rownames(fc.p.2d$greater)[sel.2d]
	# all path ids
	path.ids <- c(path.ids.h, path.ids.l, path.ids.2d)
	if(sum(regexpr("hsa[0-9]{5}",names(gs),perl=T)==rep(1,length(names(gs))))==length(names(gs)))
	{
		path.ids <- substr(path.ids, 1, 8)
	}
	## output table
	out1<-fc.p$greater[sel.h,]
	out1<-data.frame(pathway=rownames(out1),out1)
	out2<-fc.p$less[sel.l,]
	out2<-data.frame(pathway=rownames(out2),out2)
	out1d<-rbind(out1,out2)
	out2d<-fc.p.2d$greater[sel.2d,]
	out2d<-data.frame(pathway=rownames(out2d),out2d)
	write.table(out1d, file=paste(outDir,"/gage.",gs.name,".1d.txt",sep=""),sep="\t",quote=F,row.names=F)
	write.table(out2d, file=paste(outDir,"/gage.",gs.name,".2d.txt",sep=""),sep="\t",quote=F,row.names=F)
	####
	if(gs.name=="KEGG")
	{
	    # pathview, only for "KEGG"
	    picDir <- paste(outDir,"/reports/figure_gage",sep="")
	    dir.create(picDir,recursive=T, showWarnings = FALSE)
	    out.suffix="gage"
	    oriDir <- getwd()
	    setwd(picDir)
	    ## OUT/reports/figure_gage/hsa04740.gage.DESeq2.png
	    #pv.out.list <- sapply(path.ids, function(pid) pathview(kegg.dir="/lustre1/zeminz_pkuhpc/00.database/kegg/pathview",gene.data =  exp.fc, pathway.id = pid, species = "hsa", out.suffix=out.suffix))
	    pv.out.list <- sapply(path.ids, function(pid) pathview(kegg.dir="/Share/BP/zhenglt/00.database/kegg/pathview",gene.data =  exp.fc, pathway.id = pid, species = "hsa", out.suffix=out.suffix))
	    setwd(oriDir)
	    if(nrow(out1d) > 0)
	    {
		imagename <- paste0("figure_gage/",substr(out1d$pathway, 1, 8),".",out.suffix,".png")
		out1d$Image <- hwriteImage(imagename, link = imagename, table = FALSE, height=50, width=50)
	    }
	    if(nrow(out2d) > 0)
	    {
		imagename <- paste0("figure_gage/",substr(out2d$pathway, 1, 8),".",out.suffix,".png")
		out2d$Image <- hwriteImage(imagename, link = imagename, table = FALSE, height=50, width=50)
	    }
	}
	#### Report writting
	addGSetLink <- function(object, ...)
	{
		if(!is.null(gl))
		{
			object$pathway <- hwrite(as.character(object$pathway), link=gl[as.character(object$pathway)], table = FALSE)
		}
		return(object)
	}
	gageReport <- HTMLReport(shortName = paste("gage_analysis_rnaseq_",gs.name,"",sep=""), title = paste("GAGE analysis of RnaSeqData (", gs.name, ")",sep=""), reportDirectory = "reports", basePath=outDir)
	publish(hwrite("One direction (up or down)", heading=2), gageReport)
	if(nrow(out1d) > 0)
	{
		publish(out1d, gageReport, reportDir="./reports", .modifyDF=list(addGSetLink))
	}
	publish(hwrite("Two direction (up and down)", heading=2), gageReport)
	if(nrow(out2d) > 0)
	{
		publish(out2d, gageReport, reportDir="./reports", .modifyDF=list(addGSetLink))
	}
	finish(gageReport)
	list(report=gageReport)
    }
    gage.kegg=gagePathway(exp.fc,outDir,kegg.gs,"KEGG",0.1)
    gage.Canonical=gagePathway(exp.fc,outDir,c2.cp.gs$gSet,"Canonical",fdr,c2.cp.gs$gLink)
    ## return report obj for "index" web
    list(gage.kegg$report, gage.Canonical$report)
}

outputDEGene<-function(res,resSigStrict,dds,vstMat,fc.df,outDir)
{
    ## txt report
    resSig <- subset(res, padj < 0.1)
    resSigStrict <- subset(resSig, padj < opt.qval & abs(log2FoldChange) > 1)

    out_res <- as.data.frame(res)
    out_res$EntrezGeneID <- rownames(out_res)
    out_res$GeneSymbol <- entrezToXXX(out_res$EntrezGeneID)
    out_res$ENSG <- entrezToXXX(out_res$EntrezGeneID,"ENSG")
    
    if(nrow(fc.df)>0)
    {
    	out_res <- cbind(out_res, as.data.frame(vstMat)[rownames(out_res),], fc.df[rownames(out_res),c(-1,-2,-3)])
    }else
    {
    	out_res <- cbind(out_res, as.data.frame(vstMat)[rownames(out_res),])
    }
    
    out_resSig <- subset(out_res, padj < 0.1 )
    out_resSigStrict <- subset(out_resSig, padj < opt.qval & abs(log2FoldChange) > 1 )
    write.table(out_res, file=paste(outDir,"/DESeq2.txt",sep=""),sep="\t",quote=F,row.names=F)
    write.table(out_resSig, file=paste(outDir,"/DESeq2.sig.txt",sep=""),sep="\t",quote=F,row.names=F)
    write.table(out_resSigStrict, file=paste(outDir,"/DESeq2.sig.strict.txt",sep=""),sep="\t",quote=F,row.names=F)
    ## HTML report
    des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2', title = 'RNA-seq analysis of differential expression using DESeq2', reportDirectory="reports", basePath=outDir)
    if(nrow(resSigStrict)>0)
    {
	    publish(resSigStrict,des2Report, pvalueCutoff=opt.qval,lfc=1, n=Inf, annotation.db="org.Hs.eg.db", DataSet=dds, factor = colData(dds)$sampleType, reportDir="./reports")
    }
    finish(des2Report)
    des2Report
}


## resSig: result return by run.limma.from.matrixFile
outputDEGeneFromLimma <-function(res,output.prefix)
{
    ## txt report
    all.table <- res$all
    resSig <- res$sig
    #resSigStrict <- subset(resSig, adj.P.Val < 0.05 & abs(logFC) > 1 )
    write.table(all.table,paste0(output.prefix,".all.txt"),row.names = F,quote = F,sep = "\t")
    write.table(resSig,paste0(output.prefix,".sig.txt"),row.names = F,quote = F,sep = "\t")
    #write.table(resSigStrict,paste0(output.prefix,".sig.strict.txt"),row.names = F,quote = F,sep = "\t")

    suppressPackageStartupMessages(require("ReportingTools"))
    suppressPackageStartupMessages(require("hwriter"))
    ## HTML report
    deReport <- HTMLReport(shortName = 'RNA_analysis_with_Limma', title = 'RNA Analysis of differential expression using Limma', reportDirectory="reports", basePath=dirname(output.prefix))
    addEntrezLink <- function(object, ...)
    {
    	object$Entrez <- hwrite(as.character(object$Entrez), link = paste0("http://www.ncbi.nlm.nih.gov/gene/", as.character(object$Entrez)), table = FALSE)
	return(object)
    }
    if(nrow(resSig)>0)
    {
    	publish(resSig,deReport, annotation.db="org.Hs.eg.db",.modifyDF=list(addEntrezLink) )
    }
    finish(deReport)
    deReport
}

hyperGEA<-function(x,selectedIDs,universeIDs,reportName,outDir, ...)
{
    suppressPackageStartupMessages(require("ReportingTools"))
    suppressPackageStartupMessages(require("hwriter"))
    suppressPackageStartupMessages(require("GOstats"))
    suppressPackageStartupMessages(require("Category"))
    suppressPackageStartupMessages(require("KEGG.db"))
    
    addGOIDLink <- function(object, ...)
    {
    	object$GOID <- hwrite(as.character(object$GOID), link = paste0("http://amigo.geneontology.org/amigo/term/", as.character(object$GOID)), table = FALSE)
	return(object)
    }
    addKEGGIDLink <- function(object, ...)
    {
    	#http://www.genome.jp/kegg/pathway/hsa/hsa05219.html
    	object$KEGGID <- hwrite(as.character(object$KEGGID), link = paste0("http://www.genome.jp/kegg/pathway/hsa/hsa", as.character(object$KEGGID),".html"), table = FALSE)
	return(object)
    }

    aReport <- HTMLReport(shortName = paste("Hypergeometric Tests for Gene Set (",reportName,")",sep=""), title = paste("Hypergeometric Tests for Gene Set (",reportName,")",sep=""), reportDirectory="reports", basePath=outDir)

    if(length(selectedIDs)>0)
    {
	aParams <- new(x, geneIds = selectedIDs, universeGeneIds = universeIDs, annotation ="org.Hs.eg", pvalueCutoff = 0.01, testDirection = "over", ... )
	aResults <- hyperGTest(aParams)
	aResultsSummary <- summary(aResults)

	if(nrow(aResultsSummary)>0)
	{
		if(x=="GOHyperGParams")
		{
		    #http://amigo.geneontology.org/amigo/term/GO:0000070
		    names(aResultsSummary)[1]<-"GOID"
		    publish(aResultsSummary, aReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg", pvalueCutoff= 0.01, .modifyDF=list(addGOIDLink))
		}else if(x=="KEGGHyperGParams")
		{
		    publish(aResultsSummary, aReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg", pvalueCutoff= 0.01, .modifyDF=list(addKEGGIDLink))
		}
	}
    }
    finish(aReport)
    aReport
   
    # not work still
    #pfamParams <- new("PFAMHyperGParams", geneIds= selectedIDs, universeGeneIds=universeIDs, annotation="org.Hs.eg", pvalueCutoff= 0.01, testDirection="over")
    #PFAMResults <- hyperGTest(pfamParams)
    #PFAMReport <- HTMLReport(shortName = 'pfam_analysis_rnaseq', title = "PFAM analysis of RnaSeqData", reportDirectory = paste(outDir,"/reports",sep=""))
    #publish(PFAMResults, PFAMReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg",categorySize=5)
    #finish(PFAMReport)

}

pairs2 <- function (x, labels, panel = points, ..., lower.panel = panel, 
    upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
    label.pos = 0.5 + has.diag/3, line.main = 3, cex.labels = NULL, 
    font.labels = 1, row1attop = TRUE, gap = 1, log = "",xlim=NULL, ylim=NULL) 
{
    if (doText <- missing(text.panel) || is.function(text.panel)) 
        textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
            y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
        oma, ...) {
        xpd <- NA
        if (side%%2L == 1L && xl[j]) 
            xpd <- FALSE
        if (side%%2L == 0L && yl[i]) 
            xpd <- FALSE
        if (side%%2L == 1L) 
            Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for (i in seq_along(names(x))) {
            if (is.factor(x[[i]]) || is.logical(x[[i]])) 
                x[[i]] <- as.numeric(x[[i]])
            if (!is.numeric(unclass(x[[i]]))) 
                stop("non-numeric argument to 'pairs'")
        }
    }
    else if (!is.numeric(x)) 
        stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
        lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
        upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
        diag.panel <- match.fun(diag.panel)
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
        stop("only one column in the argument to 'pairs'")
    if (doText) {
        if (missing(labels)) {
            labels <- colnames(x)
            if (is.null(labels)) 
                labels <- paste("var", 1L:nc)
        }
        else if (is.null(labels)) 
            doText <- FALSE
    }
    oma <- if ("oma" %in% nmdots) 
        dots$oma
    main <- if ("main" %in% nmdots) 
        dots$main
    if (is.null(oma)) 
        oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    xl <- yl <- logical(nc)
    if (is.numeric(log)) 
        xl[log] <- yl[log] <- TRUE
    else {
        xl[] <- grepl("x", log)
        yl[] <- grepl("y", log)
    }
    for (i in if (row1attop) 
        1L:nc
    else nc:1L) for (j in 1L:nc) {
        l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", 
            ""))
        # modify here
         if (is.null(xlim) & is.null(ylim))
		localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
			type = "n", ..., log = l)
	if (is.null(xlim) & !is.null(ylim))
		localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
			type = "n", ..., log = l, ylim=ylim)
	if (!is.null(xlim) & is.null(ylim))
		localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
			type = "n", ..., log = l, xlim = xlim)
	if (!is.null(xlim) & !is.null(ylim))
		localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
			type = "n", ..., log = l, xlim = xlim, ylim=ylim)
        #localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
        #    type = "n", ..., log = l)
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
            box()
            if (i == 1 && (!(j%%2L) || !has.upper || !has.lower)) 
                localAxis(1L + 2L * row1attop, x[, j], x[, i],
                  ...)
            # modify here
	    #if (i == nc && (j%%2L || !has.upper || !has.lower)) 
            #    localAxis(3L - 2L * row1attop, x[, j], x[, i], 
            #      ...)
            #if (j == 1 && (!(i%%2L) || !has.upper || !has.lower)) 
            #    localAxis(2L, x[, j], x[, i], ...)
            if (j == nc && (i%%2L || !has.upper || !has.lower)) 
                localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if (i == j) {
                if (has.diag) 
                  localDiagPanel(as.vector(x[, i]), ...)
                if (doText) {
                  par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                  }
                  xlp <- if (xl[i]) 
                    10^0.5
                  else 0.5
                  ylp <- if (yl[j]) 
                    10^label.pos
                  else label.pos
                  text.panel(xlp, ylp, labels[i], cex = cex.labels, 
                    font = font.labels)
                }
            }
            else if (i < j) 
                localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                  i]), ...)
            else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                i]), ...)
            if (any(par("mfg") != mfg)) 
                stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots) 
            dots$font.main
        else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots) 
            dots$cex.main
        else par("cex.main")
        mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, 
            font = font.main)
    }
    invisible(NULL)
}

panel.cor <- function (x, y, digits = 2, meth = "pearson", cex.cor = 1,...)
{
	#print(meth)
	usr <- par("usr")
	on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
        r <- round(cor(x, y, method = meth), digits)
	cr <- round(cor.test(x, y, method = meth, alternative = "t")$p.value,digits)
	text(0.5, 0.65, paste("r =", r), cex = cex.cor)
	if (!cr == 0) {
		text(0.5, 0.35, paste("p =", cr), cex = cex.cor)
	}
	else text(0.5, 0.35, paste("p < 0.01"), cex = cex.cor)
}


#' group specific genes
#' 
#' @param dat.g
#' @param grps
#' @param out.prefix (string) output prefix
#' @param FDR.THRESHOLD
#' @param FC.THRESHOLD
runMultiGroupSpecificGeneTest <- function(dat.g,grps,out.prefix,mod=NULL,FDR.THRESHOLD=0.05,F.FDR.THRESHOLD=0.05,FC.THRESHOLD=1,verbose=F,n.cores=NULL,gid.mapping=NULL)
{
    suppressPackageStartupMessages(require("plyr"))
    suppressPackageStartupMessages(require("doParallel"))
    g <- unique(grps)
    dat.g <- as.matrix(dat.g)
    RhpcBLASctl::omp_set_num_threads(1)
    registerDoParallel(cores = n.cores)
    ret <- ldply(rownames(dat.g),function(v){
                  aov.out <- aov(y ~ g,data=data.frame(y=dat.g[v,],g=grps))
                  aov.out.s <- summary(aov.out)
                  t.res.f <- unlist(aov.out.s[[1]]["g",c("F value","Pr(>F)")])
                  aov.out.hsd <- TukeyHSD(aov.out)
                  hsd.name <- rownames(aov.out.hsd$g)
                  t.res.hsd <- c(aov.out.hsd$g[,"diff"],aov.out.hsd$g[,"p adj"])
                  t.res.hsd.minP <- min(aov.out.hsd$g[,"p adj"])
                  j <- which.min(aov.out.hsd$g[,"p adj"])
                  t.res.hsd.minPDiff <- if(is.na(t.res.hsd.minP)) NaN else aov.out.hsd$g[j,"diff"]
                  t.res.hsd.minPCmp <- if(is.na(t.res.hsd.minP)) NaN else rownames(aov.out.hsd$g)[j]
                  ## whether cluster specific ?
                  t.res.spe  <-  sapply(g,function(v){ all( aov.out.hsd$g[grepl(v,hsd.name,perl=T),"p adj"] < FDR.THRESHOLD )  })
                  ## wheter up across all comparison ?
                  is.up <- sapply(g,function(v){  all( aov.out.hsd$g[grepl(paste0(v,"-"),hsd.name),"diff"]>0 ) & all( aov.out.hsd$g[grepl(paste0("-",v),hsd.name),"diff"]<0 ) }) 
                  is.down <- sapply(g,function(v){  all( aov.out.hsd$g[grepl(paste0(v,"-"),hsd.name),"diff"]<0 ) & all( aov.out.hsd$g[grepl(paste0("-",v),hsd.name),"diff"]>0 ) })
                  is.clusterSpecific <- (sum(t.res.spe,na.rm = T) == 1)
                  if(is.clusterSpecific){
                    t.res.spe.lable <- names(which(t.res.spe))  
                    if(is.up[t.res.spe.lable]) { 
                        t.res.spe.direction <- "UP" 
                    }else if(is.down[t.res.spe.lable]) {
                        t.res.spe.direction <- "DOWN"
                    }else{
                        t.res.spe.direction <- "INCONSISTANT"
                    }
                  }else{
                    t.res.spe.lable <- "NA"
                    t.res.spe.direction <- "NA"
                  }
                  if(!is.null(mod) && mod=="cluster.specific") { }
                  structure(c(t.res.f,t.res.hsd,t.res.hsd.minP,t.res.hsd.minPDiff,t.res.hsd.minPCmp,
                              t.res.spe,is.clusterSpecific,t.res.spe.lable,t.res.spe.direction),
                            names=c("F","F.pvalue",paste0("HSD.diff.",hsd.name),paste0("HSD.padj.",hsd.name),
                                    "HSD.padj.min","HSD.padj.min.diff","HSD.padj.min.cmp",
                                    paste0("cluster.specific.",g),"is.clusterSpecific","cluster.lable","cluster.direction"))
                  
            },.progress = "none",.parallel=T)
    #print(str(ret))

    if(!is.null(gid.mapping)){
        cnames <- gid.mapping[rownames(dat.g)]
    }else{
        cnames <- entrezToXXX(rownames(dat.g))
    }
    f.cnames.na <- which(is.na(cnames))
    cnames[f.cnames.na] <- rownames(dat.g)[f.cnames.na]
    ret.df <- data.frame(geneID=rownames(dat.g),geneSymbol=cnames,stringsAsFactors=F)
    ret.df <- cbind(ret.df,ret)
    rownames(ret.df) <- rownames(dat.g)
    ## type conversion
    i <- 3:(ncol(ret.df)-(4+length(g)))
    ret.df[i]<-lapply(ret.df[i],as.numeric)
    i <- (ncol(ret.df)-(length(g)+2)):(ncol(ret.df)-2)
    ret.df[i]<-lapply(ret.df[i],as.logical)
    ## adjust F test's p value
    ret.df$F.adjp <- 1
    ret.df.1 <- subset(ret.df,!is.na(F.pvalue))
    ret.df.1$F.adjp <- p.adjust(ret.df.1[,"F.pvalue"],method = "BH")
    ret.df.2 <- subset(ret.df,is.na(F.pvalue))
    ret.df <- rbind(ret.df.1,ret.df.2)
    ret.df <- ret.df[order(ret.df$F.adjp,-ret.df$F,ret.df$HSD.padj.min),]
    ### select
    ret.df.sig <- subset(ret.df,F.adjp<F.FDR.THRESHOLD & HSD.padj.min<FDR.THRESHOLD & abs(HSD.padj.min.diff)>=FC.THRESHOLD)
    ### output
    write.table(ret.df.sig,file = sprintf("%s.aov.sig.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    if(verbose){
        write.table(ret.df,file = sprintf("%s.aov.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    }
    #print(str(ret.df))
    return(list(aov.out=ret.df,aov.out.sig=ret.df.sig))
}

#' calculate Hypergeometric p value
#' @param gset
#' @param gene_of_intrest
#' @param gene_universe
get.geneSet.hyper <- function(gset,gene_of_intrest,gene_universe,min.size=5,n.cores=4,verbose=F,IDMapping=NULL)
{
    suppressPackageStartupMessages(require("plyr"))
    suppressPackageStartupMessages(require("doParallel"))
    suppressPackageStartupMessages(require("magrittr"))
    registerDoParallel(cores = n.cores)
    uglist.gset <- unique(unlist(gset))
    ## update gene_universe (sequenced & in database)
    gene_universe <- intersect(uglist.gset, gene_universe)
    ## update gene_of_intrest (sequenced & in this category)
    gene_of_intrest <- intersect(uglist.gset, gene_of_intrest)
    ret <- ldply(names(gset),function(v){
                a <- intersect(gset[[v]],gene_universe)
                m <- length(a)
                total <- length(gene_universe)
                n <- total-m
                b <- intersect(a,gene_of_intrest)
                k <- length(gene_of_intrest)
                observed <- length(b)
                expected <- m*k/total
                p.value <- phyper(observed-1, m, n, k, lower.tail = F)
                o.df <- data.frame(stringsAsFactors = F)
                ## observed, expected, ratio, p value
                if(m<min.size){ 
                    o.vec <- c(observed,expected,observed/expected,k,m,n,NA,1)
                }else{
                    o.vec <- c(observed,expected,observed/expected,k,m,n,p.value,1)
                }
                o.df <- t(o.vec)
                colnames(o.df) <- c("observed","expected","ratio","gene.of.intrest","gset.size","compl.size","p.value","p.adj")
                if(verbose){
                    other.info <- c()
                    if(p.value<0.05){
                        other.info <- append(other.info,paste(b,collapse=","))
                        if(!is.null(IDMapping)){
                            other.info <- append(other.info,paste(IDMapping[b],collapse=","))
                        }
                    }else{
                        other.info <- append(other.info,c("",""))
                    }
                    o.df <- cbind(o.df,
                                   data.frame("geneID"=c(other.info[1]),stringsAsFactors = F),
                                   data.frame("geneSymbol"=c(other.info[2]),stringsAsFactors = F))
                }
                o.df
            },.progress = "none",.parallel=T)
           
    ret.df <- data.frame(geneSet=names(gset),stringsAsFactors = F)
    ret.df <- cbind(ret.df,ret)
    ret <- ret.df
    rownames(ret) <- names(gset)
    ret <- subset(ret,!is.na(p.value))
    ret$p.adj <- p.adjust(ret$p.value,method = "BH")
    return(ret[order(ret$p.adj,-ret$ratio),])
}

#' run t-test on given data(row for genes and column for samples)
#' 
#' @param dat.g1
#' @param dat.g2
#' @param out.prefix (string) output prefix
#' @param FDR.THRESHOLD
runTTest <- function(dat.g1,dat.g2,out.prefix,FDR.THRESHOLD=0.05,FC.THRESHOLD=1,verbose=F,n.cores=NULL,gid.mapping=NULL)
{
    suppressPackageStartupMessages(require("plyr"))
    suppressPackageStartupMessages(require("doParallel"))
    registerDoParallel(cores = n.cores)
    dat.g1 <- as.matrix(dat.g1)
    dat.g2 <- as.matrix(dat.g2)
    ret <- ldply(rownames(dat.g2),function(v){
           x <- dat.g1[v,]
           y <- dat.g2[v,]
           x.mean <- mean(x)
           y.mean <- mean(y)
           t.res.fc  <- y.mean-x.mean
           t.res.stat <- NA
           t.res.p <- 1
           if(t.res.fc!=0){
               t.res <- t.test(x=x,y=y)
               t.res.p <- t.res$p.value
               t.res.stat <- t.res$statistic
           }
           structure(c(t.res.p,y.mean,x.mean,t.res.fc,t.res.stat),
                     names=c("p.value","y.mean","x.mean","fc","stat"))

    },.progress = "none",.parallel=T)

    rownames(ret) <- rownames(dat.g2)
    if(!is.null(gid.mapping)){
        cnames <- gid.mapping[rownames(ret)]
    }else{
        cnames <- entrezToXXX(rownames(ret))
    }
    f.cnames.na <- which(is.na(ret))
    cnames[f.cnames.na] <- rownames(ret)[f.cnames.na]
    ret.df <- data.frame(geneID=rownames(ret),geneSymbol=cnames,stringsAsFactors=F)
    ret.df <- cbind(ret.df,ret)
    ret.df$p.adj <- 1
    ret.df.1 <- subset(ret.df,fc!=0)
    ret.df.2 <- subset(ret.df,fc==0)
    ret.df.1$p.adj <- p.adjust(ret.df.1[,"p.value"],method = "BH")
    ret.df <- rbind(ret.df.1,ret.df.2)
    ret.df <- ret.df[order(ret.df$p.adj,ret.df$p.value),]
    ret.df.sig <- subset(ret.df,p.adj<FDR.THRESHOLD & abs(fc)>=FC.THRESHOLD)
    write.table(ret.df.sig,file = sprintf("%s.ttest.sig.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    if(verbose){
        write.table(ret.df,file = sprintf("%s.ttest.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    }
    return(list(ttest.out=ret.df,ttest.out.sig=ret.df.sig))
}

#' run t-SNE analysis on given data(row for genes and column for samples), transpose the data internally
#' 
#' @param data.plot (matrix) input data (row for genes and column for samples)
#' @param out.prefix (string) output prefix
#' @param legend
#' @param col.points 
#' @param col.legend
runTSNEAnalysis <- function(in.data,out.prefix,legend,col.points,col.legend,pch=20,pch.legend=20,inPDF=TRUE,eps.clus=NULL,dims=2,k=NULL,do.dbscan=F,myseed=NULL,width.pdf=10,height.pdf=6,margin.r=1,legend.inset=-0.21,preSNE=NULL,n.cores=NULL,original.space=F,do.clustering=F,data.forDE=NULL,do.scale=T,do.pca=T,distance.metric=NULL,do.determ=F,pca.n=50,...)
{
    suppressPackageStartupMessages(require("Rtsne"))
    suppressPackageStartupMessages(require("dbscan"))
    suppressPackageStartupMessages(require("ks"))
    suppressPackageStartupMessages(require("RColorBrewer"))
    suppressPackageStartupMessages(require("factoextra"))
    suppressPackageStartupMessages(require("seriation"))
    suppressPackageStartupMessages(require("clustertend"))
    suppressPackageStartupMessages(require("cluster"))
    suppressPackageStartupMessages(require("clValid"))
    suppressPackageStartupMessages(require("fields"))
    if(is.null(k)) { k=5; }
    rowVar <- apply(in.data,1,var)
    in.data <- as.matrix(in.data[rowVar>0,,drop=F])
    if(do.scale) {
        dat.plot <- t(scale(t(in.data)))
    }else{
        dat.plot <- in.data
    }

    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }
    dat.plot <- t(dat.plot)
    f.noDup <- !duplicated(dat.plot)
    if(sum(f.noDup)!=m){
        loginfo(sprintf("number of NoDup cells: %d, total cell: %d\n",sum(f.noDup),m))
    }
    sample.used <- rownames(dat.plot)[f.noDup]
    col.points <- col.points[f.noDup]
    dat.plot <- dat.plot[f.noDup,]
    pch <- pch[f.noDup]
    if(!is.null(data.forDE)) { data.forDE <- data.forDE[,sample.used,drop=F] }

    if(is.null(myseed)){ myseed <- as.integer(Sys.time()) }
    set.seed(myseed)
    loginfo(sprintf("set.seed(%s) for tsne\n",myseed))
    ret.list <- list()
    loginfo(sprintf("dim of dat.plot, row: %d, column: %d",n,m))
    doit <- function(par.perplexity)
    {
        ###### tSNE embeding ######
        loginfo(sprintf("... begin Rtsne\n"))
        if(is.null(preSNE)) { 
            .tsne.input <- dat.plot
            .tsne.isDistance <- F
            ###.tsne.pca <- T
            .tsne.pca <- do.pca
            if(!is.null(distance.metric)){
                if(distance.metric=="spearman"){
                    .tsne.input <- as.dist(1-cor(t(dat.plot),method = "spearman"))
                    .tsne.isDistance <- T
                    .tsne.pca <- F
                }
            }
            is.TEST <- T
            if(is.TEST){
                save(.tsne.input,par.perplexity,dims,.tsne.isDistance,.tsne.pca,file=sprintf("%s.tSNE.TEST.RData",out.prefix))
            }
            if(do.determ){
                #### .tsne.input must be cell-cell similarity matrix
                diag(.tsne.input) <- 0
                .D = diag(apply(.tsne.input,MARGIN=2,FUN=sum))
                .L = .D - .tsne.input
                .eigen_L = eigen(.L)
                .U = .eigen_L$vectors
                .D = .eigen_L$values
                .U_index = seq(ncol(.U),(ncol(.U)-dims+1))
                tryCatch({
                Rtsne.res  <-  Rtsne(.tsne.input,perplexity=par.perplexity,initial_dims=pca.n,dims=dims,
                                  Y_init=.U[,.U_index],pca=F)
                },error=function(e){
                    if(grepl("Perplexity is too large",e,perl=T)){
                        cat("Perplexity is too large; try to use perplexity=5 now\n")
                        Rtsne.res <<- Rtsne(.tsne.input,perplexity=5,initial_dims=pca.n,dims=dims,
                                  Y_init=.U[,.U_index],pca=F)
                    }else{ print(e) }
                })
            }else{
                tryCatch({ 
                
                    Rtsne.res <- Rtsne(.tsne.input,perplexity=par.perplexity,initial_dims=pca.n,dims=dims,
                                                is_distance=.tsne.isDistance,pca=.tsne.pca) },
                    error=function(e){
                        if(grepl("Perplexity is too large",e,perl=T)){
                            cat("Perplexity is too large; try to use perplexity=5 now\n")
                            Rtsne.res <<- Rtsne(.tsne.input,perplexity=5,initial_dims=pca.n,dims=dims,
                                                is_distance=.tsne.isDistance, pca=.tsne.pca)
                        }else{
                            print(e)
                        }
                    })
            }
        }else{
            Rtsne.res <- preSNE[[as.character(par.perplexity)]]$Rtsne.res
        }
        rownames(Rtsne.res$Y) <- sample.used
        Rtsne.res.dist <- dist(Rtsne.res$Y)
        loginfo(sprintf("... end Rtsne\n"))
        ###### clustering based on tSNE
        if(do.dbscan){
            loginfo(sprintf("... begin dbscan\n"))
            ### optics method
            dbscan.res <- optics(Rtsne.res$Y, eps = 10,  minPts = k)
            ## to do: add this
            #par(mar=c(5,5,4,5))
            #plot(dbscan.res)
            dbscan.res.list <- list()
            dbscan.clustering.validation.df <- c()
            dbscan.diffGene.list <- list()
            ###for(.eps.for.clust in c(1.5,1.8,2.5,3,3.5,4,5,6,7,8,9,10,2))  ## for plexity=10
            if(is.null(eps.clus)) { 
                #eps.clus <- c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3,3.5,4,5)
                eps.clus <- c(0.8,0.9,1,1.1,1.2)
                #eps.clus <- c(1.2)
            }
            for(.eps.for.clust in eps.clus)
            {
                #### dbscan
                cat(sprintf("DEBUG EPS.FOR.CLUST:%4.2f\n",.eps.for.clust))
                dbscan.res <- extractDBSCAN(dbscan.res, eps_cl = .eps.for.clust)
                pdf(file=sprintf("%s.perplexity%d.dims%d.dbscan.eps%4.2f.minPts%d.pdf",
                                 out.prefix,par.perplexity,dims,.eps.for.clust,k),width=width.pdf,height=height.pdf)
                par(mar=c(5,5,4,5))
                plot(dbscan.res,main=sprintf("Reachability Plot (eps=%4.2f)",.eps.for.clust))  ## black is noise
                clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(dbscan.res$cluster))),
                                  names=as.character(unique(dbscan.res$cluster)))
                clusterColor["0"] <- "gray"
                clusterColor <- clusterColor[as.character(sort(as.numeric(names(clusterColor))))]
                ###
                par(mar=c(5,15,4,15))
                kNNdistplot(Rtsne.res$Y, k = k)
                abline(h = .eps.for.clust, lty = 2, lwd=1.5)
                #dbscan.res <- dbscan::dbscan(Rtsne.res$Y, eps = eps, minPts = k)
                loginfo(sprintf("... end dbscan\n"))
                ## cluster color
                .cluster.count <- table(dbscan.res$cluster)
                .cluster.larger1 <- as.numeric(names(.cluster.count)[.cluster.count>1])
                print(clusterColor)
                print(.cluster.count)
                ### save
                out.df <- data.frame(sample=sample.used,cluster=dbscan.res$cluster)
                write.table(out.df,file=sprintf("%s.perplexity%d.dims%d.dbscan.eps%4.2f.minPts%d.clust.txt",
                                                out.prefix,par.perplexity,dims,.eps.for.clust,k),
                            row.names = F,sep = "\t",quote = F)
                ### some plots
                if(max(.cluster.count)>1){
                    par(mar=c(5,15,4,15))
                    ####### hullplot
                    #### error: incorrect dim number (polygon(d[ch,],...))
                    tryCatch( hullplot(Rtsne.res$Y[dbscan.res$cluster %in% .cluster.larger1,],
                                       dbscan.res$cluster[dbscan.res$cluster %in% .cluster.larger1],
                                       xlab="Dim.1",ylab="Dim.2",main=sprintf("Convex Cluster Hulls (eps=%4.2f)",.eps.for.clust)),
                             error=function(e) {print(e);})
                    ####### points plot
                    .opar <- par()
                    layout(mat = matrix(c(1,2),ncol=2),widths = c(0.6,0.4))
                    par(mar=c(5,5,4,2))
                    plot(Rtsne.res$Y[,1],Rtsne.res$Y[,2],
                         col=clusterColor[as.character(dbscan.res$cluster)],
                         cex=min(1.0,0.6*7000/nrow(Rtsne.res$Y)), 
                         pch=pch,main=sprintf("BarnesHutSNE(eps=%4.2f)",.eps.for.clust),xlab="Dim1",ylab="Dim2")
                    par(mar=c(5,0,4,2))
                    plot.new()
                    legend("left",legend=sprintf("cluster%d",unique(sort(dbscan.res$cluster))),
                           fill = NULL,xpd = NA,cex=1.5,pch=pch,border =NA,
                           ncol=as.integer(ceiling(length(unique(dbscan.res$cluster))/20)),
                           col = clusterColor[as.character(unique(sort(dbscan.res$cluster)))])
                    par(.opar)
                    ####### points plot without outlier
                    .opar <- par()
                    layout(mat = matrix(c(1,2),ncol=2),widths = c(0.6,0.4))
                    par(mar=c(5,5,4,2))
                    plot(Rtsne.res$Y[dbscan.res$cluster>0,1],Rtsne.res$Y[dbscan.res$cluster>0,2],
                         col=clusterColor[as.character(dbscan.res$cluster[dbscan.res$cluster>0])],
                         cex=min(1.0,0.6*7000/nrow(Rtsne.res$Y)), 
                         pch=if(length(pch)>1) pch[dbscan.res$cluster>0] else pch,
                         main=sprintf("BarnesHutSNE(eps=%4.2f)",.eps.for.clust),xlab="Dim1",ylab="Dim2")
                    par(mar=c(5,0,4,2))
                    plot.new()
                    legend("left",legend=sprintf("cluster%d",unique(sort(dbscan.res$cluster))),
                           fill = NULL,xpd = NA,cex=1.5,pch=pch,border =NA,
                           ncol=as.integer(ceiling(length(unique(dbscan.res$cluster))/20)),
                           col = clusterColor[as.character(unique(sort(dbscan.res$cluster)))])
                    par(.opar)
                    #######
                }
                dev.off()
                dbscan.res.list[[as.character(.eps.for.clust)]] <- dbscan.res
                ### genes
                .dat.to.test <- t(dat.plot)
                .clust <- dbscan.res$cluster
                .gene.table <- my.clusterMarkerGene(if(is.null(data.forDE)) .dat.to.test[,.clust>0,drop=F] else data.forDE[,.clust>0,drop=F],
                                                     clust=.clust[.clust>0],
                                                     ann.col=col.points[.clust>0],clust.col=clusterColor,
                                                     out.prefix=sprintf("%s.perplexity%d.dims%d.dbscan.eps%4.2f.minPts%d.clust",
                                                                        out.prefix,par.perplexity,dims,.eps.for.clust,k),
                                                     n.cores=NULL,
                                                     if(is.null(data.forDE)) F else T)
                ### evaluation
                if(!original.space){
                    val.res <- my.clusteringValidation(t(Rtsne.res$Y),clust=dbscan.res$cluster,
                                                       out.prefix=sprintf("%s.perplexity%d.dims%d.dbscan.eps%4.2f.minPts%d.clust",
                                                                           out.prefix,par.perplexity,dims,.eps.for.clust,k),
                                                       n.cores=NULL,dist.obj=NULL)
                }else{
                    .clust <- dbscan.res$cluster
                    names(.clust) <- colnames(.dat.to.test)
                    .dat.scale <- t(scale(t(.dat.to.test[,.clust>0,drop=F])))
                    ###.dat.scale.dist <- as.dist(1-cor((.dat.scale),method = "spearman"))
                    val.res <- my.clusteringValidation(.dat.scale,clust=.clust[.clust>0],
                                                       out.prefix=sprintf("%s.perplexity%d.dims%d.dbscan.eps%4.2f.minPts%d.clust",
                                                                           out.prefix,par.perplexity,dims,.eps.for.clust,k),
                                                       n.cores=NULL,dist.obj=NULL)
                }
                dbscan.clustering.validation.df <- rbind(dbscan.clustering.validation.df,
                                                         structure(c(k,.eps.for.clust,
                                                                     ifelse(!is.null(val.res),val.res$c.stats$avg.silwidth,NA),
                                                                     ifelse(!is.null(val.res),val.res$c.stats$dunn,NA),
                                                                     sum(dbscan.res$cluster>0),
                                                                     length(dbscan.res$cluster),
                                                                     length(unique(dbscan.res$cluster[dbscan.res$cluster>0]))),
                                                                   names=c("minPts","eps","avg.silwidth","dunn",
                                                                           "clustered.sample","total.sample","num.clusters")))
                if(!is.null(.gene.table)){
                    dbscan.diffGene.list[[as.character(.eps.for.clust)]] <- .gene.table$aov.res$aov.out.sig
                }else{
                    dbscan.diffGene.list[[as.character(.eps.for.clust)]] <- NULL
                }
            }
            dbscan.clustering.validation.df <- as.data.frame.matrix(dbscan.clustering.validation.df)
            dbscan.clustering.validation.df$avg.silwidth.rank <- rank(dbscan.clustering.validation.df[,"avg.silwidth"],na.last = F)
            dbscan.clustering.validation.df$dunn.rank <- rank(dbscan.clustering.validation.df[,"dunn"],na.last = F)
            dbscan.clustering.validation.df$score <- apply(dbscan.clustering.validation.df[,c("avg.silwidth.rank","dunn.rank")],
                                                           1,mean)
            dbscan.clustering.validation.df.best <- subset(dbscan.clustering.validation.df[ order(dbscan.clustering.validation.df$score,decreasing = T),], clustered.sample/total.sample>0.85)
            write.table(dbscan.clustering.validation.df,
                        file=sprintf("%s.perplexity%d.dims%d.dbscan.minPts%d.clust.validation.txt",
                                     out.prefix,par.perplexity,dims,k), row.names = F,sep = "\t",quote = F)
            write.table(dbscan.clustering.validation.df.best,
                        file=sprintf("%s.perplexity%d.dims%d.dbscan.minPts%d.clust.validation.best.txt",
                                     out.prefix,par.perplexity,dims,k), row.names = F,sep = "\t",quote = F)
            #### plot
            if(nrow(dbscan.clustering.validation.df)>=3){
                pdf(sprintf("%s.perplexity%d.dims%d.dbscan.minPts%d.clust.validation.pdf",
                            out.prefix,par.perplexity,dims,k),width=8,height = 6)
                par(mar=c(3,4,4,2),cex=1.5)
                plot(dbscan.clustering.validation.df$eps, dbscan.clustering.validation.df$avg.silwidth,type="b",xlab="eps",ylab="Avg.Silhouette")
                dev.off()
            }

            ret.list[[as.character(par.perplexity)]] <<- list(Rtsne.res=Rtsne.res,dbscan.res=dbscan.res.list,
                                                              myseed=myseed,f.noDup=f.noDup,
                                                              dbscan.diffGene.list=dbscan.diffGene.list,
                                                              dbscan.clustering.validation.df=dbscan.clustering.validation.df,
                                                              dbscan.clustering.validation.df.best=dbscan.clustering.validation.df.best)
        }else if(do.clustering==T){
            loginfo(sprintf("... begin other clustering\n"))
            clmethods <- c("hierarchical","kmeans","pam")
            #clmethods <- c("hierarchical")
            #print(str(Rtsne.res$Y))
            #print(head(Rtsne.res$Y))
            ####intern <<- clValid(head(t(dat.plot),n=500), nClust = 2:6, clMethods = clmethods, validation = "internal")
            #tryCatch({
            #    intern <<- clValid(head(t(dat.plot),n=500), nClust = 2:9, clMethods = clmethods, validation = "internal")
            #    ###intern <<- clValid(t(Rtsne.res$Y), nClust = 2:9, clMethods = clmethods, validation = "internal")
            #    pdf(file=sprintf("%s.perplexity%d.dims%d.clValid.pdf",
            #                     out.prefix,par.perplexity,dims),width=width.pdf,height=height.pdf)
            #    plot(intern)
            #    dev.off()
            #}, error=function(e) {print(e);})
            ret.list[[as.character(par.perplexity)]] <<- list(Rtsne.res=Rtsne.res,myseed=myseed,f.noDup=f.noDup)
        }else {
            ret.list[[as.character(par.perplexity)]] <<- list(Rtsne.res=Rtsne.res,myseed=myseed,f.noDup=f.noDup)
        }

        ########## tsne plot ##########
        if(inPDF){
            pdf(file=sprintf("%s.perplexity%d.dims%d.pdf",out.prefix,par.perplexity,dims),width=width.pdf,height=height.pdf)
            opar <- par()
            par(mar=c(5,5,4,margin.r),cex.lab=1.5,cex.main=1.5)
        }
        if(all(pch!=20)){
            layout(mat = matrix(c(1,2),ncol=2),widths = c(0.6,0.4))
            par(mar=c(5,5,4,2))
            plot(Rtsne.res$Y[,1],Rtsne.res$Y[,2], t='n', main="BarnesHutSNE",xlab="Dim1",ylab="Dim2")
            points(Rtsne.res$Y,col=col.points,pch=pch,
                   cex=min(0.8,0.6*5000/nrow(Rtsne.res$Y)),...)
            par(mar=c(5,0,4,2))
            plot.new()
            legend("left",legend=legend,fill = NULL,xpd = NA,cex=1.2,pch=pch.legend,border =NA,col = col.legend)
        }
        ### points' map
        layout(mat = matrix(c(1,2),ncol=2),widths = c(0.6,0.4))
        par(mar=c(5,5,4,2))
        plot(Rtsne.res$Y[,1],Rtsne.res$Y[,2], t='n', main="BarnesHutSNE",xlab="Dim1",ylab="Dim2")
        points(Rtsne.res$Y,col=col.points,pch=16,
               cex=min(0.8,0.6*5000/nrow(Rtsne.res$Y)),...)
               ##cex=min(1.0,0.6*7000/nrow(Rtsne.res$Y)),...)
        par(mar=c(5,0,4,2))
        plot.new()
        legend("left",legend=legend,fill = NULL,xpd = NA,cex=1.2,pch=16,border =NA,col = col.legend)

        ### density map
        layout(mat = matrix(c(1,2),ncol=2),widths = c(0.6,0.4))
        par(mar=c(6,6,6,6))
        #plot(Rtsne.res$Y[,1],Rtsne.res$Y[,2], t='n', main="BarnesHutSNE",xlab="Dim1",ylab="Dim2")
        #points(Rtsne.res$Y,col=col.points,pch=16,
        #       cex=min(1.0,0.6*7000/nrow(Rtsne.res$Y)),...)
        .density <- kde(Rtsne.res$Y)
        .zz <- c(10,20,30,40,50,60,70,80,90)
        plot(.density,display="filled.contour2", cont=.zz,xlab="Dim1", ylab="Dim2")
        image.plot(zlim=c(0,.zz[length(.zz)]),legend.only=TRUE, col = c("transparent", rev(heat.colors(length(.zz)))),
                   axis.args=list( at=.zz, labels=sprintf("%s%%",100-.zz)), legend.width=2.5,legend.mar=4.0)
        par(mar=c(5,0,4,2))
        plot.new()
        #legend("left",legend=legend,fill = NULL,xpd = NA,cex=1.2,pch=16,border =NA,col = col.legend)

        if(inPDF){
            dev.off()
        }
        #hopkins.stat <- hopkins(Rtsne.res$Y,n = nrow(Rtsne.res$Y)-1)
        #cat(sprintf("INFO hopkins.stat:%4.4f\n",hopkins.stat))
        #png(sprintf("%s.perplexity%d.dims%d.VAT.png",out.prefix,par.perplexity,dims),500,400)
        #par(mar=c(5,4,4,2))
        #print(fviz_dist(Rtsne.res.dist, gradient = list(low = "steelblue", mid = NULL, high = "white"),show_labels = F) +
        #    theme(plot.margin = unit(c(1,2,1,2), "cm"),
        #          legend.title=element_text(size=16),
        #          legend.text=element_text(size=14)))
        #dev.off()
    }
    sapply(seq(30,30,10), function(i) { 
               tryCatch(doit(i),
                        error = function(e){ print(e); e }, 
                        finally = loginfo(sprintf("runTSNEAnalysis() finally with perplexity %d \n",i)) ) 
               })
    ##sapply(seq(30,30,5), function(i) { tryCatch(doit(i), error = function(e){ print(e); e }, finally = loginfo(sprintf("runTSNEAnalysis() finally with perplexity %d \n",i)) ) } )
    return(ret.list)
}

#' run NMF analysis on given data(row for genes and column for samples)
#' 
#' @param data.plot (matrix) input data (row for genes and column for samples)
#' @param out.prefix (string) output prefix
#' @param ann.col  column annatation of dat.plot
#' @param ann.colors  a list specify how to mapping ann.col to colors
runNMFAnalysis <- function(dat.plot,out.prefix,ann.col,ann.colors)
{
    suppressPackageStartupMessages(require("NMF"))

    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }
    doNMF <- function(dat.plot)
    {
        nmf.res <- nmf(dat.plot, c(2:10,12,14,16), nrun = 100, .opt = "vp8", seed = 123456)
        pdf(file=sprintf("%s.pdf",out.prefix),width=12,height=10)
        print(plot(nmf.res))
        ###nmf.options(grid.patch=TRUE)
        invisible(sapply(seq(length(nmf.res$fit)), function(i,x){ 
                         consensusmap(x[[i]],
                                      annCol = ann.col,
                                      annColors = ann.colors,
                                      fontsize=20) 
                    }, nmf.res$fit))
        invisible(sapply(seq(length(nmf.res$fit)), function(i,x){ 
                         basismap(x[[i]],main = paste0("Metagenes (Rank=",i+1,")"),fontsize=20) 
                    }, nmf.res$fit))
        invisible(sapply(seq(length(nmf.res$fit)), function(i,x){ 
                         coefmap(x[[i]],main = paste0("Metagene contributions in each sample (Rank=",i+1,")"),
                                 annCol = ann.col,
                                 annColors = ann.colors,
                                 fontsize=20)
                    }, nmf.res$fit))
        dev.off()
        save(nmf.res,file = sprintf("%s.RData",out.prefix))
    }
    tryCatch(doNMF(dat.plot), 
             error = function(e) e, 
             finally = loginfo(sprintf("runNMFAnalysis()'s finally")) )
}


#' simple clustering
#' @importFrom cluster
#' @param in.data row genes, column samples
runSimpleClusteringAnalysis <- function(in.data,out.prefix,sampleType,colSet,do.scale=T,k.batch=2:10,B=100,
                                        col.points=NULL,legend,col.legend,myseed=NULL,data.forDE=NULL,method="kmeans",
                                        val.space="tSNE",gid.mapping=NULL,select.ncp=FALSE,transform="none",
                                        opt.rho=NULL,opt.delta=NULL)
{
    suppressPackageStartupMessages(require("cluster"))
    suppressPackageStartupMessages(require("RColorBrewer"))
    suppressPackageStartupMessages(require("factoextra"))
    suppressPackageStartupMessages(require("FactoMineR"))
    suppressPackageStartupMessages(require("densityClust"))
    suppressPackageStartupMessages(require("Rtsne"))
    file.obj <- sprintf("%s.simple.obj.RData",out.prefix)
    if(file.exists(file.obj)){
        lenv <- loadToEnv(file.obj)
        return(lenv[["ret.obj"]])
    }

    if(do.scale) {
        dat.plot <- t(scale(t(in.data)))
    }else{
        dat.plot <- in.data
    }
    dat.plot <- t(dat.plot)
    f.noDup <- !duplicated(dat.plot)
    sample.used <- rownames(dat.plot)[f.noDup]
    col.points <- col.points[f.noDup]
    dat.plot <- dat.plot[f.noDup,]
    in.data <- in.data[,sample.used,drop=F]
    sampleType <- sampleType[f.noDup]
    if(!is.null(data.forDE)) { data.forDE <- data.forDE[,sample.used,drop=F] }

    tsne.res.old <- runTSNEAnalysis(in.data,sprintf("%s.ori",out.prefix),
                                     col.points = col.points,
                                     legend=legend,
                                     col.legend=col.legend,
                                     myseed=myseed,do.dbscan = F,do.scale = do.scale)

    res.list <- list()
    clustering.validation.df <- c()
    diffGene.list <- list()
    tsne.res.de.list <- list()

    ### columns for genes(variables) and rows for samples
    dat.transformed <- dat.plot
    if(transform=="pca"){
        ##pca.res <- PCA(dat.plot,ncp=min(50,nrow(dat.plot),ncol(dat.plot)),graph=F,scale.unit=T)
        pca.res <- runPCAAnalysis(t(dat.plot),sprintf("%s.PCA",out.prefix),sampleType,colSet,ntop=NULL,verbose=FALSE,do.clust=FALSE,do.scale=TRUE,myseed=myseed,gid.mapping=gid.mapping)
        test.pca.res <<- pca.res
        eigengap <- sapply(seq_len(length(pca.res$eig$eigenvalue)-1),function(i){ pca.res$eig$eigenvalue[i]-pca.res$eig$eigenvalue[i+1] })
        eigengap.max <- which.max(eigengap)
        loginfo(sprintf("max eigengap: %d ( %4.4f)\n",eigengap.max,eigengap[eigengap.max]))
        if(select.ncp){ 
            pca.ncp <- max(eigengap.max)
        }else{
            pca.ncp <- ncol(pca.res$ind$coord)
        }
        loginfo(sprintf("number of components to use: %d",pca.ncp))
        dat.transformed <- pca.res$ind$coord[,1:pca.ncp,drop=F]
    }else if(transform=="pca_tsne"){
        dat.transformed <- tsne.res.old$`30`$Rtsne.res$Y
    }else if(transform=="iCor_tsne"){
        library("lsa")
        #.iCor.dat <- cosine(t(dat.plot))
        #.iCor.dat <- stats::dist(dat.plot,upper = T,diag = T)
        #.iCor.dat.mtx <- as.matrix(.iCor.dat)
        #.iCor.iCor.dat <- cor(.iCor.dat.mtx,method="spearman")
        #save(col.points,legend,col.legend,myseed,dat.plot,.iCor.dat,.iCor.dat.mtx,.iCor.iCor.dat,file = sprintf("%s.forTEST.00.RData",out.prefix))
        #.iCor.dat <- .iCor.iCor.dat

        .iCor.dat <- cor(t(dat.plot),method="pearson")
        #.iCor.dat <- cor(t(dat.plot),method="spearman")
        #.iCor.dat <- cor(in.data,method="spearman")
        #.iCor.dat <- 1-cor(t(dat.plot),method="spearman")

        .iCor.dat.tsne <- runTSNEAnalysis(.iCor.dat,sprintf("%s.iCor",out.prefix),
                                     col.points = col.points,
                                     legend=legend,
                                     col.legend=col.legend,
                                     myseed=myseed,do.dbscan = F,do.scale = F,do.pca=F,do.determ=T)

        dat.transformed <- .iCor.dat.tsne$`30`$Rtsne$Y
        rownames(dat.transformed) <- rownames(.iCor.dat)
        save(dat.transformed,col.points,legend,col.legend,myseed,dat.plot,.iCor.dat,file = sprintf("%s.forTEST.RData",out.prefix))
###        tryCatch({ .iCor.dat.tsne <- Rtsne(.iCor.dat,perplexity=30,dims=2) },
###            error=function(e){ if(grepl("Perplexity is too large",e,perl=T)){
###                        cat("Perplexity is too large; try to use perplexity=5 now\n")
###                        .iCor.dat.tsne <- Rtsne(.iCor.dat,perplexity=5,dims=2)
###                    }else{
###                        print(e)
###                    }
###                })
###        dat.transformed <- .iCor.dat.tsne$Y
###        rownames(dat.transformed) <- rownames(.iCor.dat)
    }
    if(is.null(dat.transformed)){ return(NULL) }
    test.dat.transformed <<- dat.transformed
    clust.res <- NULL
    for(k in k.batch){
        if(method=="kmeans"){
            clust.res <- kmeans(dat.transformed,k,iter.max=1000,nstart=50)
        }else if(method=="hclust"){
            clust.res <- eclust(dat.transformed, "hclust", k = k, method = "complete", graph = FALSE)
        }else if(method=="density"){
            if(is.null(clust.res)){
                dist.obj <- dist(dat.transformed)
                density.clust.res <- densityClust::densityClust(dist.obj, gaussian = T)
                if(is.null(opt.rho)){ opt.rho <- quantile(density.clust.res$rho, probs = 0.90) }
                if(is.null(opt.delta)){ opt.delta <- quantile(density.clust.res$delta, probs = 0.95) }
                set.seed(myseed)
                pdf(sprintf("%s.decision.pdf",out.prefix),width = 5,height = 5)
                density.clust.res <- densityClust::findClusters(density.clust.res, rho = opt.rho, delta = opt.delta,plot=T)
                abline(v=opt.rho,lty=2)
                abline(h=opt.delta,lty=2)
                plot(sort(density.clust.res$rho * density.clust.res$delta,decreasing = T),ylab="rho * delta")
                ###
                .density <- kde(dat.transformed)
                par(mar=c(5,4,5,6))
                .zz <- c(10,20,30,40,50,60,70,80,90)
                plot(.density,display="filled.contour2", cont=.zz,xlab="Dim1", ylab="Dim2")
                image.plot(zlim=c(0,.zz[length(.zz)]),legend.only=TRUE, col = c("transparent", rev(heat.colors(length(.zz)))),
                           axis.args=list( at=.zz, labels=sprintf("%s%%",100-.zz)), legend.width=2.0,legend.mar=4.5)
                plot(.density,display="filled.contour2", cont=.zz,xlab="Dim1", ylab="Dim2")
                points(dat.transformed[names(density.clust.res$peaks),,drop=F],pch=3,cex=2,col="black")
                image.plot(zlim=c(0,.zz[length(.zz)]),legend.only=TRUE, col = c("transparent", rev(heat.colors(length(.zz)))),
                           axis.args=list( at=.zz, labels=sprintf("%s%%",100-.zz)), legend.width=2.0,legend.mar=4.5)
                dev.off()
                clust.res <- list("cluster"=density.clust.res$clusters)
                names(clust.res$cluster) <- names(density.clust.res$rho)
                print(str(clust.res$cluster))
            }
        }
        if(method=="density" && length(res.list)>0 ){ break  }
        loginfo(sprintf("... k=%d ",k))
        ### save
        out.df <- data.frame(sample=colnames(in.data),cluster=clust.res$cluster,stringsAsFactors = F)
        write.table(out.df,file=sprintf("%s.k%d.clust.txt",out.prefix,k),
                    row.names = F,sep = "\t",quote = F)

        clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(clust.res$cluster))),
                                  names=sort(as.character(unique(clust.res$cluster))))
        ### genes
        .dat.to.test <- in.data
        if(!is.null(data.forDE)){
            .dat.to.test <- data.forDE
        }
        .clust <- clust.res$cluster
        print("TEST: XXXXX")
        save(.dat.to.test,.clust,col.points,clusterColor,sampleType,colSet,data.forDE,file="XXX.test.RData")

        .gene.table <- my.clusterMarkerGene(.dat.to.test,
                                            clust=.clust,
                                            ann.col=col.points,
                                            clust.col=clusterColor,
                                            out.prefix=sprintf("%s.k%d.clust", out.prefix,k),
                                            n.cores=NULL,
                                            original.labels=F,
                                            sampleType = sampleType,sampleTypeColSet = colSet,gid.mapping = gid.mapping)

                                            ####original.labels=if(is.null(data.forDE)) F else T,
        print("TEST: YYYYYY")
        ##### visualization by tSNE
        ## a reference tSNE
        tsne.res.cls <- runTSNEAnalysis(in.data,sprintf("%s.k%d.viz",out.prefix,k),
                                        col.points = clusterColor[as.character(.clust)],
                                        legend=sprintf("C%s",names(clusterColor)),
                                        col.legend=clusterColor,preSNE = tsne.res.old,
                                        myseed=tsne.res.old$`30`$myseed,do.dbscan = F,do.scale = do.scale)
        ## tSNE using DE genes
        if(!is.null(.gene.table) && nrow(.gene.table$aov.res$aov.out.sig)>30){
            tsne.res.cls.de <- runTSNEAnalysis(.dat.to.test[rownames(.gene.table$aov.res$aov.out.sig),],
                                               sprintf("%s.k%d.de.viz",out.prefix,k),
                            col.points = clusterColor[as.character(.clust)],
                            legend=sprintf("C%s",names(clusterColor)),
                            col.legend=clusterColor,preSNE = NULL,
                            myseed=tsne.res.old$`30`$myseed,do.dbscan = F,do.scale = do.scale)
            tsne.res.de.list[[as.character(k)]] <- tsne.res.cls.de
        }
        ### evaluation
        .dat.to.val <- in.data
        if(val.space=="tSNE"){
            .dat.to.val <- t(tsne.res.old$`30`$Rtsne.res$Y)
        }else if(val.space=="pca" && transform=="pca"){
            .dat.to.val <- t(pca.res$ind$coord)
        }
        val.res <- my.clusteringValidation(.dat.to.val,clust=clust.res$cluster,
                                           out.prefix=sprintf("%s.k%d.clust",out.prefix,k),
                                           n.cores=NULL,dist.obj=NULL)
        ### gather result
        res.list[[as.character(k)]] <- clust.res
        clustering.validation.df <- rbind(clustering.validation.df,
                                             structure(c(k,NA,
                                                         ifelse(!is.null(val.res),val.res$c.stats$avg.silwidth,NA),
                                                         ifelse(!is.null(val.res),val.res$c.stats$dunn,NA),
                                                         sum(clust.res$cluster>0),
                                                         length(clust.res$cluster),
                                                         length(unique(clust.res$cluster))),
                                                       names=c("k","note","avg.silwidth","dunn",
                                                               "clustered.sample","total.sample","num.clusters")))
        if(!is.null(.gene.table)){
            diffGene.list[[as.character(k)]] <- .gene.table$aov.res$aov.out.sig
        }else{
            diffGene.list[[as.character(k)]] <- NULL
        }
    }
    clustering.validation.df <- as.data.frame.matrix(clustering.validation.df)
    print(clustering.validation.df)
    
    clustering.validation.df$avg.silwidth.rank <- rank(clustering.validation.df[,"avg.silwidth"],na.last = F)
    clustering.validation.df$dunn.rank <- rank(clustering.validation.df[,"dunn"],na.last = F)
    clustering.validation.df$score <- apply(clustering.validation.df[,c("avg.silwidth.rank","dunn.rank")],1,mean)
    clustering.validation.df.best <- subset(clustering.validation.df[ order(clustering.validation.df$avg.silwidth,decreasing = T),], clustered.sample/total.sample>0.85)
    ###clustering.validation.df.best <- subset(clustering.validation.df[ order(clustering.validation.df$score,decreasing = T),], clustered.sample/total.sample>0.85)
    write.table(clustering.validation.df, file=sprintf("%s.val.txt",out.prefix), row.names = F,sep = "\t",quote = F)
    write.table(clustering.validation.df.best, file=sprintf("%s.val.best.txt",out.prefix), row.names = F,sep = "\t",quote = F)

    #### plot
    if(nrow(clustering.validation.df)>=3){
        pdf(sprintf("%s.val.pdf",out.prefix),width=8,height = 6)
        par(mar=c(3,4,4,2),cex=1.5)
        plot(clustering.validation.df$k, clustering.validation.df$avg.silwidth,type="b",xlab="k",ylab="Avg.Silhouette")
        dev.off()
    }
    ####
    ret.obj <- list(res.list=res.list,
                clustering.validation.df=clustering.validation.df,
                clustering.validation.df.best=clustering.validation.df.best,
                diffGene.list=diffGene.list,
                tsne.res.de.list=tsne.res.de.list,
                tsne.res.old=tsne.res.old,
                sample.used=sample.used)
    save(ret.obj,file=file.obj)
    return(ret.obj)

    ### gap ??

    #print(table(sampleType,k.res$cluster))
    #gskmn <- clusGap(dat.plot, FUN = kmeans, nstart = 50, K.max = 30, B = B)
        
    #pdf(file=sprintf("%s.clusGap.pdf",out.prefix),width=10,height=8)
    #par(mar=c(5,5,4,8),cex.lab=1.5)
    #plot(gskmn)
    #dev.off()

}

plot.tsne.points <- function(x,out.file,tsne.col.points,col.tsne.legend,tsne.legend,pch,nclusters=NULL,
                             peak=NULL,main="BarnesHutSNE",x.dim=c(1,2))
{
    pdf(out.file,width=10,height=6)
    tryCatch({
    opar <- par()
    par(mar=c(5,5,4,1),cex.lab=1.5,cex.main=1.5)
    layout(mat = matrix(c(1,2),ncol=2),widths = c(0.6,0.4))
    par(mar=c(5,5,4,2))
    #print(head(x))
    plot(x[,x.dim[1]],x[,x.dim[2]], col=tsne.col.points, cex=auto.point.size(nrow(x)),
         pch=pch,main=main,xlab=sprintf("Dim%d",x.dim[1]),ylab=sprintf("Dim%d",x.dim[2]))
    if(!is.null(peak)){
        points(x[peak,,drop=F],pch=3,cex=2,col="black")
    }
    par(mar=c(5,0,4,2))
    plot.new()
    legend("left",tsne.legend,
           fill = NULL,xpd = NA,cex=1.5,pch=pch,border =NA,
           ncol=as.integer(ceiling(length(tsne.legend)/20)),
           col = col.tsne.legend)
           ##ncol=as.integer(ceiling(nclusters/20)),
    par(opar)
    },error=function(e){print(e);e})
    dev.off()
}

#' SC3 clustering
#' @importFrom SC3
#' @param in.data row genes, column samples
runSC3Analysis <- function(in.data,out.prefix,sampleType,colSet,do.log.scale=FALSE,myseed=NULL,n.cores=4,gid.mapping=NULL,data.forDE=NULL,clust.k=6,frac.clustered=0.85,update.tsne=FALSE,rm.lowSil=FALSE)
{
    suppressPackageStartupMessages(library("R.utils"))
    file.sc3.obj <- sprintf("%s.sc3.obj.RData",out.prefix)
    if(file.exists(file.sc3.obj)){
        lenv <- loadToEnv(file.sc3.obj)
        return(lenv[["ret.obj"]])
    }

    suppressPackageStartupMessages(require("SC3"))
	suppressPackageStartupMessages(require("RColorBrewer"))
    ### hijack the original code
	if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/SC3/shiny-funcs.R"))
	{
		source("/Share/BP/zhenglt/02.pipeline/cancer/lib/SC3/shiny-funcs.R")
		source("/Share/BP/zhenglt/02.pipeline/cancer/lib/SC3/iwanthue.R")
	}else{
		source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/SC3/shiny-funcs.R")
		source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/SC3/iwanthue.R")
	}
    dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

    if(!is.null(gid.mapping)){
        cnames <- gid.mapping[rownames(in.data)]
    }else{
        cnames <- entrezToXXX(rownames(in.data))
    }
    cnames.na <- which(is.na(cnames))
    cnames[cnames.na] <- rownames(in.data)[cnames.na]
    rownames(in.data) <- cnames

    col.points <- colSet[sampleType]
    ### check samples
    .c.dat.plot <- t(in.data)
    f.noDup <- !duplicated(.c.dat.plot)
    sample.used <- rownames(.c.dat.plot)[f.noDup]
    col.points <- col.points[f.noDup]
    sampleType <- sampleType[f.noDup]
    .c.dat.plot <- .c.dat.plot[f.noDup,]
    in.data <- in.data[,sample.used,drop=F]
    if(!is.null(data.forDE)) { data.forDE <- data.forDE[,sample.used,drop=F] }
    
    colSet <- colSet[names(colSet) %in% unique(sampleType)]
    tsne.res.old  <- NULL
    tsne.res.old <- runTSNEAnalysis(in.data,sprintf("%s.ori",out.prefix),
                                     col.points = colSet[sampleType],
                                     legend=names(colSet),
                                     col.legend=colSet,
                                     myseed=myseed,do.dbscan = F,do.scale = F)
    if(is.null(tsne.res.old)){ return(NULL) }

    ### helper functions
    output.cluster.labels <- function(values,input.param,kk)
    {
        clusts <- cutree(values$hc, kk)
        cell.order <- order.dendrogram(as.dendrogram(values$hc))
        d <- input.param$dataset
        colnames(d) <- clusts
        d <- d[ , cell.order]
        values$original.labels <- input.param$cell.names[cell.order]
        ## prepare a consensus matrix for downloading
        #tmp <- data.frame(values$consensus[cell.order, cell.order])
        #colnames(tmp) <- values$original.labels
        #rownames(tmp) <- seq(length=nrow(tmp))
        #values$consensus.download <- tmp

        values$new.labels <- reindex_clusters(colnames(d))
        colnames(d) <- values$new.labels
        values$dataset <- d
        values$col.gaps <- which(diff(as.numeric(colnames(d))) != 0)
        values$cell.labels <- data.frame(new.labels = as.numeric(values$new.labels), 
                                        original.labels = values$original.labels, 
                                        stringsAsFactors = FALSE)
        rownames(values$cell.labels) <- values$cell.labels$original.labels
        return(values)
    }

    plot.marker.genes <- function(out.prefix,values,input.param,input.threhold,kk,col.args.clusterAndCellType)
    {
        suppressPackageStartupMessages(require("ROCR"))
        d <- values$dataset
        col.gaps <- values$col.gaps

        mark.res <- get_marker_genes(d,
                                     as.numeric(colnames(d)),
                                     as.numeric(input.threhold$auroc.threshold),
                                     as.numeric(input.threhold$p.val.mark))
        if(is.null(mark.res)) {
            cat(sprintf("Unable to find significant marker genes from obtained clusters (k=%s).\n",kk))
            values$marker.genes <- data.frame(new.labels=c(),gene=c(),AUROC=c(),p.value=c(),stringsAsFactors = F)
            return(values)
        }
        ##colnames(mark.res) <- c("AUC","clusts","p.value")
        d.param <- mark_gene_heatmap_param(mark.res, unique(colnames(d)))
        #sc3.res$mark <- TRUE
        values$marker.genes <- data.frame(new.labels = as.numeric(mark.res[,2]),
                                        gene = rownames(mark.res),
                                        AUROC = as.numeric(mark.res[,1]),
                                        p.value = as.numeric(mark.res[,3]),
                                        stringsAsFactors = FALSE)
        write.table(values$marker.genes,sprintf("%s.markerGene.txt",out.prefix),sep = "\t",row.names = F,quote = F)

        dat.plot <- d[rownames(d.param$mark.res.plot), , drop = FALSE]
        colnames(dat.plot) <- as.character(values$hc$order)
        pheatmap::pheatmap(dat.plot,
                         show_colnames = FALSE,
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = d.param$row.ann,
                         annotation_col = col.args.clusterAndCellType$col.ann.df, 
                         annotation_colors = col.args.clusterAndCellType$ann.colors,
                         annotation_names_row = FALSE,
                         gaps_row = d.param$row.gaps, gaps_col = col.gaps,
                         annotation_legend = FALSE,
                         fontsize = 20, fontsize_col = 12*55/ncol(d),fontsize_row = 250/max(nrow(d.param$mark.res.plot),25),
                         cellheight = 10)
        return(values)
    }

    output.outlier <- function(out.prefix,values,input.param,kk,sampleType)
    {
        suppressPackageStartupMessages(require("colorspace"))
        suppressPackageStartupMessages(require("ggplot2"))
        suppressPackageStartupMessages(require("rrcov"))
        # compute outlier cells
        #sc3.res$outl <- FALSE
        d <- values$dataset
        col.gaps <- values$col.gaps
        outl.res <- outl_cells_main(d, input.param$chisq.quantile)
        t <- as.data.frame(outl.res)
        colnames(t)[1] <- "outl"
        t$Cluster <- names(outl.res)
        t$Cells <- 1:dim(t)[1]
        t$Cluster <-  factor(t$Cluster, levels = unique(as.character(sort(as.numeric(t$Cluster)))))

        cols <- iwanthue(length(unique(t$Cluster)))
        Cells <- outl <- Cluster <- NULL

        #sc3.res$outl <- TRUE
        values$cells.outliers <- data.frame(new.labels = as.numeric(names(outl.res)),
                                          original.labels = values$original.labels,
                                          MCD.dist = outl.res,
                                          cell.type = sampleType[values$hc$order],
                                          stringsAsFactors = FALSE)
        rownames(values$cells.outliers) <- values$original.labels
        write.table(values$cells.outliers,sprintf("%s.outlier.txt",out.prefix),sep = "\t",row.names = F,quote = F)

        print(ggplot(t, aes(x = t$Cells, y = t$outl, fill = t$Cluster, color = t$Cluster)) +
            geom_bar(stat = "identity") +
            geom_point() +
            scale_fill_manual(values = cols) +
            scale_color_manual(values = cols) +
            guides(color = FALSE, fill = FALSE) +
            labs(y = "Outlier score") +
            theme_bw())
        ###ggsave(filename=sprintf("%s.outlier.pdf",out.prefix),width=8,height=6)
        return(values)
    }
    plot.DE.genes <- function(out.prefix,values,input.param,kk,sampleType,col.args.clusterAndCellType,FDR.THRESHOLD=0.05,FC.THRESHOLD=1,n.cores=n.cores)
    {
        dat.g <- values$dataset
        d.original.labels <- values$original.labels
        d.new.labels <- values$new.labels
        colnames(dat.g) <- d.original.labels
        sampleType <- sampleType[values$hc$order]
        names(sampleType) <- d.original.labels
        grps.list <- c()
        cmp.grp <- sprintf("Grp%02d",d.new.labels)
        for(j in unique(d.new.labels)){ 
            grps.list[[sprintf("Grp%02d",j)]] <- d.original.labels[d.new.labels==j]
        }
        annDF <- data.frame(cluster=structure(cmp.grp,names=colnames(dat.g)),stringsAsFactors = F)
        annDF[is.na(annDF)] <- "NA"
        annCol <- list(cluster=structure(auto.colSet(length(grps.list),name = "Dark2"),names=unique(cmp.grp)))
        if(length(unique(d.new.labels))>2){
            mgeneTest.out <- runMultiGroupSpecificGeneTest(if(is.null(data.forDE)) dat.g else data.forDE,
                                                           cmp.grp,
                                                           sprintf("%s.FDR%g",out.prefix,100*FDR.THRESHOLD),
                                                           FDR.THRESHOLD=FDR.THRESHOLD,
                                                           FC.THRESHOLD=FC.THRESHOLD, verbose=T,n.cores = n.cores,gid.mapping=gid.mapping)
            
            do.heatmap.plot <- function(dat.g,g.list,extra="",kkk=1)
            {
                dat.g.sig <- dat.g[g.list,,drop=F]
                dat.tmp <- c()
                for(i in seq_along(grps.list))
                {
                    dat.tmp.1 <- dat.g.sig[,grps.list[[i]],drop=F]
                    if(ncol(dat.tmp.1)>=2){
                        dat.tmp <- cbind(dat.tmp,dat.tmp.1[,hclust(dist(t(dat.tmp.1)))$order,drop=F])
                    }else{
                        dat.tmp <- cbind(dat.tmp,dat.tmp.1)
                    }
                }
                dat.g.sig <- dat.tmp
                #annDF <- data.frame(cluster=structure(cmp.grp,names=colnames(dat.g.sig)),stringsAsFactors = F)
                #annDF[is.na(annDF)] <- "NA"
                #annCol <- list(cluster=structure(auto.colSet(length(grps.list),name = "Dark2"),names=unique(cmp.grp)))
                runHierarchicalClusteringAnalysis(dat.g.sig,mytitle = sprintf("K=%d",kk),
                                              pdf.width=18,pdf.height=15,do.clustering.col=F,
                                              sprintf("%s.aov.sig.noClusteringCol.FDR%g%s",
                                                      out.prefix,100*FDR.THRESHOLD,extra), 
                                              sampleType=sampleType[colnames(dat.g.sig)], 
                                              colSet=colSet,
                                              ann.extra.df = annDF,
                                              ann.extra.df.col = annCol,
                                              ann.bar.height = 0.8,
                                              k.row=kkk,clonotype.col=NULL,ntop=NULL, 
                                              row.names.original=TRUE,
                                              complexHeatmap.use=TRUE,verbose=FALSE,gid.mapping=gid.mapping)
            }
            ## aov genes
            g.list <- as.character(rownames(mgeneTest.out$aov.out.sig))
            if(length(g.list)>=3){ 
                do.heatmap.plot(if(is.null(data.forDE)) dat.g else data.forDE,g.list,extra="")
                if(length(g.list)>30){ do.heatmap.plot(if(is.null(data.forDE)) dat.g else data.forDE,head(g.list,n=30),extra=".top30") }
            }
            ## cluster specific genes
            g.list <- as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction!="INCONSISTANT" )))
            if(length(g.list)>=3){ do.heatmap.plot(if(is.null(data.forDE)) dat.g else data.forDE,g.list,extra=".cluster.specific", kk=1) }
            g.list <- c()
            for(t.grp in unique(cmp.grp)){
                g.list <- append(g.list, head(as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction!="INCONSISTANT" & cluster.lable==t.grp ))),n=10))
            }
            ###cat(sprintf("#### TEST kk=%d ####",kk))
            ###print(g.list)
            ###print(str(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction!="INCONSISTANT" & cluster.lable==cmp.grp[1])))
            if(length(g.list)>=3){ do.heatmap.plot(if(is.null(data.forDE)) dat.g else data.forDE,g.list,extra=".cluster.specific.top10", kk=1) }
            ## cluster specific genes & up-regulated
            g.list <- as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction=="UP")))
            if(length(g.list)>=3){ do.heatmap.plot(if(is.null(data.forDE)) dat.g else data.forDE,g.list,extra=".cluster.specific.UP", 
                                                   kk=length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))) }
            g.list <- c()
            for(t.grp in unique(cmp.grp)){
                g.list <- append(g.list, head(as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction=="UP" & cluster.lable==t.grp))),n=10))
            }
            if(length(g.list)>=3){ do.heatmap.plot(if(is.null(data.forDE)) dat.g else data.forDE,
                                                   g.list,extra=".cluster.specific.UP.top10", 
                                                   kk=length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))) }

        }else if(length(unique(d.new.labels))==2){
            if(length(grps.list[[1]])<3 || length(grps.list[[2]])<3)
            {
                loginfo(sprintf("Two few samples: ncol(dat.g1)=%d, ncol(dat.g2)=%d",length(grps.list[[1]]),length(grps.list[[1]])))
            }else{
                ttest.out <<- runTTest(if(is.null(data.forDE)) dat.g[,grps.list[[1]]] else data.forDE[,grps.list[[1]]],
                                       if(is.null(data.forDE)) dat.g[,grps.list[[2]]] else data.forDE[,grps.list[[2]]],
                                       sprintf("%s.%s.FDR%g",out.prefix,"ttest",100*FDR.THRESHOLD),
                                       FDR.THRESHOLD,FC.THRESHOLD,verbose=T,n.cores = n.cores,gid.mapping=gid.mapping)
                do.heatmap.plot <- function(dat.g,g.list,extra="",kkk=1){
                    if(length(g.list)>=3){
                            dat.g.1 <- dat.g[g.list,grps.list[[1]],drop=F]
                            dat.g.2 <- dat.g[g.list,grps.list[[2]],drop=F]
                            dat.g.sig <- cbind(dat.g.1[,hclust(dist(t(dat.g.1)))$order,drop=F], 
                                           dat.g.2[,hclust(dist(t(dat.g.2)))$order,drop=F])
                            runHierarchicalClusteringAnalysis(dat.g.sig,mytitle = sprintf("K=%d",kk),
                                                              pdf.width=18,pdf.height=15,
                                                              sprintf("%s.%s.sig.FDR%g.noClusteringCol%s",
                                                                      out.prefix,"ttest",100*FDR.THRESHOLD,extra), 
                                                              do.clustering.col=F,
                                                              sampleType=sampleType[colnames(dat.g.sig)], 
                                                              colSet=colSet,
                                                              ann.extra.df = annDF,
                                                              ann.extra.df.col = annCol,
                                                              ann.bar.height = 0.8,
                                                              k.row=length(unique(ttest.out$ttest.out.sig$fc>0)),
                                                              clonotype.col=NULL,ntop=NULL,
                                                              row.names.original=FALSE,
                                                              complexHeatmap.use=TRUE,verbose=FALSE,gid.mapping=gid.mapping)
                            dat.g.sig.df <- data.frame(geneSymbol=rownames(dat.g.sig))
                            dat.g.sig.df <- cbind(dat.g.sig.df,dat.g.sig)
                            write.table(dat.g.sig.df, 
                                        file = sprintf("%s.%s.sig.FDR%g.noClusteringCol%s.txt",
                                                       out.prefix,"ttest",100*FDR.THRESHOLD,extra), 
                                        quote = F,sep = "\t",row.names = F,col.names = T)
                    }
                } 
                ## for genes diff significantly
                g.list.1 <- as.character(rownames(ttest.out$ttest.out.sig))
                do.heatmap.plot(if(is.null(data.forDE)) dat.g else data.forDE,g.list.1,extra="")
                if(length(g.list.1)>30){ do.heatmap.plot(if(is.null(data.forDE)) dat.g else data.forDE,head(g.list.1,n=30),extra=".top30") }
                ## 
                g.list.1 <- c()
                g.list.1 <- append(g.list.1,head(as.character(rownames(subset(ttest.out$ttest.out.sig,x.mean > y.mean))),n=10))
                g.list.1 <- append(g.list.1,head(as.character(rownames(subset(ttest.out$ttest.out.sig,y.mean > x.mean))),n=10))
                if(length(g.list.1)>3){ do.heatmap.plot(if(is.null(data.forDE)) dat.g else data.forDE,g.list.1,extra=".ClusterTop10") }
            }
        }
    }
    ### old version(SC3_0.99.25)
    k.max <- min(clust.k,ncol(in.data)-1)
    sc3(in.data, ks = 2:k.max, cell.filter = FALSE,gene.filter = FALSE,log.scale = do.log.scale,interactivity=FALSE,svm.num.cells=25000,use.max.cores=FALSE,n.cores=n.cores)
    ### latest version
    ###sc3(in.data, ks = 2:10, cell.filter = FALSE,gene.filter = FALSE,log.scale = do.log.scale,interactivity=FALSE,svm.num.cells=2500,n.cores=n.cores)
    ###input.threhold=list(auroc.threshold=0.85,p.val.mark=0.05)
    input.threhold=list(auroc.threshold=0.80,p.val.mark=0.05)
    ###save.image(file = sprintf("%s.tmp.RData",out.prefix))
    ret.list=list()
    clustering.validation.df <- c()
    clustering.validation.df.best <- c()
    #for(kk in 2:16)
    for(kk in 2:k.max)
    {
        cluster.toRet <- NULL

        cat(sprintf("kk==%s\n",kk))
        ## extreac result
        sc3.cons<-with(sc3.interactive.arg, cons.table[cons.table[,1]=="spearman" & cons.table[,2]=="spectral" & cons.table[,3]==kk,4][[1]])
        sc3.res<-list(consensus=sc3.cons[[1]],labels=sc3.cons[[2]],hc=sc3.cons[[3]],silh=sc3.cons[[4]])

        ## output cluster lables  
        sc3.res<-output.cluster.labels(sc3.res,sc3.interactive.arg,kk)
        ## plot silhouette
        ##pdf(sprintf("%s.k%s.pdf",out.prefix,kk),width = 10,height = 8)
        ##png(sprintf("%s.k%s.silhouette.png",out.prefix,kk),width = 1000,height = 800)
        pdf(sprintf("%s.k%s.silhouette.pdf",out.prefix,kk),width = 10,height = 8)
        plot(sc3.res$silh)
        dev.off()
        ## plot consensus
        col.ann.df <- data.frame(CellType=sampleType,Cluster=as.character(sc3.res$cell.labels[colnames(in.data),"new.labels"]))
        rownames(col.ann.df) <- colnames(sc3.res$consensus)
        #col.ann.df <- col.ann.df[sc3.res$hc$order,,drop=F]
        #ann.colors <- list(CellType=colSet,Cluster=as.character(seq_len(kk)))
        if(kk<=12){
            ann.colors <- list(CellType=colSet,Cluster=structure(brewer.pal(max(kk,3),"Set3"), .Names=as.character(seq_len(max(kk,3)))))
        }else{
            ann.colors <- list(CellType=colSet,Cluster=structure(colorRampPalette(brewer.pal(12,"Set3"))(kk), .Names=as.character(seq_len(max(kk,3)))))
        }
        png(sprintf("%s.k%s.consensus.png",out.prefix,kk),width = 1000,height = 800)
        pheatmap::pheatmap(sc3.res$consensus,
                           cluster_rows = sc3.res$hc, cluster_cols = sc3.res$hc,
                           cutree_rows = kk, cutree_cols = kk,
                           annotation_col = col.ann.df, annotation_colors = ann.colors,
                           fontsize = 20, fontsize_col = 12*55/max(ncol(in.data),25),
                           labels_col = colnames(in.data),
                           drop_levels = TRUE,
                           show_rownames = FALSE, show_colnames = TRUE)
        dev.off()
        ## plot marker genes
        png(sprintf("%s.k%s.marker.png",out.prefix,kk),width = 1000,height = 800)
        sc3.res<-plot.marker.genes(sprintf("%s.k%s",out.prefix,kk),sc3.res,sc3.interactive.arg,input.threhold,kk,list(col.ann.df=col.ann.df,ann.colors=ann.colors))
        dev.off()
        ## outerlier
        png(sprintf("%s.k%s.outlier.png",out.prefix,kk),width = 1000,height = 800)
        sc3.res<-output.outlier(sprintf("%s.k%s",out.prefix,kk),sc3.res,sc3.interactive.arg,kk,sampleType)
        dev.off()
        
        #tryCatch(
        #    plot.DE.genes(sprintf("%s.k%s",out.prefix,kk),sc3.res,sc3.interactive.arg,kk,sampleType,
        #              list(col.ann.df=col.ann.df,ann.colors=ann.colors),
        #              FDR.THRESHOLD=0.05,FC.THRESHOLD=1,n.cores=8),  
        #    error = function(e) e, 
        #    finally = loginfo(sprintf("plot.DE() finally")) )

        ### cluster labels
        sc3.res$cluster <- structure(sc3.res$new.labels,names=sc3.res$original.labels)
        sc3.res$cluster <- sc3.res$cluster[colnames(in.data)]
        clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(sc3.res$cluster))),
                                  names=sort(as.character(unique(sc3.res$cluster))))
        clusterColor["0"] <- "gray"
        ### first-evaluation
        tryCatch({
        ### genes
        .clust <- sc3.res$cluster
        .dat.to.test <- in.data
        .gene.table <- my.clusterMarkerGene(if(is.null(data.forDE)) .dat.to.test[,.clust>0,drop=F] else data.forDE[,.clust>0,drop=F],
                                            clust=.clust[.clust>0],
                                            ann.col=col.points[.clust>0],
                                            clust.col=clusterColor,
                                            out.prefix=sprintf("%s.k%s.clust", out.prefix,kk),
                                            n.cores=n.cores,
                                            original.labels = T,sampleType = sampleType,sampleTypeColSet = colSet)
                                            ###if(is.null(data.forDE)) F else T)
        },error=function(e){print(e);e})

        if(update.tsne){
            if(is.null(data.forDE)){
                in.data.update <- in.data[rownames(.gene.table$gene.table),,drop=F]
            }else{
                in.data.update <- data.forDE[rownames(.gene.table$gene.table),,drop=F]
            }
            if(nrow(in.data.update)>10){
                tsne.res.update <- NULL
                tsne.res.update <- runTSNEAnalysis(in.data.update,sprintf("%s.update.k%d",out.prefix,kk),
                                         col.points = colSet[sampleType],
                                         legend=names(colSet),
                                         col.legend=colSet,
                                         myseed=myseed,do.dbscan = F,do.scale = F)
                if(is.null(tsne.res.update)){
                    tsne.Y.forVal <- tsne.res.old$`30`$Rtsne.res$Y
                    cat("tsne.res.update is NULL\n")
                }else{
                    tsne.Y.forVal <- tsne.res.update$`30`$Rtsne.res$Y
                    tryCatch({
                                ###tsne.col.points=clusterColor[as.character(sc3.res$cluster)],
                    plot.tsne.points(tsne.Y.forVal, 
                                sprintf("%s.k%s.tsne.updated.all.pdf",out.prefix,kk),
                                tsne.col.points=clusterColor[as.character(sc3.res$cluster[rownames(tsne.Y.forVal)])],
                                col.tsne.legend=clusterColor,
                                tsne.legend=sprintf("C%s",names(clusterColor)),
                                pch=20,nclusters=length(clusterColor))
                    },error=function(e){ print(sprintf("ERROR: %s.k%s.tsne.updated.all.pdf",out.prefix,kk)); print(str(tsne.Y.forVal));print(e);e})
                }
            }else{
                tsne.Y.forVal <- tsne.res.old$`30`$Rtsne.res$Y
            }
        }else{
            tsne.Y.forVal <- tsne.res.old$`30`$Rtsne.res$Y
            
        }
        tryCatch({
        plot.tsne.points(tsne.res.old$`30`$Rtsne.res$Y, 
                        sprintf("%s.k%s.tsne.all.pdf",out.prefix,kk),
                        tsne.col.points=clusterColor[as.character(sc3.res$cluster)],
                        col.tsne.legend=clusterColor,
                        tsne.legend=sprintf("C%s",names(clusterColor)),
                        pch=20,nclusters=length(clusterColor))
        },error=function(e){ print(str(tsne.res.old$`30`$Rtsne.res$Y));print(e);e})

        tryCatch({
        val.res <- my.clusteringValidation(t(tsne.Y.forVal),clust=sc3.res$cluster[rownames(tsne.Y.forVal)],
                                           out.prefix=sprintf("%s.k%d.clust.00",out.prefix,kk),
                                           n.cores=n.cores,dist.obj=NULL)
        ##sil.kk.df <- data.frame(sample=colnames(in.data),stringsAsFactors = F,
        sil.kk.df <- data.frame(sample=rownames(tsne.Y.forVal),stringsAsFactors = F,
                                cluster.SC3=val.res$sil[,1],
                                neighbor=val.res$sil[,2],
                                sil_width=val.res$sil[,3])
        sil.kk.df$cluster <- sil.kk.df$cluster.SC3
        #### outlier by sil
        sil.kk.df$cluster[sil.kk.df$sil_width<0] <- 0
        rownames(sil.kk.df) <- sil.kk.df$sample
        write.table(sil.kk.df,file = sprintf("%s.k%d.clust.sil.txt",out.prefix,kk),row.names = F,quote = F,sep = "\t")
        
        ### return clusters with low sil
        cluster.toRet <- sc3.res$cluster
        #print("TEST")
        #print(str(cluster.toRet))
        ##sc3.res$cluster[sil.kk.df$sil_width<0] <- 0
        sc3.res$cluster[rownames(subset(sil.kk.df,sil_width<0))] <- 0
        ### return clusters without low sil
        if(rm.lowSil){
            cluster.toRet <- sc3.res$cluster
        }
        #if(kk==3) { test.val.res <<- val.res }

        },error=function(e){print(e);e})
        
        ###### vis by tSNE
        .clust <- sc3.res$cluster
        ### without outlier
        tryCatch({
        plot.tsne.points(tsne.Y.forVal, 
                        sprintf("%s.k%s.tsne.grayC0.pdf",out.prefix,kk),
                        tsne.col.points=clusterColor[as.character(sc3.res$cluster[rownames(tsne.Y.forVal)])],
                        col.tsne.legend=clusterColor,
                        tsne.legend=sprintf("C%s",names(clusterColor)),
                        pch=20,nclusters=length(clusterColor))
        .f.plotTSNE <- .clust[rownames(tsne.Y.forVal)]>0
        plot.tsne.points(tsne.Y.forVal[.f.plotTSNE,,drop=F], 
                        sprintf("%s.k%s.tsne.noC0.pdf",out.prefix,kk),
                        tsne.col.points=clusterColor[as.character(sc3.res$cluster[rownames(tsne.Y.forVal[.f.plotTSNE,,drop=F])])],
                        col.tsne.legend=clusterColor,
                        tsne.legend=sprintf("C%s",names(clusterColor)),
                        pch=20,nclusters=length(clusterColor))
        },error=function(e){print(str(tsne.Y.forVal));print(e);e})

        ### re-evaluation
        tryCatch({
        val.res <- my.clusteringValidation(t(tsne.Y.forVal)[,.f.plotTSNE,drop=F],clust=sc3.res$cluster[rownames(tsne.Y.forVal)[.f.plotTSNE]],
                                           out.prefix=sprintf("%s.k%d.clust",out.prefix,kk),
                                           n.cores=n.cores,dist.obj=NULL)
        clustering.validation.df <- rbind(clustering.validation.df,
                                             structure(c(kk,NA,
                                                         ifelse(!is.null(val.res),val.res$c.stats$avg.silwidth,NA),
                                                         ifelse(!is.null(val.res),val.res$c.stats$dunn,NA),
                                                         sum(sc3.res$cluster>0),
                                                         length(sc3.res$cluster),
                                                         length(unique(sc3.res$cluster))),
                                                       names=c("k","note","avg.silwidth","dunn",
                                                               "clustered.sample","total.sample","num.clusters")))
        },error=function(e){print(e);e})
        ###
        if(!is.null(cluster.toRet)){
            sc3.res$cluster <- cluster.toRet
            ret.list[[kk]] <- sc3.res
        }
    }
    if(!is.null(clustering.validation.df)){
        print(str(clustering.validation.df))
        clustering.validation.df <- as.data.frame.matrix(clustering.validation.df)
        clustering.validation.df$avg.silwidth.rank <- rank(clustering.validation.df[,"avg.silwidth"],na.last = F)
        clustering.validation.df$dunn.rank <- rank(clustering.validation.df[,"dunn"],na.last = F)
        #clustering.validation.df$score <- apply(clustering.validation.df[,c("avg.silwidth.rank","dunn.rank")],1,mean)
        clustering.validation.df$score <- clustering.validation.df[,"avg.silwidth"]
        print(clustering.validation.df)
        clustering.validation.df.best <- subset(clustering.validation.df[ order(clustering.validation.df$score,decreasing = T),], clustered.sample/total.sample>frac.clustered)
        write.table(clustering.validation.df, file=sprintf("%s.val.txt",out.prefix), row.names = F,sep = "\t",quote = F)
        write.table(clustering.validation.df.best, file=sprintf("%s.val.best.txt",out.prefix), row.names = F,sep = "\t",quote = F)
    }else{
        cat("clustering.validation.df is NULL\n")
        clustering.validation.df <- data.frame()
        clustering.validation.df.best <- data.frame()
    }

    #### plot
    if(nrow(clustering.validation.df)>=3){
        pdf(sprintf("%s.val.pdf",out.prefix),width=8,height = 6)
        par(mar=c(3,4,4,2),cex=1.5)
        plot(clustering.validation.df$k, clustering.validation.df$avg.silwidth,type="b",xlab="k",ylab="Avg.Silhouette")
        dev.off()
    }
    ret.obj <- list(ret.list=ret.list,
                clustering.validation.df=clustering.validation.df,
                clustering.validation.df.best=clustering.validation.df.best,
                sample.used=sample.used,
                tsne.res=tsne.res.old)
    save(ret.obj,file=file.sc3.obj)
    return(ret.obj)
}

#' run PCA analysis on given data(row for genes and column for samples)
#' 
#' @param data.plot (matrix) input data (row for genes and column for samples)
#' @param out.prefix (string) output prefix
#' @param sampleType  vector of sample type
#' @param colSet  named vector of colors, names is the sample type, values is the mapped colors
#' @param ntop  select the ntop genes with the largest variance
#' @param verbose
runPCAAnalysis <- function(dat.plot,out.prefix,sampleType,colSet,ntop=NULL,verbose=FALSE,do.clust=FALSE,do.scale=TRUE,myseed=NULL,gid.mapping=NULL,do.dimdesc=F,...)
{
    #require("DESeq2")
    suppressPackageStartupMessages(require("ggplot2"))
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("factoextra"))
    suppressPackageStartupMessages(require("FactoMineR"))
    rowVar <- apply(dat.plot,1,var)
    dat.plot <- as.matrix(dat.plot[rowVar>0,,drop=F])
    rowVar <- apply(dat.plot,1,var)
    print(dim(dat.plot))
    print(length(sampleType))
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }
    if(!is.null(ntop)) { 
        select <- order(rowVar,decreasing = TRUE)[seq_len(min(ntop, length(rowVar)))]
        print(str(select))
        print(summary(select))
        print(str(dat.plot))
        dat.plot <- dat.plot[select,]
        n <- nrow(dat.plot)
    }
    dat.plot.t <- t(dat.plot)
    if(!is.null(gid.mapping)){
        cnames <- gid.mapping[colnames(dat.plot.t)]
    }else{
        cnames <- entrezToXXX(colnames(dat.plot.t))
    }
    cnames.na <- which(is.na(cnames))
    cnames[cnames.na] <- colnames(dat.plot.t)[cnames.na]
    colnames(dat.plot.t) <- cnames
    if(n<3) { loginfo(paste0("Too few genes: n=",n)); return(NULL) }
    if(m<3) { loginfo(paste0("Too few samples: n=",m)); return(NULL)
    }
    #### PCA
    if(is.null(myseed)){ myseed <- as.integer(Sys.time()) }
    set.seed(myseed)
    loginfo(sprintf("set.seed(%s) for PCA\n",as.character(myseed)))

    #### prcomp
    ####pca <- prcomp(dat.plot.t,scale=do.scale)
    ####percentVar <- pca$sdev^2/sum(pca$sdev^2)
    ####pca.eig <- (pca$sdev)^2
    ####variance.in.percent <- pca.eig*100/sum(pca.eig)
    ####cumvar.in.percent <- cumsum(variance.in.percent)
    ####pca.eig.df <- data.frame(eig = pca.eig, variance = variance.in.percent, cumvariance = cumvar.in.percent)

    ### FactoMineR
    pca <- PCA(dat.plot.t,ncp=min(50,n,m),graph=F,scale.unit=T)
    loginfo(sprintf("PCA() done.\n"))
    res.desc <- NULL
    if(do.dimdesc){
        res.desc <- dimdesc(pca, axes = c(1:min(10,n,m-1)))
        loginfo(sprintf("dimdesc() done.\n"))
        invisible(sapply(names(res.desc),function(x){
                         oo <- data.frame(gene=rownames(res.desc[[x]][["quanti"]]))
                         oo <- cbind(oo,res.desc[[x]][["quanti"]])
                         write.table(oo,sprintf("%s.variable.association.%s",out.prefix,x),
                                        row.names = F,col.names = T,sep="\t",quote = F)
                 }))
    }

    pca.eig.df <- pca$eig

    pca.var <- get_pca_var(pca)
    pca.ind <- get_pca_ind(pca)

    d.contri <- data.frame(geneSymbol=rownames(pca.var$contrib))
    d.contri <- cbind(d.contri,pca.var$contrib)
    d.ind.coord <- data.frame(sampleType=sampleType, sampleName = colnames(dat.plot))
    d.ind.coord <- cbind(d.ind.coord,pca$ind$coord)
    if(verbose){
        ozfile <- gzfile(sprintf("%s.contribution.txt.gz",out.prefix),open = "w")
        write.table(d.contri,file = ozfile,sep = "\t",col.names = T,row.names = F,quote = F)
        close(ozfile)
        
        ozfile <- gzfile(sprintf("%s.txt.gz",out.prefix),open = "w")
        write.table(d,file = ozfile,sep = "\t",col.names = T,row.names = F,quote = F)
        close(ozfile)
    }
   
    my.plot.score <- function(pc.x,pc.y,pc.x.lab,pc.y.lab,pc.x.pve,pc.y.pve,e.x.lab=NULL,e.y.lab=NULL){
        layout(mat = matrix(c(1,2),ncol=2),widths = c(0.6,0.4))
        opar <- par(mar=c(5,5,4,2),cex.lab=1.5)
        plot(x=NULL,y=NULL, xlim = range(pc.x), ylim = range(pc.y), type = "n", cex.main=1.5,cex=0.8,
             xlab=sprintf("%s: %s%% variance ",ifelse(is.null(e.x.lab),pc.x.lab,e.x.lab), pc.x.pve), 
             ylab=sprintf("%s: %s%% variance ",ifelse(is.null(e.y.lab),pc.y.lab,e.y.lab), pc.y.pve)) 
        invisible(lapply(levels(d.ind.coord$sampleType),function(x){ points(subset(d.ind.coord,sampleType==x,select=pc.x.lab)[,1],
                                                                  subset(d.ind.coord,sampleType==x,select=pc.y.lab)[,1],
                                                                  col=colSet[x],pch=20,cex=auto.point.size(length(pc.x))) }))
        par(mar=c(5,0,4,2))
        plot.new()
        legend("left",legend=names(colSet),fill = NULL,xpd = NA,cex=min(1.2,1.5*7/length(names(colSet))),pch=20,border =NA,col = colSet)
        par(opar)
    }
    ### loading plot
    my.plot.loading <- function(ld.x,ld.y,ld.x.lab,ld.y.lab){
        a <- seq(0, 2*pi, length = 100)
        plot( cos(a), sin(a), type = 'l', col="gray", xlab = ld.x.lab,  ylab = ld.y.lab)
        abline(h = 0, v = 0, lty = 2)
        arrows(0, 0, ld.x, ld.y, length = 0.1, angle = 15, code = 2)
        text(ld.x,ld.y, labels=names(ld.x), col="#0000FF99", cex = 1, adj=1)
        ####arrows(0, 0, pca.var$coord[, 1], pca.var$coord[, 2], length = 0.1, angle = 15, code = 2)
        ####text(pca.var$coord[, 1], pca.var$coord[, 2], labels=rownames(pca.var$coord), cex = 1, adj=1)
    }

    nDim.plot <- min(30,nrow(pca.eig.df))
    save(pca,pca.eig.df,colSet,sampleType,d.ind.coord,out.prefix,file=sprintf("%s.tmp.RData",out.prefix))
    ##source("/Share/BP/zhenglt/02.pipeline/cancer/lib/my.getContrib.R")
    pdf(file=sprintf("%s.contribInd.pdf",out.prefix),width=8,height=5)
    print(fviz_contrib(pca, choice = "ind", axes = 1, top = nDim.plot)+labs(title = "Contribution of cells to PC1")+theme(plot.margin = margin(20,40,20,20), plot.title = element_text(size = rel(2)), axis.title.y = element_text(size = rel(2)), axis.text.y=element_text(size = rel(1.8)), axis.text.x=element_text(size = rel(1.1),vjust=0.85)))
    print(fviz_contrib(pca, choice = "ind", axes = 2, top = nDim.plot)+labs(title = "Contribution of cells to PC2")+theme(plot.margin = margin(20,40,20,20), plot.title = element_text(size = rel(2)), axis.title.y = element_text(size = rel(2)), axis.text.y=element_text(size = rel(1.8)), axis.text.x=element_text(size = rel(1.1),vjust=0.85)) )
    print(fviz_contrib(pca, choice = "ind", axes=1:2, top = nDim.plot)+labs(title = "Contribution of cells to PC1 and PC2")+theme(plot.margin = margin(20,40,20,20), plot.title = element_text(size = rel(2)), axis.title.y = element_text(size = rel(2)), axis.text.y=element_text(size = rel(1.8))))
    dev.off()

    pdf(file=sprintf("%s.contribVar.pdf",out.prefix),width=8,height=5)
    print(fviz_contrib(pca, choice = "var", axes = 1, top = nDim.plot)+labs(title = "Contribution of genes to PC1")+theme(plot.margin = margin(20,40,20,20), plot.title = element_text(size = rel(2)), axis.title.y = element_text(size = rel(2)), axis.text.y=element_text(size = rel(1.8)), axis.text.x=element_text(size = rel(1.1),vjust=0.85)))
    print(fviz_contrib(pca, choice = "var", axes = 2, top = nDim.plot)+labs(title = "Contribution of genes to PC2")+theme(plot.margin = margin(20,40,20,20), plot.title = element_text(size = rel(2)), axis.title.y = element_text(size = rel(2)), axis.text.y=element_text(size = rel(1.8)), axis.text.x=element_text(size = rel(1.1),vjust=0.85)) )
    print(fviz_contrib(pca, choice = "var", axes=1:2, top = nDim.plot)+ labs(title = "Contribution of genes to PC1 and PC2")+ theme(plot.margin = margin(20,40,20,20), plot.title = element_text(size = rel(2)), axis.title.y = element_text(size = rel(2)), axis.text.y=element_text(size = rel(1.8))))
    dev.off()
    
    #my.plot.loading(pca$rotation[,1],pca$rotation[,2],"PC1","PC2")
    #my.plot.loading(pca$rotation[,1],pca$rotation[,3],"PC1","PC3")
    #my.plot.loading(pca$rotation[,2],pca$rotation[,3],"PC2","PC3")

    ### score plot
    pdf(file=sprintf("%s.score.PC1-PC2.pdf",out.prefix),width=10,height=6)
    my.plot.score(d.ind.coord$Dim.1,d.ind.coord$Dim.2,"Dim.1","Dim.2",round(pca.eig.df[1,2],2),round(pca.eig.df[2,2],2),e.x.lab="PC1",e.y.lab="PC2")
    dev.off()
    pdf(file=sprintf("%s.score.PC1-PC3.pdf",out.prefix),width=10,height=6)
    my.plot.score(d.ind.coord$Dim.1,d.ind.coord$Dim.3,"Dim.1","Dim.3",round(pca.eig.df[1,2],2),round(pca.eig.df[3,2],2),e.x.lab="PC1",e.y.lab="PC3")
    dev.off()
    pdf(file=sprintf("%s.score.PC2-PC3.pdf",out.prefix),width=10,height=6)
    my.plot.score(d.ind.coord$Dim.2,d.ind.coord$Dim.3,"Dim.2","Dim.3",round(pca.eig.df[2,2],2),round(pca.eig.df[3,2],2),e.x.lab="PC2",e.y.lab="PC3")
    dev.off()
    pdf(file=sprintf("%s.score.PC1-PC5.pdf",out.prefix),width=10,height=6)
    my.plot.score(d.ind.coord$Dim.1,d.ind.coord$Dim.5,"Dim.1","Dim.5",round(pca.eig.df[1,2],2),round(pca.eig.df[5,2],2),e.x.lab="PC1",e.y.lab="PC5")
    dev.off()

    pdf(file=sprintf("%s.pdf",out.prefix),width=6,height=5)
    par(mar=c(5,5,4,4),cex.lab=1.2)
    ### scree plot
    xat <- barplot(pca.eig.df[1:nDim.plot, 2], names.arg=1:nDim.plot, main = "Variances", xlab = "Principal Components", ylab = "Percentage of variances", col ="steelblue")
    abline(h = 1*100/sum(pca.eig.df[,2]), lty = 2, col = "red", lwd=2)
    barplot(pca.eig.df[1:nDim.plot, 3], names.arg=1:nDim.plot, main = "Variances", xlab = "Principal Components", ylab = "Cumulative Variances", col ="steelblue")
    barplot(pca.eig.df[1:nDim.plot, 1], names.arg=1:nDim.plot, main = "Eigenvalue", xlab = "Principal Components", ylab = "Eigenvalues", col ="steelblue")

    #print(fviz_screeplot(pca, ncp=nDim.plot))
    print(fviz_pca_var(pca, axes=c(1,2), col.var="contrib",geom=c("point","text"),gradient.cols=c("white","blue","red"),ggtheme = theme_minimal()))
    print(fviz_contrib(pca, choice = "var", axes = 1, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 2, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 1:2, top = nDim.plot))
    print(fviz_contrib(pca, choice = "ind", axes = 1, top = nDim.plot))
    print(fviz_contrib(pca, choice = "ind", axes = 2, top = nDim.plot))
    print(fviz_contrib(pca, choice = "ind", axes = 1:2, top = nDim.plot))
    #print(fviz_pca_var(pca, col.var="contrib") + scale_color_gradient2(low="white", mid="blue", high="red", midpoint=50) + theme_minimal())
    ### low sample size lead to low dimensions (<6)
    tryCatch({
        #### extra components
        print(fviz_pca_var(pca, axes=c(3,4), col.var="contrib",geom=c("point","text"),gradient.cols=c("white","blue","red"),ggtheme = theme_minimal()))
        print(fviz_pca_var(pca, axes=c(5,6), col.var="contrib",geom=c("point","text"),gradient.cols=c("white","blue","red"),ggtheme = theme_minimal()))
        print(fviz_pca_var(pca, axes=c(3,6), col.var="contrib",geom=c("point","text"),gradient.cols=c("white","blue","red"),ggtheme = theme_minimal()))
        print(fviz_pca_var(pca, axes=c(4,5), col.var="contrib",geom=c("point","text"),gradient.cols=c("white","blue","red"),ggtheme = theme_minimal()))
        print(fviz_pca_var(pca, axes=c(1,5), col.var="contrib",geom=c("point","text"),gradient.cols=c("white","blue","red"),ggtheme = theme_minimal()))
        ###print(fviz_pca_ind(pca, axes=c(3,6),col.ind=colSet[sampleType],label="none"))
        print(fviz_contrib(pca, choice = "var", axes = 3, top = nDim.plot))
        print(fviz_contrib(pca, choice = "var", axes = 4, top = nDim.plot))
        print(fviz_contrib(pca, choice = "var", axes = 5, top = nDim.plot))
        print(fviz_contrib(pca, choice = "var", axes = 6, top = nDim.plot))
        print(fviz_contrib(pca, choice = "var", axes = c(3,4), top = nDim.plot))
        print(fviz_contrib(pca, choice = "var", axes = c(5,6), top = nDim.plot))
        print(fviz_contrib(pca, choice = "var", axes = c(3,6), top = nDim.plot))
        print(fviz_contrib(pca, choice = "var", axes = c(4,5), top = nDim.plot))
        print(fviz_contrib(pca, choice = "var", axes = c(1,5), top = nDim.plot))
    },error = function(e){ print(e); e })
    #print(fviz_pca_ind(pca, label = "none", col.ind = patientcolors ,habillage = group, addEllipses = T) + theme_minimal())
    tryCatch({
        sampleType.factor <- factor(sampleType)
        print(fviz_pca_ind(pca, label = "none", habillage=sampleType.factor,addEllipses = F,palette = colSet[levels(sampleType.factor)]) + theme_minimal())
    }, error = function(e){ print(e); e }, finally = loginfo(sprintf("runPCAAnalysis() finally")) )
    dev.off()
    return(pca)
}

#' run hierachical clustering analysis on given data(row for genes and column for samples)
#' 
#' @param data.plot (matrix) input data (row for genes and column for samples)
#' @param out.prefix (string) output prefix
#' @param sampleType  vector of sample type
#' @param colSet  named vector of colors, names is the sample type, values is the mapped colors
#' @param k   number of clusters of the samples
#' @param clonotype.col   list for clonotype 
#' @param ntop  select the ntop genes with the largest variance
#' @param complexHeatmap.use  whether to use complexHeatmap
#' @param verbose
runHierarchicalClusteringAnalysis <- function(dat.plot,out.prefix,sampleType=NULL,colSet=NULL,k.row=1,k.col=1,
                                              clonotype.col=NULL,patient.col.list=NULL,ntop=NULL,
                                              complexHeatmap.use=FALSE,verbose=FALSE,gid.mapping=NULL,
                                              row.names.original=FALSE,
                                              pdf.width=16,pdf.height=15,ann.extra.df=NULL,ann.extra.df.col=NULL,ann.bar.height=1.5,
                                              do.scale=TRUE,z.lo=-2.5,z.hi=2.5,z.step=1,z.title="Exp", annotation_legend_param=list(),
                                              do.cuttree=F,ha.col=NULL,
                                              do.clustering.row=T,do.clustering.col=T,clustering.distance="spearman",clustering.method="complete",mytitle="",save.obj=F,...)
{
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("ComplexHeatmap"))
    suppressPackageStartupMessages(require("circlize"))
    suppressPackageStartupMessages(require("gridBase"))
    suppressPackageStartupMessages(require("dendextend"))
	suppressPackageStartupMessages(require("RColorBrewer"))
   
    if(is.null(colSet) && !is.null(sampleType)){
        sampleType.names <- unique(sampleType)
        colSet <- structure(auto.colSet(length(sampleType.names),name="Accent"),names=sampleType.names)
    }
    dat.plot <- as.matrix(dat.plot)
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }
    if(!is.null(ntop)) { 
        rowVar <- apply(dat.plot,1,var)
        select <- order(rowVar,decreasing = TRUE)[seq_len(min(ntop, length(rowVar)))]
        dat.plot <- dat.plot[select,]
        n <- nrow(dat.plot)
    }
    if(!is.null(gid.mapping)){
        cnames <- gid.mapping[rownames(dat.plot)]
    }else{
        cnames <- entrezToXXX(rownames(dat.plot))
    }
    cnames.na <- which(is.na(cnames))
    cnames[cnames.na] <- rownames(dat.plot)[cnames.na]
    if(!row.names.original){ rownames(dat.plot) <- cnames }

    ###rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
    if(verbose)
    {
        ##dat.plot.export <- data.frame(geneID=rownames(dat.plot),geneSymbol=if(!row.names.original) entrezToXXX(rownames(dat.plot)) else rownames(dat.plot) )
        dat.plot.export <- data.frame(geneID=rownames(dat.plot),geneSymbol=cnames)
        dat.plot.export <- cbind(dat.plot.export,dat.plot)
        write.table(dat.plot.export,sprintf("%s.dat.txt",out.prefix),sep="\t",col.names = T,row.names = F,quote = F)
    }
    if(is.null(complexHeatmap.use))
    {
        patientcolors <- colSet[as.character(sampleType)]
        pdf(sprintf("%s.pdf",out.prefix),width=12,height=12)
        ####heatmap.2(dat.plot,col=bluered(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", keysize=1.2, trace="none", margin=c(15, 20), main="Most variable genes",cexRow=min(1.8,55/n),cexCol=min(1.8,55/m),distfun=function(x){ as.dist(1-cor(t(x),method = "spearman")) } )
        heatmap.2(dat.plot,col=bluered(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", keysize=1.2, trace="none", margin=c(15, 20),cexRow=min(1.8,55/n),cexCol=min(1.8,55/m),...)
        legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5)
        dev.off()
    }else
    {
        require("ComplexHeatmap")
        require("circlize")
        require("gridBase")
        require("dendextend")
        require("moduleColor")
        require("dynamicTreeCut")
        
        dat.plot.unscale <- dat.plot
        ###print(dat.plot[1:4,1:8])
        #### scale by row
        if(do.scale)
        {
            rowM <- rowMeans(dat.plot, na.rm = T)
            rowSD <- apply(dat.plot, 1, sd, na.rm = T)
            dat.plot <- sweep(dat.plot, 1, rowM)
            dat.plot <- sweep(dat.plot, 1, rowSD, "/")
            dat.plot[dat.plot < z.lo] <- z.lo
            dat.plot[dat.plot > z.hi] <- z.hi
            ###print(dat.plot[1:4,1:8])
        }else{
            tmp.var <- pretty(abs(dat.plot),n=8)
            #tmp.var <- pretty((dat.plot),n=8)
            ##tmp.fffff <- abs(dat.plot)
            ##tmp.fffff[tmp.fffff>12] <- 12
            ##tmp.var <- pretty(tmp.fffff,n=8)
            
            z.lo <- tmp.var[1]
            z.hi <- tmp.var[length(tmp.var)]
            z.step <- tmp.var[2]-tmp.var[1]
        }
        branch.row <- FALSE
        branch.col <- FALSE
        obj.hclust.col <- NULL
        obj.hclust.row <- NULL
        if(do.clustering.col){
            if(clustering.distance=="spearman" || clustering.distance=="pearson"){
                tryCatch({
                        obj.hclust.col <- hclust(as.dist(1-cor(dat.plot.unscale,method=clustering.distance)),method=clustering.method)
                        branch.col <- color_branches(as.dendrogram(obj.hclust.col),k=k.col)
                    },
                    error = function(e) { cat("using spearman or pearson as distance failed; will try to fall back to use euler distance ... \n"); e }
                    ###error = function(e) { tmp.dat <<- dat.plot.unscale; print(clustering.method); print(k.col); e }
                    )
            }
            if(is.logical(branch.col) && !branch.col){
                obj.hclust.col <- hclust(dist(t(dat.plot.unscale)),method=clustering.method)
                dat.plot.unscale <- dat.plot.unscale[,obj.hclust.col$order]
                dat.plot <- dat.plot[,obj.hclust.col$order]
                sampleType <- sampleType[obj.hclust.col$order]
                #tryCatch({
                #    obj.hclust.col.dend <- as.dendrogram(obj.hclust.col)
                #    branch.col <- color_branches(obj.hclust.col.dend,k=k.col)
                #},error=function(e){ e })
            }
        }
        if(do.clustering.row){
            if(clustering.distance=="spearman" || clustering.distance=="pearson"){
                tryCatch({
                        obj.hclust.row <- hclust(as.dist(1-cor(t(dat.plot.unscale),method=clustering.distance)),method=clustering.method)
                        branch.row <- color_branches(as.dendrogram(obj.hclust.row),k=k.row)
                    },
                    error = function(e) { cat("using spearman or pearson as distance failed; will try to fall back to use euler distance ... \n"); e }
                    )
            }
            if(is.logical(branch.row) && !branch.row){
                obj.hclust.row <- hclust(dist(dat.plot.unscale),method=clustering.method)
                branch.row <- color_branches(as.dendrogram(obj.hclust.row),k=k.row)
            }
        }
        dend.cutree <- NULL
        if(do.cuttree && !is.null(obj.hclust.col)){
            ### cut tree
            tryCatch({
                dend.cutree <- cutree(obj.hclust.col, k=2:10, order_clusters_as_data = T)
                colnames(dend.cutree) <- sprintf("K=%s",colnames(dend.cutree))
                ### output clusters
                dend.cutree.df <- data.frame(sampleID=rownames(dend.cutree))
                dend.cutree.df <- cbind(dend.cutree.df,dend.cutree)
                write.table(dend.cutree.df,sprintf("%s.cutree.txt",out.prefix),sep = "\t",row.names = F,quote = F)
                ### make plot
                pdf(sprintf("%s.cutree.pdf",out.prefix),width=10,height=12)
                par(mar=c(5,4,4,2))
                layout(matrix(c(1,2),nrow = 2),heights = c(0.6,0.4))
                plot(obj.hclust.col,sub="",xlab="",hang=-1,cex=1.0*50/max(m,32))
                par(mar=c(2,4,0,2))
                colSet.cls <- auto.colSet(10)
                col.cls <- t(apply(dend.cutree,1,function(x){ colSet.cls[x] }))
                plotHclustColors(obj.hclust.col, colors=col.cls, cex.rowLabels = 1.1)
                dev.off()
            },error=function(e) e)
        }
        pdf(sprintf("%s.pdf",out.prefix),width=pdf.width,height=pdf.height)
        #par(mar=c(5,14,4,4))
        par(mar=c(4,12,4,4))
        plot.new()
        title(main = mytitle,cex.main=2)
        #legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)
        ### Integrating Grid Graphics Output with Base Graphics Output
        vps <- baseViewports()
        pushViewport(vps$inner, vps$figure, vps$plot)
        
        if(!is.null(sampleType)){
            annDF <- data.frame(sampleType=sampleType)
            annColList <- list(sampleType=colSet)
        }else{
            annDF <- data.frame()
            annColList <- list()
        }
        if(!is.null(ann.extra.df))
        {
            if(ncol(annDF)==0) { 
                annDF <- ann.extra.df 
            }else{ 
                annDF <- cbind(annDF,ann.extra.df) 
            }
            annColList <- c(annColList,ann.extra.df.col)
        }
        if(!is.null(clonotype.col))
        {
            annDF$clonotype=clonotype.col$ctypeII[colnames(dat.plot)]
            annColList$clonotype=clonotype.col$colII
        }
        g.show.legend <- T
        if(do.cuttree && !is.null(dend.cutree)){
            annDF <- cbind(annDF,dend.cutree)
            for(i in seq_len(ncol(dend.cutree))){
                tmp.list <- list(structure(colSet.cls,names=seq_along(colSet.cls)))
                names(tmp.list) <- colnames(dend.cutree)[i]
                annColList <- c(annColList,tmp.list)
                ### list(TReg=list(color <- bar="continuous",legend <- width=unit(2, "cm"),legend <- height=unit(4, "cm"))
                annotation_legend_param <- c(annotation_legend_param,list())
            }
            ann.bar.height <- max(0.5,ann.bar.height/ncol(annDF))
            g.show.legend <- F
            legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)
        }
        #print(str(annDF))
        #print(str(annColList))

        if(save.obj){
            save(dat.plot,z.title,z.lo,z.hi,z.step,do.scale,branch.col,branch.row,
                 pdf.width,pdf.height,mytitle,annDF,annColList,annotation_legend_param,ann.bar.height,
                 file=sprintf("%s.RData",out.prefix))
        }
        print(annColList)
        if(is.null(ha.col) && ( !is.null(sampleType) || !is.null(ann.extra.df))){
            ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = g.show.legend, annotation_legend_param = annotation_legend_param)
        }
        #else{
        #    ha.col <- NULL
        #}
        top_annotation_height <- unit(ann.bar.height * ncol(annDF), "cm")
        ###col = colorRamp2(c(0, 0.609, 1, 10), c("darkblue", "darkblue", "yellow", "red")),
        ####
                    ###col = colorRamp2(seq(-bk.range[2],bk.range[2],length=5), bluered(5),space="LAB"),
                    ##column_names_gp = gpar(fontsize = 12*55/max(m,32)),row_names_gp = gpar(fontsize = 10*55/max(n,32)),
        bk.range <- quantile(abs(dat.plot),probs=c(0.01,1),na.rm=T)
        ht <- Heatmap(dat.plot, name=z.title,
                    ##col = colorRamp2(seq(-bk.range[2],bk.range[2],length=100), 
                    col = colorRamp2(seq(z.lo,z.hi,length=100), 
                                     colorRampPalette(rev(brewer.pal(n = 7, name = ifelse(do.scale,"RdBu","RdYlBu"))))(100), space="LAB"),
                                     ####colorRampPalette(rev(brewer.pal(n = 7, name =  "RdBu")))(100), space="LAB"),
                    column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"),
                    column_names_gp = gpar(fontsize = 12*28/max(m,32)),row_names_gp = gpar(fontsize = 10*28/max(n,32)),
                    show_heatmap_legend = T, row_names_max_width = unit(10,"cm"),
                    top_annotation_height = top_annotation_height,
                    cluster_columns = branch.col,
                    cluster_rows = branch.row,
                    heatmap_legend_param = list(grid_width = unit(0.8, "cm"), 
                                                grid_height = unit(0.8, "cm"), 
                                                at = seq(z.lo,z.hi,z.step),
                                                title_gp = gpar(fontsize = 14, fontface = "bold"),
                                                label_gp = gpar(fontsize = 12), color_bar = "continuous"),
                    top_annotation = ha.col)
        ComplexHeatmap::draw(ht, newpage= FALSE)
        ##for(i in seq_along(names(annColList))){
        if(!is.null(ha.col)){
            for(i in seq_along(names(ha.col@anno_list))){
                decorate_annotation(names(ha.col@anno_list)[i], 
                                    {grid.text(names(ha.col@anno_list)[i], unit(-4, "mm"),gp=gpar(fontsize=14),just = "right")})
            }
        }
        dev.off()
    }
}

qcAndViz<-function(vsd,vstMat,designM,outDir,extra="",sfilter=NULL, gfilter=NULL,pair=FALSE)
{
    if(!is.null(sfilter))
    {
	    vstMat <- vstMat[,sfilter]
	    designM <- designM[sfilter,]
    }
    if(!is.null(gfilter))
    {
	    vstMat <- vstMat[rownames(vstMat) %in% gfilter[!is.na(gfilter)],]
    }
	vstMat <- vstMat[, rownames(designM)]
#    ## visualization
#    pdf(paste(outDir,"/DESeq2.plot.QA.pdf",sep=""))
#    ## MA plot
#    plotMA(res, main="DESeq2", ylim=c(-2,2))
#    ## dispersion plot
#    plotDispEsts(x)
#    dev.off()
    ## PCA
    suppressPackageStartupMessages(require("DESeq2"))
    suppressPackageStartupMessages(require("ggplot2"))
    suppressPackageStartupMessages(require("gplots"))
    cat(sprintf("%s\t---- PCA plot\n", Sys.time()))
    if(pair)
    { 
        plotPCA.dat <- plotPCA(vsd, intgroup=c("sampleType", "patient"),returnData=T)
	    percentVar <- round(100 * attr(plotPCA.dat, "percentVar"))
            ggplot(data = plotPCA.dat, aes_string(x = "PC1", y = "PC2", color = "group")) + geom_point(size = 3) + 
		    xlab(paste0("PC1: ", round(percentVar[1]), "% variance")) + 
		    ylab(paste0("PC2: ", round(percentVar[2]), "% variance"))
	    ggsave(filename=paste(outDir,"/DESeq2.PCA",extra,".pdf",sep=""),width=8,height=6)
	    write.table(plotPCA.dat,paste(outDir,"/DESeq2.PCA",extra,".txt",sep=""),col.names = T,row.names = F,quote = F)

    }else
    {
        plotPCA.dat <- plotPCA(vsd, intgroup=c("sampleType"),returnData=T) 
	    percentVar <- round(100 * attr(plotPCA.dat, "percentVar"))
            ggplot(data = plotPCA.dat, aes_string(x = "PC1", y = "PC2", color = "group")) + geom_point(size = 3) + 
		    xlab(paste0("PC1: ", round(percentVar[1]), "% variance")) + 
		    ylab(paste0("PC2: ", round(percentVar[2]), "% variance"))
	    ggsave(filename=paste(outDir,"/DESeq2.PCA",extra,".pdf",sep=""),width=8,height=6)
	    write.table(plotPCA.dat,paste(outDir,"/DESeq2.PCA",extra,".txt",sep=""),col.names = T,row.names = F,quote = F)

    }
    ## heatmap by most variable genes
    cat(sprintf("%s\t---- clustering variable genes\n", Sys.time()))

    if(length(levels(designM$sampleType))<=9)
    {
	    suppressPackageStartupMessages(require("RColorBrewer"))
	    colSet <- brewer.pal(9,"Set1")
	    patientcolors <- colSet[as.numeric(designM$sampleType)]
    }else
    {
        colSet <- rainbow(length(levels(designM$sampleType)))
	    patientcolors <- colSet[as.numeric(designM$sampleType)]
    }
    #print(patientcolors)
    designM$color <- patientcolors
    legendDat <- unique(designM[,c("sampleType","color")])
    
    rowVar <- apply(vstMat,1,var)
    for(n in c(50,100,150,200,250,300,350,400,450,500))
    {
        ###hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
        select <- order(rowVar,decreasing = TRUE)[seq_len(min(n, length(rowVar)))]
        dat.plot <-vstMat[select,]
        
        dat.plot.export <- data.frame(geneID=rownames(dat.plot),geneSymbol=entrezToXXX(rownames(dat.plot)))
        dat.plot.export <- cbind(dat.plot.export,dat.plot)
	    write.table(dat.plot.export,paste(outDir,"/DESeq2.QC.Cluster.Var.n",n,extra,".txt",sep=""),sep="\t",col.names = T,row.names = F,quote = F)
        
        nn <- min(n, length(rowVar))
        m <- ncol(dat.plot)
        rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
        pdf(paste(outDir,"/DESeq2.QC.Cluster.Var.n",n,extra,".pdf",sep=""),width=12,height=12)
        heatmap.2(dat.plot,col=bluered(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", keysize=1.2, trace="none", margin=c(15, 20), main="Most variable genes",cexRow=min(1.8,55/nn),cexCol=min(1.8,55/m))
        legend("topright",legend=legendDat$sampleType,fill=legendDat$color,border=legendDat$color,cex=1.5)
        dev.off()
    }
    cat(sprintf("%s\t---- qcAndViz finished\n", Sys.time()))

}

qcAndVizMat<-function(vstMat,designM,outDir,colSet=NULL,intgroup=c("sampleType"),ntop=500,extra="",sfilter=NULL, gfilter=NULL,complexHeatmap.use=NULL,clonotype.col=NULL,runNMF=FALSE,...)
{
    dir.create(outDir,showWarnings = F,recursive = T)
    if(!is.null(sfilter))
    {
	    vstMat <- vstMat[, colnames(vstMat) %in% sfilter, drop = FALSE]
	    designM <- designM[rownames(designM) %in% sfilter,, drop = FALSE]
    }
    if(!is.null(gfilter))
    {
	    vstMat <- vstMat[rownames(vstMat) %in% gfilter[!is.na(gfilter)],, drop = FALSE]
    }
    ###print(dim(designM))
    ###print(dim(vstMat))
	vstMat <- vstMat[, rownames(designM), drop = FALSE]
    ###print(str(vstMat))
    ###return(NULL)

    if(is.null(colSet))
    {
        nLevel <- length(levels(designM$sampleType))
        if(nLevel<=9)
        {
            suppressPackageStartupMessages(require("RColorBrewer"))
            colSet <- brewer.pal(nLevel,"Set1")
            names(colSet) <- levels(designM$sampleType)
        }else
        {
            colSet <- auto.colSet(n=nLevel)
            names(colSet) <- levels(designM$sampleType)
        }
    }else
    {
      colSet <- colSet[names(colSet) %in% levels(designM$sampleType)]
    }
    suppressPackageStartupMessages(require("ggplot2"))
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("factoextra"))
    rowVar <- apply(vstMat,1,var)
    f.var <- rowVar > 0
    vstMat <- vstMat[f.var,]
    rowVar <- rowVar[f.var]
    ###print(dim(vstMat))
    #rowVar <- apply(vstMat,1,function(x){ var(x)/(mean(x)^2) } )
    if(is.null(ntop))
    {
        ntop=length(rowVar)
    }
    select <- order(rowVar, decreasing = TRUE)[seq_len(min(ntop, length(rowVar)))]
    
    ## PCA
    runPCAAnalysis(vstMat[select, ],sprintf("%s/PCA%s",outDir,extra),
                   designM$sampleType,colSet,ntop=NULL,main=sub("^.","",extra))
    ## t-SNE
    runTSNEAnalysis(vstMat[select, ],sprintf("%s/tsne%s",outDir,extra),
                    names(colSet),colSet[as.character(designM$sampleType)],colSet)
    ### hclustering
    for(n in unique(c(50,100,150,200,250,500,ntop)))
    #for(n in c(50,100,150,200,250,300,350,400,450,500,1000,2000,ntop))
    {
        if(!is.null(clonotype.col))
        {
            runHierarchicalClusteringAnalysis(vstMat,
                                              sprintf("%s/QC.Cluster.Var.n%s%s.share",outDir,n,extra),
                                              designM$sampleType,colSet,
                                              clonotype.col=clonotype.col[["share"]],
                                              ntop=n,
                                              complexHeatmap.use=complexHeatmap.use,
                                              verbose=TRUE,main="Variable Genes")
            runHierarchicalClusteringAnalysis(vstMat,
                                              sprintf("%s/QC.Cluster.Var.n%s%s.strict",outDir,n,extra),
                                              designM$sampleType,colSet,
                                              clonotype.col=clonotype.col[["strict"]],
                                              ntop=n,
                                              complexHeatmap.use=complexHeatmap.use,
                                              verbose=TRUE,main="Variable Genes")
        }else{
            runHierarchicalClusteringAnalysis(vstMat,
                                              sprintf("%s/QC.Cluster.Var.n%s%s",outDir,n,extra),
                                              designM$sampleType,colSet,
                                              clonotype.col=clonotype.col,
                                              ntop=n,
                                              complexHeatmap.use=complexHeatmap.use,
                                              verbose=TRUE,main="Variable Genes")
        }
        ## sample distance
#        if(!file.exists(paste(outDir,"/DESeq2.QC.SampleDist",extra,".pdf")))
#        {
#            distsRL <- dist(t(vstMat))
#            mat <- as.matrix(distsRL)
#            pdf(paste(outDir,"/DESeq2.QC.SampleDist",extra,".pdf",sep=""),width=12,height=12)
#            heatmap.2(mat, RowSideColors=patientcolors, ColSideColors=patientcolors, density.info="none", keysize=1.2, trace="none", margin=c(15, 15), main="Sample Distance",cexRow=1.8,cexCol=1.8)
#            dev.off()
#        }
    }
    # NMF
    if(runNMF)
    {
        runNMFAnalysis(vstMat[select, ],sprintf("%s/NMF%s",outDir,extra),
                       designM[,"sampleType",drop=F],list(sampleType=colSet))
    }
    loginfo("qcAndVizMat() finished")
}


#read.ExpData <- function(exp.file,design.file) {
read.ExpData <- function(exp.file,design=NULL) {
  if( grepl(".RData$",tab.file,perl=T) )
  {
      suppressPackageStartupMessages(library("R.utils"))
      lenv <- loadToEnv(exp.file)
      #vsd.blind <- lenv[["vsd.blind"]]
      in.table <- lenv[["vstMat.blind"]]

  }else{
      in.table <- read.table(exp.file,sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
      in.table <- as.matrix(in.table[,-1])
  }
  #design<-read.table(design.file,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
  if(is.null(design))
  {
      in.table
  }else
  {
      in.table[,rownames(design)]
  }
}

read.SampleTypeColor <- function(in.file) {
    if(is.null(in.file) || !file.exists(in.file)) { return(NULL); }
    in.table <- read.table(in.file,header = T,check.names = F,stringsAsFactors = F,sep="\t",comment.char="!")
    ret <- in.table$color
    names(ret) <- in.table$sampleType
    return(ret)
}

patientColorListFromMyDesign <- function(d) {
    ctype <- d$patient
    names(ctype) <- rownames(d)
    c.palette.fct <- as.factor(d$patient)
    c.palette <- seq_along(levels(c.palette.fct))+1
    names(c.palette) <- levels(c.palette.fct)
    list(ctype=ctype,col=c.palette)
}

read.clonotype <- function(in.file,ctype.col) {
    suppressPackageStartupMessages(require("RColorBrewer"))
    if(is.null(in.file) || !file.exists(in.file)) { return(NULL) }
    tryCatch({ in.table <- read.table(in.file,row.names = "Cell_Name",header = T,check.names = F,stringsAsFactors = F,sep="\t",comment.char="!",quote = "") },
                error=function(e){
                    in.table <<- read.table(in.file,row.names = "sample",header = T,check.names = F,stringsAsFactors = F,sep="\t",comment.char="!",quote = "")
                })
    ctype <- in.table[,ctype.col]
    ctype <- strsplit(x = ctype,split = ":",perl = T)
    ###ctype <- sapply(ctype,function(x){ if(as.integer(x[2])>1) { x[1] } else { "NoClonal" }  })
    ctypeI <- sapply(ctype,function(x){ if(as.integer(x[2])>1) { "Clonal" } else { "NoClonal" }  })
    ctypeII <- sapply(ctype,function(x){ n=as.integer(x[2]); if(n>4) n=4; sprintf("N%d",n) })
    ctypeIII <- sapply(ctype,function(x){ n=as.integer(x[2]); if(n==1) { sprintf("unique") } else { sprintf("%s:%d",x[1],n) } })
    names(ctypeI) <- rownames(in.table)
    names(ctypeII) <- rownames(in.table)
    names(ctypeIII) <- rownames(in.table)
    #ctypeI <- as.factor(ctypeI)
    #ctypeII <- as.factor(ctypeII)
    c.type.level.I <- unique(ctypeI)
    c.type.level.II <- unique(ctypeII)
    c.type.level.III <- unique(ctypeIII)
    ##c.palette <- rainbow(length(c.type.level))
    #names(c.palette) <- c.type.level
    ##c.palette[names(c.palette)=="NoClonal"]="gray60"
    ##c.palette[names(c.palette)=="Clonal"]="#7FC97F"
    c.palette.I <- structure(auto.colSet(n = length(c.type.level.I) ,name = "Paired"),names=c.type.level.I)
    ###c.palette.II <- structure(auto.colSet(n = length(c.type.level.II) ,name = "Paired"),names=c.type.level.II)
    c.palette.II <- structure(brewer.pal(11,"BrBG")[c(7,8,9,10)],names=c("N1","N2","N3","N4"))
    c.palette.III <- structure(auto.colSet(n = length(c.type.level.III)-1,name = "Paired"),names=setdiff(c.type.level.III,"unique"))
    c.palette.I <- c(c.palette.I, structure(c("gray"), names = c("NA")))
    ##c.palette.II <- c(c.palette.II, structure(c("gray","gray"), names = c("NA","N0")))
    c.palette.II <- c(c.palette.II, structure(c("gray"), names = c("NA")))
    c.palette.III <- c(c.palette.III, structure(c("gray","black"), names = c("NA","unique")))
    #c.color <- c.palette.I[ctypeI]
    #names(c.color) <- rownames(in.table)
    list(itable=in.table,ctype=ctypeI,col=c.palette.I,ctypeII=ctypeII,colII=c.palette.II,ctypeIII=ctypeIII,colIII=c.palette.III)
}

runTopGOAnalysis <- function( goiIDs, universeIDs , mapping.db="org.Hs.eg.db" )
{
    require("topGO")
    sapply( c( "MF", "BP", "CC" ), function( ont ) { 
           alg <- factor( as.integer( universeIDs %in% goiIDs ) )
           print(str(alg))
           names(alg) <- goiIDs
           tgd <- new( "topGOdata", ontology=ont, allGenes = alg, nodeSize=5, annot=annFUN.org, mapping=mapping.db)
           resultTopGO <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
           GenTable( tgd, resultTopGO, topNodes=15 ) 
       }, simplify=FALSE )
}

require("methods")
require("DESeq2")
require("statmod")
## class definition
SCDenoise <- setClass("SCDenoise", slots = c(raw.endo = "matrix", normalized.endo = "matrix", final.endo = "matrix", size.factor.endo = "vector",
                                             raw.ERCC = "matrix", normalized.ERCC = "matrix", final.ERCC = "matrix", size.factor.ERCC = "vector",
                                             withERCC = "logical"))
setValidity("SCDenoise",
            function(object) {
              msg <- NULL
              if ( ! is.matrix(object@raw.endo) ){
                msg <- c(msg, "input data must be data.frame")
              }else if ( nrow(object@raw.endo) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@raw.endo) < 2 ){
                msg <- c(msg, "input data must have more than one column")
              }
              if (is.null(msg)) TRUE
              else msg
            }
            )
setMethod("initialize",
          signature = "SCDenoise",
          definition = function(.Object, expdata, ignore.ERCC=FALSE){
            ###.Object@normalized.endo <- expdata
            ###.Object@final.endo <- expdata
            .Object@raw.ERCC <- expdata[grepl(pattern = "^ERCC-",x = rownames(expdata),perl = T),]
            .Object@raw.endo <- expdata[!grepl(pattern = "^ERCC-",x = rownames(expdata),perl = T),]
            if(nrow(.Object@raw.ERCC)==0 || ignore.ERCC==TRUE) { 
                .Object@withERCC <- FALSE 
            }else { 
                .Object@withERCC <- TRUE 
            }
            validObject(.Object)
            return(.Object)
          }
          )

setGeneric("SCDenoise.normalize", function(object,useERCCSizeFactor=FALSE) standardGeneric("SCDenoise.normalize"))
setMethod("SCDenoise.normalize",
          signature = "SCDenoise",
          definition = function(object,useERCCSizeFactor=FALSE) {
              endo.sf <- estimateSizeFactorsForMatrix(object@raw.endo)
              object@size.factor.endo <- endo.sf
              if(object@withERCC==FALSE) { 
                  useERCCSizeFactor <- FALSE 
                  ######object@size.factor.ERCC <- NULL
                  object@normalized.ERCC <- object@raw.ERCC
              }else
              {
                  try({
                      ERCC.sf <- estimateSizeFactorsForMatrix(object@raw.ERCC)
                      object@size.factor.ERCC <- ERCC.sf
                      object@normalized.ERCC <- t( t(object@raw.ERCC) / ERCC.sf )
                  })
              }
              if(useERCCSizeFactor) {
                  object@normalized.endo <- t( t(object@raw.endo) / ERCC.sf )
              }else{
                  object@normalized.endo <- t( t(object@raw.endo) / endo.sf )
              }
              return(object)
          }
          )

setGeneric("SCDenoise.fitTechnicalNoise", function(object, use_ERCC = TRUE,fit_type = "counts",plot=TRUE, fit_opts=NULL) standardGeneric("SCDenoise.fitTechnicalNoise"))
setMethod("SCDenoise.fitTechnicalNoise",
          signature = "SCDenoise",
          definition = function(object, use_ERCC = TRUE,fit_type = "counts",plot=TRUE, fit_opts=NULL) {
                nCountsEndo <- object@normalized.endo
                nCountsERCC <- object@normalized.ERCC
                #### check for parameters
                if(object@withERCC==FALSE && use_ERCC==TRUE) { 
                    cat("No ERCC available, use endo-genes to do fitting (set use_ERCC to FALSE)\n")
                    use_ERCC <- FALSE 
                }
                ##if( use_ERCC==FALSE &&  (fit_type %in% c('counts','logvar'))){
                ##    warning("Without ERCCs 'fit_type' 'log' is recommedned")
                ##}
                if((fit_type %in% c('counts', 'log','logvar'))==F){stop("'fit_type' needs to be 'counts', 'log' or 'logvar'")}
                ##if(fit_type=="counts" & use_ERCC==FALSE){
                ##   print("Without ERCCs fit needs to be perfromed in log-space")
                ##   use_ERCC = FALSE  
                ##}
 
                if(use_ERCC==TRUE){
                    if(fit_type=="counts"){
                        meansEndo <- rowMeans( nCountsEndo )
                        varsEndo <- rowVars( nCountsEndo )
                        cv2Endo <- varsEndo / meansEndo^2

                        meansERCC <- rowMeans( nCountsERCC )
                        varsERCC <- rowVars( nCountsERCC )
                        cv2ERCC <- varsERCC / meansERCC^2
      
                        #Do fitting of technical noise
                        if(!is.null(fit_opts)){
                            if("mincv2" %in% names(fit_opts)){mincv2 = fit_opts$mincv2}else{mincv2=.3}
                            if("quan" %in% names(fit_opts)){quan = fit_opts$quan}else{quan=0.8}
                        }else{
                            mincv2 = 0.3
                            quan = 0.8
                        }
      
                        #normalised counts (with size factor)
                        minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > mincv2 ) ], quan ) )
                        useForFitA <- meansERCC >= minMeanForFitA
                        fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ), cv2ERCC[useForFitA] )
     
                        ###xi <- mean( 1 / object@size.factor.ERCC )
                        ###a0 <- unname( fitA$coefficients["a0"] )
                        ###a1 <- unname( fitA$coefficients["a1tilde"] - xi ) 
                        #4. Transform to log-space and propagate error
                        eps=1
                        LogNcountsEndo=log10(nCountsEndo+eps)
                        dLogNcountsEndo=1/((meansEndo+eps)*log(10))
                        var_techEndo=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansEndo)*meansEndo^2
                        LogVar_techEndo=(dLogNcountsEndo*sqrt(var_techEndo))^2 #error propagation 
      
                        if(plot==TRUE){
                            #plot fit
                            par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                            plot( meansERCC, cv2ERCC, log="xy", col=1+2*useForFitA, pch=19, xlab = 'Means', ylab = expression("CV" ^ 2))
                            xg <- 10^seq( -3, 5, length.out=100 )
                            lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='blue' )
                            segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
                            meansERCC[useForFitA], fitA$fitted.values, col="gray" )
                            legend('bottomleft',c('Genes used for fit', 'Fit technical noise'),pch=c(19, NA),lty =c(NA,1),col=c('green','blue'),cex=0.8)
                            title(expression("CV" ^ 2 * " ~ Mean (using ERCC)"))
                            
                            #plot fot with all genes
                            plot( meansEndo, cv2Endo, log="xy", col=1, xlab = 'Means', ylab = expression("CV" ^ 2), pch=20, cex=0.3)
                            points(meansERCC, cv2ERCC, col='#0000FFF0', pch=15, cex=0.6)
                            xg <- 10^seq( -3, 5, length.out=100 )
                            lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='#0000FFF0',lwd=2 )

                            # Plot quantile lines around the fit
                            df <- ncol(nCountsEndo) - 1
                            lines( xg, ( coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg ) * qchisq( .975, df ) / df, col="#0000FFF0", lwd=2, lty="dashed" )
                            lines( xg, ( coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg ) * qchisq( .025, df ) / df, col="#0000FFF0", lwd=2, lty="dashed" )

                            legend('bottomleft',c('Endogenous genes','ERCCs', 'Fit technical noise'),pch=c(20,15, NA),lty =c(NA,NA,1),col=c('black','blue','blue'),cex=0.8)        
                            title(expression("CV" ^ 2 * " ~ Mean (using ERCC)"))
                            ###par(mfrow=c(1,1))
                        }
                        res = list()
                        res$fit = fitA
                        res$techNoiseLog = LogVar_techEndo
                    }else{#with ERCCs in log space
                        if(fit_type=="log"){
                            LCountsEndo <- log10(nCountsEndo+1)
                            LmeansEndo <- rowMeans( LCountsEndo )
                            LvarsEndo <- rowVars( LCountsEndo )
                            Lcv2Endo <- LvarsEndo / LmeansEndo^2
                            
                            LCountsERCC = log10(nCountsERCC+1)
                            LmeansERCC <- rowMeans( LCountsERCC )
                            LvarsERCC <- rowVars( LCountsERCC )
                            Lcv2ERCC <- LvarsERCC / LmeansERCC^2
                            
                            if(!is.null(fit_opts)){
                                if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=2}
                            }else{
                                minmean = .5
                            }
                            LogNcountsList=list()
                            useForFitL=LmeansERCC>minmean
                            LogNcountsList$mean=LmeansERCC[useForFitL]
                            LogNcountsList$cv2=Lcv2ERCC[useForFitL]
                            fit_loglin=nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=20,k=1))
                            LogVar_techEndo_logfit <- coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*LmeansEndo)*LmeansEndo^2
                            
                            if(plot==TRUE){
                                plot( LmeansEndo, Lcv2Endo, log="y", col=1,ylim=c(1e-3,1e2),xlab='meansLogEndo',ylab='cv2LogEndo')
                                xg <- seq( 0, 5.5, length.out=100 )
                                lines( xg, coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*xg ),lwd=2,col='green' )
                                points(LmeansERCC, Lcv2ERCC,col='blue',pch=15,cex=1.1)
                                legend('topright',c('Endo','ERCC'),pch=c(1,1,15),col=c('black','blue'))
                            }
                            res = list()
                            res$fit = fit_loglin
                            res$techNoiseLog = LogVar_techEndo_logfit
                        }else{#with ERCCs fit variance in log space with loess
                            LCountsEndo <- log10(nCountsEndo+1)
                            LmeansEndo <- rowMeans( LCountsEndo )
                            LvarsEndo <- rowVars( LCountsEndo )
                            Lcv2Endo <- LvarsEndo / LmeansEndo^2
                            
                            LCountsERCC = log10(nCountsERCC+1)
                            LmeansERCC <- rowMeans( LCountsERCC )
                            LvarsERCC <- rowVars( LCountsERCC )
                            
                            if("span" %in% names(fit_opts)){span = fit_opts$span}else{span=0.8}
                            if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=0.5}
                            
                            useForFitA <- LmeansERCC >= minmean
                            fit_var2 = loess(LvarsERCC[useForFitA] ~ LmeansERCC[useForFitA], span=span, control=loess.control(surface="direct"))
                            xg <- seq( 0, 5.5, length.out=100 )
                            Var_techEndo_logfit_loess <-  predict(fit_var2, xg)
                            
                            minVar_ERCC = min(LvarsERCC[LmeansERCC>3])
        
                            if(any(xg>3 & (Var_techEndo_logfit_loess<0.8*minVar_ERCC))){
                                idx_1 = which(xg>3 & (Var_techEndo_logfit_loess<0.8*minVar_ERCC))[1]
                                idx_end = length(Var_techEndo_logfit_loess)
                                Var_techEndo_logfit_loess[idx_1:idx_end] = 0.8*minVar_ERCC        
                            }
                            
                            if(plot==TRUE){
                                plot( LmeansEndo, LvarsEndo, col=1,ylim=c(1e-3,150.5),log="y",xlab='meansLogEndo',ylab='VarLogEndo')
                                points(LmeansERCC, LvarsERCC,col='blue',pch=15,cex=1.1)
                                lines(xg, Var_techEndo_logfit_loess,lwd=3,col='blue',lty=1)  
                                legend('topright',c('Endo. genes','ERCC', 'Tech. noise fit'),pch=c(1,15,NA), lty = c(NA,NA,1),col=c('black','blue', 'blue'))
                            }
                            
                            #use model for endogenous genes
                            xg=LmeansEndo
                            Var_techEndo_logfit_loess <-  predict(fit_var2, xg)      
                            
                            if(any(xg>3 & Var_techEndo_logfit_loess<0.8*minVar_ERCC)){
                            idx_1 = which(xg>3 & Var_techEndo_logfit_loess<0.8*minVar_ERCC)[1]
                            idx_end = length(Var_techEndo_logfit_loess)
                            Var_techEndo_logfit_loess[idx_1:idx_end] = 0.8*minVar_ERCC       
                            }          
                            
                            res = list()
                            res$fit = fit_var2
                            res$techNoiseLog = Var_techEndo_logfit_loess
                        }
                    }
                }else{#no ERCCs available
                    if(fit_type=="log"){
                        LCountsEndo <- log10(nCountsEndo+1)
                        LmeansEndo <- rowMeans( LCountsEndo )
                        LvarsEndo <- rowVars( LCountsEndo )
                        Lcv2Endo <- LvarsEndo / LmeansEndo^2
                        
                        if(!is.null(fit_opts)){
                            if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{fit_opts$minmean=0.3}
                            if("offset" %in% names(fit_opts)){offset = fit_opts$offset}else{fit_opts$offset=1}
                        }else{
                            fit_opts$minmean = 0.3
                            fit_opts$offset=1
                        }
                        
                        LogNcountsList = list()
                        useForFitL = LmeansEndo>fit_opts$minmean
                        LogNcountsList$mean = LmeansEndo[useForFitL]
                        LogNcountsList$cv2 = Lcv2Endo[useForFitL]
                        fit_loglin = nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=10,k=2))
                        fit_loglin$opts = fit_opts      
                        LogVar_techEndo_logfit <- fit_opts$offset* coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*LmeansEndo)*LmeansEndo^2
                        
                        if(plot==TRUE){
                            par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                            plot( LmeansEndo, Lcv2Endo, log="y", col=1,ylim=c(1e-3,1e2),
                                 xlab='means',ylab=expression("CV" ^ 2),
                                 main=expression("CV" ^ 2 * " ~ Mean (using Endo in log space)"),pch=20,cex=0.3)
                            xg <- seq( 0, 5.5, length.out=100 )
                            lines( xg, fit_opts$offset*coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*xg ),
                                  lwd=2,col='green' )
                        }
                        res = list()
                        res$fit = fit_loglin
                        res$techNoiseLog = LogVar_techEndo_logfit
                    }
                    if(fit_type=='counts'){
                        meansEndo <- rowMeans( nCountsEndo )
                        varsEndo <- rowVars( nCountsEndo )
                        cv2Endo <- varsEndo / meansEndo^2
                        
                        #Do fitting of technical noise
                        if(!is.null(fit_opts)){
                            if("mincv2" %in% names(fit_opts)){mincv2 = fit_opts$mincv2}else{mincv2=.3}
                            if("quan" %in% names(fit_opts)){quan = fit_opts$quan}else{quan=0.8}
                        }else{
                            mincv2 = 0.3
                            quan=0.8
                        }
                        
                        #normalised counts (with size factor)
                        minMeanForFitA <- unname( quantile( meansEndo[ which( cv2Endo > mincv2 ) ], quan ) )
                        useForFitA <- meansEndo >= minMeanForFitA
                        fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansEndo[useForFitA] ),
                        cv2Endo[useForFitA] )
                        
                        #4. Transform to log-space and propagate error
                        eps=1
                        LogNcountsEndo=log10(nCountsEndo+eps)
                        dLogNcountsEndo=1/((meansEndo+eps)*log(10))
                        var_techEndo=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansEndo)*meansEndo^2
                        LogVar_techEndo=(dLogNcountsEndo*sqrt(var_techEndo))^2 #error propagation 
                        
                        if(plot==TRUE){
                            #plot fit 
                            isTmp <- TRUE
                            if(isTmp){
                                par(cex.lab=2.0,cex.main=2.0,cex.axis=1.8,mar=c(5,6,4,2))
                                ###plot( meansEndo, cv2Endo, log="xy", col=1+2*useForFitA, pch=20, cex=0.3, xlab = '', ylab = "",xaxt="s",yaxt="s",ylim=c(0.1,100))
                                plot( meansEndo, cv2Endo, log="xy", col=1+2*useForFitA, pch=20, cex=0.3, xlab = '', ylab = "",xaxt="s",yaxt="s")
                                xg <- 10^seq( -3, 5, length.out=100 )
                                lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='#0000FFF0' )
                                #print(pretty(meansEndo,n=10))
                                #print(pretty(cv2Endo,n=10))
                                #axis(1,c(0.01,1,1e2,1e4),labels = F)
                                #axis(2,c(0.1,0.5,1,5,10,50,100),labels = F)
                                legend('bottomleft',c('Genes used for fit', 'Fit baseline variation'),
                                       pch=c(20, NA),lty =c(NA,1),col=c('green','blue'),cex=1.8)
                                title(expression("CV" ^ 2 * " ~ Mean (using endogeneous genes)"))

                            }else{       
                                par(cex.lab=2.0,cex.main=2.0,cex.axis=1.8,mar=c(5,6,4,2))
                                plot( meansEndo, cv2Endo, log="xy", col=1+2*useForFitA, pch=20, cex=0.3, xlab = 'Means', ylab = expression("CV" ^ 2))
                                xg <- 10^seq( -3, 5, length.out=100 )
                                lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='#0000FFF0' )
                                legend('bottomleft',c('Genes used for fit', 'Fit baseline variation'),
                                       pch=c(20, NA),lty =c(NA,1),col=c('green','blue'),cex=1.8)
                                title(expression("CV" ^ 2 * " ~ Mean (using endogeneous genes)"))
                            }
                        }
                        res = list()
                        res$fit = fitA
                        res$techNoiseLog = LogVar_techEndo
                    }
                    if(fit_type=='logvar'){
                        LCountsEndo <- log10(nCountsEndo+1)
                        LmeansEndo <- rowMeans( LCountsEndo )
                        LvarsEndo <- rowVars( LCountsEndo )
                        Lcv2Endo <- LvarsEndo / LmeansEndo^2
                        
                        if("span" %in% names(fit_opts)){span = fit_opts$span}else{span=0.8}
                        if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=0.5}
                        
                        useForFitA <- LmeansEndo >= minmean
                        fit_var2 = loess(LvarsEndo[useForFitA] ~ LmeansEndo[useForFitA], span=span, control=loess.control(surface="direct"))
                        xg <- seq( 0, 5.5, length.out=100 )
                        Var_techEndo_logfit_loess <-  predict(fit_var2, xg)
                        
                        minVar_ERCC = min(LvarsEndo[LmeansEndo>3])
                        
                        if(any(xg>2.5 & (Var_techEndo_logfit_loess<0.6*minVar_ERCC))){
                            idx_1 = which(xg>2.5 & (Var_techEndo_logfit_loess<0.6*minVar_ERCC))[1]
                            idx_end = length(Var_techEndo_logfit_loess)
                            Var_techEndo_logfit_loess[idx_1:idx_end] = 0.6*minVar_ERCC        
                        }
                        
                        if(plot==TRUE){
                            par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                            plot( LmeansEndo, LvarsEndo, col=1,ylim=c(1e-3,150.5),log="y",xlab='means',ylab='Var',main=expression("Var~ Mean (using Endo in log space)"),pch=20,cex=0.3)
                            lines(xg, Var_techEndo_logfit_loess,lwd=3,col='blue',lty=1)  
                            legend('topright',c('Endo. genes', 'Tech. noise fit'),pch=c(1,NA), lty = c(NA,1),col=c('black', 'blue'))
                        }
                        
                        #use model for endogenous genes
                        xg=LmeansEndo
                        Var_techEndo_logfit_loess <-  predict(fit_var2, xg)      
                        
                        if(any(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_ERCC)){
                            idx_1 = which(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_ERCC)[1]
                            idx_end = length(Var_techEndo_logfit_loess)
                            Var_techEndo_logfit_loess[idx_1:idx_end] = 0.6*minVar_ERCC       
                        }  
                        
                        res = list()
                        res$fit = fit_var2
                        res$techNoiseLog = Var_techEndo_logfit_loess
                    }
                }
                res$fit_opts = fit_opts
                res  
            })

setGeneric("SCDenoise.getVariableGenes", function(object, fit, method = "fit", threshold = 0.1, fit_type=NULL, plot=T, fitB=NULL) standardGeneric("SCDenoise.getVariableGenes"))
setMethod("SCDenoise.getVariableGenes",
            signature = "SCDenoise",
            definition = function(object, fit, method = "fit", threshold = 0.1, fit_type=NULL, plot=T, fitB=NULL) {
                nCountsEndo <- object@normalized.endo
                nCountsERCC <- object@normalized.ERCC
                sfEndo <- object@size.factor.endo
                sfERCC <- object@size.factor.ERCC
                res.list <- list()
                if(!(method %in% c("fdr","fit"))){ stop("'method' needs to be either 'fdr' or 'fit'") }
                if(is.null(fit_type)){
                  print("No 'fit_type' specified. Trying to guess its from parameter names")
                  if("a0" %in% names(coefficients(fit)) & "a1tilde" %in% names(coefficients(fit))){fit_type="counts"}else{
                    if("a" %in% names(coefficients(fit)) & "k" %in% names(coefficients(fit))){fit_type="log"}else{
                      if(is.call(fit$call)){fit_type="logvar"}
                    }
                  }
                  print(paste("Assuming 'fit_type' is ","'",fit_type,"'",sep=""))
                }
                if(is.null(fit_type)){stop("Couldn't guess fit_type. Please specify it or run the fitTechnicalNoise function to obtain the fit")}
                if(!(fit_type %in% c("counts","log", "logvar")) & !is.null(fit_type)){ stop("'fit_type' needs to be either 'fdr' or 'fit'") }
                if(method=='fdr' & fit_type!="counts"){stop("method='fdr', can only be used with fit_type 'counts'")}
                if(method=='fdr' & object@withERCC==FALSE) { method <- "fit" }
                if(method=='fdr' & (is.null(sfERCC) | is.null(sfEndo))){stop("Please specify sfERCC and sfEndo when using method='fdr'")}
                if(method=='fdr'){
                  meansEndo <- rowMeans( nCountsEndo )
                  varsEndo <- rowVars( nCountsEndo )
                  cv2Endo <- varsEndo / meansEndo^2
                   
                  meansERCC <- rowMeans( nCountsERCC )
                  varsERCC <- rowVars( nCountsERCC )
                  cv2ERCC <- varsERCC / meansERCC^2
              
                  minBiolDisp <- .5^2
                    
                  xi <- mean( 1 / sfERCC )
                  ###a0 <- unname( fit$coefficients["a0"] )
                  ###a1 <- unname( fit$coefficients["a1tilde"] - xi )
                  m <- ncol(nCountsEndo)
                  psia1thetaA <- mean( 1 / sfEndo ) + ( coefficients(fit)["a1tilde"] - xi ) * mean( sfERCC / sfEndo )
                  cv2thA <- coefficients(fit)["a0"] + minBiolDisp + coefficients(fit)["a0"] * minBiolDisp
                  testDenomA <- ( meansEndo * psia1thetaA + meansEndo^2 * cv2thA ) / ( 1 + cv2thA/m )
                  pA <- 1 - pchisq( varsEndo * (m-1) / testDenomA, m-1 )
                  padjA <- p.adjust( pA, "BH" )
                  #print( table( padjA < .1 ))
                  is_het =  padjA < threshold
                  is_het[is.na(is_het)] = FALSE
                  res.list[["padjA"]] <- padjA
                  res.list[["residual"]] <- cv2Endo - (coefficients(fit)[[1]] + coefficients(fit)[[2]]/meansEndo)
                  res.list[["x"]] <- meansEndo
                  res.list[["y"]] <- cv2Endo
                  
                  if(plot==TRUE){
                    par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                    plot( meansEndo, cv2Endo, log="xy", col=1+is_het,ylim=c(0.1,250), xlab='Means', ylab = expression("CV" ^ 2), pch=20, cex=0.3)
                    xg <- 10^seq( -3, 5, length.out=100 )
                    lines( xg, coefficients(fit)[1] + coefficients(fit)[2]/xg,lwd=2,col='#0000FFF0' )
                    try({
                        points( meansERCC, cv2ERCC, col="blue", pch=15, cex=0.6 )
                        # Plot quantile lines around the fit
                        df <- ncol(nCountsEndo) - 1
                        lines( xg, ( coefficients(fit)["a0"] + coefficients(fit)["a1tilde"]/xg ) * qchisq( .975, df ) / df, col="#0000FFF0", lwd=2, lty="dashed" )
                        lines( xg, ( coefficients(fit)["a0"] + coefficients(fit)["a1tilde"]/xg ) * qchisq( .025, df ) / df, col="#0000FFF0", lwd=2, lty="dashed" )
                        lines( xg, psia1thetaA/xg + coefficients(fit)["a0"] + minBiolDisp, lty="dashed", col="#FF0000F0", lwd=2 )
                    })
                    if(!is.null(fitB))
                    {
                        print(rbind(byERCC=coefficients(fit),byEndo=coefficients(fitB)))
                        lines( xg, coefficients(fitB)["a0"] + coefficients(fitB)["a1tilde"]/xg, col='#00FF00F0', lwd=2 )
                        legend('bottomleft',c('Endo. genes','Var. genes','ERCCs',"Fit","Fit(Endo)"),
                               pch=c(20,20,15,NA,NA),lty = c(NA,NA,NA,1,1),
                               col=c('black','red','blue', '#0000FFF0', '#00FF00F0'),cex=1.2)
                    }else
                    {
                        legend('bottomleft',c('Endo. genes','Var. genes','ERCCs',"Fit"),
                               pch=c(20,20,15,NA),lty = c(NA,NA,NA,1),
                               col=c('black','red','blue', '#0000FFF0'),cex=1.2)
                    }
                    title(expression("Variable genes by Mean ~ CV" ^ 2 * " (using ERCCs)"))
                  }
                }
                if(method=='fit' & fit_type=='log')
                {
                  LCountsEndo <- log10(nCountsEndo+1)
                  LmeansEndo <- rowMeans( LCountsEndo )
                  Lcv2Endo = rowVars(LCountsEndo)/LmeansEndo^2
                  is_het = (fit$opts$offset * coefficients(fit)["a"] *10^(-coefficients(fit)["k"]*LmeansEndo) < Lcv2Endo) &  LmeansEndo>fit$opts$minmean 
                  res.list[["residual"]] <- Lcv2Endo - (fit$opts$offset * coefficients(fit)["a"] *10^(-coefficients(fit)["k"]*LmeansEndo))
                  res.list[["x"]] <- LmeansEndo
                  res.list[["y"]] <- Lcv2Endo
                  
                  if(plot==TRUE){
                    par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                    plot( LmeansEndo, Lcv2Endo, log="y", col=1+is_het,pch=20,cex=0.3,
                         xlab='means',ylab=expression("CV" ^2),main=expression("CV"^2*" ~ Means (log transformed)"))
                    xg <- seq( 0, 5.5, length.out=100 )
                    lines( xg, fit$opts$offset * coefficients(fit)[1] *10^(-coefficients(fit)[2]*xg ),lwd=2,col='green' )
                    legend('bottomright',c('Endo. genes','Var. genes',"Fit"),
                           pch=c(20,20,NA),lty = c(NA,NA,1),col=c('black','red', 'blue'),cex=1.2)
                  }
                }
                if(method=='fit' & fit_type=='counts')
                {
                  meansEndo <- rowMeans( nCountsEndo )
                  varsEndo <- rowVars( nCountsEndo )
                  cv2Endo <- varsEndo/meansEndo^2
                  is_het = (coefficients(fit)[[1]] + coefficients(fit)[[2]]/meansEndo) < cv2Endo #&  meansEndo>2
                  res.list[["residual"]] <- cv2Endo - (coefficients(fit)[[1]] + coefficients(fit)[[2]]/meansEndo)
                  res.list[["x"]] <- meansEndo
                  res.list[["y"]] <- cv2Endo
                  
                  if(plot==TRUE){
                    par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                    plot( meansEndo, cv2Endo, log="xy", col=1+is_het,ylim=c(0.1,95), xlab='Means', ylab=expression("CV" ^ 2), pch=20, cex=0.3)
                    xg <- 10^seq( -3, 5, length.out=100 )
                    lines( xg, coefficients(fit)[1] + coefficients(fit)[2]/xg,lwd=2,col='#00FF00F0' )
                    legend('bottomright',c('Endo. genes','Var. genes',"Fit"),
                           pch=c(20,20,NA),lty = c(NA,NA,1),col=c('black','red', '#00FF00F0'),cex=1.2)   
                    title(expression("Variable genes by Mean ~ CV" ^ 2 * " (using endogenous genes)"))
                  }
                }
                if(method=='fit' & fit_type=='logvar')
                {
                  LCountsEndo <- log10(nCountsEndo+1)
                  LmeansEndo <- rowMeans( LCountsEndo )
                  LVarsEndo <- rowVars( LCountsEndo )
                  
                  xg = LmeansEndo
                  
                  Var_techEndo_logfit_loess =  predict(fit, LmeansEndo)
                  
                  minVar_Endo = min(LVarsEndo[LmeansEndo>2.5])
                  
                  if(any(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)){
                    idx = which(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)
                    Var_techEndo_logfit_loess[idx] = 0.6*minVar_Endo       
                  }      
                  
                  is_het = (Var_techEndo_logfit_loess < LVarsEndo) &  LmeansEndo>0.3
                  res.list[["residual"]] <- LVarsEndo - Var_techEndo_logfit_loess
                  res.list[["x"]] <- LmeansEndo
                  res.list[["y"]] <- LVarsEndo
                  
                  if(plot==TRUE){

                    par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                    plot( LmeansEndo, LVarsEndo, log="y", col=1+is_het,pch=20,cex=0.3,
                         xlab='means',ylab='Vars',main=expression("Vars ~ Means (log transformed)"))
                    xg <- seq( 0, 5.5, length.out=100 )
                    Var_techEndo_logfit_loess =  predict(fit, xg)
                    if(any(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)){
                      idx_1 = which(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)[1]
                      idx_end = length(Var_techEndo_logfit_loess)
                      Var_techEndo_logfit_loess[idx_1:idx_end] = 0.6*minVar_Endo       
                    }
                    lines( xg, Var_techEndo_logfit_loess,lwd=2,col='green' )
                    legend('bottomright',c('Endo. genes','Var. genes',"Fit"),
                           pch=c(20,20,NA),lty = c(NA,NA,1),col=c('black','red', 'green'),cex=1.2)
                  }
                }
                is_het[is.na(is_het)] = FALSE
                res.list[["is_het"]] <- is_het
                ##res.list[["geneID"]] <- rownames(nCountsEndo)
                ###is_het
                res.list
            })


#  “memo sort” method by B. Arman Aksoy (https://gist.github.com/armish/564a65ab874a770e2c26)
#' 
#' @param M (matrix) input data
memoSort <- function(M) {
    geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
    scoreCol <- function(x) {
        score <- 0;
        for(i in 1:length(x)) { 
            if(x[i]) {
                score <- score + 2^(length(x)-i);
            }
        }
        return(score);
    }
    scores <- apply(M[geneOrder, ], 2, scoreCol);
    sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
    return(M[geneOrder, sampleOrder]);
}

auto.colSet <- function(n=2,name="Set1"){
    suppressPackageStartupMessages(require("RColorBrewer"))
    if(n<=8){
         ret <- brewer.pal(max(n,3),name)[seq_len(n)]
    }else{
        ret <- colorRampPalette(brewer.pal(12,"Paired"))(n)
    }
    return(ret)
}

run.GSVA <- function(exp.dat,gset,gsva.method="gsva",bool.rnaseq=TRUE,ncores=4,min.sz=1,ssgsea.norm=T)
{
    suppressPackageStartupMessages(require("GSVA"))
    #suppressPackageStartupMessages(require("doParallel"))
    #registerDoParallel(cores = ncores)
    gsva.out <- gsva(exp.dat,gset,method=gsva.method,rnaseq=bool.rnaseq,mx.diff=T,parallel.sz=ncores,min.sz=min.sz,ssgsea.norm=ssgsea.norm)
    return(gsva.out)
}

calInfiltrationScore <- function(x,use.scale=F){
    suppressPackageStartupMessages(require("circlize"))
    suppressPackageStartupMessages(require("RColorBrewer"))
    x.scale <- t(apply(x,1,scale))
    colnames(x.scale) <- colnames(x)
    TSet <- c("T_cells","CD8_T_cells","T_helper_cells","Tcm","Tem","Th1_cells","Th2_cells","Th17_cells","TReg")
    OSet <- c(c("Macrophages","DC","B_cells","Eosinophils","Mast_cells","Neutrophils","NK_cells"),TSet)
    inf.score <- data.frame(TIS=apply(x.scale[TSet,],2,mean),
                            OIS=apply(x.scale[OSet,],2,mean),
                            CD8_Treg_Ratio=if(use.scale) apply(x.scale[c("CD8_T_cells","TReg"),],2,function(x){x[1]-x[2]}) else 
                                apply(x[c("CD8_T_cells","TReg"),],2,function(x){x[1]-x[2]}),
                            Th17_Th2_Ratio=if(use.scale) apply(x.scale[c("Th17_cells","Th2_cells"),],2,function(x){x[1]-x[2]}) else 
                                apply(x[c("Th17_cells","Th2_cells"),],2,function(x){x[1]-x[2]})
                            )
    inf.score$TIS.raw <- inf.score$TIS
    inf.score$OIS.raw <- inf.score$OIS
    q.tis <- quantile(inf.score$TIS)
    q.ois <- quantile(inf.score$OIS)
    print(q.tis)
    print(q.ois)
    #print(head(inf.score$TIS))
    #print(head(inf.score$OIS))
    for(i in 1:4){ inf.score$TIS[ q.tis[i] <= inf.score$TIS.raw & inf.score$TIS.raw <= q.tis[i+1] ] <- i }
    for(i in 1:4){ inf.score$OIS[ q.ois[i] <= inf.score$OIS.raw & inf.score$OIS.raw <= q.ois[i+1] ] <- i }
    #print(head(inf.score$TIS))
    #print(head(inf.score$OIS))
    #inf.score$TIS[inf.score$TIS > 3] <- 3
    #inf.score$TIS[inf.score$TIS < -3] <- -3
    #inf.score$OIS[inf.score$OIS > 3] <- 3
    #inf.score$OIS[inf.score$OIS < -3] <- -3
    i.col <- colorRamp2(1:4,colorRampPalette((brewer.pal(11,"BrBG")[c(7,8,9,10)]))(4), space="LAB")
    return(list(score=inf.score,
                col=list(TIS=i.col, OIS=i.col),
                x.scaled=x.scale
               ))
}

### code copy from https://github.com/mylesmharrison/colorRampPaletteAlpha
# addalpha()
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}
# colorRampPaletteAlpha()
colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
    # Create the color ramp normally
    cr <- colorRampPalette(colors, interpolate=interpolate)(n)
    # Find the alpha channel
    a <- col2rgb(colors, alpha=T)[4,]
    # Interpolate
    if (interpolate=='linear') { 
        l <- approx(a, n=n)
    } else { 
        l <- spline(a, n=n) 
    }
    l$y[l$y > 255] <- 255 # Clamp if spline is > 255
    cr <- addalpha(cr, l$y/255.0)
    return(cr)
}

getEntrezFromGO <- function(term,direct=T)
{ 
    if (require(org.Hs.eg.db)) {
        if(direct){ xxGO <- AnnotationDbi::as.list(org.Hs.egGO2EG) }
        else{ xxGO <- AnnotationDbi::as.list(org.Hs.egGO2ALLEGS) }
    } else { stop("Install org.Hs.eg.db package for retrieving gene lists from GO") }
    ret.EG <- unlist(xxGO[term])
    ret.EG
}

do.aov.1way <- function(dat.ana,xname,yname,out.prefix,ylab.txt="y"){
    aov.out <- aov(as.formula(sprintf("%s ~ %s",yname,xname)),data=dat.ana)
    aov.out.s <- summary(aov.out)
    aov.out.F <- unlist(aov.out.s[[1]][xname,c("F value","Pr(>F)")])
    aov.out.hsd <- TukeyHSD(aov.out)
    ## plot
    pdf(sprintf("%s.aov.pdf",out.prefix),width=6,height=4)
    p <- ggplot(dat.ana, aes_string(xname, yname)) + 
            geom_dotplot(binaxis='y', stackdir='center') + 
                stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                             geom="errorbar", color="red", width=0.2) +
                stat_summary(fun.y=mean, geom="point", color="red") + labs(y=ylab.txt) + 
                theme_bw(base_size = 12) + 
                theme(axis.text.x = element_text(angle = ifelse(length(table(dat.ana[,xname]))>3,30,0), hjust = 1), 
                      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    print(p)
    psig=as.numeric(apply(aov.out.hsd[[xname]][,2:3,drop=F],1,prod)>=0)+1
    op=par(mar=c(5,14,4,2))
    plot(aov.out.hsd,col=psig,yaxt="n")
    for(j in 1:length(psig)){ 
        axis(2,at=j,labels=rownames(aov.out.hsd[[xname]])[length(psig)-j+1], las=1,cex.axis=.8,col.axis=psig[length(psig)-j+1])
    }
    par(op)
    dev.off()
}

#' @importFrom ROCR prediction performance
#' @importFrom stats aggregate wilcox.test
## code from SC3
my.getAUC <- function(gene, labels,use.rank=T)
{
    suppressPackageStartupMessages(library("ROCR"))
    if(use.rank){
        score <- rank(gene)
    }else{
        score <- gene
    }
    # Get average score for each cluster
    ms <- aggregate(score ~ labels, FUN = mean)
    # Get cluster with highest average score
    posgroup <- ms[ms$score == max(ms$score), ]$labels
    # Return negatives if there is a tie for cluster with highest average score
    # (by definition this is not cluster specific)
    if(length(posgroup) > 1) {
        return (c(-1,-1,1))
    }
    # Create 1/0 vector of truths for predictions, cluster with highest
    # average score vs everything else
    truth <- as.numeric(labels == posgroup)
    #Make predictions & get auc using RCOR package.
    pred <- prediction(score,truth)
    val <- unlist(performance(pred,"auc")@y.values)
    pval <- suppressWarnings(wilcox.test(score[truth == 1],
                                         score[truth == 0])$p.value)
    return(c(val,posgroup,pval))
}

expressedFraction <- function(exp.bin,group,ncores=8){
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("doParallel"))
    registerDoParallel(cores = ncores)
    out.res <- ldply(rownames(exp.bin),function(v){
                .res <- aggregate(exp.bin[v,],by=list(group),FUN=function(x){ sum(x==1)/length(x) })
                structure(.res[,2],names=.res[,1])
			},.progress = "none",.parallel=T)
    rownames(out.res) <- rownames(exp.bin)
    colnames(out.res) <- sprintf("HiFrac.%s",colnames(out.res))
    out.res <- as.matrix(out.res)
    return(out.res)
}

expressedFraction.HiExpressorMean <- function(exp.bin,exp.tpm,group,ncores=8){
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("doParallel"))
    registerDoParallel(cores = ncores)
    exp.bin[exp.bin<1] <- 0
    .exp <- exp.bin*exp.tpm
    ###out.res <- ldply(rownames(.exp)[1000:1010],function(v){
    out.res <- ldply(rownames(.exp),function(v){
                .n <- aggregate(exp.bin[v,],by=list(group),FUN=function(x){ sum(x==1) })
                .res <- aggregate(.exp[v,],by=list(group),FUN=function(x){ sum(x) })
                structure(.res[,2]/.n[,2],names=.res[,1])
			},.progress = "none",.parallel=T)
    rownames(out.res) <- rownames(exp.bin)
    colnames(out.res) <- sprintf("AvgHiExpr.%s",colnames(out.res))
    out.res <- as.matrix(out.res)
    return(out.res)
}


#' get marker genes given data and cluster assignment
#' @importFrom plyr
#' @importFrom doParallel
#' @param dat.to.test row genes, column samples
#' @param clust
#' @param out.prefix
#' @param n.cores
my.clusterMarkerGene <- function(dat.to.test,clust,out.prefix,n.cores=NULL,ann.col=NULL,
                                 clust.col=NULL,original.labels=F,
                                 sampleType=NULL,sampleTypeColSet=NULL,
                                 exp.bin=NULL,gid.mapping=NULL,exp.tpm=NULL,F.FDR.THRESHOLD=0.05,verbose=F,minCell=5)
{
    suppressPackageStartupMessages(require("plyr"))
    suppressPackageStartupMessages(require("dplyr"))
    suppressPackageStartupMessages(require("doParallel"))
    suppressPackageStartupMessages(require("factoextra"))
    ##if(length(unique(clust))<2 || is.null(dat.to.test) || !all(table(clust) > 3)){
    if(length(unique(clust))<2 || is.null(dat.to.test) || !all(table(clust) > 1)){
        cat("WARN: clusters<2 or no data.to.test provided or not all clusters have more than 3 samples\n")
        return(NULL)
    }
    #### marker genes by AUC
    registerDoParallel(cores = n.cores)

    if(!is.null(exp.bin)){
        .f.gid <- intersect(rownames(exp.bin),rownames(dat.to.test))
        exp.frac <- expressedFraction(exp.bin[.f.gid,,drop=F],clust,n.cores)
        print(str(exp.frac))
        print(head(exp.frac))
#        f.bin <- apply(exp.frac,1,function(x){ nE <- sum(x>0.05); nNE <- sum(x<0.95); return(nE > 0 && nNE > 0) })
#        dat.to.test <- dat.to.test[f.bin,]
        f.notAllZero <- apply(dat.to.test,1,function(x){ nE <- sum(x>0); return( nE > minCell & nE/length(x) > 0.01 ) })
        dat.to.test <- dat.to.test[f.notAllZero,]
    }else{
        f.notAllZero <- apply(dat.to.test,1,function(x){ nE <- sum(x>0); return( nE > minCell & nE/length(x) > 0.01 ) })
        dat.to.test <- dat.to.test[f.notAllZero,]
    }

    .gene.table <- ldply(rownames(dat.to.test),function(x){
        my.getAUC(dat.to.test[x,],clust,use.rank = F)
    },.progress = "none",.parallel=T)
    colnames(.gene.table) <- c("AUC","cluster","score.p.value")
    print("str(.gene.table)")
    print(str(.gene.table))
    if(is.character(.gene.table$AUC)){ .gene.table$AUC <- as.numeric(.gene.table$AUC) }
    if(is.character(.gene.table$score.p.value)){ .gene.table$score.p.value <- as.numeric(.gene.table$score.p.value) }
    if(is.numeric(.gene.table$cluster)){ .gene.table$cluster <- sprintf("C%s",.gene.table$cluster) }
    ##.gene.table$score.q.value <- p.adjust(.gene.table$score.p.value,method = "BH")
                                    ####geneSymbol=if(original.labels) rownames(dat.to.test) else entrezToXXX(rownames(dat.to.test)),
    .geneID2GeneName <- function(x){
        if(!original.labels){
            if(!is.null(gid.mapping)){
                .gene.name <- gid.mapping[x]
            }else{
                .gene.name <- entrezToXXX(x)
            }
        }else{
            return(x)
        }
    }
    .gene.table$score.q.value <- 1
    .gene.table <- cbind(data.frame(geneID=rownames(dat.to.test),
                                    geneSymbol=.geneID2GeneName(rownames(dat.to.test)),
                                    stringsAsFactors = F),
                         .gene.table)
    
    test.gene.table.11 <<- .gene.table
   
    #### diff genes
    .aov.res <- runMultiGroupSpecificGeneTest(dat.to.test,
                                              grps = if(is.numeric(clust)) sprintf("C%s",clust) else clust,
                                              out.prefix,mod=NULL,
                                              FDR.THRESHOLD=0.05,F.FDR.THRESHOLD=F.FDR.THRESHOLD,FC.THRESHOLD=1,
                                              verbose=F, n.cores=n.cores,gid.mapping=gid.mapping)

    ####.gene.table <- left_join(x = .gene.table,y = .aov.res$aov.out)
    ##.gene.table <- inner_join(x = .gene.table,y = .aov.res$aov.out.sig,by="geneID")
    if(verbose){
        .gene.table <- inner_join(x = .gene.table,y = .aov.res$aov.out)
    }else{
        .gene.table <- inner_join(x = .gene.table,y = .aov.res$aov.out.sig)
    }
    ### average expression of each cluster
    ###avg.exp <- ldply(rownames(dat.to.test[as.character(.gene.table$geneID),,drop=F]),function(v){
    avg.exp <- ldply(as.character(.gene.table$geneID),function(v){
                .res <- aggregate(dat.to.test[v,],by=list(clust),FUN=mean)
                .res.2 <- aggregate(dat.to.test[v,],by=list(clust),FUN=sd)
                structure(c(.res[,2],.res.2[,2]),names=c(sprintf("avg.%s",.res[,1]),sprintf("sd.%s",.res[,1])))
			},.progress = "none",.parallel=T)
    ##rownames(avg.exp) <- rownames(dat.to.test[as.character(.gene.table$geneID),,drop=F])
    rownames(avg.exp) <- as.character(.gene.table$geneID)
    ####colnames(avg.exp) <- sprintf("avg.%s",colnames(avg.exp))
    avg.exp.df <- data.frame(geneID=rownames(avg.exp),
                             geneSymbol=.geneID2GeneName(rownames(avg.exp)),
                             stringsAsFactors = F)
    avg.exp.df <- cbind(avg.exp.df,avg.exp)
    .gene.table <- inner_join(x=.gene.table,y=avg.exp.df)
    ### binarized expression fraction
    if(!is.null(exp.bin)){
        #.f.gid <- intersect(rownames(exp.bin),as.character(.gene.table$geneID))
        #exp.frac <- expressedFraction(exp.bin[.f.gid,,drop=F],clust,n.cores)
        exp.frac.df <- data.frame(geneID=rownames(exp.frac), 
                                  geneSymbol=.geneID2GeneName(rownames(exp.frac)),
                                  stringsAsFactors = F)
        exp.frac.df <- cbind(exp.frac.df,exp.frac)
        .gene.table <- inner_join(x = .gene.table,y = exp.frac.df)
    }
    if(!is.null(exp.tpm)){
        avg.exp.hiExpr <- expressedFraction.HiExpressorMean(exp.bin,exp.tpm,clust,ncores=8)
        avg.exp.hiExpr.df <- data.frame(geneID=rownames(avg.exp.hiExpr),
                                        geneSymbol=.geneID2GeneName(rownames(avg.exp.hiExpr)),
                                        stringsAsFactors = F)
        avg.exp.hiExpr.df <- cbind(avg.exp.hiExpr.df,avg.exp.hiExpr)
        .gene.table <- inner_join(x = .gene.table,y = avg.exp.hiExpr.df)
    }
    rownames(.gene.table) <- .gene.table$geneID
    if(verbose){
        .gene.table <- .gene.table[order(.gene.table$F,decreasing = T),]
        write.table(.gene.table, file = sprintf("%s.geneTable.all.txt",out.prefix),
                    row.names = F,quote = F,sep = "\t")
    }
    f.isSig <- intersect(rownames(.aov.res$aov.out.sig),rownames(.gene.table))
    .gene.table <- .gene.table[f.isSig,]
    .gene.table$score.q.value <- p.adjust(.gene.table$score.p.value,method = "BH")
    order.gene <- order(.gene.table$cluster,-.gene.table$AUC)
    .gene.table <- .gene.table[order.gene,,drop=F]

    #test.gene.table <<- .gene.table
    #test.dat.to.test <<- dat.to.test
    #test.clust <<- clust
    #test.aov.res <<- .aov.res

    #rownames(.gene.table) <- .gene.table$geneID
    ### save .txt file
    write.table(.gene.table, file = sprintf("%s.markerGene.all.txt",out.prefix),
                row.names = F,quote = F,sep = "\t")
    write.table(subset(.gene.table,score.q.value<0.01),
                file = sprintf("%s.markerGene.q01.txt",out.prefix),
                row.names = F,quote = F,sep = "\t")
    ### plot heatmap
    gene.for.heatmap <- ldply(unique(subset(.gene.table,score.q.value<0.01 & AUC>0.7 & cluster>0,select="cluster",drop=T)),
                              function(x){ head(subset(.gene.table,cluster==x),n=10) })
    names(clust) <- colnames(dat.to.test)
    #### cluster color
    if(is.null(clust.col)){
        .colSet <- auto.colSet(length(unique(sort(clust[clust>0]))))
        names(.colSet) <- sprintf("C%s",unique(sort(clust[clust>0])))
    }else{
        .colSet <- clust.col
        if(is.numeric(clust)){
            names(.colSet) <- sprintf("C%s",names(.colSet)) 
        }
    }

    nSamples <- sum(clust>0)
    names(ann.col) <- colnames(dat.to.test)
    s.heatmap <- names(sort(clust[clust>0]))
    ### sampleType color
    if(!is.null(sampleTypeColSet)){
        annC <- sampleType
        names(annC) <- colnames(dat.to.test)
        annCSet <- sampleTypeColSet
    }else{
        annC <- ann.col
        names(annC) <- colnames(dat.to.test)
        annCSet=structure(unique(annC),names=unique(annC))
    }

    print(.colSet)
    print(str(annC))
    print(str(annCSet))
    #print(str(annC[s.heatmap]))
    runHierarchicalClusteringAnalysis(dat.to.test[gene.for.heatmap$geneID,s.heatmap,drop=F],
                                      mytitle = "", pdf.width=if(nSamples>1000) 14 else 12,pdf.height=8,
                                      sprintf("%s.markerGene",out.prefix),
                                      do.clustering.col=F, 
                                      do.clustering.row=F, 
                                      sampleType=annC[s.heatmap],
                                      colSet=annCSet,
                                      ann.extra.df=data.frame(cluster=if(is.numeric(clust)) sprintf("C%s",clust[s.heatmap]) else clust[s.heatmap]),
                                      ann.extra.df.col=list(cluster=.colSet),
                                      row.names.original=original.labels, 
                                      ann.bar.height=0.8,
                                      gid.mapping=gid.mapping,
                                      complexHeatmap.use=TRUE,verbose=FALSE,do.scale=T)

    return(list(gene.table=.gene.table,aov.res=.aov.res))
}

#' clustering validation given data and cluster assignment
#' @importFrom plyr
#' @importFrom doParallel
#' @param dat.to.test row genes, column samples
#' @param clust
#' @param out.prefix
#' @param n.cores
my.clusteringValidation <- function(dat.to.test,clust,out.prefix,n.cores=NULL,dist.obj=NULL)
{
    suppressPackageStartupMessages(require("plyr"))
    suppressPackageStartupMessages(require("doParallel"))
    suppressPackageStartupMessages(require("factoextra"))
    suppressPackageStartupMessages(require("cluster"))
    suppressPackageStartupMessages(require("fpc"))
    
    names(clust) <- colnames(dat.to.test)
    clust.noZero <- clust[clust>0]
    
    if(length(unique(clust.noZero))<2){
        return(NULL)
    }

    if(is.null(dist.obj)){
        dist.obj <- dist(t(dat.to.test[,names(clust.noZero)]))
    }
    #print("TEST:XXX")
    #print(str(clust.noZero))
    #print(str(dist.obj))
    ### silhouette by cluster package
    if(is.numeric(clust)){
        sil <- silhouette(clust.noZero, dist.obj)
    }else{
        sil <- silhouette(as.numeric(factor(clust.noZero)), dist.obj)
    }
    #print("TEST:YY")
    ###print(head(sil[,1:3]))
    ###print(str(sil))
    pdf(sprintf("%s.silhouette.pdf",out.prefix),width=8,height=8)
    plot(sil, main ="Silhouette plot")
    dev.off()
    ### other stats
    loginfo("...... (begin) cal some cluster.stats ...\n")
    if(is.numeric(clust)){
        c.stats <- cluster.stats(dist.obj, clust.noZero,sepindex=F,wgap=F)
    }else{
        c.stats <- cluster.stats(dist.obj, as.numeric(factor(clust.noZero)),sepindex=F,wgap=F)
    }
    loginfo("...... (end) cal some cluster.stats ...\n")
    return(list(c.stats=c.stats,sil=sil))
}

processInput <- function(designFile,cellTypeColorFile,inputFile,args.notFilter,geneFile,args.center,args.log,args.norm.exprs=F,args.measure=NULL)
{
    suppressPackageStartupMessages(library("R.utils"))
    suppressPackageStartupMessages(library("scater"))
    ### designFile
    ##myDesign<-read.table(designFile,header=T,check.names=F,stringsAsFactors = F)
    myDesign<-read.delim(designFile,header=T,check.names=F,stringsAsFactors = F)
    rownames(myDesign) <- myDesign$sample
    ### cellTypeColorFile
    ####sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
    if(!is.null(cellTypeColorFile) && file.exists(cellTypeColorFile)){
        sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
        sampleTypeColor <- sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[,"sampleType"]))]
    }else{
        sampleTypeColor <- auto.colSet(n=length(unique(myDesign[,"sampleType"])),name="Paired")
        names(sampleTypeColor) <- unique(myDesign[,"sampleType"])
    }
    if("patient" %in% colnames(myDesign)){
        patient.col.list <- patientColorListFromMyDesign(myDesign)
    }else{
        patient.col.list <- NULL
    }

    ### input expression file
    if(grepl("\\.scran\\..+RData$",inputFile,perl=T) || grepl("\\.scran\\.RData$",inputFile,perl=T)){
        lenv <- loadToEnv(inputFile)
        if(!is.null(args.measure) && args.measure=="tpm"){
            Y <- tpm(lenv[["sce.norm"]])
        }else if(!is.null(args.measure) && args.measure=="counts"){
            Y <- counts(lenv[["sce.norm"]])
        }else if(args.norm.exprs || (!is.null(args.measure) && args.measure=="norm_exprs")){
            Y <- norm_exprs(lenv[["sce.norm"]])
            cat("range of norm_exprs Y:\n")
            print(range(Y))
            args.notFilter <- T
            args.log <- F
            args.center <- F
        }else if(!args.norm.exprs || (!is.null(args.measure) && args.measure=="centered_norm_exprs")){
            Y <- assay(lenv[["sce.norm"]],"centered_norm_exprs")
            args.notFilter <- T
            args.log <- F
            args.center <- F
        }else{
            ##Y <- exprs(lenv[["sce.norm"]])
            Y <- assay(lenv[["sce.norm"]],"centered_norm_exprs")
            args.notFilter <- T
            args.log <- F
            args.center <- F
        }
        ##g.GNAME <- fData(lenv[["sce.norm"]])[,"geneSymbol"]
        g.GNAME <- rowData(lenv[["sce.norm"]])[,"geneSymbol"]
        names(g.GNAME) <- rownames(Y)
    }else if(grepl("RData$",inputFile,perl=T)){
        lenv <- loadToEnv(inputFile)
        if("Y" %in% names(lenv)){
            Y <- lenv[["Y"]]
        }else if("exp.bin" %in% names(lenv)){
            Y <- lenv[["exp.bin"]]
        }else{
            stop("can not find expression data. processInput()")
        }
        g.GNAME <- lenv[["g.GNAME"]]
        if(is.null(g.GNAME)){
            g.GNAME <- entrezToXXX(rownames(Y))
            names(g.GNAME) <- rownames(Y)
        }
        ###obj.scdn <- lenv[["obj.scdn"]]
        ###Y <- obj.scdn@normalized.endo
    }else{
        in.table <- read.table(inputFile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
        rownames(in.table) <- in.table[,1]
        Y <- as.matrix(in.table[,c(-1,-2)])
        #g.GNAME <- entrezToXXX(rownames(Y))
        g.GNAME <- in.table[,2]
        names(g.GNAME) <- rownames(Y)
    }
    f.gname.na <- is.na(g.GNAME)
    g.GNAME[f.gname.na] <- names(g.GNAME)[f.gname.na]
    f.gname.dup <- duplicated(g.GNAME)
    g.GNAME[f.gname.dup] <- names(g.GNAME)[f.gname.dup]

    sname <- intersect(rownames(myDesign),colnames(Y))
    myDesign <- myDesign[sname,,drop=F]
    Y <- Y[,sname,drop=F]
    #### other options
    if(!args.notFilter){
        f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
        Y <- Y[f,]
    }
    if(!is.null(geneFile) && file.exists(geneFile)){
        geneTable <- read.table(geneFile,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
        Y <- Y[unique(as.character(geneTable$geneID)),]
    }

    if(args.log) { Y <- log2(Y+1) }
    if(args.center){
        Y.new <- c()
        for(pp in unique(myDesign$patient)){
            Y.block <- t(scale(t(Y[,subset(myDesign,patient==pp,select="sample",drop=T)]),center = T,scale = F))
            Y.new <- cbind(Y.new,Y.block)
            print(apply(Y.block[1:4,],1,mean))
            print(apply(Y.block[1:4,],1,sd))
        }
        Y <- Y.new
        cat("test: all.equal(colnames(Y),rownames(myDesign))\n")
        print(all.equal(colnames(Y),rownames(myDesign)))
        Y <- Y[,sname,drop=F]
        cat("test: all.equal(colnames(Y),rownames(myDesign))\n")
        print(all.equal(colnames(Y),rownames(myDesign)))
    }
    #if(args$onlyY){
    #    save(Y,g.GNAME,file=sprintf("%s/%s.Y.RData",out.dir,sample.id))
    #    loginfo("Y data saved.")
    #    q()
    #}
    #myDesign <<- myDesign
    #sampleTypeColor <<- sampleTypeColor
    #patient.col.list <<- patient.col.list
    #Y <<- Y
    #g.GNAME <<- g.GNAME
    return(list(myDesign=myDesign,
                sampleTypeColor=sampleTypeColor,
                patient.col.list=patient.col.list,
                Y=Y,
                g.GNAME=g.GNAME))
}

#' dynamicTreeCut
#' @importFrom dynamicTreeCut
#' @importFrom dendextend
#' @param dat.to.sort row genes, column samples
run.quickClustBydynamicTreeCut <- function(dat.to.sort,k=4,useDyn=T,col.sort=F,clustering.distance = "spearman", clustering.method = "complete")
{
    suppressPackageStartupMessages(library("dynamicTreeCut"))
    suppressPackageStartupMessages(library("dendextend"))
    my.dist <- as.dist(1-cor(dat.to.sort,method=clustering.distance))
    my.hclust <- hclust(my.dist,method=clustering.method)
    if(useDyn){
        my.hclust.cutDyn <- cutreeDynamic(my.hclust, distM=as.matrix(my.dist), verbose=0)
        my.clusters <- unname(my.hclust.cutDyn)
    }else{
        ### cut tree
        dend.cutree <- cutree(my.hclust, k=k, order_clusters_as_data = T)
        my.clusters <- dend.cutree ### dropped if kk is scalar
    }
    ## save result to txt file
    list(my.clusters=structure(my.clusters,names=colnames(dat.to.sort)),my.hclust=my.hclust)
}

### zero.as.low: if True, zero counts as low expression (0); else as "-1"
binarizedExp <- function(x,ofile=NULL,G=NULL,e.TH=NULL,e.name="Exp",verbose=F, draw.CI=T, zero.as.low=T,...)
{
  require(mclust)
  ### bin.Exp == -1: drop out
  o.df <- data.frame(sample=names(x),bin.Exp=-1,stringsAsFactors = F)
  rownames(o.df) <- o.df$sample
  #f<-is.finite(x) & x>0
  f<-is.finite(x)
  if(sum(f)<3){
      colnames(o.df) <- c("sample",e.name)
      if(verbose){
        return(list(x_mix=NULL,o.df=o.df))
      }else{
        return(structure(o.df[,e.name],names=rownames(o.df)))
      }
  }

  x<-x[f]
  x_mix<-densityMclust(x,G=G,modelNames=c("E","V"))
  x_mix_summary<-summary(x_mix)
  
  if(verbose && !is.null(ofile)){
	  print(x_mix_summary)
      if(grepl(".png$",ofile,perl = T)){
          png(ofile,width=800,height=600)
	      old_par<-par(no.readonly=T)
	      layout(matrix(c(1,1,1,1,1,1,2,3,4), 3, 3, byrow = TRUE))
	      a_par<-par(cex.axis=2,cex.lab=2,cex.main=1.8,mar=c(5,6,4,2)+0.1)
      }else{
          pdf(ofile,width=8,height=6)
      }
	  plot(x_mix,what="density",data=x,breaks=50,col="darkgreen",lwd=2,main="",...)
	  abline(v=x_mix_summary$mean,lty=2)
	  if(!is.null(e.TH)){
		abline(v=e.TH,lty=2,col="red")
	  }
	  
	  for(i in 1:x_mix_summary$G)
	  {
		i_mean<-x_mix_summary$mean[i]
		i_sd<-sqrt(x_mix_summary$variance[i])
		i_pro<-x_mix_summary$pro[i]
		#i_sd<-RC_mix_summary$variance[i]
		d<-qnorm(c(0.0013,0.9987),i_mean,i_sd)
		e<-i_pro*dnorm(i_mean,i_mean,i_sd)
		lines(seq(d[1],d[2],by=0.01),i_pro*dnorm(seq(d[1],d[2],by=0.01),i_mean,i_sd),col="orange",lwd=2)
        if(draw.CI){ rect(d[1],0,d[2],e+0.02,col=NA,border="blue") }
		#rect(d[1],0,d[2],e+0.02,col=rgb(0,0,0.8,0.2),border=NA)
	  }
	  plot(x_mix,data=x,breaks=20,col="darkgreen",lwd=2,what="BIC")
	  densityMclust.diagnostic(x_mix,type = "cdf",cex.lab=1.5)
	  densityMclust.diagnostic(x_mix,type = "qq")
	  dev.off()
  }
  
  o.df[names(x_mix$classification),"bin.Exp"] <- x_mix$classification 
  colnames(o.df) <- c("sample",e.name)
  if(!is.null(G) && x_mix_summary$G==3){
	  ### determin which classes 'not-expressed' and which classes 'expressed'	
	  i_mean<-x_mix_summary$mean
	  i_sd<-sqrt(x_mix_summary$variance)
	  ci95.1 <- qnorm(c(0.0013,0.9987),i_mean[1],i_sd[1])
	  ci95.1.overlap <- pnorm(ci95.1[2],i_mean[2],i_sd[2])
	  if(ci95.1.overlap>0.3333){
		  C.min <-  2
	  }else{
		  C.min <- 1
	  }
	  ci95.2 <- qnorm(c(0.0013,0.9987),i_mean[2],i_sd[2])
	  ci95.3 <- qnorm(c(0.0013,0.9987),i_mean[3],i_sd[3])
	  ci95.3.overlap <- pnorm(ci95.3[1],i_mean[2],i_sd[2],lower.tail = F)
	  if(ci95.3.overlap>0.3333){
		  C.max <- 2
	  }else{
		  C.max <- 3
	  }
#	  ci95.2.overlap <- pnorm(ci95.2[2],i_mean[3],i_sd[3])
#	  if(ci95.2.overlap>0.3333){
#	  	  C.max <- 2
#	  }else{
#	  	  C.max <- 3
#	  }
      if(zero.as.low){
	    f.low <- o.df[,e.name] <= C.min
      }else{
	    f.low <- o.df[,e.name] <= C.min & o.df[,e.name] != -1
      }
	  f.hi <-  o.df[,e.name] >= C.max
	  f.mid <- (o.df[,e.name] > C.min) & (o.df[,e.name] < C.max)
	  o.df[f.low,e.name] <- 0
	  o.df[f.mid,e.name] <- 0.5
	  o.df[f.hi,e.name] <- 1
	  #f.G <- o.df[,e.name] < C.max
	  #o.df[f.G,e.name] <- 0
	  #o.df[!f.G,e.name] <- 1
  }else if(x_mix_summary$G==2){
      if(zero.as.low){
	    f.low <- o.df[,e.name] <= 1
      }else{
	    f.low <- o.df[,e.name] <= 1 & o.df[,e.name] != -1
      }
	  f.hi <-  o.df[,e.name] == 2
	  o.df[f.low,e.name] <- 0
	  o.df[f.hi,e.name] <- 1
  }
  if(verbose){
	return(list(x_mix=x_mix,o.df=o.df))
  }else{
	return(structure(o.df[,e.name],names=rownames(o.df)))
  }
}

auto.point.size <- function(n){
    if(n<=100){
        return(1.2)
    }else if(n>=5000){
        return(0.6)
    }else{
        return(-0.6*n/4900+1.212002)
    }
}


#' @param d row samples, must have rowname
plot.mark.onTSNE <- function(x,mark.name,d,out.file,x.dim=c(1,2),colSet=NULL)
{
    suppressPackageStartupMessages(require("RColorBrewer"))
    ss <- intersect(rownames(x),rownames(d))
    x <- x[ss,,drop=F]
    d <- d[ss,,drop=F]
    if(class(d[,mark.name])=="character" || class(d[,mark.name])=="factor" || class(d[,mark.name])=="logical"){
        nn <- sort(unique(as.character((d[,mark.name]))))
        if(is.null(colSet)){
            colSet <- structure(auto.colSet(n = length(nn),name = "Dark2"),
                                names=nn)
        }
        print("XXXX")
        print(nn)
        print(colSet)
        plot.tsne.points(x,out.file,tsne.col.points=colSet[d[rownames(x),mark.name]],
                         col.tsne.legend=colSet,
                         tsne.legend=names(colSet),pch=20,nclusters=length(nn),peak=NULL,main=mark.name,x.dim=x.dim)
    }else if(class(d[,mark.name])=="integer" || class(d[,mark.name])=="numeric"){
        x <- x[,x.dim,drop=F]
        s.f <- order(d[,mark.name],decreasing=F)
        Y.level <- pretty(d[,mark.name],n=10)
        ppalette <- brewer.pal(9,"YlGnBu")
        Y.color <- as.character(cut(d[s.f,mark.name],
                                      breaks=quantile(c(Y.level[1],Y.level[length(Y.level)]), seq(0,1,0.01)),
                                      labels=addalpha(colorRampPalette(ppalette)(100),alpha = 1)))
        pdf(out.file,width=9,height=8)
        par(mar=c(5,5,4,9),cex.lab=1.8,cex.main=2.5,cex.axis=1.5)
        plot(x[s.f,1],x[s.f,2],
             t='n',pch=20,col="lightgray", main=sprintf("%s",mark.name),
             xlab=sprintf("Dim%d",x.dim[1]),
             ylab=sprintf("Dim%d",x.dim[2]))
        points(x[s.f,1],x[s.f,2],col=Y.color,pch=20,cex=auto.point.size(nrow(x)))
        image.plot(zlim=c(Y.level[1],Y.level[length(Y.level)]),
                   legend.only=TRUE,col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=8.0)
        dev.off()
    }
}

plot.gene.onTSNE <- function(Y,dat.map,gene.to.show,out.prefix,args.scale=F){
    #gene.to.show <- read.table(args.geneOnPCA,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
    #gene.to.show <- structure(as.character(gene.to.show[,1]),names=gene.to.show[,2])
    #dir.create(sprintf("%s.perGene.PCA",out.prefix),showWarnings = F,recursive = T)
    dir.create(sprintf("%s.perGene.tSNE",out.prefix),showWarnings = F,recursive = T)
    for(i in seq_along(gene.to.show)){
        gid <- gene.to.show[i]
        gname <- names(gene.to.show)[i]
        if(!(gid %in% rownames(Y))){ 
            cat(sprintf("Gene not in data: %s (%s)\n",gid, gname))
            next
        }
        i.Y <- Y[gid,]
        s.f <- order(i.Y,decreasing=F)
        Y.level <- pretty(Y[gid,],n=10)
        ppalette <- brewer.pal(9,"YlOrRd")
#        if(args.scale){
#            i.Y <- scale(i.Y)
#            #Y.level <- pretty(i.Y,n=10)
#            Y.level <- c(-2.5,2.5)
#            ppalette <- rev(brewer.pal(9,"RdBu"))
#            print(gname)
#            print(sprintf("mean: %4.4f, sd: %4.4f", mean(i.Y), sd(i.Y)))
#            print(summary(i.Y))
#            i.Y[i.Y>=2.5] <- 2.5
#            i.Y[i.Y<=-2.5] <- -2.5
#        }
#        gid.color <- as.character(cut(i.Y[s.f],
#                                      breaks=quantile(c(Y.level[1],Y.level[length(Y.level)]), seq(0,1,0.01)),
#                                      labels=addalpha(colorRampPalette(ppalette)(100),alpha = 1)))
#        ####### PCA #######
#        pdf(sprintf("%s.perGene.PCA/PCA.%s.pdf",out.prefix,gname),width = 9,height = 8)
#        par(mar=c(5,5,4,8),cex.lab=1.5,cex.main=3.0,cex.axis=1.5)
#        plot(pca.res.sampleType$ind$coord[,"Dim.1"],pca.res.sampleType$ind$coord[,"Dim.2"],
#             t='n',pch=16,col="lightgray", main=sprintf("%s",gname),xlab="Dim1",ylab="Dim2")
#        points(pca.res.sampleType$ind$coord[s.f,"Dim.1"],pca.res.sampleType$ind$coord[s.f,"Dim.2"],col=gid.color,pch=16)
#        image.plot(zlim=c(Y.level[1],Y.level[length(Y.level)]),legend.only=TRUE,
#                   col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=5.0)
#        dev.off()
        ####### tSNE #######
        slist.tsneUsed <- rownames(dat.map)
        i.Y <- Y[gid,slist.tsneUsed]
        s.f <- order(i.Y,decreasing=F)
        if(args.scale){
            i.Y <- scale(i.Y)
            print(gname)
            print(sprintf("mean: %4.4f, sd: %4.4f", mean(i.Y), sd(i.Y)))
            print(summary(i.Y))
            i.Y[i.Y>=2.5] <- 2.5
            i.Y[i.Y<=-2.5] <- -2.5
            #Y.level <- pretty(i.Y,n=10)
            Y.level <- c(-2.5,2.5)
            ppalette <- rev(brewer.pal(9,"RdBu"))
        }
        gid.color <- as.character(cut(i.Y[s.f],
                                      breaks=quantile(c(Y.level[1],Y.level[length(Y.level)]), seq(0,1,0.01)),
                                      labels=addalpha(colorRampPalette(ppalette)(100),alpha = 1)))
        
        pdf(sprintf("%s.perGene.tSNE/tSNE.%s.pdf",out.prefix,gname),width = 9,height = 8)
        par(mar=c(5,5,4,8),cex.lab=1.5,cex.main=3.0,cex.axis=1.5)
        plot(dat.map,
             t='n',pch=16,col="lightgray", main=sprintf("%s",gname),xlab="Dim1",ylab="Dim2")
        points(dat.map[s.f,],col=gid.color,pch=16)
        image.plot(zlim=c(Y.level[1],Y.level[length(Y.level)]),legend.only=TRUE,
                   col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=5.0)
        dev.off()
    }
}


#' @param Y row samples, column genes
my.rmBatchEffect  <-  function(Y, batch, keep.intercept=T)
{
  #X <- model.matrix(~batch)
  batch <- as.factor(batch)
  contrasts(batch) <- contr.sum(levels(batch))
  X <- model.matrix(~batch)[, , drop = FALSE]
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  if(keep.intercept){
    Y.corrected <- Y-X[,-1,drop=F] %*% beta[-1,,drop=F]
  }else{
    Y.corrected <- Y-X %*% beta
  }
  return(Y.corrected)
}

#' @param Y row genes, column samples
#' @param d row samples, have column "patient" and "sample"
my.centerData <- function(Y, d, do.scale=F)
{
    Y.new <- c()
    for(pp in unique(d$patient)){
        Y.block <- t(scale(t(Y[,subset(d,patient==pp,select="sample",drop=T)]),center = T,scale = do.scale))
        Y.new <- cbind(Y.new,Y.block)
        print(apply(Y.block[1:4,],1,mean))
        print(apply(Y.block[1:4,],1,sd))
    }
    return(Y.new)
}

#' @param dat.to.plot row samples, column variables
plot.venn.as.heatmap <- function(dat.to.plot,out.prefix,do.sort=T,link.list=NULL,colSet=NULL)
{
    require("grid")
    require("gridBase")
    require("ComplexHeatmap")
    require("RColorBrewer")
    pdf(sprintf("%s.venn.heatmap.pdf",out.prefix),width=5,height=8)
    par(mar=c(2,4,4,2))
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    ha.link <- NULL
    if(!is.null(link.list)){
        ann.link.at <- match(link.list,rownames(dat.to.plot))
        ann.link.label <- link.list
        ha.link <- rowAnnotation(link = row_anno_link(at = ann.link.at, labels = ann.link.label,labels_gp=gpar(fontsize=8)), 
                                 width = unit(1, "cm") + max_text_width(ann.link.label))
    }
    m <- ncol(dat.to.plot)
    n <- nrow(dat.to.plot)
    for(i in seq_len(m)){
        dat.to.plot[,i] <- as.numeric(dat.to.plot[,i])
    }
    if(do.sort){
        hclust.res <- hclust(dist(dat.to.plot))
        dat.to.plot <- dat.to.plot[hclust.res$order,]
    }
    print(str(dat.to.plot))
    print(head(dat.to.plot))
    if(is.null(colSet)){
        colSet <- structure(c("gray","brown1"), names = c("0", "1"))
    }
    ht <- Heatmap(dat.to.plot,"Called",
                col=colSet,  
                column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"),
                column_names_gp = gpar(fontsize = 12*28/max(m,32)),row_names_gp = gpar(fontsize = 10*28/max(n,32)),
                show_row_names=T,
                show_heatmap_legend = T, row_names_max_width = unit(10,"cm"),
                cluster_columns = F,
                cluster_rows = F)
    if(!is.null(ha.link)){
        ComplexHeatmap::draw(ht+ha.link, newpage= FALSE)
    }else{
        ComplexHeatmap::draw(ht, newpage= FALSE)
    }
    dev.off()
}

cal.signatureScore.Fred <- function(obj,GSx=c("CD3D","CD3E","TRDC","TRGC1","TRGC2"),
									GSy=c("CD8A","CD8B"),
									col.name="Score.gammaDeltaT",th.score=0.35,out.prefix=NULL)
{
	if(!all(GSx %in% rowData(obj)$display.name) || !all(GSy %in% rowData(obj)$display.name)){
		warning(sprintf("not all genes in GSx and GSy in the obj !!\n"))
		return(obj)
	}

	f.gene <- which(rowData(obj)$display.name %in% GSx)
	GSx <- rowData(obj)$display.name[f.gene]
	f.gene <- which(rowData(obj)$display.name %in% GSy)
	GSy <- rowData(obj)$display.name[f.gene]
	count.tot <- colSums(assay(obj,"counts"))
	GSx.tot <- colSums(assay(obj,"counts")[names(GSx),])
	GSy.tot <- colSums(assay(obj,"counts")[names(GSy),])
	GSx.score <- GSx.tot/count.tot
	GSy.score <- -GSy.tot/count.tot
	GSx.score <- (GSx.score-min(GSx.score))/(max(GSx.score)-min(GSx.score))
	GSy.score <- (GSy.score-min(GSy.score))/(max(GSy.score)-min(GSy.score))
	comb.score <- GSx.score * GSy.score
	colData(obj)[[col.name]] <- comb.score
	colData(obj)[[sprintf("%s.pred.th",col.name)]] <- (comb.score > th.score)
	print(table(obj$stype,obj[[sprintf("%s.pred.th",col.name)]]))
	bin.tb <- sscClust::classify.outlier(comb.score,out.prefix=out.prefix,e.name="Score",th.score=th.score)
	all(bin.tb$sample==colnames(obj))
	colData(obj)[[sprintf("%s.pred.auto",col.name)]] <-(bin.tb$score.cls.tb$classification==1)
	print(table(obj$stype,obj[[sprintf("%s.pred.auto",col.name)]]))
	freq.tb <- unclass(table(obj$stype,obj[[sprintf("%s.pred.auto",col.name)]]))
	freq.tb <- 100*sweep(freq.tb,1,rowSums(freq.tb),"/")

	### subtypes: delta 1 and delta 2
	gene.GC1 <- c("TRGC1")
	gene.GC2 <- c("TRGC2")
	f.gene <- which(rowData(obj)$display.name %in% gene.GC1)
	gene.GC1 <- rowData(obj)$display.name[f.gene]
	f.gene <- which(rowData(obj)$display.name %in% gene.GC2)
	gene.GC2 <- rowData(obj)$display.name[f.gene]
	exp.diff <- assay(obj,"norm_exprs")[names(gene.GC1),] - assay(obj,"norm_exprs")[names(gene.GC2),]
	f.isGammDeltaT <- obj[[sprintf("%s.pred.auto",col.name)]]==T
	colData(obj)[[sprintf("%s.auto.subtype",col.name)]] <- "NotGammaDeltaT"
	colData(obj)[[sprintf("%s.auto.subtype",col.name)]][exp.diff==0 & f.isGammDeltaT] <- "undetermined"
	colData(obj)[[sprintf("%s.auto.subtype",col.name)]][exp.diff>0 & f.isGammDeltaT] <- "delta2"
	colData(obj)[[sprintf("%s.auto.subtype",col.name)]][exp.diff<0 & f.isGammDeltaT] <- "delta1"
	colData(obj)[[sprintf("%s.auto.subtype",col.name)]] <- factor(colData(obj)[[sprintf("%s.auto.subtype",col.name)]],
																  levels=c("delta1","delta2","undetermined","NotGammaDeltaT"))
	table(obj[[sprintf("%s.pred.auto",col.name)]],obj[[sprintf("%s.auto.subtype",col.name)]])

	return(obj)
}

inSilico.TGammaDelta <- function(obj,out.prefix=NULL,assay.name="norm_exprs",vis.v=c(0.25,0.5),
                           Th.CD3=0.25,Th.DC=0.25,Th.GC1=0.25,Th.GC2=0.25)
{
    library("data.table")
    library("ggplot2")
    gene.to.test <- c("CD3D","CD3G","TRDC","TRGC1","TRGC2")
    f.gene <- which(rowData(obj)$display.name %in% gene.to.test)
    gene.to.test <- structure(rowData(obj)$display.name[f.gene],names=rownames(obj)[f.gene])
    gene.to.test <- gene.to.test[order(gene.to.test)]

    dat.plot <- as.data.frame(t(as.matrix(assay(obj,assay.name)[names(gene.to.test),])))
    setDT(dat.plot,keep.rownames=T)
    colnames(dat.plot)[-1] <- gene.to.test
    dat.plot[,"CD3"] <- apply(dat.plot[,c("CD3D","CD3G")],1,mean)
    dat.plot.melt <- melt((dat.plot),id.vars="rn")
    colnames(dat.plot.melt) <- c("cell","gene","norm_exprs")

    p <- ggplot(dat.plot.melt, aes(norm_exprs, fill = gene, colour = gene)) +
        geom_density(alpha = 0.1) +
		theme_bw() +
        geom_vline(xintercept = vis.v,linetype=2) +
        facet_wrap(~gene,ncol=3,scales="free_y")
    ggsave(sprintf("%s.inSiliso.marker.density.pdf",out.prefix),width=7,height=3.5)

    obj$TGammaDelta <- "unknown"
    obj$TGammaDelta[dat.plot[["CD3"]] > Th.CD3 & dat.plot[["TRDC"]] > Th.DC ] <- "undetermined"
    obj$TGammaDelta[dat.plot[["CD3"]] > Th.CD3 & dat.plot[["TRDC"]] > Th.DC & ( dat.plot[["TRGC1"]] > Th.GC1 & dat.plot[["TRGC2"]] < Th.GC2 ) ] <- "delta2"
    obj$TGammaDelta[dat.plot[["CD3"]] > Th.CD3 & dat.plot[["TRDC"]] > Th.DC & ( dat.plot[["TRGC1"]] < Th.GC1 & dat.plot[["TRGC2"]] > Th.GC2 ) ] <- "delta1"
    table(obj$TGammaDelta)

    #bin.tb <- sscClust::classify.outlier(dat.plot[["TRDC"]],out.prefix=out.prefix,e.name="TRDC",th.score=0.25)

    ##write.table(colData(obj),sprintf("%s.cellInfo.txt",out.prefix),row.names=F,sep="\t",quote=F)
    return(obj)
}

inSilico.TCell <- function(sce, out.prefix, assay.name="norm_exprs",vis.v=c(0.25,0.5),
                           Th.CD3=0.25,Th.CD8=0.5,Th.CD4=0.5,Th.TH=0.25,Th.TR=0.25,do.zscore=F,
			   Th.DC=0.25,Th.GC1=0.25,Th.GC2=0.25,
			   do.rescue=T)
{
    #### in silico classification
    library("data.table")
    library("ggplot2")
    library("ggpubr")
    ##gene.to.test <- c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G","CD40LG","FOXP3","IL2RA")
    gene.to.test <- c("CD4","CD8A","CD8B","CD3D","CD3G","CD40LG","FOXP3","IL2RA",
		      "TRDC","TRGC1","TRGC2")
    f.gene <- which(rowData(sce)$display.name %in% gene.to.test)
    if(do.zscore){
	sce.z <- sce[f.gene,]
	sce.z <- ssc.scale(sce.z,gene.symbol=gene.to.test,assay.name=assay.name,
			   adjB="batchV",do.scale=T)
	gene.to.test <- structure(rowData(sce.z)$display.name,names=rownames(sce.z))
	gene.to.test <- gene.to.test[order(gene.to.test)]
	dat.plot <- as.data.frame(t(as.matrix(assay(sce.z,sprintf("%s.scale",assay.name))[names(gene.to.test),])))
	colnames(dat.plot) <- gene.to.test
    }else{
	gene.to.test <- structure(rowData(sce)$display.name[f.gene],names=rownames(sce)[f.gene])
	gene.to.test <- gene.to.test[order(gene.to.test)]
	dat.plot <- as.data.frame(t(as.matrix(assay(sce,assay.name)[names(gene.to.test),])))
	colnames(dat.plot) <- gene.to.test
    }

    ##dat.plot[,"CD3"] <- apply(dat.plot[,c("CD3D","CD3E","CD3G")],1,mean)
    dat.plot[,"CD3"] <- apply(dat.plot[,c("CD3D","CD3G")],1,mean)
    dat.plot[,"CD8"] <- apply(dat.plot[,c("CD8A","CD8B")],1,mean)
    dat.plot[,"TH"] <- apply(dat.plot[,c("CD4","CD40LG")],1,mean)
    dat.plot[,"TR"] <- apply(dat.plot[,c("CD4","FOXP3")],1,mean)
    dat.plot.melt <- melt(as.matrix(dat.plot))
    colnames(dat.plot.melt) <- c("cell","gene","norm_exprs")

    #p <- ggplot(dat.plot.melt, aes(norm_exprs, fill = gene, colour = gene)) +
    #    geom_density(alpha = 0.1) +
    p <- ggdensity(dat.plot.melt,x="norm_exprs",fill="gene",color="gene",
		   alpha=0.1) +
        geom_vline(xintercept = vis.v,linetype=2) +
        facet_wrap(~gene,ncol=3,scales="free_y") +
	theme_bw() +
	theme(legend.position="none")
    ggsave(sprintf("%s.inSiliso.marker.density.pdf",out.prefix),width=7,height=6)

    sce$stype <- "unknown"
    #sce$stype[dat.plot[,"CD3"] < 0.25] <- "noCD3"
    sce$stype[dat.plot[,"CD3"] > Th.CD3 &
	      dat.plot[,"CD8"] > Th.CD8 & dat.plot[,"CD4"] < Th.CD4] <- "CD8"
    sce$stype[dat.plot[,"CD3"] > Th.CD3 &
	      dat.plot[,"CD8"] > Th.CD8 & dat.plot[,"CD4"] > Th.CD4] <- "DP"
    sce$stype[dat.plot[,"CD3"] > Th.CD3 &
	      dat.plot[,"CD8"] < Th.CD8 & dat.plot[,"CD4"] > Th.CD4] <- "CD4"
    sce$stype[dat.plot[,"CD3"] > Th.CD3 &
	      dat.plot[,"CD8"] < Th.CD8 & dat.plot[,"CD4"] < Th.CD4] <- "DN"
    sce$stype.strict <- sce$stype

    #####
    has.gd.gene <- F
    if(all(c("TRDC","TRGC1","TRGC2") %in% colnames(dat.plot))){
	x.gdType <- rep("unknown",ncol(sce))
	x.gdType[ dat.plot[,"CD3"] > Th.CD3 & dat.plot[,"TRDC"] > Th.DC  ] <- "undetermined"
	x.gdType[ dat.plot[,"CD3"] > Th.CD3 & dat.plot[,"TRDC"] > Th.DC &
		 (dat.plot[,"TRGC1"] > Th.GC1 & dat.plot[,"TRGC2"] < Th.GC2 )  ] <- "delta2"
	x.gdType[ dat.plot[,"CD3"] > Th.CD3 & dat.plot[,"TRDC"] > Th.DC &
		 (dat.plot[,"TRGC1"] < Th.GC1 & dat.plot[,"TRGC2"] > Th.GC2 )  ] <- "delta1"
	sce$gdType <- x.gdType
	has.gd.gene <- T
    }else{
	warning(sprintf("There are no gd genes found in the data\n"))
    }
    #####

    f.cell <- sce$stype=="DN" & (dat.plot[,"TH"] > Th.TH | dat.plot[,"TR"] > Th.TR)
    ###f.cell <- sce$stype=="DN" & dat.plot[,"CD4"] > 0 & (dat.plot[,"TH"] > Th.TH | dat.plot[,"TR"] > Th.TR)
    #print(table(sce$stype))
    cat(sprintf("numbe of cells can be rescued: %d\n",sum(f.cell)))
    if(do.rescue){
	### todo: add requirement: DN and not gamma delta T ?
	sce$stype[f.cell] <- "CD4"
	sce$stype.rescue <- sce$stype
    }

    ####write.table(colData(sce),sprintf("%s.cellInfo.txt",out.prefix),row.names=F,sep="\t",quote=F)
    return(sce)
}

my.plot.volcano <- function(dat.plot,out.prefix,g.f=NULL,col.x="logFC",col.y="adj.P.Val",
							sample.id="",g.GNAME=NULL,my.xlim=c(-7,7))
{
	require("ggplot2")
	require("ggrepel")
	dat.plot <- dat.plot[order(dat.plot[,col.y],-dat.plot[,col.x]),]
	rownames(dat.plot) <- dat.plot$geneID
	print(head(dat.plot[,1:5],n=30))
	if(is.null(g.f)){
		g.f <- head(dat.plot$geneID,n=30)
	}
	dat.plot$Significant <- abs(dat.plot[,col.x])>1 & dat.plot[,col.y]<0.05
	dat.plot$logP <- -log10(dat.plot[,col.y])
	dat.plot$logP[dat.plot$logP>=300] <- 300
	y.pretty <- pretty(dat.plot$logP)
	y.range <- y.pretty[c(1,length(y.pretty))]
	n.up <- sum(dat.plot$Significant & dat.plot[,col.x]>0)
	n.down <- sum(dat.plot$Significant & dat.plot[,col.x]<0)
	#annotations <- data.frame(xpos= -Inf,ypos=Inf,hjustvar=-0.2,vjustvar=1.2,
	#						  annotateText=sprintf("Up: %d\nDown: %d",n.up,n.down))
	annotations <- data.frame(xpos= c(-Inf,-Inf),ypos=c(Inf,Inf),
							  hjustvar=c(-0.2,-0.15),vjustvar=c(2.0,3.4),
							  annotateText=c(sprintf("Up: %d",n.up),sprintf("Down: %d",n.down)))
	my.plot.volcano <- ggplot(dat.plot, aes_string(x = col.x, y = "logP")) + 
		geom_point(aes(color = Significant),size=0.5,data = subset(dat.plot,Significant==F)) + 
		geom_point(aes(color = Significant),size=1.0,data = subset(dat.plot,Significant==T)) + 
		scale_color_manual(values = c("grey", "red","#1B9E77", "#D95F02", "#7570B3")) +  
		labs(title = sample.id) + xlab("log2FoldChange") + ylab("-log10(pvalue)") + 
		#coord_cartesian(ylim = c(0, y.range[2]),xlim=c(-10,10)) +
		coord_cartesian(ylim = c(0, y.range[2]),xlim=my.xlim) +
		theme_bw(base_size = 12) + 
		theme(legend.position = "bottom", 
			  title=element_text(size=18),
			  axis.text=element_text(size=14), 
			  axis.title=element_text(size=16)) + 
		geom_text_repel(data = dat.plot[g.f,], 
						#aes(label = g.GNAME[rownames(dat.plot[g.f,])]),
						aes(label = if(is.null(g.GNAME)) g.f else g.GNAME[g.f] ),
						size = 4, box.padding = unit(0.35, "lines"),fontface=3, 
						point.padding = unit(0.3, "lines")) +
		geom_text(data=annotations,
				  aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))
	ggsave(sprintf("%s.volcano.pdf",out.prefix),width=6,height=6)
}


findCorrelatedGeneToAnchor.clust <- function(sce,assay.name="norm_exprs",a.gene.vec,out.prefix,
                                          TH.score=0.17,n.downsample=NULL,gene.ignore=NULL,
                                          seed=1234,nthreads=8,K=20,dyn.cut=T)
{
    require("dynamicTreeCut")
    require("gplots")
    dir.create(dirname(out.prefix),F,T)
    set.seed(seed)
    if(!is.null(gene.ignore)){
        sce <- sce[!(rownames(sce) %in% gene.ignore),]
    }
    if(!is.null(n.downsample)){
        f.cells <- sample(ncol(sce),min(n.downsample,ncol(sce)))
    }else{
        f.cells <- seq_len(ncol(sce))
    }
    f.gene.anchor <- match(a.gene.vec,rowData(sce)$display.name)
    
    #anchor.mtx <- assay(sce[f.gene.anchor,f.cells],"counts")
    #other.mtx <- assay(sce[,f.cells],"counts")
    anchor.mtx <- assay(sce[f.gene.anchor,f.cells],assay.name)
    f.zero <- rowSums(assay(sce[,f.cells],"counts"))>0
    print(summary(f.zero))
    other.mtx <- assay(sce[f.zero,f.cells],assay.name)

    cor.gene <-  sscClust:::cor.BLAS(anchor.mtx, y = other.mtx, method = "spearman", nthreads = nthreads)
    print(cor.gene[,1:3])
    this.score <- apply(cor.gene,2,max)
    score.table <- data.table(geneID=names(this.score),
                              geneSymbol=rowData(sce)[names(this.score),"display.name"],
                              score=this.score)
    ##this.score <- apply(cor.gene,2,median)

    f.gene.cor <- which(score.table$score>TH.score)

    cor.mtx <- sscClust:::cor.BLAS(other.mtx[f.gene.cor,], method = "spearman", nthreads = 8)
    cor.mtx[1:5,1:3]

    dist.obj <- as.dist(1-cor.mtx)
    hclust.obj <- hclust(dist.obj,method="ward.D2")
    #hclust.obj <- hclust(dist.obj,method="complete")
    if(dyn.cut){
        cluster.obj <- dynamicTreeCut::cutreeDynamic(hclust.obj, distM=as.matrix(dist.obj), verbose=0)
    }else{    
        cluster.obj <- cutree(hclust.obj,K)
    }
    cluster.df <- data.table(geneID=colnames(cor.mtx),
                             geneSymbol=rowData(sce)[colnames(cor.mtx),"display.name"],
                             cluster=unname(cluster.obj))
                             ##cluster=unname(hclust.cutDyn.obj))

    ### average correlation of the modules
    plot.cor.matix <- function(X,cls){
        gname <- colnames(X)
        gname <- rowData(sce)[gname,"display.name"]
        colnames(X) <- gname
        rownames(X) <- gname
          
        dat.plot <- X
        hc <- hclust(as.dist((1-dat.plot)/2))
        dat.plot <-dat.plot[hc$order, hc$order]
        diag(dat.plot) <- NA
        v.lim <- c(-0.4,0.4)
        dat.plot[ dat.plot > v.lim[2] ]  <- v.lim[2]
        dat.plot[ dat.plot < v.lim[1] ]  <- v.lim[1]
        dat.plot <- data.table::melt(dat.plot)
        p <- ggplot(data = dat.plot, aes(Var2, Var1, fill = value))+
            geom_tile(color = NA,width=1.0,height=1.0)+xlab("")+ylab("")+
            scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                 midpoint = 0, limit = v.lim, space = "Lab",
                                 name="Correlation") +
            #theme_bw()+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
            coord_fixed()
        ggsave(sprintf("%s.cluster.%d.png",out.prefix,cls),width=6,height=6)
    }

    average.upper.tri <- function(X,cls){
        f <- upper.tri(X)
        return(median(X[f]))
    }
    cluster.avg.cor <- cluster.df[,.(N=.N,avg.cor=average.upper.tri(cor.mtx[.SD$geneID,.SD$geneID])),
                                  by="cluster"]
    for(cls in cluster.avg.cor$cluster){
        gid <- cluster.df[cluster==cls,][["geneID"]]
        plot.cor.matix(cor.mtx[gid,gid],cls)
    }
    cluster.df <- cluster.df[cluster.avg.cor,,on="cluster"][order(-avg.cor,cluster,geneSymbol),]
    return(cluster.df)
    
    #anchor.cluster <- cluster.df[geneID %in% rownames(anchor.mtx),]
    #cluster.flt.df <- cluster.df[cluster %in% anchor.cluster[["cluster"]] & avg.cor>TH.score,
    #                             ][order(-avg.cor,cluster),]
    #write.table(cluster.flt.df,sprintf("%s.cutree.glist",out.prefix),row.names=F,sep="\t",quote=F)
    #return(cluster.flt.df)
}


my.FindTransferAnchors <- function (reference, query, reference.assay = NULL, query.assay = NULL,
    reduction = "pcaproject", project.query = FALSE, features = NULL,
    npcs = 30, l2.norm = TRUE, dims = 1:30, k.anchor = 5, k.filter = 200,
    k.score = 30, max.features = 200, eps = 0, approx.pca = TRUE,
    verbose = TRUE)
{
    if (length(x = reference) > 1 | length(x = query) > 1) {
        stop("We currently only support transfer between a single query and reference")
    }
    if (!reduction %in% c("pcaproject", "cca", "pcaqueryproject")) {
        stop("Please select either pcaproject, cca, or pcaqueryproject for the reduction parameter.")
    }
    query <- RenameCells(object = query, new.names = paste0(Cells(object = query),
        "_", "query"))
    reference <- RenameCells(object = reference, new.names = paste0(Cells(object = reference),
        "_", "reference"))
    features <- if(!is.null(features)) features else  VariableFeatures(object = reference)
    reference.assay <- if(!is.null(reference.assay)) reference.assay else DefaultAssay(object = reference)
    query.assay <- if(!is.null(query.assay)) query.assay else  DefaultAssay(object = query)
    DefaultAssay(object = reference) <- reference.assay
    DefaultAssay(object = query) <- query.assay
#    print(str(features))

    if (reduction == "cca") {
#        reference <- ScaleData(object = reference, features = features,
#                               do.scale=F,vars.to.regress = c("batchV"),
#                               verbose = T)
#        query <- ScaleData(object = query, features = features,
#                               do.scale=F,vars.to.regress = c("patient"),
#                               verbose = T)

        exp.r.mat <- sscClust:::simple.removeBatchEffect(GetAssayData(reference,"data")[features,],
                                    batch = reference$batchV)
        exp.q.mat <- sscClust:::simple.removeBatchEffect(GetAssayData(query,"data")[features,],
                                    batch = query$patient)
        reference <- SetAssayData(reference,"scale.data",exp.r.mat)
        query <- SetAssayData(query,"scale.data",exp.q.mat)

        combined.ob <- RunCCA(object1 = reference, object2 = query,
            features = features, num.cc = max(dims), renormalize = FALSE,
            rescale = FALSE, verbose = verbose)
    }
    if (l2.norm) {
        combined.ob <- L2Dim(object = combined.ob, reduction = reduction)
        reduction <- paste0(reduction, ".l2")
    }
    anchors <- Seurat:::FindAnchors(object.pair = combined.ob, assay = c(reference.assay,
        query.assay), cells1 = colnames(x = reference), cells2 = colnames(x = query),
        reduction = reduction, dims = dims, k.anchor = k.anchor,
        k.filter = k.filter, k.score = k.score, max.features = max.features,
        eps = eps, verbose = verbose)
    anchor.set <- new(Class = "AnchorSet", object.list = list(combined.ob),
        reference.cells = colnames(x = reference), query.cells = colnames(x = query),
        anchors = anchors, anchor.features = features)
    return(anchor.set)
}

fill.contamination <- function(seu,out.prefix,g.name="plasmaB",g.test=c("CD79A", "JCHAIN", "SDC1"),
							   score.t=1,vis.v=c(0.25,0.5,1))
{
	require("ggplot2")
	require("data.table")
	## plasma contamination
	B.score.tb <- data.table(cell=colnames(seu),
				 B.score= colMeans(GetAssayData(seu,"data")[g.test,,drop=F]))
	B.score.tb <- cbind(B.score.tb,as.matrix(t(GetAssayData(seu,"data")[g.test,,drop=F])))
	seu[[sprintf("%s.score",g.name)]] <- B.score.tb$B.score
	seu[[sprintf("%s.class",g.name)]] <- seu[[sprintf("%s.score",g.name)]] > score.t

	B.score.tb.melt <- melt(B.score.tb,id="cell")
	colnames(B.score.tb.melt) <- c("cell","gene","norm_exprs")

	p <- ggplot(B.score.tb.melt, aes(norm_exprs, fill = gene, colour = gene)) +
			geom_density(alpha = 0.1) +
			geom_vline(xintercept = vis.v,linetype=2) +
			facet_wrap(~gene,ncol=2,scales="free_y")
	ggsave(sprintf("%s.marker.%s.density.pdf",out.prefix,g.name),width=8,height=8)
	return(seu)
}



