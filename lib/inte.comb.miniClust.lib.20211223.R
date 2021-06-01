library("sscClust")
library("Seurat")
library("harmony") 
library("tictoc")
library("plyr")
#library("dplyr")
#library("tibble")
library("doParallel")
library("Matrix")
library("data.table")
library("R.utils")
library("gplots")
library("ggpubr")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("limma")
library("lisi")
library("RColorBrewer")
library("ComplexHeatmap")
#library("metacell")

source("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/PanC.T/lib/plot.func.R")

### special genes

g.GeneOnUmap.list <- list(
			  "Epi"=c("EPCAM","CLDN4","CLDN7","KRT18","KRT8","KRT19","CAPS","WFDC2","CD24"),
			  "Endo"=c("PECAM1","VWF","CDH5","PLVAP","CLDN5","SLCO2A1","RAMP2","EGFL7","EMCN"),
			  "Fibro"=c("COL1A1","COL1A2","COL6A1","COL6A2","COL6A3","DCN","FAP","THY1","LUM"),
			  "myofibroblast"=c("RGS5", "NDUFA4L2", "HIGD1B", "BGN", "COX4I2", "LHFP", "SPARC", "MYL9", "GJA4"),
			  "SMC.01"=c("MYH11", "ACTA2", "TAGLN", "CNN1", "MYL9", "DES", "TPM2", "MYLK", "PLN"),
			  "SMC.02"=c("MYH11", "TAGLN", "ACTA2", "MYL9", "PLN", "TPM2", "CNN1", "CSRP1", "ADIRF"),
			  "DC.01"=c("LILRA4","LAMP3", "CLEC9A","CLEC10A","CD1C", "FCER1A", "IDO1","CCL19", "CD274"),
			  "DC.02"=c("ITGAX","IRF4", "IRF8", "LTB", "CCL21","IL3RA","NRP1","CLE4C","RSG2", "THBD"),
			  "Monocyte"=c("CD14","VCAN", "S100A9","FCGR3A","FCEG1G","TIMP1",
				       "LST1","ITGAM", "HLA-DQA1","PTPRC"),
			  "Macrophase.00"=c("CD68","CD14","FCGR3A","THBS1","MARCO", "C1QA","APOE",
					    "SEPP1","ITGAM", "HLA-DQA1","CSF3R","CD3D"),
			  "Macrophase.01"=c("CD68","CD14","FCGR3A","THBS1","MARCO", "C1QA","APOE",
					    "SEPP1","RNASE1","CD163","FCN1","VCAN","LST1","LILRB2"),
			  "Macrophase.02"=c("HLA-DQA1", "CST3", "LYZ","CCR5","TFRC","CD52", "ISGF6","S100A4","S100A9",
					    "FUT4", "CD33","IFI30","FCGR1A","GPNMB"),
			  "Macrophase.03"=c("HLA-DQA1", "HLA-DRB1", "HLA-DRA", "ISGF6","IFI35", "","IFI30","FCGR1A","LILRA4"),
			  "Neutrophils"=c("FUT4","CD3E","TLR2","CD14","CD33","FCGR1A", "CSF3R","IFITM2","S100A8"),
			  "NK"=c("FCGR3A", "NCR1","NCAM1","CD3E", "NKG7", "KLRD1","KIR2DL3","KIR3DL2","KLRB1"),
			  "T.01"=c("CD4", "CD8A", "CD8B", "CD2", "CD3D", "CD3G", "NKG7", "GNLY", "FGFBP2"),
			  "T.02"=c("CX3CR1", "GZMK", "ZNF683", "IL7R","GZMB", "CD40LG","FOXP3","CCR7","MKI67"),
			  "T.03"=c("TRGC1","TRGC2","TRBC1","TRBC2","TRDC","TRAC","FOXP3","IL2RA","CD40LG"),
			  "T.04"=c("CCR7","GZMK","LAYN","HAVCR2","CXCL13","CX3CR1","ZNF683","S1PR1","RGS1"),
			  "T.05"=c("CXCR6","SLC4A10","CCL5","CD160","CD6","TMIGD2","MAL","SELL","KIR2DL4"),
			  "T.06"=c("ADAM19","GZMA","CXCR6","HAVCR2","CXCR5","BCL6","RORC","IL23R","IFNG"),
			  "T.07"=c("SELL","MAL","TCF7","CTLA4","IL10","IL2RA","ANXA1","ANXA2","TNFRSF9"),
			  "T.08"=c("CCL4L2","CCL4","XCL1","XCL2","CRTAM","RGCC","LMNA","FCER1G","KRT81"),
			  "T.09"=c("MALAT1","RPS27","RPL32","PTPRC","TXNIP","HBB","HBA2", "HBA1","MCL1","TSC22D3",
				   "ASH1L","CLDN15", "ZNF117", "NF1", "ZNF546","HLA-A","HLA-B","HLA-C"),
			  "T.10"=c("TRAF3IP3","LCP2","GIMAP7","CCL4","MALAT1","MT-ND4","TC2N","MT-CYB","TXNIP"),
			  "T.11"=c("SLC25A5","ARHGDIB","ARPC1B","CORO1A","NDUFA12","ATP5C1","CAPZB","SUB1", "TXNIP"),
			  "T.12"=c("RUNX2","ANKRD28","TMIGD2","CTSH","GAB3","IL23R","FURIN","BHLHE40","RORC",
				   "TNFSF13B","CCL20","ABCB1","CCR6"),
			  "T.13"=c("TOX","PDCD1","MYO1E","MYO7A","GEM","ACP5","ETV1","AKAP5","LAG3"),
			  "T.14"=c("CISH","MYL12A","TBC1D10C","HCST","GIMAP1","RAC2","EVL","CKLF","GZMH"),
			  "T.15"=c("PFN1","LGALS1","PPP1CA","ACTB","SH3BGRL3","EMP3","NDUFB7","S100A4","S100A6"),
			  "T.16"=c("CORO1B", "ACTB", "SIT1", "CPNE2", "CENPW", "NME1", "NME2","FOXP3","CD40LG"),
			  "T.17"=c("PASK","TCF7","CCR7","LEAF1", "SESN3", "RNASET2", "TXK", "SELL", "ACTN1", "ADSL", "LDHB"),
			  "T.18"=c("PARK7","FKBP1A","ATP5G3","RHOA","PSME2","VDAC1","CD63","GAPDH","ENO1"),
			  "T.19"=c("KLRG1","TNFRSF1B","VPS37B","SLC7A5","CREM","ZNF331","DUSP4","RGS1","SRGN"),
			  "T.20"=c("ITM2C","GLIPR2","GZMH","S100A10","CLU","VIM","CD52","HLA-DPA1","TMSB4X"),
			  "T.21"=c("EOMES","KLRG1","CD74","CXCR5","CXCR4","CD28","CD27","SH2D1A","TBX21"),
			  "cellCycle"=c("TYMSOS","PCLAF","TYMS","ZWINT","MCM2","MIR3917","UHRF1","TK1","MKI67"),
			  "IFN"=c("ISG15","ISG20","IFIT1","IFIT2","IRF7","OAS1","OAS2","STAT1"),
			  "stress"=c("DNAJB1","DUSP1","FOS","JUN","HSP90AA1","HSP90AB1","HSPA1A","HSPA6"),
			  "other.01"=c("PTPRC", "RGS1", "PDCD1", "ITGAE","HAVCR2", "CXCL13","TNFRSF9","CTLA4","LAYN"),
			  "other.02"=c("CX3CR1","GZMK","IFNG", "CXCL13","HAVCR2","CXCR5","FOXP3","ENTPD1","IFIT1"),
			  "other.03"=c("IL7R","CCR7","TCF7","ZFP36","TNF","CAPG","CTSH","CCL5","CREM"),
			  "other.04"=c("IL7R","GZMK","CX3CR1","ZNF683","SLC4A10",
				       "CXCL13","TYROBP","LAYN","IFIT1", "TCF7","MAL","ZFP36"),
			  "other.05"=c("CX3CR1","GZMK","IFNG", "CXCL13","HAVCR2","CXCR5","FOXP3",
				       "TNFRSF9","IFIT1","IL7R","TCF7","ZFP36","TNF","CTSH","CCL5","CREM"),
			  "other.06"=c("CCR7","SELL","LEF1","TCF7","ANXA1","PASK","ADSL","ZFP36",
				       "IL7R","CREM","TNF","CCL5","GZMA","CAPG","GZMK","CX3CR1",
				       "RORC","CXCR5","CXCL13","CD40LG","FOXP3","RTKN2","S1PR1","RGS1","TNFRSF9",
				       "OAS1","IFIT1","NME1","ACTB1")
					)


HVG.From.GeneRankTb <- function(gene.rank.tb,n.common=1000,n.specific=1000,th.rank=0.1)
{
    gene.rank.tb <- gene.rank.tb[order(median.F.rank),]
    gene.rank.tb[,rank:=seq_len(nrow(gene.rank.tb))]
    gene.rank.mat <- as.matrix(gene.rank.tb[,-c("geneID","median.F.rank","rank")])
    rownames(gene.rank.mat) <- gene.rank.tb$geneID
    gene.rank.mat[1:4,1:5]
    gene.rank.comp.mat <- 1 - gene.rank.mat
    
    ###
    f.common <- gene.rank.tb[rank <= n.common,][["geneID"]]
    gene.rank.tail.mat <- gene.rank.mat[!(rownames(gene.rank.mat) %in% f.common),]
    f.specific.tb <- as.data.table(ldply(seq_len(nrow(gene.rank.tail.mat)),function(i){
			    x <- gene.rank.tail.mat[i,]	
			    ff <- x < th.rank
			    nSpe <- sum(ff)
			    out.tb <- data.table(geneID=rownames(gene.rank.tail.mat)[i],
						 nDataSets=nSpe,
						 fDataSets=nSpe/length(x),
						 medianRankSpeGene=if(nSpe > 0) median(x[ff]) else 1)
			    out.tb[,hasSpeGene:= (nSpe>=3 & fDataSets > 0.1) ]
					}))
    f.specific.tb <- f.specific.tb[hasSpeGene==T,][order(medianRankSpeGene),]
    f.specific.tb[,rank:=seq_len(nrow(f.specific.tb))]
    print(f.specific.tb[geneID %in% c("IL17A","IL17F","IL23R","RORC","IL26","IL22","IL21","CD4","CD8A","CD8B"),])
    f.specific <- f.specific.tb[rank <= n.specific,][["geneID"]]

    ####
#    sscVis:::plotMatrix.simple(gene.rank.tail.mat[f.specific,],
#			       out.prefix=sprintf("%s.gene.rank.datasetSpecific",out.prefix),
#			       show.dendrogram=T,clust.row=T,clust.column=T,
#			       pdf.width = 8,pdf.height = 10,
#			       exp.name="1-percRank")
    ret.tb <- rbind(data.table(geneID=f.common,category="common"),
		    data.table(geneID=f.specific,category="specific"))
    return(ret.tb)
    ##return(c(f.common,f.specific))
}


#### x: row for sample, column for variable
#### implementation in Seurat (FindClusters) can only detect singletons :-(
run.Leiden <- function(dat.pca,SNN.k=20,myseed=123456,...)
{
	####x <- Embeddings(seu.merged,"harmony")[,1:15]
	snn.gr <- scran::buildSNNGraph(t(dat.pca), k=SNN.k,d=NA)
	clust.res <- leiden::leiden(snn.gr,seed=myseed,...)
	return(clust.res)
}

mergeDataFromFileTable <- function(exp.list.table,gene.de.common,seu.list,sce.list,
								   res.hi,method.clustering="louvain",contamination.vec=NULL,cor.MALAT1=F,cor.DIG=T)
{
    ret.list <- llply(seq_len(nrow(exp.list.table)),function(i){
							data.id <- exp.list.table$data.id[i]
							##### seu and sce --> seu.x
							seu <- seu.list[[i]]
							sce <- sce.list[[i]]
							cat(sprintf("all(colnames(seu)==colnames(sce)) in %s: %s\n",data.id,
										all(colnames(seu)==colnames(sce))))
							dat.x <- assay(sce,"norm_exprs")
							if(!"seu.id" %in% colnames(rowData(sce))){
								rowData(sce)[,"seu.id"] <- gsub("_","-",rowData(sce)[,"display.name"])
							}
							rownames(dat.x) <- rowData(sce)[,"seu.id"]
							#### pad zero
							{
							    gene.pad <- setdiff(gene.de.common,rownames(dat.x))
							    cat(sprintf("For dataset %s, those genes will be imputed by 0: %s\n",
										    data.id, 
										    paste(gene.pad,collapse=", ")))
							    if(length(gene.pad)>0){
								    value.pad <-  0
								    mtx.pad <- matrix(rep(value.pad,length(gene.pad)*ncol(dat.x)),nrow=length(gene.pad))
								    rownames(mtx.pad) <- gene.pad
								    dat.x <- rbind(dat.x,mtx.pad)
							    }
							}

						    
							score.MALAT1 <- assay(sce,"norm_exprs")[rowData(sce)[,"display.name"]=="MALAT1",]

							if(data.id %in% c("HCC.YaoHe10X","HCC.YaoHeSS2") ){
							    gene.MALAT1.vec <- c("MALAT1-ENSG00000251562")
							}else{
							    gene.MALAT1.vec <- c("MALAT1")
							}
							#seu <- AddModuleScore(seu, features=list("score.MALAT1"=gene.MALAT1.vec), name="score.MALAT1",
							#					  pool = NULL, nbin = 24, ctrl = 100)

							sce <- NULL
							if(!is.null(contamination.vec)){
							    f.cell.cont <- sprintf("%s.%s",data.id,colnames(dat.x)) %in% contamination.vec
							    cat(sprintf("Number of contaminated cells (%s): %d\n",data.id,sum(f.cell.cont)))
							    dat.x <- dat.x[,!f.cell.cont]
							}
							seu.x <- CreateSeuratObject(dat.x,project="panC", meta.data=seu[[]][colnames(dat.x),])
							### seu$score.MALAT11, two "1"
							#score.MALAT1 <- seu$score.MALAT11
							seu.x$score.MALAT1 <- score.MALAT1
							seu <- NULL
							gc()
							##### regression
							meta.extra.tb <- data.frame(cellID.uniq=sprintf("%s.%s",data.id,colnames(seu.x)),
														S.Score=seu.x$S.Score,
														G2M.Score=seu.x$G2M.Score,
														Phase=seu.x$Phase,
														DIG.Score1=seu.x$DIG.Score1,
														score.MALAT1=seu.x$score.MALAT1,
														percent.mito=NA,
														stringsAsFactors=F)
							adj.cov <- c("S.Score","G2M.Score")
							if(cor.DIG){
							    adj.cov <- c(adj.cov,"DIG.Score1")
							}
							if(cor.MALAT1){
							    adj.cov <- c(adj.cov,"score.MALAT1")
							}
							if("percent.mito" %in% colnames(seu.x[[]])){
							    adj.cov <- c("percent.mito",adj.cov)
							    meta.extra.tb$percent.mito <- seu.x$percent.mito
							}
							nBatch <- length(table(seu.x$batchV))
							###adj.cov <- c()
							if(nBatch>1){
							    adj.cov <- c("batchV",adj.cov)
							}
							cat(sprintf("adj.cov (%s): %s\n",data.id,paste(adj.cov,collapse=", ")))
							
							seu.x <- ScaleData(seu.x,do.scale=T,features=gene.de.common,
											   vars.to.regress = adj.cov,verbose=F)
							###### Seurat high resolution ####
							seu.x <- RunPCA(seu.x,features=rownames(seu.x),npcs= 15,verbose = FALSE)
							if(method.clustering=="leiden"){
								clust.x <- run.Leiden(Embeddings(seu.x,"pca"),SNN.k=20,myseed=123456,resolution_parameter=res.hi)
								seu.x@meta.data[[sprintf("RNA_snn_res.%d",res.hi)]] <- clust.x
							}else if(method.clustering=="louvain"){
								seu.x <- FindNeighbors(seu.x, reduction = "pca",
													   #k=if(ncol(seu.x)<500) 5 else 10,
													   #k=if(ncol(seu.x)<500) 5 else 5,
													   k=if(ncol(seu.x)<500) 10 else 10,
													   #k=if(ncol(seu.x)<500) 10 else 20,
													   dims = 1:15,nn.eps=0,force.recalc=T)
								res.hi <- if(ncol(seu.x)<500){ 25 } else res.hi
								seu.x <- FindClusters(seu.x,resolution =res.hi, algorithm=1)
							}

							###### single cell #####
							sce.x <- ssc.build(GetAssayData(seu.x,"scale.data"))
							print(all(colnames(sce.x)==sprintf("%s",colnames(seu.x))))
							sce.x$ClusterID <- sprintf("%s.C%04d",data.id,
										   as.integer(as.character(seu.x@meta.data[[sprintf("RNA_snn_res.%d",res.hi)]])))
							dat.avg <- ssc.average.cell(sce.x,column="ClusterID",ret.type="data.mtx")
							colnames(dat.avg) <- sprintf("%s", (colnames(dat.avg)))
							meta.extra.tb$miniCluster <- sce.x$ClusterID
							return(list("dat.avg"=dat.avg[gene.de.common,],
										"meta.extra.tb"=meta.extra.tb,
										"clust"=structure(sce.x$ClusterID,
												  names=sprintf("%s.%s",data.id,colnames(sce.x)))))
	    },.parallel=ncores)
    return(ret.list)
}

mergeSCEDataFromFileTable <- function(exp.list.table,gene.common,sce.list,group.vec,ncores=6,contamination.vec=NULL)
{
    ret.list <- llply(seq_len(nrow(exp.list.table)),function(i){
							data.id <- exp.list.table$data.id[i]
							##### seu and sce --> seu.x
							sce <- sce.list[[i]]
							colnames(sce) <- sprintf("%s.%s",data.id,colnames(sce))
							if(!is.null(contamination.vec)){
								f.cell.cont <- colnames(sce) %in% contamination.vec
								cat(sprintf("Number of contaminated cells (%s, mergeSCEDataFromFileTable): %d\n",
											data.id,sum(f.cell.cont)))
								sce <- sce[,!f.cell.cont]
							}
							sce$miniCluster <- group.vec[colnames(sce)]
							sce <- ssc.scale(sce,gene.symbol=gene.common,assay.name="norm_exprs",
											 adjB="batchV",do.scale=T)
							dat.avg <- ssc.average.cell(sce,"norm_exprs.scale",
										    column="miniCluster",ret.type="data.mtx")
							rownames(dat.avg) <- rowData(sce)$display.name
							#### pad zero
							{
							    gene.pad <- setdiff(gene.common,rownames(dat.avg))
							    cat(sprintf("For dataset %s, %d genes will be imputed by 0\n",
										    data.id, 
										    length(gene.pad)))
							    if(length(gene.pad)>0){
								    value.pad <-  0
								    mtx.pad <- matrix(rep(value.pad,length(gene.pad)*ncol(dat.avg)),nrow=length(gene.pad))
								    rownames(mtx.pad) <- gene.pad
								    dat.avg <- rbind(dat.avg,mtx.pad)
							    }
							}
							return(dat.avg[gene.common,])
	    },.parallel=ncores)
    return(ret.list)
}

run.inte.metaClust <- function(exp.list.table,
			       out.prefix,
			       gene.exclude.file,
			       ncores=12,npc=15,TH.gene.occ=1,
			       res.hi=50,method.clustering="louvain",
			       contamination.vec=NULL,cor.MALAT1=F,cor.DIG=T)
{
	RhpcBLASctl::omp_set_num_threads(1)
	doParallel::registerDoParallel(cores = ncores)

	env.misc <- loadToEnv(gene.exclude.file)
	cat(sprintf("black gene list:\n"))
	env.misc$all.gene.ignore.df %>% head %>% print

	cat(sprintf("contamination.vec:\n"))
	print(str(contamination.vec))

	loginfo("read data ...")
	seu.list <- llply(seq_len(nrow(exp.list.table)),function(i){
						  seu <- readRDS(exp.list.table$seufile[i])
						  data.id <- exp.list.table$data.id[i]
						  seu$dataset <- data.id
						  seu 
	},.parallel=T)
	names(seu.list) <- exp.list.table$data.id

	sce.list <- llply(seq_len(nrow(exp.list.table)),function(i){
	  sce <- readRDS(exp.list.table$scefile[i])
	  sce$dataset <- exp.list.table$data.id[i]
	  if(!("norm_exprs" %in% assayNames(sce)))
	  {
	      assay(sce,"norm_exprs") <- assay(sce,"log2TPM")
	  }
	  if(exp.list.table$data.id[i] %in% c("HCC.YaoHe10X","HCC.YaoHeSS2") ){
	      f.gene.tmp <- which(rowData(sce)[,"display.name"]=="MALAT1-ENSG00000251562")
	      rowData(sce)[f.gene.tmp,"display.name"] <- "MALAT1"
	  }
	  sce
	  ##assay.name <- "norm_exprs"
	},.parallel=T)
	names(sce.list) <- exp.list.table$data.id

	#### get genes used for PCA. input from  sce.list and gene.de.list
	{

	    gene.occ <- table(unlist(sapply(sce.list,function(x){ unname(rowData(x)[,"display.name"]) })))
	    gene.occ <- gene.occ/length(sce.list)
	    gene.occ <- sort(gene.occ,decreasing=T)
	    gene.common <- names(gene.occ)[gene.occ >= TH.gene.occ]
	    gene.common.all <- names(gene.occ)[gene.occ >= 0.85]
	    loginfo(sprintf("total %d common genes obtain ",length(gene.common)))

	    l_ply(names(sce.list),function(x){
			      ##cat(sprintf("gene %s in dataset %s: %s\n","CX3CR1",x,"CX3CR1" %in% rowData(sce.list[[x]])[,"display.name"]))
			      #cat(sprintf("gene %s in dataset %s: %s\n","MALAT1",x,"MALAT1" %in% rowData(sce.list[[x]])[,"display.name"]))
			      #cat(sprintf("gene %s in dataset %s: %s\n","IL17A",x,"IL17A" %in% rowData(sce.list[[x]])[,"display.name"]))
			      cat(sprintf("gene %s in dataset %s: %s\n","CCR8",x,"CCR8" %in% rowData(sce.list[[x]])[,"display.name"]))
	    },.parallel=T)

	    gene.de.list <- list()
	    for(i in seq_len(nrow(exp.list.table))){
		id.d <- exp.list.table$data.id[i]
		ifile <- exp.list.table$defile[i]
		de.out <- readRDS(ifile)
		gene.de.list[[id.d]] <- de.out$all
		gene.de.list[[id.d]]$geneID <- gene.de.list[[id.d]]$geneSymbol
	    }
	    names(gene.de.list) <- exp.list.table$data.id

	    gene.rank.tb <- as.data.table(ldply(names(gene.de.list),function(x){
						    ret.tb <- unique(gene.de.list[[x]][,c("geneID","F.rank")])
						    ret.tb$dataset.id <- x
						    ret.tb <- ret.tb[geneID %in% gene.common,]
						    return(ret.tb) }))
	    gene.rank.tb <- dcast(gene.rank.tb,geneID~dataset.id,value.var="F.rank",fill=1)
	    gene.rank.tb$median.F.rank <- rowMedians(as.matrix(gene.rank.tb[,-c("geneID"),with=F]))
	    gene.rank.tb <- gene.rank.tb[order(median.F.rank),]
	    ##gene.rank.tb <- gene.rank.tb[geneID %in% gene.common,]
	    #rowData(sce.pb)$median.F.rank <- gene.rank.tb[["median.F.rank"]][match(rownames(sce.pb),gene.rank.tb$geneID)]

	    #### ASH1L MALAT1  CLDN15  ZNF117 NF1 ZNF546 HLA-A HLA-B HLA-C HLA-E PDIA3 IL2RG PTPRC B2M 
	    f.gene.blackList <- (gene.rank.tb$geneID %in% env.misc$all.gene.ignore.df$geneSymbol) |
							    grepl("^RP[LS]",gene.rank.tb$geneID,perl=T) |
							    gene.rank.tb$geneID=="MALAT1"
	    
	    #### select genes
	    #gene.de.common <- head(gene.rank.tb[!f.gene.blackList,][["geneID"]],n=1500)
	    saveRDS(gene.rank.tb[!f.gene.blackList,],sprintf("%s.gene.rank.tb.flt.rds",out.prefix))
	    gene.de.common.tmp.tb <- HVG.From.GeneRankTb(gene.rank.tb[!f.gene.blackList,],
						  n.common=1500,n.specific=1000,th.rank=0.1)
	    write.table(gene.de.common.tmp.tb,
			file=sprintf("%s.gene.de.common.try00.1500.1000.tb",out.prefix),
			row.names=F,sep="\t",quote=F)
	    gene.de.common.tmp.tb <- HVG.From.GeneRankTb(gene.rank.tb[!f.gene.blackList,],
						  ##n.common=1500,n.specific=0,th.rank=0.1)
						  ##n.common=750,n.specific=750,th.rank=0.1)
						  n.common=1000,n.specific=1000,th.rank=0.1)
	    write.table(gene.de.common.tmp.tb,
			file=sprintf("%s.gene.de.common.1000.1000.tb",out.prefix),
			row.names=F,sep="\t",quote=F)

	    gene.de.common <- gene.de.common.tmp.tb$geneID

	    cat(sprintf("RP gene in gene.de.common:\n"))
	    print(gene.de.common[grepl("^RP[LS]",gene.de.common,perl=T)])
		
	}

	#### cellInfo
	{
	    seu.res <- 2
	    meta.tb <- ldply(seq_len(nrow(exp.list.table)),function(i){
				#m.tb <- fread(exp.list.table$cfile[i])
				m.tb <- as.data.table(colData(sce.list[[i]]))
				m.tb$dataset <- exp.list.table$data.id[i]
				a.res <- ifelse(exp.list.table$platform[i]=="SmartSeq2",
								sprintf("RNA_snn_res.%s",seu.res),
								sprintf("SCT_snn_res.%s",seu.res))
				m.tb$ClusterID <- sprintf("%s.%s",m.tb$dataset,m.tb[[a.res]])
				m.tb$dataset <- exp.list.table$data.id[i]
				o.tb <- m.tb[,c("patient","cellID","libraryID","cancerType","loc",
						"batchV","TCR","dataset","ClusterID")]
				##o.tb <- m.tb[,.(NCell=.N),by=c("cancerType","dataset","ClusterID")]
				return(o.tb)
			    })
	    rownames(meta.tb) <- sprintf("%s.%s",meta.tb$dataset,meta.tb$cellID)
	    meta.tb$dataset.tech <- gsub("^.+?\\.","",meta.tb$dataset)
	}

	loginfo("merge data ...")
	dat.merged.list <- mergeDataFromFileTable(exp.list.table,gene.de.common,seu.list,sce.list,res.hi,
						  contamination.vec=contamination.vec,cor.MALAT1=cor.MALAT1,cor.DIG=cor.DIG)
	####dat.merged.list <- mergeDataFromFileTable(exp.list.table[1:2],gene.de.common,seu.list[1:2],sce.list[1:2],res.hi,contamination.vec=contamination.vec)
	
	dat.merged.mtx <- do.call(cbind,llply(dat.merged.list,function(x){ x[["dat.avg"]] }))
	rownames(dat.merged.mtx) <- gsub("_","-",rownames(dat.merged.mtx))
	dim(dat.merged.mtx)
	dat.merged.mtx[1:4,1:5]

	meta.extra.tb <- ldply(dat.merged.list,function(x){ x[["meta.extra.tb"]]})
	rownames(meta.extra.tb) <- meta.extra.tb$cellID.uniq
	f.cell <- intersect(rownames(meta.tb),rownames(meta.extra.tb))
	meta.tb <- cbind(meta.tb[f.cell,],meta.extra.tb[f.cell,])
	cellID2MiniClust.vec <- do.call(c,llply(dat.merged.list,function(x){ x[["clust"]] }))
	###meta.tb$miniCluster.00 <- cellID2MiniClust.vec[rownames(meta.tb)]
	meta.tb <- as.data.table(meta.tb)

	############# miniCluster size
        dat.plot <- meta.tb[,.N,by=c("dataset","miniCluster")]
	p <- ggboxplot(dat.plot,x="dataset",y="N") +
			theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
	ggsave(file=sprintf("%s.miniCluster.size.png",out.prefix),width=10,height=5)

	######################

	#seu.list <- NULL
	#sce.list <- NULL
	#gc()

	loginfo("begin Seurat pipeline ...")
	##### seurat object
	seu.merged <- CreateSeuratObject((dat.merged.mtx),project="panC")
	##seu.merged$dataset.tech <- gsub("^.+?\\.","",seu.merged$dataset)
	seu.merged <- SetAssayData(seu.merged,"scale.data",dat.merged.mtx)
	seu.merged$dataset <- gsub("\\.C\\d\\d\\d\\d$","",colnames(seu.merged))
	seu.merged$dataset.tech <- gsub("^.+?\\.","",seu.merged$dataset)

	seu.merged <- RunPCA(seu.merged,features = rownames(seu.merged), npcs = npc, verbose = FALSE)
	seu.merged <- RunUMAP(seu.merged,reduction="pca",dims=1:npc)
	seu.merged <- RunHarmony(seu.merged, c("dataset"))
	seu.merged <- RunUMAP(seu.merged,reduction = "harmony",reduction.name = "harmony.umap", dims = 1:npc)
	#seu.merged <- RunTSNE(seu.merged,reduction ="pca",tsne.method="FIt-SNE",dims=1:npc,reduction.name="tsne")
	#seu.merged <- RunTSNE(seu.merged,reduction ="harmony",tsne.method="FIt-SNE",dims=1:npc,reduction.name="harmony.tsne")

	loginfo("begin Rtsne ...")
	seu.merged <- RunTSNE(seu.merged,reduction ="pca",tsne.method="Rtsne",dims=1:npc,reduction.name="tsne.Rtsne")
	seu.merged <- RunTSNE(seu.merged,reduction ="harmony",tsne.method="Rtsne",dims=1:npc,reduction.name="harmony.tsne.Rtsne")
	loginfo("end Rtsne ...")

	resolution.vec <- seq(0.1,5,0.1)
	###resolution.vec <- seq(2.5,3.2,0.1)
	#resolution.vec <- seq(3.3,4.0,0.1)
	if(method.clustering=="louvain") {
		seu.merged <- FindNeighbors(seu.merged,reduction = "harmony", dims = 1:npc)
		seu.merged <- FindClusters(seu.merged,resolution = resolution.vec)
		seu.merged <- FindNeighbors(seu.merged,reduction = "pca", dims = 1:npc,graph.name="RNA_pca_snn")
		seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,graph.name="RNA_pca_snn")
	}else if(method.clustering=="leiden") {

		loginfo("begin leiden (harmony) ...")
		clust.res.harmony <- do.call(cbind,llply(resolution.vec,function(x){
							   run.Leiden(Embeddings(seu.merged,"harmony")[,1:npc],SNN.k=20,
										  resolution_parameter=x,
										  myseed=123456)
								   },.parallel=T))
		colnames(clust.res.harmony) <- sprintf("RNA_snn_res.%s",resolution.vec)
		rownames(clust.res.harmony) <- colnames(seu.merged)
		loginfo("end leiden ...")

		loginfo("begin leiden (pca) ...")
		clust.res.pca <- do.call(cbind,llply(resolution.vec,function(x){
							   run.Leiden(Embeddings(seu.merged,"pca")[,1:npc],SNN.k=20,
										  resolution_parameter=x,
										  myseed=123456)
								   },.parallel=T))
		colnames(clust.res.pca) <- sprintf("RNA_pca_snn_res.%s",resolution.vec)
		rownames(clust.res.pca) <- colnames(seu.merged)
		loginfo("end leiden ...")
		seu.merged@meta.data <- cbind(seu.merged@meta.data,clust.res.harmony,clust.res.pca)
		###seu.merged@meta.data <- cbind(seu.merged@meta.data,clust.res.harmony)
	}

	#for(a.res in resolution.vec){
	#	seu.merged@meta.data[[sprintf("RNA_snn_res.%s",a.res)]] <- NULL
	#	colData(sce.merged)[[sprintf("RNA_snn_res.%s",a.res)]] <- NULL
	#}

	#sce.merged <- ssc.build((dat.merged.mtx))
	##### sce.merged contain all common genes
	save(gene.common.all,exp.list.table,cellID2MiniClust.vec,contamination.vec,
	     file=sprintf("%s.data.debug.mergeSCEDataFromFileTable.RData",out.prefix))
	sce.merged <- ssc.build(do.call(cbind, mergeSCEDataFromFileTable(exp.list.table, gene.common.all,sce.list,
									 cellID2MiniClust.vec,
									 ncores=1,
									 contamination.vec=contamination.vec)))

	###colData(sce.merged) <- DataFrame(meta.tb[colnames(sce.merged),])
	rowData(sce.merged)$gene.de.common <- rowData(sce.merged)$display.name %in% gene.de.common

	f.cell <- intersect(colnames(seu.merged),colnames(sce.merged))
	seu.merged <- seu.merged[,f.cell]
	sce.merged <- sce.merged[,colnames(seu.merged)]

	all(colnames(seu.merged)==colnames(sce.merged))
	sce.merged$dataset <- seu.merged$dataset
	sce.merged$dataset.tech <- seu.merged$dataset.tech
	reducedDim(sce.merged,"pca") <- Embeddings(seu.merged,"pca")
	reducedDim(sce.merged,"umap") <- Embeddings(seu.merged,"umap")
	reducedDim(sce.merged,"harmony") <- Embeddings(seu.merged,"harmony")
	reducedDim(sce.merged,"harmony.umap") <- Embeddings(seu.merged,"harmony.umap")
	#reducedDim(sce.merged,"tsne") <- Embeddings(seu.merged,"tsne")
	#reducedDim(sce.merged,"harmony.tsne") <- Embeddings(seu.merged,"harmony.tsne")
	reducedDim(sce.merged,"tsne.Rtsne") <- Embeddings(seu.merged,"tsne.Rtsne")
	reducedDim(sce.merged,"harmony.tsne.Rtsne") <- Embeddings(seu.merged,"harmony.tsne.Rtsne")

	for(a.res in resolution.vec){
		sce.merged[[sprintf("RNA_pca_snn_res.%s",a.res)]] <- seu.merged@meta.data[,sprintf("RNA_pca_snn_res.%s",a.res)]
		sce.merged[[sprintf("RNA_snn_res.%s",a.res)]] <- seu.merged@meta.data[,sprintf("RNA_snn_res.%s",a.res)]
	}

	loginfo("save result ...")
	saveRDS(seu.merged,file=sprintf("%s.seu.merged.rds",out.prefix))
	saveRDS(sce.merged,file=sprintf("%s.sce.merged.rds",out.prefix))
	saveRDS(meta.tb,file=sprintf("%s.meta.tb.rds",out.prefix))
	
	return(list("seu.merged"=seu.merged,"sce.merged"=sce.merged,"meta.tb"=meta.tb))
}


run.inte.metaClust.post <- function(de.limma.tb,out.prefix,ncore=8,
				    min.ncells=30,min.ncellsStudy=200,
				    gset.list=NULL, de.mode="multiAsTwo",column.exp="meanScale",
				    gene.used=NULL,colSet=list())
{
    ######################
    ### cellID libraryID cancerType loc patient batchV TCR majorCluster stype dataset
    #### load all

    #### de genes
    {
	sce.list <- list()
	gene.common <- c()
	gene.de.list <- list()
	for(i in seq_len(nrow(de.limma.tb)))
	{
	    id.d <- de.limma.tb$data.id[i]
	    dfile <- de.limma.tb$dfile[i]
	    de.out <- readRDS(dfile)
	    #### number of cells in each ClusterID
	    #de.out$all[1:2,]
	    if(de.mode=="multi"){
		    ncells.vec <- unlist(unique(as.data.table(de.out$all[,grep("^length\\.",
									       colnames(de.out$all),
									       value=T),with=F])))
		    names(ncells.vec) <- gsub("^length\\.","",names(ncells.vec))
		    ncell.df <- data.frame(ClusterID=names(ncells.vec),
						       ncells=unname(ncells.vec))
	    }else if(de.mode=="multiAsTwo"){
		    ##unique(as.data.table(de.out$all[,c("cluster",grep("^length\\.",colnames(de.out$all),value=T)),with=F]))
		    ncell.df <- unique(as.data.table(de.out$all[,c("cluster","length._case"),with=F]))
		    colnames(ncell.df) <- c("ClusterID","ncells")
		    ncells.vec <- structure(ncell.df$ncells,names=ncell.df$ClusterID)
	    }
	    
	    #### information from de.out$fit
	    geneID.mapping.tb <- unique(de.out$all[,c("geneID","geneSymbol"),with=F])
	    geneID.mapping.vec <- structure(geneID.mapping.tb$geneSymbol,names=geneID.mapping.tb$geneID)
	    es.tb <- as.data.table(ldply(names(ncells.vec),function(group.id){
						     fit <- de.out$fit[[group.id]]
						     n1i <- ncells.vec[group.id]
						     if(de.mode=="multi"){
							    n2i <- mean(ncells.vec[setdiff(names(ncells.vec),group.id)])
							    ES <- effectsize(fit$t,((n1i*n2i)/(n1i+n2i)),(fit$df.prior+fit$df.residual))
							    #### the p value report by limma is two side, cannot use it directly
							    ###### pLimma=2*(pt(abs(fit$t),df=(fit$df.prior+fit$df.residual),lower.tail=F))
							    ###### all( fit$p.value == pLimma[,1])
							    #### convert one-sided P-values to z
							    #zp <- qnorm(pt(fit$t,df=(fit$df.prior+fit$df.residual),lower.tail=F)[,1])
							    zp <- qnorm(pt(-(fit$t),df=(fit$df.prior+fit$df.residual))[,1])
							    #zp <- qnorm(fit$p.value*0.5)
						     }else if(de.mode=="multiAsTwo"){
							    n2i <- sum(ncells.vec[setdiff(names(ncells.vec),group.id)])
							    ES <- effectsize(fit$t[,"II",drop=F],((n1i*n2i)/(n1i+n2i)),(fit$df.prior+fit$df.residual))
							    zp <- qnorm(pt(-(fit$t[,"II",drop=F]),df=(fit$df.prior+fit$df.residual))[,1])
						     }

						     print(all(rownames(ES)==names(zp)))
						     out.tb <- data.table(geneSymbol=geneID.mapping.vec[rownames(ES)],cluster=group.id)
						     out.tb <- cbind(out.tb,ES,zp)
						     return(out.tb)
					    },.parallel=T))

	    #### information from de.out$all
	    gene.de.list[[id.d]] <- de.out$all
	    gene.de.list[[id.d]]$geneID <- gene.de.list[[id.d]]$geneSymbol
	    #if("SNR" %in% colnames(de.out$all)){
	    #	mnames <- c("logFC","t","P.Value","adj.P.Val","meanExp","meanScale","SNR")
	    #}else{
	    #	mnames <- c("logFC","t","P.Value","adj.P.Val","meanExp","meanScale")
	    #}
	    mnames <- c("logFC","t","P.Value","adj.P.Val","meanExp","meanScale","freq._case","freq._control")

	    ### significance
	    dat.d.sig <- de.out$all[,c("geneSymbol","cluster",mnames),with=F]
	    f.sig <- sprintf("%s.%s",dat.d.sig$geneSymbol,dat.d.sig$cluster) %in% 
		    sprintf("%s.%s",de.out$sig[logFC>0,][["geneSymbol"]],de.out$sig[logFC>0,][["cluster"]])
	    dat.d.sig$sig <- 0
	    dat.d.sig$sig[f.sig] <- 1

	    ####
	    dat.d.ext <- merge(dat.d.sig,es.tb,by=c("geneSymbol","cluster"))
	    mnames.ext <- c("d","vard","dprime","vardprime","zp")
	    ####### construct sce
	    dat.d <- llply(c(mnames,"sig",mnames.ext),function(x){
						    dd <- dcast(dat.d.ext,geneSymbol~cluster,value.var=x)
						    dd.mtx <- as.matrix(dd[,-1])
						    rownames(dd.mtx) <- dd[[1]]
						    if(!is.null(gene.used)){
							    gene.pad <- setdiff(gene.used,rownames(dd.mtx))
							    if(length(gene.pad)>0){
								    value.pad <- if(x %in% c("logFC","t","meanExp","meanScale",
														     "freq._case","freq._control",
														     "sig",
														     "d","vard","dprime","vardprime")) 0 else 1
								    mtx.pad <- matrix(rep(value.pad,length(gene.pad)*ncol(dd.mtx)),nrow=length(gene.pad))
								    rownames(mtx.pad) <- gene.pad
								    dd.mtx <- rbind(dd.mtx,mtx.pad)
							    }
							    dd.mtx[gene.used,]
						    }
						    return(dd.mtx)
				       })
	    loginfo("building sce object ...")
	    names(dat.d) <- c(mnames,"sig",mnames.ext)
	    sce.obj <- ssc.build(dat.d[[1]])
	    for(x in c(mnames,"sig",mnames.ext)){
		    assay(sce.obj,x) <- dat.d[[x]]
	    }
	    assay(sce.obj,"exprs") <- assay(sce.obj,column.exp)
	    sce.obj$ClusterID <- colnames(sce.obj)
	    sce.obj$nCellsStudy <- sum(ncells.vec)
	    sce.obj$nCellsCluster <- ncells.vec[sce.obj$ClusterID]
	    ###sce.obj <- sce.obj[,sce.obj$nCellsCluster >= min.ncells]
	    sce.obj <- sce.obj[,sce.obj$nCellsCluster >= min.ncells & sce.obj$nCellsStudy >= min.ncellsStudy]

	    ################
	    sce.list[[id.d]] <- sce.obj
	    if(length(gene.common)==0)
		    ##gene.common <- rownames(sce.list[[id.d]])
		    gene.common <- rowData(sce.list[[id.d]])$display.name
	    else
		    gene.common <- intersect(gene.common,rowData(sce.list[[id.d]])$display.name)
	    cat(sprintf("%d, %s\n",i,de.limma.tb$data.id[i]))
	}
    }

    #### check
    l_ply(names(sce.list),function(aid){
      #print(table(sce.list[[aid]]$ClusterID))
      cat(sprintf("CXCR5 in dataset (%s): %s\n",aid,"CXCR5" %in% rownames(sce.list[[aid]])))
      #rowData(sce.list[[aid]])$gene.common <- rownames(sce.list[[aid]]) %in% gene.common
      #table(rowData(sce.list[[aid]])$gene.common) %>% print
    })

    l_ply(names(gene.de.list),function(aid){
	cat(sprintf("%s, %d genes\n",aid,length(unique(gene.de.list[[aid]]$geneID))))
    })

    #####################################
    sce.pb <- integrate.by.avg(sce.list,sprintf("%s",out.prefix),
				     assay.name=column.exp, is.avg=T,n.downsample=NULL,
				     ncores=ncore,gene.de.list=gene.de.list, de.thres=1000,
				     do.clustering=F,
				     avg.by="ClusterID")

    #### some meta info
    mm <- regexec(".+\\.(CD[48]\\.c\\d+.+)$",colnames(sce.pb),perl=T)
    sce.pb$meta.cluster <- sapply(regmatches(colnames(sce.pb),mm),"[",2)
    sce.pb$cancerType <- sce.pb$dataset.id
    sce.pb$dataset.id <- NULL
    sce.pb$dataset.old <- gsub("\\.CD[48]\\..+$","",sce.pb$ClusterID,perl=T)
    sce.pb$dataset <- sce.pb$dataset.old
    dataset.mapping <- structure(de.limma.tb$platform,names=de.limma.tb$data.id)
    sce.pb$tech <- ifelse(dataset.mapping[sce.pb$dataset.old]=="SmartSeq2","SmartSeq2","Droplet")
    sce.pb$weight.tech <- ifelse(sce.pb$tech=="SmartSeq2",1,0.5)
    #### HCC and LIHC
    sce.pb <- changeSomeNames(sce.pb)
    metadata(sce.pb)$ssc$colSet <- colSet
    sce.pb <- resetSig(sce.pb)
    

    ##### core signature genes
    {
	################# sig genes table
	gene.desc.top <- sscClust:::rank.de.gene(sce.pb,group="meta.cluster",
						 weight.adj="weight.tech", group.2nd="cancerType")

	gene.desc.top[,sig:=comb.padj < 0.01 & comb.ES>0.15]
	gene.desc.top[,sig.cate:="notSig"]
	gene.desc.top[sig==T ,sig.cate:="Tier3"]
	gene.desc.top[sig==T & comb.ES>0.5,sig.cate:="Tier1"]
	gene.desc.top[sig==T & comb.ES<0.5 & comb.padj<1e-20,sig.cate:="Tier2"]

	if(!is.null(gset.list)){
	    for(gset in names(gset.list)){
		    gene.desc.top[[sprintf("geneSet.%s",gset)]] <- gene.desc.top$geneSymbol %in% gset.list[[gset]]
		    rowData(sce.pb)[[sprintf("geneSet.%s",gset)]] <- rownames(sce.pb) %in% gset.list[[gset]]
	    }
	}

	metadata(sce.pb)$ssc$gene.desc.top <- gene.desc.top
	print(head(metadata(sce.pb)$ssc$gene.desc.top,n=4))
	############
	saveRDS(gene.desc.top,file=sprintf("%s.gene.desc.tb.rds",out.prefix))
	saveRDS(sce.pb,file=sprintf("%s.sce.pb.rds",out.prefix))
	##sce.pb <- readRDS(sprintf("%s.sce.pb.rds",out.prefix))
	##gene.desc.top <- readRDS(file=sprintf("%s.gene.desc.tb.rds",out.prefix))

    }

    #### top core genes by median.rank
    {
	ntop <- 50
	sig.prevelance <- 0.4
	gene.top.mrank.tb <- gene.desc.top[median.rank<0.01 & freq.sig >= 0.15,][order(meta.cluster),]
	gene.top.mrank.tb <- gene.top.mrank.tb[,head(.SD,n=ntop),by=c("meta.cluster")]
	if(!is.null(out.prefix)){
		write.table(gene.top.mrank.tb,sprintf("%s.gene.top.mrank.n%d.prev%4.2f.mrank%4.2f.txt",
											  out.prefix,ntop,0.15,0.01),
					row.names=F,sep="\t",quote=F)
	}
	gene.top.mrank.plot.tb <- gene.top.mrank.tb[freq.sig >=sig.prevelance, ][!duplicated(geneID),]

	##### heatmap
	ssc.plot.heatmap(sce.pb,out.prefix=sprintf("%s.gene.top%d.mrank.prev%4.2f.mrank%4.2f",
						   out.prefix,ntop,sig.prevelance,0.01),
			 columns="meta.cluster",columns.order="meta.cluster",colSet=colSet,
			 gene.desc=gene.top.mrank.plot.tb,
			 row.split=gene.top.mrank.plot.tb$Group,
			 column.split=sce.pb$meta.cluster,
			 #border = TRUE,
			 #row_title_rot = 0,
			 par.heatmap=list(row_title_rot=0,border=T,
					  row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
					  column_title_gp = gpar(fontsize = 0)),
			 mytitle=sprintf("Core Genes (top%d.prev%4.2f)",ntop,sig.prevelance),
			 pdf.width=20,pdf.height=10,do.scale=F,
			 palette.name="RdBu",
			 z.lo=-0.5,z.hi=0.5,z.step=0.25,
			 do.clustering.row=F,
			 do.clustering.col=T)

    }
    
    #### top core genes by comb.ES
    {
	ntop <- 50
	es.top <- 0.25
	sig.prevelance <- 0.4
	gene.top.ES.tb <- gene.desc.top[sig==T,][order(meta.cluster,-comb.ES),]
	gene.top.ES.tb <- gene.top.ES.tb[,head(.SD,n=ntop),by=c("meta.cluster")]
	if(!is.null(out.prefix)){
		write.table(gene.top.ES.tb,sprintf("%s.gene.top.ES.n%d.txt",out.prefix,ntop),
					row.names=F,sep="\t",quote=F)
	}
	gene.top.ES.plot.tb <- gene.top.ES.tb[comb.ES >=es.top & freq.sig >= sig.prevelance, ][!duplicated(geneID),]

	sigGeneHeatmap(out.prefix=sprintf("%s.gene.top%d.ES.prev%4.2f.ES%4.2f",
									  out.prefix,ntop,sig.prevelance,es.top),
				   gene.desc.top,sce.pb,gene.to.show.tb=gene.top.ES.plot.tb,
				   colSet=colSet,
				   mytitle=sprintf("Core Genes (top%d.prev%4.2f)", ntop,sig.prevelance))

	cat(sprintf(".... to plot sigGeneHeatmap(..universal.TF..)\n"))
	#### universal TF
	{
		g.tb <- gene.desc.top[sig==T & freq.sig > 0.5 & geneSet.TF==T,][order(meta.cluster,-comb.ES),]
		write.table(g.tb,sprintf("%s.universal.TF.txt",out.prefix),
						row.names=F,sep="\t",quote=F)
		g.tb <- g.tb[!duplicated(geneID),]

		sigGeneHeatmap(out.prefix=sprintf("%s.universal.TF.ht",out.prefix),
				   gene.desc.top,sce.pb,gene.to.show.tb=g.tb,
				   colSet=colSet, mytitle="Universal TF")
	}

    }
    return(sce.pb)
}

make.geneTableLong <- function(obj,out.prefix,th.adj.P=0.01,th.dprime=0.15)
{
    out.list <- (llply(unique(sort(obj[["meta.cluster"]])),function(x){
	    dat.long <- ssc.toLongTable(obj[,obj$meta.cluster==x],gene.id=NULL,
								    assay.name=c("meanExp","meanScale","logFC",
											     "zp","dprime","vardprime",
											     "P.Value","adj.P.Val","sig","freq._case"),
								    col.idx=colnames(colData(obj)))
	    dat.collapse <- sscVis::collapseEffectSizeLong(dat.long,mode.collapse="comb",
										       group.2nd="cancerType",
										       th.adj.P=th.adj.P,th.dprime=th.dprime)
	    dat.collapse$meta.cluster <- x
	    dat.collapse <- dat.collapse[order(dprime),]
	    return(list("dat.long"=dat.long,"dat.long.collapsed"=dat.collapse))
    }))
    out.long <- as.data.table(ldply(out.list,function(x){ x$dat.long  }))
    out.long.collapsed <- as.data.table(ldply(out.list,function(x){ x$dat.long.collapsed  }))
    saveRDS(out.long,sprintf("%s.geneTableLong.rds",out.prefix))
    saveRDS(out.long.collapsed,sprintf("%s.geneTableLong.collapsed.rds",out.prefix))
}



run.clusterProfiler <- function(out.prefix, gene.oi, gene.bg,
								my.title="",db.name="GO",es=NULL,semData=NULL,
								pdf.width=c(9,7,7,18),gset.hl=NULL,...)
{
	require("clusterProfiler")
	require("org.Hs.eg.db")
	require("GOSemSim")

	ret.enrich <- NULL
	if(db.name=="GO")
	{
	    ret.enrich.raw <- clusterProfiler::enrichGO(gene.oi,universe=gene.bg,...)
	    ret.enrich <- clusterProfiler::simplify(clusterProfiler::gofilter(ret.enrich.raw,c(4,5)),
						    cutoff=0.6, by="pvalue", select_fun=min, measure="Jiang",
						    semData=semData)
	}else if(db.name=="wikipathways")
	{
	    wpgmtfile <- "/lustre1/zeminz_pkuhpc/00.database/wikipathways/wikipathways-20200210-gmt-Homo_sapiens.gmt"
	    wp2gene <- read.gmt(wpgmtfile)
	    wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
	    wpid.id = bitr(wp2gene$gene, toType="SYMBOL", fromType="ENTREZID", OrgDb="org.Hs.eg.db")
	    wpid.id.vec <- structure(wpid.id$SYMBOL,names=wpid.id$ENTREZID)
	    wp2gene$symbol <- wpid.id.vec[wp2gene$gene]
	    wpid2gene <- wp2gene %>% dplyr::filter(!is.na(symbol)) %>% dplyr::select(wpid, symbol) #TERM2GENE
	    wpid2name <- wp2gene %>% dplyr::filter(!is.na(symbol)) %>% dplyr::select(wpid, name) #TERM2NAME

	    ret.enrich <- enricher(gene.oi, universe=gene.bg,
						       qvalueCutoff=0.1,
						       TERM2GENE=wpid2gene,TERM2NAME=wpid2name)
	}else if(db.name=="KEGG")
	{
	    eg.oi = bitr(gene.oi, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	    eg.bg = bitr(gene.bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	    ret.enrich <- clusterProfiler::enrichKEGG(eg.oi$ENTREZID, organism = 'hsa',
							     universe=eg.bg$ENTREZID,
							     qvalueCutoff=0.1)
	    ret.enrich <- setReadable(ret.enrich,OrgDb=org.Hs.eg.db,keyType="ENTREZID")

	}else if(db.name=="MKEGG")
	{
	    eg.oi = bitr(gene.oi, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	    eg.bg = bitr(gene.bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	    ret.enrich <- enrichMKEGG(eg.oi$ENTREZID, organism = 'hsa',
							      universe=eg.bg$ENTREZID,
							      qvalueCutoff=0.1)
	    ret.enrich <- setReadable(ret.enrich,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
	}
	##write.table(ret.enrich@result,file=sprintf("%s.enrich.txt",out.prefix),row.names=F,sep="\t",quote=F)
	write.table(cbind(data.table(aid=my.title),as.data.frame(ret.enrich)),
				file=sprintf("%s.enrich.txt",out.prefix),row.names=F,sep="\t",quote=F)

	if(nrow(as.data.frame(ret.enrich))==0){
	    return(list("plot"=NULL,
			"enrichResult"=ret.enrich))
	}

	p1 <- dotplot(ret.enrich, showCategory=30) +
		ggtitle(sprintf("%s",my.title)) +
		scale_colour_distiller(palette = "Reds",breaks=c(0,0.01,0.05,0.1),limits=c(0,0.1),na.value="lightgray")
	ggsave(sprintf("%s.dotplot.pdf",out.prefix),width=pdf.width[1],height=7,useDingbats=F)

	p2 <- emapplot(ret.enrich, pie_scale=1.5,layout="kk") +
		scale_colour_distiller(palette = "Reds",breaks=c(0,0.01,0.05,0.1),limits=c(0,0.1),na.value="lightgray")
	ggsave(sprintf("%s.emaplot.pdf",out.prefix),width=pdf.width[2],height=7,useDingbats=F)


	if(!is.null(es))
	{
	    if(!is.null(gset.hl)){
		    ret.enrich.hl <- clusterProfiler.dplyr::filter(ret.enrich,Description %in% gset.hl)
		    rownames(ret.enrich.hl@result) <- ret.enrich.hl@result$ID
	    }else{
		    ret.enrich.hl <- ret.enrich
	    }
	    p3 <- cnetplot(ret.enrich.hl, foldChange=es) +
		    scale_colour_gradient(name = "ES", low="blue",high="red",
							      na.value = "#E5C494",
							      limits=c(0,2),breaks=c(0.15,0.5,1,2),
							      guide = guide_colorbar(title = "ES"))
	    ggsave(sprintf("%s.cnet.pdf",out.prefix),width=pdf.width[3],height=7,useDingbats=F)

	    options(stringsAsFactors=T)
	    p4 <- clusterProfiler::heatplot(ret.enrich, foldChange=es) +
			    theme(axis.text.x=element_text(size=3)) +
			    geom_tile(aes_(fill = ~foldChange), color = NA) +
			    scale_fill_continuous(low="blue", high="red", name = "ES",
								      breaks=c(0.15,0.5,1,2),limits=c(0,2),
								      na.value="lightgray",
								      guide = guide_colorbar(title = "ES"))
	    ggsave(sprintf("%s.heatmap.pdf",out.prefix),width=pdf.width[4],height=6)
	    options(stringsAsFactors=F)
	}
	return(list("plot"=list("p1"=p1,"p2"=p2,"p3"=p3,"p4"=p4),
				"enrichResult"=ret.enrich))
}



run.gene.analysis.CmpCancerType <- function(gene.desc.top,out.prefix,mcls="CD8.c12.Tex.CXCL13",
					    gset.cmp=c("sig","universal","diverse","loFreqHiES"),
					    plot.type="png",vector.friendly=T,my.parallel=T,
					    th.geneRatio=c(0,0.1),th.qvalue=0.1,th.showCategory=30,font.size.cmp=9,
					    pdf.width.dotplot=NULL,verbose=0,gene.long.collapsed.tb=NULL)
{
	library("org.Hs.eg.db")
	library("clusterProfiler")
	library("GOSemSim")
	dir.create(dirname(out.prefix),F,T)

	dat.plot <- gene.desc.top[meta.cluster==mcls,]

	makeFig.ES.Freq <- function(dat.plot,out.prefix,
				    col.x="dprime.max.group.2nd",col.y="N.sig.group.2nd",col.N.total="N.total.group.2nd",
				    my.xlab="Maximum of Effect Size",my.ylab="Significant in # Cancer Types",my.title="",
				    limits.x=NULL,th.N.sig=0,vline=c(0.15,0.25))
	{
		dat.plot[[col.x]][ dat.plot[[col.x]] > 2 ] <- 2
		### filter out genes significant in less than th.N.sig cancer types
		dat.plot <- dat.plot[ dat.plot[[col.y]] >= th.N.sig,]
		
		if(is.null(limits.x)){
			limits.x <- pretty(dat.plot[[col.x]])
			limits.x <- limits.x[c(1,length(limits.x))]
		}
		print(limits.x)
		nbins <- max(dat.plot[[col.N.total]])
		th.y <- ceiling(nbins/2)
		th.x <- vline[1]
		### cell numbers
		text.cellNumber <- c(sum(dat.plot[[col.x]] <= th.x & dat.plot[[col.y]] >= th.y),
							 sum(dat.plot[[col.x]] > th.x & dat.plot[[col.y]] >= th.y),
							 sum(dat.plot[[col.x]] > th.x & dat.plot[[col.y]] < th.y),
							 sum(dat.plot[[col.x]] <= th.x & dat.plot[[col.y]] < th.y))
		cat("text.cellNumber:\n")
		print(text.cellNumber)
		print(th.x)
		print(th.y)
		print(limits.x)
		###breaks.x=seq(-0.45,4.5,0.05)
		breaks.x=seq(limits.x[1],limits.x[2],0.05)
		breaks.y <- 0:(nbins+1)
		x.tick.pos.text <- unique(c(0,vline,1,2))
		x.tick.pos <- ((x.tick.pos.text-limits.x[1])/0.05+1)
		### filter out genes significant in less than th.N.sig cancer types
		#dat.plot <- dat.plot[ dat.plot[[col.y]] >= th.N.sig,]

		dat.plot[["x.idx"]] <- dat.plot[[col.x]]/0.05
		dat.plot[["bin_x"]] <- cut(dat.plot[[col.x]],breaks=breaks.x,include.lowest=T)
		dat.plot[["bin_y"]] <- cut(dat.plot[[col.y]],breaks=breaks.y,include.lowest=T,right=F)
		dat.plot %>% dplyr::count(bin_x, bin_y) %>%
			 tidyr::complete(bin_x, bin_y, fill = list(n = 0)) -> dat.plot.tmp
		
		p <- ggplot(dat.plot.tmp,aes(bin_x, bin_y, fill = log10(n+1))) + geom_tile() +
				theme_minimal() +
				theme(axis.text.x=element_text(angle=45,hjust=1,size=3)) +
				scale_fill_continuous(type = "viridis") +
				scale_x_discrete(expand=c(0,0)) +
				scale_y_discrete(expand=c(0,0))
		##ggsave(sprintf("%s.test.01.png",out.prefix),width=8,height=7)

		tmpfilename <- sprintf("%s.tmp.%s.CmpCancerType.png",out.prefix,"a")
		ggsave(filename=tmpfilename,
			   plot = p + theme_void() + theme(legend.position = "none",
							   axis.line.x = element_blank(), axis.line.y = element_blank(),
							   axis.title.x = element_blank(), axis.title.y = element_blank(),
							   axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
							   plot.title = element_blank()),
						 width=7,height=6)
		pbuild.params <- ggplot_build(plot = p)$layout$panel_params[[1]]
		range.values <- c( pbuild.params$x.range, pbuild.params$y.range)
		img <- png::readPNG(source = tmpfilename)
		blank <- ggplot(data = dat.plot,mapping = aes_string(x="x.idx",y=col.y)) +
					geom_blank() +
					theme_classic() +
					labs(title=my.title,x=my.xlab,y=my.ylab) +
					theme(plot.title = element_text(hjust = 0.5)) +
					coord_cartesian(xlim = range.values[1:2],
									#ylim = range.values[3:4],
									ylim = c(0,nbins+1),
									expand = F)
		blank <- blank + annotation_raster(raster = img,
						   xmin = range.values[1], xmax = range.values[2],
						   ymin = 0, ymax = nbins+1)
		blank <- blank + geom_vline(xintercept=((vline-breaks.x[1])/0.05+1)-0.5,
							   linetype="dashed",alpha=0.8,color="lightgray") +
					geom_hline(yintercept=c(ceiling(nbins/2)),linetype="dashed",alpha=0.8,color="lightgray") +
					geom_text(data=data.table(x=c(1,(th.x-limits.x[1])/0.05+1,
								      (th.x-limits.x[1])/0.05+1, 1),
								  y=c(nbins,nbins,floor(nbins/2),floor(nbins/2))),
						  aes(x=x,y=y,label=text.cellNumber),
						  hjust=0,
						  color="white") +
					scale_x_continuous(breaks=x.tick.pos-0.5,labels=x.tick.pos.text) +
					scale_y_continuous(breaks=seq(1,nbins+1,5)-0.5,labels=seq(1,nbins+1,5)-1)
		legend.blank <- cowplot::get_legend(p)
		pp <- cowplot::plot_grid(blank, legend.blank, rel_widths = c(4, 1))
		##ggsave(sprintf("%s.test.03.pdf",out.prefix),width=5,height=4)
		ggsave(sprintf("%s.pdf",out.prefix),width=5,height=4)
		
		file.remove(tmpfilename)

	}

	###makeFig.ES.Freq(dat.plot,sprintf("%s.maxES.freq.NSig0",out.prefix),my.title=mcls,th.N.sig=0)
	makeFig.ES.Freq(dat.plot,sprintf("%s.maxES.freq.NSig1",out.prefix),my.title=mcls,th.N.sig=1,vline=0.25)
	makeFig.ES.Freq(dat.plot,sprintf("%s.combES.freq",out.prefix),my.title=mcls,
					my.xlab="Combined Effect Size",col.x="comb.ES",vline=0.15)

	##return(NULL)

	{
#		p <- ggplot(dat.plot,aes(comb.ES,dprime.max.group.2nd)) +
#				geom_point() +
#				geom_hline(yintercept=c(0.15,0.25),linetype="dashed",color="lightgray") +
#				geom_vline(xintercept=c(0.15,0.25),linetype="dashed",color="lightgray") +
#				stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
#				scale_fill_viridis_c()+
#				theme_pubr()
#		ggsave(sprintf("%s.ES.comb.max.pdf",out.prefix),width=5,height=4)
	}


	f.sig <- with(dat.plot,sig==T)
	f.universal <- with(dat.plot,freq.sig.group.2nd > 0.5)
	f.diverse <- with(dat.plot,N.sig.group.2nd > 0 & freq.sig.group.2nd < 0.5 )
	#f.diverse <- with(dat.plot,N.sig.group.2nd > 1 & freq.sig.group.2nd < 0.5 )
	f.loFreqHiES <- with(dat.plot,freq.sig.group.2nd > 0 & freq.sig.group.2nd < 0.5 &
					  ##dprime.max.group.2nd > 0.5 & dprime.max.adj.P.group.2nd < 0.01)
					  dprime.max.group.2nd > 0.25 & dprime.max.adj.P.group.2nd < 0.01)

	write.table(dat.plot[f.sig,],
				file=sprintf("%s.CmpCancerType.%s.sig.txt",out.prefix,mcls),
				sep="\t",quote=F,row.names=F)
	write.table(dat.plot[f.universal,],
				file=sprintf("%s.CmpCancerType.%s.universal.txt",out.prefix,mcls),
				sep="\t",quote=F,row.names=F)
	write.table(dat.plot[f.diverse,],
				file=sprintf("%s.CmpCancerType.%s.diverse.txt",out.prefix,mcls),
				sep="\t",quote=F,row.names=F)
	write.table(dat.plot[f.loFreqHiES,],
				file=sprintf("%s.CmpCancerType.%s.loFreqHiES.txt",out.prefix,mcls),
				sep="\t",quote=F,row.names=F)

	out.prefix.plot <- sprintf("%s/clusterProfiler/%s",dirname(out.prefix),basename(out.prefix))
	dir.create(dirname(out.prefix.plot),F,T)

	gene.list.test <- list("sig"=dat.plot[f.sig,][["geneID"]],
						   "universal"=dat.plot[f.universal,][["geneID"]],
						   "diverse"=dat.plot[f.diverse,][["geneID"]],
						   "loFreqHiES"=dat.plot[f.loFreqHiES,][["geneID"]])
	gene.list.bg <- dat.plot$geneID

	hsGO <- llply(c("CC","MF","BP"),function(x){ GOSemSim::godata('org.Hs.eg.db', ont=x) },
				  .parallel=my.parallel,
				  .paropts=list(.packages=c("GOSemSim","clusterProfiler","AnnotationDbi"),
								.export="org.Hs.eg.db"))
	names(hsGO) <- c("CC","MF","BP")
	cat(sprintf("good here\n"))

	es <- list()
	es[["sig"]] <- structure(dat.plot[f.sig,][["comb.ES"]],
					names=dat.plot[f.sig,][["geneID"]])
	es[["sig"]][es[["sig"]]>2] <- 2
	es[["universal"]] <- structure(dat.plot[f.universal,][["comb.ES"]],
					names=dat.plot[f.universal,][["geneID"]])
	es[["universal"]][es[["diverse"]]>2] <- 2
	es[["diverse"]] <- structure(dat.plot[f.diverse,][["dprime.max.group.2nd"]],
					names=dat.plot[f.diverse,][["geneID"]])
	es[["diverse"]][es[["diverse"]]>2] <- 2
	es[["loFreqHiES"]] <- structure(dat.plot[f.loFreqHiES,][["dprime.max.group.2nd"]],
					names=dat.plot[f.loFreqHiES,][["geneID"]])
	es[["loFreqHiES"]][es[["loFreqHiES"]]>2] <- 2

	if(verbose==1)
	{
	    out.enrich.go <- list()
	    for(gene.cate in c("sig","universal","diverse","loFreqHiES"))
	    {
		out.enrich.go[[gene.cate]] <- llply(c("CC","MF","BP"),function(x){
				run.clusterProfiler(sprintf("%s.%s.GO.%s",out.prefix.plot,gene.cate,x),
									gene.list.test[[gene.cate]],
									gene.list.bg,
									my.title=sprintf("%s (%s)",x,gene.cate),
									db.name="GO",qvalueCutoff=0.1,
									OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont=x,semData=hsGO[[x]],
									es=es[[gene.cate]],
									pdf.width=c(9,7,7,18))
					},.parallel=my.parallel,
				.paropts=list(.packages=c("clusterProfiler","org.Hs.eg.db","GOSemSim","AnnotationDbi")))
		names(out.enrich.go[[gene.cate]]) <- c("CC","MF","BP")
	    }
	    saveRDS(out.enrich.go,file=sprintf("%s.out.enrich.GO.rds",out.prefix.plot))
	    cat(sprintf("enrich.GO done.\n"))

	    out.enrich.pathway <- list()
	    for(gene.cate in c("sig","universal","diverse","loFreqHiES"))
	    {
		out.enrich.pathway[[gene.cate]] <- llply(c("wikipathways","KEGG","MKEGG"),function(x){
					run.clusterProfiler(sprintf("%s.%s.%s",out.prefix.plot,gene.cate,x),
									gene.list.test[[gene.cate]],
									gene.list.bg,
									my.title=sprintf("%s (%s)",x,gene.cate),
									db.name=x,es=es[[gene.cate]],
									pdf.width=c(9,7,7,18),gset.hl=NULL)
					},.parallel=my.parallel)
		names(out.enrich.pathway[[gene.cate]]) <- c("wikipathways","KEGG","MKEGG")
	    }
	    saveRDS(out.enrich.pathway,file=sprintf("%s.out.enrich.pathway.rds",out.prefix.plot))
	}

	######### GO comparison #######
	out.cmp.go <- llply(c("CC","MF","BP"),function(x){
		ck <- compareCluster(geneCluster = gene.list.test[gset.cmp], fun = "enrichGO",
							 keyType="SYMBOL",ont=x,universe=gene.list.bg,
							 qvalueCutoff=th.qvalue, OrgDb='org.Hs.eg.db')
		ck.slim <- clusterProfiler::simplify(gofilter(ck,c(4,5)),cutoff=0.6, by="pvalue", select_fun=min,
							measure="Jiang",semData=hsGO[[x]])
		ck.slim.fortify <- fortify(ck.slim, showCategory=30, by="geneRatio",
								   includeAll=TRUE, split=NULL)
		###write.table(ck.slim@compareClusterResult,file=sprintf("%s.cmp.GO.%s.txt",out.prefix.plot,x),
		###			row.names=F,sep="\t",quote=F)
		#####
		ck.tmp <- compareCluster(geneCluster = gene.list.test[gset.cmp], fun = "enrichGO",
							 keyType="SYMBOL",ont=x,universe=gene.list.bg,
							 qvalueCutoff=Inf,pvalueCutoff=Inf, OrgDb='org.Hs.eg.db')
		enrich.ref <- clusterProfiler::enrichGO(gene.list.test$sig,OrgDb=org.Hs.eg.db,
							universe=gene.list.bg,keyType="SYMBOL",ont=x,
							pvalueCutoff=Inf,qvalueCutoff=Inf)
		enrich.ref.tmp <- clusterProfiler.dplyr::filter(enrich.ref,
								as.character(Description) %in%
								    as.character(ck.slim.fortify$Description))
		enrich.ref.tmp <- clusterProfiler::simplify(gofilter(enrich.ref.tmp,c(4,5)),
							    cutoff=0.6,by="pvalue",select_fun=min,
							    measure="Jiang",semData=hsGO[[x]])
		ck.slim.plot <- clusterProfiler.dplyr::filter(ck.tmp,as.character(Description) %in%
								as.character(enrich.ref.tmp@result$Description) &
							      p.adjust < 0.05 & qvalue<th.qvalue)
		write.table(ck.slim.plot@compareClusterResult,file=sprintf("%s.cmp.GO.%s.txt",out.prefix.plot,x),
					row.names=F,sep="\t",quote=F)
		p1 <- dotplot(ck.slim.plot, showCategory=th.showCategory,font.size=font.size.cmp) +
			ggtitle(x) +
			scale_size(breaks=c(0.025,0.05,0.075),range=c(0.2,6),labels=c("2.5%","5%","7.5%"),
					   limits=th.geneRatio) +
			scale_colour_distiller(palette = "Reds",breaks=c(0,0.01,0.05,0.1),limits=c(0,0.1),na.value="lightgray")
		ggsave(sprintf("%s.cmp.GO.%s.dotplot.pdf",out.prefix.plot,x),
			   width=if(!is.null(pdf.width.dotplot)) pdf.width.dotplot else { if(length(gset.cmp)>3) 7 else 9 },
			   height=7,useDingbats=F)
		####
		#p1 <- dotplot(ck.slim, showCategory=30,font.size=8) +
		#	ggtitle(x) +
		#	scale_colour_distiller(palette = "Reds",breaks=c(0,0.01,0.05,0.1),limits=c(0,0.1),na.value="lightgray")
		#ggsave(sprintf("%s.cmp.GO.%s.dotplot.pdf",out.prefix.plot,x),width=9,height=7)

		p2 <- emapplot(ck.slim,showCategory=th.showCategory,pie="count", pie_scale=1.0, layout="kk")
		ggsave(sprintf("%s.cmp.GO.%s.emaplot.pdf",out.prefix.plot,x),width=9,height=7)

		return(list("plot"=list("p1"=p1,"p2"=p2),
				"cmpResult"=ck.slim.plot))
	},.parallel=my.parallel)
	names(out.cmp.go) <- c("CC","MF","BP")
	saveRDS(out.cmp.go,file=sprintf("%s.out.cmp.GO.rds",out.prefix.plot))
	pp <- plot_grid(out.cmp.go$BP$plot$p1+theme(legend.position = "none"),
					out.cmp.go$MF$plot$p1+theme(legend.position = "none"),
					out.cmp.go$CC$plot$p1+theme(legend.position = "none"),
					get_legend(out.cmp.go$BP$plot$p1),
					#rel_widths=c(1,1,1,0.3),
					align="hv",nrow=2)
	ggsave(sprintf("%s.cmp.GO.merge.dotplot.pdf",out.prefix.plot),width=13,height=14,useDingbats=F)
	pp <- plot_grid(out.cmp.go$MF$plot$p1+theme(legend.position = "none"),
					out.cmp.go$CC$plot$p1+theme(legend.position = "none"),
					get_legend(out.cmp.go$BP$plot$p1),
					rel_widths=c(1,1,0.3),
					align="hv",nrow=1)
	ggsave(sprintf("%s.cmp.GO.merge.n2.dotplot.pdf",out.prefix.plot),width=15,height=7,useDingbats=F)
	############################



}


run.gene.specificity.byPattern <- function(gene.desc.top,gene.long.collapsed.tb,
										   out.prefix,
										   mcls="CD8.c12.Tex.CXCL13",
										   gene.highlight=NULL,gene.highlight.freq=0.01,
										   col.freq.sig="freq.sig.group.2nd",
										   col.ES.max="dprime.max.group.2nd",
										   pdf.width=6.0, pdf.height=8.0,my.parallel=F,verbose=1)
{
	gene.loFreqHiES.tb <- gene.desc.top[ gene.desc.top[["meta.cluster"]]==mcls &
									 gene.desc.top[[col.freq.sig]] < 0.5 &
									 gene.desc.top[[col.freq.sig]] > 0 &
									 gene.desc.top[[col.ES.max]] > 0.25 ,]
	#####  cancer type specific gene
	data.ana <- gene.long.collapsed.tb[meta.cluster==mcls,]
	###
	f.gene.00 <- data.ana[,.(aa=sum(.SD$freq._case)),by=c("meta.cluster","geneID")][aa==0,][["geneID"]]
	data.ana <- data.ana[!(geneID %in% f.gene.00),]
	###
	data.ana.debug <<- data.ana
	data.ana.max <- data.ana[,{
							    ### pattern by ES
								ES.v.order <- order(.SD$dprime,decreasing=T)
								ES.i.max <- ES.v.order[1]
								ES.i.2nd <- ES.v.order[2]
								ES.pattern <- rep(0,.N)
								ES.pattern[ES.i.max] <- 1
								res.cor.ES <- cor.test(ES.pattern,.SD$dprime,method="pearson")
								### pattern by frequency
								freq.v.order <- order(.SD$freq._case,decreasing=T)
								freq.i.max <- freq.v.order[1]
								freq.i.2nd <- freq.v.order[2]
								freq.pattern <- rep(0,.N)
								freq.pattern[freq.i.max] <- 1
								res.JSD <- 1-philentropy::JSD(rbind(.SD$freq._case/sum(.SD$freq._case),
																	freq.pattern),unit = "log")^0.5
								res.cor.freq <- cor.test(freq.pattern,.SD$freq._case,method="pearson")
								.(
								  sig.freq=sum(.SD$sig)/.N,
								  ES.max=.SD$dprime[ES.i.max],
								  ES.2nd=.SD$dprime[ES.i.2nd],
								  ES.cancerType=.SD$cancerType[ES.i.max],
								  ES.max.freq=.SD$freq._case[ES.i.max],
								  ES.max.sig=.SD$sig[ES.i.max],
								  ES.max.adj.P.Val=.SD$adj.P.Val[ES.i.max],
								  ##ES.diff=.SD$dprime[i.max]-.SD$dprime[ES.i.2nd],
								  ES.specificity=res.cor.ES$estimate,
								  ES.specificity.p=res.cor.ES$p.value,
								  freq.max=.SD$freq._case[freq.i.max],
								  freq.2nd=.SD$freq._case[freq.i.2nd],
								  freq.cancerType=.SD$cancerType[freq.i.max],
								  freq.max.ES=.SD$dprime[freq.i.max],
								  freq.max.sig=.SD$sig[freq.i.max],
								  freq.max.adj.P.Val=.SD$adj.P.Val[freq.i.max],
								  freq.specificity=res.cor.freq$estimate,
								  freq.specificity.p=res.cor.freq$p.value,
								  JSD.specificity=res.JSD
								  )},by=c("meta.cluster","geneID")][order(-ES.max), ]
	data.ana.max[,ES.specificity.FDR:=p.adjust(ES.specificity.p,method="BH")]
	data.ana.max[,freq.specificity.FDR:=p.adjust(freq.specificity.p,method="BH")]
	
	print("ha?")
	#### using ES.specificity
	###data.ana.max.flt <- data.ana.max[ES.max.sig==1 & sig.freq<0.5 & ES.specificity.FDR<0.05 & ES.max>0.25 & ES.max.freq > 0.01,
	data.ana.max.flt <- data.ana.max[ES.max.sig==1 & sig.freq<0.5 & ES.specificity.FDR<0.05 & ES.max>0.25,
									 ][order(ES.cancerType,-ES.specificity),]
	print(data.ana.max.flt[,table(ES.cancerType)])
	#### using JSD.specificity
	#data.ana.max[,quantile(JSD.specificity,c(0.75,0.9,0.95,0.99))]
	#th.JSD <- quantile(data.ana.max$JSD.specificity,0.75)
	#data.ana.max.flt <- data.ana.max[JSD.specificity > th.JSD &
	#								 freq.max.ES > 0.15 &
	#								 0 < sig.freq & sig.freq<0.5 & freq.max > 0.01,
	#								 ][order(freq.cancerType,-JSD.specificity),]
	#print(data.ana.max.flt[,table(freq.cancerType)])

	#### plot heatmap
	data.plot <- dcast(gene.long.collapsed.tb[meta.cluster==mcls & geneID %in% data.ana.max.flt$geneID,],
					   geneID~cancerType,value.var="dprime")
	data.plot.mtx <- as.matrix(data.plot[,-1])
	rownames(data.plot.mtx) <- data.plot[[1]]
	data.plot.mtx <- data.plot.mtx[data.ana.max.flt$geneID,]

	data.ana.max.flt.debug <<- data.ana.max.flt
	print(gene.highlight.freq)
	print("good")

	if(is.null(gene.highlight)){
		gene.highlight <- data.ana.max.flt[,.(geneID=head(.SD$geneID,n=max(1,as.integer(.N*gene.highlight.freq)))),
										   by="ES.cancerType"][["geneID"]]
	}
	print("good too")
	#gene.highlight <- c("RORA","TSHZ1","IL21","KLRD1","CCL20","MKI67","FOXP3")
	gene.highlight.at <- match(gene.highlight,rownames(data.plot.mtx))
	gene.highlight.at <- gene.highlight.at[!is.na(gene.highlight.at)]
	ha <- rowAnnotation(geneH = anno_mark(at = gene.highlight.at,labels = rownames(data.plot.mtx)[gene.highlight.at]))
	sscVis::plotMatrix.simple(data.plot.mtx,
								  out.prefix=sprintf("%s.specificity.ES",out.prefix),
								  mytitle=sprintf("%s",mcls),
								  palatte=rev(brewer.pal(n = 7,name = "RdBu")),
								  ##palatte=viridis::magma(7),
								  z.lo=-0.5,z.hi=0.5,
								  par.legend=list(at = seq(-0.5,0.5,0.25)),
								  pdf.width=pdf.width,
								  pdf.height=pdf.height,
								  clust.row=F,
								  par.heatmap=list(cex.column=0.8,cex.row=0,
													right_annotation =  ha,
													row_gap = unit(0, "mm"),row_title_rot=0,
												    row_split=data.ana.max.flt$ES.cancerType,border=T
												   ),
								  exp.name="ES")

	if(verbose==0){return(NULL)}

	###
	gene.loFreqHiES.tb[,cancerTypeSpecific:=geneID %in% data.ana.max.flt$geneID]

	#### GO enrich (no specific, each cancer type's )
	gene.list.test <- list()
	gene.list.test[["NoSpe"]] <- structure(gene.loFreqHiES.tb[cancerTypeSpecific==F,][[col.ES.max]],
										   names=gene.loFreqHiES.tb[cancerTypeSpecific==F,][["geneID"]])
	gene.list.test[["Spe"]] <- structure(gene.loFreqHiES.tb[cancerTypeSpecific==T,][[col.ES.max]],
										 names=gene.loFreqHiES.tb[cancerTypeSpecific==T,][["geneID"]])
	for(xtype in sort(unique(data.ana.max.flt$ES.cancerType))){
		gene.x <- structure(data.ana.max.flt[ES.cancerType==xtype,][["ES.max"]],
							names=data.ana.max.flt[ES.cancerType==xtype,][["geneID"]])
		if(length(gene.x)>10){
			gene.list.test[[xtype]] <- gene.x
		}
	}
	gene.list.bg <- unique(gene.desc.top[meta.cluster==mcls,][["geneID"]])
	####gene.list.bg <- gene.loFreqHiES.tb[["geneID"]]

	hsGO <- llply(c("CC","MF","BP"),function(x){ GOSemSim::godata('org.Hs.eg.db', ont=x) },
				  .parallel=my.parallel,
				  .paropts=list(.packages=c("GOSemSim","clusterProfiler","AnnotationDbi"),
								.export="org.Hs.eg.db"))
	names(hsGO) <- c("CC","MF","BP")
	cat(sprintf("good here\n"))


	out.prefix.plot <- sprintf("%s/clusterProfiler/%s",dirname(out.prefix),basename(out.prefix))
	dir.create(dirname(out.prefix.plot),F,T)
	
	out.enrich.go <- list()
	for(gene.cate in names(gene.list.test))
	{
	    out.enrich.go[[gene.cate]] <- llply(c("CC","MF","BP"),function(x){
			    es.x <- gene.list.test[[gene.cate]]
			    es.x[es.x > 2 ] <- 2
			    run.clusterProfiler(sprintf("%s.%s.loFreqHiES.GO.%s",out.prefix.plot,gene.cate,x),
								    names(gene.list.test[[gene.cate]]),
								    gene.list.bg,
								    my.title=sprintf("%s (%s)",x,gene.cate),
								    db.name="GO",qvalueCutoff=0.1,
								    OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont=x,semData=hsGO[[x]],
								    es=es.x,
								    pdf.width=c(9,7,7,18))
				    },.parallel=my.parallel,
			    .paropts=list(.packages=c("clusterProfiler","org.Hs.eg.db","GOSemSim","AnnotationDbi")))
	    names(out.enrich.go[[gene.cate]]) <- c("CC","MF","BP")
	}
	saveRDS(out.enrich.go,file=sprintf("%s.out.enrich.GO.loFreqHiES.rds",out.prefix.plot))
	#out.enrich.go <- readRDS(file=sprintf("%s.out.enrich.GO.loFreqHiES.rds",out.prefix.plot))
	cat(sprintf("enrich.GO done.\n"))

	####

	gset.list <- fetechGeneSetList(gmode=1)
	gset.tb <- as.data.table(ldply(names(gset.list),function(x){
									   data.table(term=x,gene=gset.list[[x]])
									}))
	res.enricher <- list()
	res.enricher[["NoSpe"]] <- enricher(names(gene.list.test[["NoSpe"]]), TERM2GENE = gset.tb,qvalueCutoff=0.1)
	res.enricher[["Spe"]] <- enricher(names(gene.list.test[["Spe"]]), TERM2GENE = gset.tb,qvalueCutoff=0.1)

	data.save <- list("gene.loFreqHiES.tb"=gene.loFreqHiES.tb,
					  "gene.cancerTypeSpecific.tb"=data.ana.max.flt,
					  "res.enricher"=res.enricher,
					  "gene.list.test"=gene.list.test)
	####
	return(data.save)
}

run.GRN.cor <- function(sce,out.prefix,do.parallel=F)
{

	tic("cor.BLAS (pearson)")
	dat.cor.p <- sscVis:::cor.BLAS(assay(sce))
	toc()
	tic("cor.BLAS (spearman)")
	dat.cor.s <- sscVis:::cor.BLAS(assay(sce),method="spearman")
	toc()

	## determine the cut off automatically
	dat.cor.uniq.tb <- rbind(data.table(corr=dat.cor.p[upper.tri(dat.cor.p)],cate="pearson"),
							 data.table(corr=dat.cor.s[upper.tri(dat.cor.s)],cate="spearman"))
	cutoff.list <- list()
	cutoff.list[["pearson"]] <- dat.cor.uniq.tb[cate=="pearson",quantile(corr,c(0.99,0.995,0.998,0.999))]
	cutoff.list[["spearman"]] <- dat.cor.uniq.tb[cate=="spearman",quantile(corr,c(0.99,0.995,0.998,0.999))]
	print(cutoff.list)

	p.list <- llply(names(cutoff.list),function(cor.type){
				p <- ggdensity(dat.cor.uniq.tb[cate==cor.type,],x="corr") +
					coord_cartesian(xlim = c(-1,1), expand = TRUE) +
					geom_vline(xintercept=c(0.35,0.3,0.25,0.2,-0.35,-0.3,-0.25,-0.2),
							   color="green", linetype="dashed", size=1)
					p <- p + geom_vline(xintercept=c(cutoff.list[[cor.type]],-cutoff.list[[cor.type]]),
										color="blue", linetype="dashed", size=1)
							 },.parallel=do.parallel)
	p <- cowplot::plot_grid(plotlist=p.list,ncol = 1,align = "hv")
	ggsave(sprintf("%s.corr.pdf",out.prefix),width=5,height=6)

	ret.list <- list("dat.cor.p"=dat.cor.p,
					 "dat.cor.s"=dat.cor.s)
	return(ret.list)

}


run.GRN.p <- function(dat.list,gene.desc.tb,out.prefix,gset.list=fetechGeneSetList(),
					   gene.selection=NULL,th.es=0.25,
					   g.layout=NULL,
					   th.cor=0.35,
					   th.compsize=NULL,
					   gene.highlight=NULL,
					   gene.cluster=NULL,
					   par.node.text=list(nudge_x=0.5,nudge_y=0.5,
										  repel=T,force=1,max.iter=1000,box.padding=0.15,
										  size=1.5),
					   par.leiden=list(resolution_parameter=1),
					   do.polish=F,plot.type="pdf",verbose=1,
					   pdf.net.width=c(6.5,4.3),
					   pdf.net.height=c(3,3))
{
	library(scales)
	library(ggraph)
	library(igraph)

	dat.cor.p <- dat.list$dat.cor.p
	dat.cor.s <- dat.list$dat.cor.s

	######### code from BigSCale
	network=((dat.cor.p>th.cor & dat.cor.s>th.cor) | (dat.cor.p<(-th.cor) & dat.cor.s<(-th.cor)))
	gene.names <- rownames(network)
	diag(network)=FALSE
	network[is.na(network)]=FALSE
	degree=Rfast::rowsums(network)
	if(!is.null(gene.selection)){
		to.include=which(degree>0 & (gene.names %in% gene.selection))
	}else{
		to.include=which(degree>0)
	}

	G=igraph::graph_from_adjacency_matrix(adjmatrix = network[to.include,to.include],mode = 'undirected')
	G=igraph::set_vertex_attr(graph = G,name = "name", value = gene.names[to.include])
	rm(network)
	gc()

	print('Calculating the final score ...')
	#Df=(Ds+Dp)/float::fl(2)
	Df=(dat.cor.p+dat.cor.s)/float::fl(2)

	print(sprintf('Inferred the raw regulatory network: %g nodes and %g edges (ratio E/N)=%f',length(igraph::V(G)),
				  length(igraph::E(G)),length(igraph::E(G))/length(igraph::V(G))))

	if(do.polish==T){
		G=bigSCale:::polish.graph(G)
	}

	comp=igraph::components(G)
	#small.comp=which(comp$csize<0.01*sum(comp$csize))
	#small.comp=which(comp$csize<2)
	if(is.null(th.compsize)){
		small.comp=which(comp$csize<max(comp$csize))
	}else{
		small.comp=which(comp$csize<th.compsize)
	}
	to.remove=which(is.element(comp$membership,small.comp))
	G=igraph::delete_vertices(G, to.remove)
	comp=igraph::components(G)
	cat('\n')
	print(sprintf('Final network after GO filtering: %g nodes and %g edges (ratio E/N)=%f and %g components',length(igraph::V(G)),length(igraph::E(G)),length(igraph::E(G))/length(igraph::V(G)), comp$no))

	print('Computing the centralities')
	Betweenness=igraph::betweenness(graph = G,directed=FALSE,normalized = TRUE)
	Degree=igraph::degree(graph = G)
	PAGErank=igraph::page_rank(graph = G,directed = FALSE)$vector
	Closeness=igraph::closeness(graph = G,normalized = TRUE)

	#### check scale-free topology
	dat.sf <- (table(Degree))
	dat.sf.tb <- data.table(k=as.integer(names(dat.sf)),
							count=as.vector(dat.sf))
	dat.sf.tb[,freq:=count/sum(count)]
	dat.sf.tb[,log10_k:=log10(k)]
	dat.sf.tb[,log10_pk:=log10(freq)]
	summary(lm(log10_pk~log10_k,dat.sf.tb))

	p <- ggscatter(dat.sf.tb,x="log10_k",y="log10_pk",color="steelblue") +
		xlab(expression(log[10](k))) + ylab(expression(log[10](p(k)))) +
		geom_smooth(method='lm') +
		stat_cor(aes(label = paste(..rr.label..,..r.label.., ..p.label.., sep = "~`,`~")))
    ggsave(sprintf("%s.scaleFree.pdf",out.prefix),width=4,height=3)

	G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.Degree',value = Degree)
	G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.PAGErank',value = PAGErank)
	G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.Closeness',value = Closeness)
	G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.Betweenness',value = Betweenness)

	ret.list <- list(graph=G,
					 correlations=Df,
					 correlations.sub=Df[names(Degree),names(Degree)],
					 centrality=cbind(data.table(geneID=names(Degree),
												 comp=comp$membership[names(Degree)]),
									  data.table(cbind(Degree,Betweenness,Closeness,PAGErank))),
					 cutoff.p=th.cor)

	if(!is.null(gset.list)){
		for(gset in names(gset.list)){
			ret.list$centrality[[sprintf("geneSet.%s",gset)]] <- ret.list$centrality$geneID %in% gset.list[[gset]]
		}
	}

	if(verbose==0){
		return(ret.list)
	}

	##### heatmap
	dat.ht <- ret.list$correlations.sub
	## dat.ht <- grn.man.list$correlations.sub

	####
	if(is.null(gene.cluster)){
		res.avai <- c(0.1,0.2,0.3,0.4,0.5,1,1.5,2)
		clust.res.list <- llply(res.avai,function(x){
			clust.res <- do.call(leiden::leiden,c(list(object=G,seed=12345,
													   resolution_parameter=x)))
			clust.res <- sprintf("G%02d",clust.res)
			names(clust.res) <- vertex_attr(G,"name")
			clust.res
						 })
		names(clust.res.list) <- sprintf("leiden.res.%s",res.avai)
		res.louvain <- igraph::cluster_louvain(G)
		clust.res.list[["louvain"]] <- structure(sprintf("G%02d",res.louvain$membership),
												 names=res.louvain$names)
	}else{
		clust.res.list <- list()
		clust.res.list[["geneModule"]] <- gene.cluster
	}
	#print(str(clust.res.list))
	##### 

#	sce.gene <- ssc.build(dat.ht)
#	rowData(sce.gene)$gene.use <- T
#	sce.gene <- ssc.reduceDim(sce.gene,method="pca",method.vgene="gene.use",
#							  pca.npc=10)
#	sce.gene <- ssc.reduceDim(sce.gene,method="umap",method.vgene="gene.use",ncore=1)
#	g.cluster.label.leiden <- run.Leiden(reducedDim(sce.gene,"pca"),SNN.k=5,myseed=123456,resolution_parameter=2)
#	sce.gene$cluster.leiden <- sprintf("G%02d",g.cluster.label.leiden)
#	p <- ssc.plot.tsne(sce.gene,columns=c("cluster.leiden"),reduced.name="umap")
#	ggsave(file=sprintf("%s.gene.cluster.leiden.pdf",out.prefix),width=6,height=4.5)
	#####

	#ret.list$centrality
	##all(ret.list$centrality$geneID==rownames(sce.gene))
	#ret.list$centrality$geneCluster <- clust.res[ret.list$centrality$geneID]
	#ret.list$centrality$geneSig <- ret.list$centrality$geneID %in% gene.desc.tb[sig==T & meta.cluster=="CD8.c12.Tex.CXCL13",][["geneID"]]

	#all(vertex_attr(G,"name")==rownames(sce.gene))
	##gene.desc.sig.tb <- gene.desc.tb[sig==T & geneSet.TF,][order(-comb.ES),]
	gene.desc.sig.tb <- gene.desc.tb[sig==T,][order(-comb.ES),]
	gene.meta.cluster <- sapply(vertex_attr(G,"name"),function(x){
									tt <- match(x,gene.desc.sig.tb$geneID)[1]
									if(is.na(tt)){
										"NA"
									}else{
										gene.desc.sig.tb$meta.cluster[tt]
									}
							  })

	#G=igraph::set.vertex.attribute(graph = G,name = 'geneCluster',value = clust.res)
	for(a.res in names(clust.res.list)){
		G=igraph::set.vertex.attribute(graph = G,name = a.res,value = clust.res.list[[a.res]])
		ret.list$centrality[[a.res]] <- clust.res.list[[a.res]][ret.list$centrality$geneID]
	}
	###G=igraph::set.vertex.attribute(graph = G,name = 'geneCluster',value = sce.gene$cluster.leiden)
	#G=igraph::set.vertex.attribute(graph = G,name = 'geneSigTex',value = rownames(sce.gene) %in% gene.desc.tb[sig==T & meta.cluster=="CD8.c12.Tex.CXCL13",][["geneID"]])
	G=igraph::set.vertex.attribute(graph = G,name = 'geneSet.TF',value = vertex_attr(G,"name") %in% gset.list$TF)
	G=igraph::set.vertex.attribute(graph = G,name = 'geneMetaCluster',value = gene.meta.cluster[vertex_attr(G,"name")])

	ret.list$graph <- G
	bigSCale::toCytoscape(G,file.name=sprintf("%s.cytoscape.json",out.prefix))
	write.table(ret.list$centrality,file=sprintf("%s.geneInfo.txt",out.prefix),sep="\t",quote=F,row.names=F)

	if(is.null(gene.highlight)){
		gene.highlight <- gene.meta.cluster[gene.meta.cluster!="NA"]
		gene.highlight.es <- sapply(names(gene.highlight),function(x){
				   gene.desc.tb[meta.cluster==gene.highlight[x] & geneSymbol==x,][["comb.ES"]]
							  })
		gene.highlight <- names(gene.highlight.es)[gene.highlight.es>th.es]
		#gene.highlight <- gene.desc.tb[sig==T,][,.N,by="geneSymbol"][N==1,][["geneSymbol"]]
		#gene.highlight <- intersect(gene.highlight,rownames(sce.gene))
	}
	print(str(gene.highlight))

	## 'mds','drl', 'gem','kk'
	if(is.null(g.layout)){
		g.layout <- c('dh','fr','graphopt', 'lgl')
	}
	print(g.layout)
	l_ply(g.layout,function(opt.layout){
		set.seed(11234)
		p <- ggraph(G,layout=opt.layout) +
			geom_edge_link(alpha=0.8,colour="lightgray",edge_width=0.1) +
			#geom_node_point(aes(size=bigSCale.PAGErank,color=geneCluster),alpha=0.8) +
			#geom_node_point(aes(size=bigSCale.PAGErank,color=geneSet.TF & geneSigTex),alpha=0.8) +
			#scale_colour_manual(values=c("TRUE"="steelblue","FALSE"="gray")) +
			geom_node_point(aes_string(
								#size="bigSCale.PAGErank",
								#size="bigSCale.Betweenness",
								shape="geneSet.TF",
								color="geneMetaCluster"),alpha=0.8,size=2) +
								##fill=geneMetaCluster),colour="none", alpha=0.8) +
			#scale_fill_manual(values=c(colSet$meta.cluster,c("NA"="gray"))) +
			scale_color_manual(values=c(colSet$meta.cluster,c("NA"="gray"))) +
			#scale_shape_manual(values=c("TRUE"=21,"FALSE"=22)) +
			scale_shape_manual(values=c("TRUE"=16,"FALSE"=15)) +
			guides(color=guide_legend(ncol=2))+
			theme_void()
		p <- p + do.call(geom_node_text,c(list(mapping=aes_string(label="name",color="geneMetaCluster"),
											   data=p$data[p$data$name %in% gene.highlight,],
											   show.legend=F),
										  par.node.text))

		ggsave(sprintf("%s.ggraph.%s.pdf",out.prefix,opt.layout),
			   width=pdf.net.width[1],
			   height=pdf.net.height[1])

		l_ply(names(clust.res.list),function(a.res){
			gene.cluster.avai <- sort(unique(clust.res.list[[a.res]]))
			set.seed(11234)
			p2 <- ggraph(G,layout=opt.layout) +
					geom_edge_link(alpha=0.8,colour="lightgray",edge_width=0.1) +
					geom_node_point(aes_string(
										shape="geneSet.TF",
										color=a.res),alpha=0.8,size=2) +
					scale_color_manual(values=structure(sscVis:::auto.colSet(length(gene.cluster.avai)),
														names=gene.cluster.avai)) +
					scale_shape_manual(values=c("TRUE"=16,"FALSE"=15)) +
					guides(color=guide_legend(ncol=2))+
					theme_void()
			p2 <- p2 + do.call(geom_node_text,c(list(mapping=aes_string(label="name"),
													 data=p2$data[p2$data$name %in% gene.highlight,],
													 show.legend=F),
												par.node.text))
			ggsave(sprintf("%s.ggraph.%s.%s.pdf",out.prefix,opt.layout,a.res),
				   width=pdf.net.width[2],
				   height=pdf.net.height[2])
							  },.parallel=F)

		},.parallel=T)

	#####

	

	#ret.list$p <- p
	######

	#g.cluster <- sscVis:::run.cutreeDynamic(dat.ht,method.hclust="ward.D2",method.distance="spearman",deepSplit=4)
	g.cluster <- sscVis:::run.cutree(dat.ht,method.hclust="ward.D2",method.distance="spearman",k=15)
	names(g.cluster$cluster) <- colnames(dat.ht)

	ann.bar.height <- 0.4
	gene.ann.tb <- data.frame(isTF=colnames(dat.ht) %in% gset.list$TF,
							  isSigTex=colnames(dat.ht) %in% gene.desc.tb[sig==T & meta.cluster=="CD8.c12.Tex.CXCL13",][["geneID"]],
							  stringsAsFactors=F)
	rownames(gene.ann.tb) <- colnames(dat.ht)
	ha.col <- ComplexHeatmap::HeatmapAnnotation(df = gene.ann.tb,
												col = list(isTF=c("TRUE"="steelblue","FALSE"="gray"),
														   isSigTex=c("TRUE"="steelblue","FALSE"="gray")),
												##col = structure(c("blue","gray"), names=c("TRUE","FALSE")),
                                    show_legend = T, annotation_height = unit(rep(ann.bar.height,ncol(gene.ann.tb)), "cm"))
	sscVis::plotMatrix.simple(dat.ht,out.prefix=sprintf("%s.ht.weight",out.prefix),
								  show.dendrogram=T,
								  #do.clust=T,
								  #clust.column="cutreeDynamic",
								  #clust.row="cutreeDynamic",
								  do.clust=T,
								  clust.column=g.cluster$branch,
								  clust.row=g.cluster$branch,
								  column_split=g.cluster$cluster,
								  row_split=g.cluster$cluster,
								  column_gap = unit(0, "mm"),row_gap = unit(0, "mm"),border=T,
								  exp.name="Corr",
								  z.lo=-0.5,z.hi=0.5,palatte=rev(brewer.pal(n = 7,name = "RdBu")),
								  par.legend=list(at = seq(-0.5,0.5,0.25)), top_annotation = ha.col,
								  pdf.width=8,pdf.height=8,fig.type=plot.type)

	#####

	return(ret.list)

}

run.WGCNA <- function()
{
	library("WGCNA")
	## The following setting is important,
	options(stringsAsFactors = FALSE);
	enableWGCNAThreads()

	RhpcBLASctl::omp_set_num_threads(1)

	datExpr <- t(assay(sce.merged[,f.use]))
	#### filtr out some genes
	g.var <- colVars(datExpr,na.rm=T)
	f.gene <- is.na(g.var) | (g.var < quantile(g.var,c(0.25),na.rm=T) )
	datExpr <- datExpr[,!f.gene]

	## Choose a set of soft-thresholding powers
	powers = c(c(1:10), seq(from = 12, to=20, by=2))

	## Call the network topology analysis function
	tic("pickSoftThreshold")
	sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
							###corFnc=cor.new,corOptions = list())
	toc()

	pdf(sprintf("%s.WGCNA.sft.pdf",out.prefix),width=12,heigh=6)
	## Plot the results:
	#sizeGrWindow(9, 5)
	par(mfrow = c(1,2));
	cex1 = 0.9;
	## Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		 xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
		 main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		 labels=powers,cex=cex1,col="red");
	# this line corresponds to using an R^2 cut-off of h
	abline(h=0.90,col="red")
	# Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
		 xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
		 main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	dev.off()

	softPower = 4;
	#softPower = 3;
	#adjacency = adjacency(datExpr, power = softPower);

	tic("blockwiseModules(...)")
	net = blockwiseModules(datExpr, power = softPower,
						   maxBlockSize = 20000,
						   TOMType = "signed",
						   minModuleSize = 20, deepSplit = 4,
						   reassignThreshold = 0,
						   mergeCutHeight = 0.05,
						   numericLabels = TRUE, pamRespectsDendro = FALSE,
						   saveTOMs = TRUE,
						   saveTOMFileBase = sprintf("CD8TOM.%d",softPower), verbose = 3)
	toc()

	# open a graphics window
	pdf(sprintf("%s.WGCNA.gene.dend.softPower%d.pdf",out.prefix,softPower),width=8,height=6)
	#sizeGrWindow(12, 9)
	# Convert labels to colors for plotting
	mergedColors = labels2colors(net$colors)
	# Plot the dendrogram and the module colors underneath
	plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
						"Module colors", dendroLabels = FALSE, hang = 0.03,
						addGuide = TRUE, guideHang = 0.05)
	dev.off()

	moduleLabels = net$colors
	moduleColors = labels2colors(net$colors)
	MEs = net$MEs;
	geneTree = net$dendrograms[[1]];
	save(MEs, moduleLabels, moduleColors, geneTree,
		 file = sprintf("%s.WGCNA.RData",out.prefix))

	table(net$colors)
	table(moduleColors)

#	lname <- load(sprintf("CD8TOM.%d-block.1.RData",softPower))
#	TOM <- as.matrix(TOM)
#	colnames(TOM) <- names(net$colors)[net$goodGenes]
#	rownames(TOM) <- names(net$colors)[net$goodGenes]

#
#	{
#
#			adj.mtx = abs(sscVis:::cor.BLAS(t(datExpr)))^softPower
#			adj.mtx[adj.mtx>1] <- 1
#			TOM = TOMsimilarity(adj.mtx);
#			dissTOM = 1-TOM
#
#			# Call the hierarchical clustering function
#			#method.hclust <- "ward.D2"
#			method.hclust <- "average"
#			geneTree = fastcluster::hclust(as.dist(dissTOM), method = method.hclust);
#
#			pdf(sprintf("%s.WGCNA.gene.dend.01.pdf",out.prefix),width=8,height=6)
#			# Plot the resulting clustering tree (dendrogram)
#			#sizeGrWindow(12,9)
#			plot(geneTree, xlab="", sub="", main= "Gene clustering on TOM-based dissimilarity",
#				 labels = FALSE, hang = 0.04);
#			dev.off()
#
#			# We like large modules, so we set the minimum module size relatively high:
#			minModuleSize = 20;
#			# Module identification using dynamic tree cut:
#			dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
#										#method = "tree",deepSplit=T,
#										deepSplit = 2,
#										pamRespectsDendro = T,
#										minClusterSize = minModuleSize);
#			table(dynamicMods)
#
#			# Convert numeric lables into colors
#			dynamicColors = labels2colors(dynamicMods)
#			table(dynamicColors)
#
#			# Plot the dendrogram and colors underneath
#			pdf(sprintf("%s.WGCNA.gene.dend.02.pdf",out.prefix),width=8,height=6)
#			#sizeGrWindow(8,6)
#			plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
#								dendroLabels = FALSE, hang = 0.03,
#								addGuide = TRUE, guideHang = 0.05,
#								main = "Gene dendrogram and module colors")
#			dev.off()
#
#			# Calculate eigengenes
#			MEList = moduleEigengenes(datExpr, colors = dynamicColors,softPower=softPower)
#			MEs = MEList$eigengenes
#			# Calculate dissimilarity of module eigengenes
#			MEDiss = 1-cor(MEs);
#			# Cluster module eigengenes
#			##METree = hclust(as.dist(MEDiss), method = "average");
#			METree = hclust(as.dist(MEDiss), method = method.hclust);
#			# Plot the result
#
#			pdf(sprintf("%s.WGCNA.gene.dend.03.pdf",out.prefix),width=8,height=6)
#			#sizeGrWindow(7, 6)
#			plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
#			MEDissThres = 0.05
#			# Plot the cut line into the dendrogram
#			abline(h=MEDissThres, col = "red")
#			# Call an automatic merging function
#			merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
#			# The merged module colors
#			mergedColors = merge$colors;
#			# Eigengenes of the new merged modules:
#			mergedMEs = merge$newMEs;
#			dev.off()
#
#			pdf(sprintf("%s.WGCNA.gene.dend.04.pdf",out.prefix),width=8,height=6)
#			#sizeGrWindow(12, 9)
#			plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
#								c("Dynamic Tree Cut", "Merged dynamic"),
#								dendroLabels = FALSE, hang = 0.03,
#								addGuide = TRUE, guideHang = 0.05)
#			dev.off()
#
#			# Rename to moduleColors
#			moduleColors = mergedColors
#			# Construct numerical labels corresponding to the colors
#			colorOrder = c("grey", standardColors(50));
#			moduleLabels = match(moduleColors, colorOrder)-1;
#			MEs = mergedMEs;
#
#	}
#
	# Save module colors and labels for use in subsequent parts
	#save(MEs, moduleLabels, moduleColors, geneTree,
	#	 file = sprintf("%s.WGCNA.moduleColors.RData",out.prefix))

}


run.nichenet <- function(gene.oi,gene.bg,out.prefix,
						 comb.width=12,comb.height=7,comb.rel.height=c(10,2),comb.rel.width=c(0.1,0.9),
						 lr.width=5,lr.height=7,n.top=20,pearson.max=0.2,
						 ligands_all=NULL,targets_all=NULL,es.df=NULL,
						 do.eval=T,eval.k.fold=5,eval.n.round=10)
{
	library("nichenetr")
	library("tidyverse")

	weighted_networks <- readRDS("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Nichnet/weighted_networks.rds")
	ligand_target_matrix <- readRDS("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Nichnet/ligand_target_matrix.rds")
	lr_network <- readRDS("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Nichnet/lr_network.rds")
	str(ligand_target_matrix)
	lr_network_expressed <- lr_network %>% filter(to %in% gene.bg)
	potential_ligands  <-  lr_network_expressed %>% pull(from) %>% unique()
	ligand_activities <- predict_ligand_activities(geneset = gene.oi,
						       background_expressed_genes = gene.bg,
						       ligand_target_matrix = ligand_target_matrix,
						       potential_ligands = potential_ligands)
	ligand_activities %>% arrange(-pearson)
	best_upstream_ligands <- ligand_activities %>% top_n(n.top, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
	###
	p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) +
		geom_histogram(color="black", fill="darkorange")  +
		# geom_density(alpha=.1, fill="orange") +
		geom_vline(aes(xintercept=min(ligand_activities %>% top_n(n.top, pearson) %>% pull(pearson))),
				   color="red", linetype="dashed", size=1) +
		labs(x="ligand activity (PCC)", y = "# ligands") +
		theme_classic()
	ggsave(sprintf("%s.dist.pearson.pdf",out.prefix),width=5,height=4)
	###
	active_ligand_target_links_df <- best_upstream_ligands %>%
	    lapply(get_weighted_ligand_target_links,geneset = gene.oi,
		   ligand_target_matrix=ligand_target_matrix, n=250) %>%
	    bind_rows()
	active_ligand_target_links_df <- active_ligand_target_links_df %>% filter(!is.na(weight))
	active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df=active_ligand_target_links_df,
									  ligand_target_matrix = ligand_target_matrix,
									  cutoff = 0.25)
	order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
	order_targets = active_ligand_target_links_df$target %>% unique()
	order_targets = intersect(order_targets,rownames(active_ligand_target_links))
	active_ligand_target_links.debug <<- active_ligand_target_links
	order_targets.debug <<- order_targets
	order_ligands.debug <<- order_ligands
	vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
	colnames(vis_ligand_target) <- make.names(colnames(vis_ligand_target))
	gene.unexplained <- setdiff(gene.oi,colnames(active_ligand_target_links))
	vis_ligand_target.debug <<- vis_ligand_target
	vis_ligand_target[vis_ligand_target>0.01] <- 0.01

	p_ligand_target_network = vis_ligand_target %>%
				    make_heatmap_ggplot("Prioritized ligands","genes in receiver cells",
							color = "purple",legend_position = "top",
							x_axis_position = "top",
							legend_title = "Regulatory potential") +
					scale_fill_gradient2(low = "whitesmoke",  high = "purple",
							     limits=c(0,0.01), breaks = c(0,0.005,0.01)) +
					theme(axis.text.x = element_text(face = "italic"),legend.key.width=unit(2,"cm"))
	#ggsave(file=sprintf("%s.ligand_target_network.tmp.pdf",out.prefix),width=8,height=5)

	####
	ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>%
				    magrittr::set_rownames(ligand_activities$test_ligand)

	vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
	vis_ligand_pearson.debug <<- vis_ligand_pearson
	vis_ligand_pearson[vis_ligand_pearson > pearson.max ] <- pearson.max
	p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity",
							color = "darkorange",legend_position = "top",
							x_axis_position = "top",
							legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") +
				scale_fill_gradient2(low = "whitesmoke",  high = "darkorange", limits=c(0,pearson.max))
	p_ligand_pearson <- p_ligand_pearson + theme(legend.key.width=unit(2,"cm"))
	#if(!is.null(keywidth)){
	#	p_ligand_pearson <- p_ligand_pearson + theme(legend.key.width=keywidth)
	#}
	#ggsave(file=sprintf("%s.ligand_pearson.pdf",out.prefix),width=5,height=5)

	{
	####
	figures_without_legend = plot_grid(p_ligand_pearson +
					   theme(legend.position = "none", axis.ticks = element_blank()) +
					   theme(axis.title.x = element_text()),
				       p_ligand_target_network + theme(legend.position = "none",
								       axis.text.y=element_blank(),
								       axis.title.y = element_blank(),
								       axis.ticks = element_blank()),
					   #### donnot let align contain "v" here
					   align = "h", nrow = 1,
					   rel_widths = comb.rel.width,
					   rel_heights = c(nrow(vis_ligand_pearson),
							   nrow(vis_target_tumor_expression_scaled) + 3))

	legends = plot_grid(as_ggplot(get_legend(p_ligand_pearson)),
						as_ggplot(get_legend(p_ligand_target_network)),
						nrow = 2,
						align = "h")
	pp <- plot_grid(figures_without_legend, legends, rel_heights = comb.rel.height, nrow = 2, align = "hv")
	ggsave(file=sprintf("%s.comb.pdf",out.prefix),width=comb.width,height=comb.height)
	}

	#### ligand-receptor network
	{
	    # get the ligand-receptor network of the top-ranked ligands
	    lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% gene.bg) %>% distinct(from,to)
	    best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

	    # get the weights of the ligand-receptor interactions as used in the NicheNet model
	    lr_network_top_df = weighted_networks$lr_sig %>%
							    filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

	    # convert to a matrix
	    lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
	    lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

	    # perform hierarchical clustering to order the ligands and receptors
	    dist_receptors = dist(lr_network_top_matrix, method = "binary")
	    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
	    order_receptors = hclust_receptors$labels[hclust_receptors$order]

	    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
	    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
	    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

	    vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
	    p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>%
									    make_heatmap_ggplot("Prioritized ligands","Receptors expressed",
														    color = "mediumvioletred",
														    x_axis_position = "top",
														    legend_title = "Prior interaction potential")
	    ggsave(file=sprintf("%s.lr.matrix.pdf",out.prefix),width=lr.width,height=lr.height)

	}

	if(!is.null(ligands_all)){
		if(is.null(targets_all)){
			#targets_all <- active_ligand_target_links_df %>% filter(ligand %in% ligands_all) %>% pull(target)
			targets_all <- as.data.frame(t(active_ligand_target_links)) %>% rownames_to_column("ligand") %>%
				gather(target,weight,-ligand) %>% 
				filter(ligand %in% ligands_all & weight>0) %>% pull(target) %>% unique()
		}
		ligand_tf_matrix <- readRDS("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Nichnet/ligand_tf_matrix.rds")
		sig_network <- readRDS("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Nichnet/signaling_network.rds")
		gr_network <- readRDS("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Nichnet/gr_network.rds")
		active_signaling_network <- get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix,
															 ligands_all = ligands_all,
															 targets_all = targets_all,
															 weighted_networks = weighted_networks)

		# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
		active_signaling_network_min_max = active_signaling_network
		active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
		active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

		graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max,
														  ligands_all = ligands_all,
														  targets_all = targets_all,
														  sig_color = "indianred", gr_color = "steelblue")

		###gg <- DiagrammeR::render_graph(graph_min_max, layout = "tree")
		###DiagrammeR::export_graph(graph_min_max,file_name=sprintf("%s.test.png",out.prefix),file_type="png",width=7,height=7)
		data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,
														   lr_network = lr_network,
														   sig_network = sig_network,
														   gr_network = gr_network)
		head(data_source_network) 
		####
		{
			write_output = TRUE # change to TRUE for writing output

			# weighted networks ('import network' in Cytoscape)
			if(write_output){
				bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"),
						  active_signaling_network$gr %>% mutate(layer = "regulatory")) %>%
					write_tsv(sprintf("%s.weighted_signaling_network.txt",out.prefix))
			}

			# networks with information of supporting data sources ('import network' in Cytoscape)
			if(write_output){
				data_source_network %>% write_tsv(sprintf("%s.data_source_network.txt",out.prefix))
			}

			# Node annotation table ('import table' in Cytoscape)
			specific_annotation_tbl = bind_rows(
												tibble(gene = ligands_all, annotation = "ligand"),
												tibble(gene = targets_all, annotation = "target"),
												tibble(gene = c(data_source_network$from, data_source_network$to) %>%
														unique() %>% setdiff(c(targets_all,ligands_all)) %>%
														intersect(lr_network$to %>% unique()), annotation = "receptor"),
												tibble(gene = c(data_source_network$from, data_source_network$to) %>%
													   unique() %>% setdiff(c(targets_all,ligands_all)) %>%
													   intersect(gr_network$from %>% unique()) %>%
													   setdiff(c(data_source_network$from, data_source_network$to) %>%
															   unique() %>% intersect(lr_network$to %>% unique())),
													   annotation = "transcriptional regulator"))
			non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>%
												 unique() %>% setdiff(specific_annotation_tbl$gene),
											 annotation = "signaling mediator")

			node.annotation.tbl <- bind_rows(specific_annotation_tbl,non_specific_annotation_tbl)
			if(!is.null(es.df)){
				node.annotation.tbl <- left_join(node.annotation.tbl,es.df,by="gene")
			}
			if(write_output){
				 write_tsv(node.annotation.tbl,sprintf("%s.annotation_table.txt",out.prefix))
			}

		}

	}
	
	### assess how well top-ranked ligands can predict a gene set of interest
	eval.out.list <- NULL
	if(do.eval)
	{
		# change rounds and folds here, to two rounds to reduce time: normally: do multiple rounds
		#eval.k.fold = 5
		#eval.n.round = 2
		eval_gene_predictions_top20_list = seq(eval.n.round) %>% lapply(assess_rf_class_probabilities,
															 folds = eval.k.fold, 
															 geneset = gene.oi,
															 background_expressed_genes = gene.bg,
															 ligands_oi = best_upstream_ligands,
															 ligand_target_matrix = ligand_target_matrix)
		target_prediction_performances_cv = eval_gene_predictions_top20_list %>%
													lapply(classification_evaluation_continuous_pred_wrapper) %>%
													bind_rows() %>%
													mutate(round=seq(1:nrow(.)))

		target_prediction_performances_discrete_cv = eval_gene_predictions_top20_list %>%
													lapply(calculate_fraction_top_predicted, quantile_cutoff = 0.95) %>%
													bind_rows() %>% ungroup() %>%
													mutate(round=rep(1:length(eval_gene_predictions_top20_list), each = 2))
		target_prediction_performances_discrete_fisher = eval_gene_predictions_top20_list %>%
															lapply(calculate_fraction_top_predicted_fisher,
																   quantile_cutoff = 0.95)

		top_predicted_genes = seq(length(eval_gene_predictions_top20_list)) %>%
								lapply(get_top_predicted_genes,eval_gene_predictions_top20_list) %>%
								reduce(full_join, by = c("gene","true_target"))
		top_predicted_genes %>% filter(true_target)
		f.gene.consensus <- apply(top_predicted_genes[,grep("^predicted_top_target_round",
															colnames(top_predicted_genes),value=T)],1,
								  function(x){ mean(x,na.rm=T) >= 0.5 })
		gene.target.missing <- top_predicted_genes[f.gene.consensus,] %>% pull(gene) %>% setdiff(gene.oi,.)
		eval.out.list <- list(target_prediction_performances_cv=target_prediction_performances_cv,
							  target_prediction_performances_discrete_cv=target_prediction_performances_discrete_cv,
							  top_predicted_genes=top_predicted_genes,
							  gene.target.missing=gene.target.missing,
							  avg.auroc=target_prediction_performances_cv$auroc %>% mean(),
							  avg.aupr=target_prediction_performances_cv$aupr %>% mean(),
							  avg.pearson=target_prediction_performances_cv$pearson %>% mean(),
							  avg.true_target.fraction_positive_predicted= target_prediction_performances_discrete_cv %>%
								  filter(true_target) %>% .$fraction_positive_predicted %>% mean(),
							  avg.false_target.fraction_positive_predicted= target_prediction_performances_discrete_cv %>%
								  filter(!true_target) %>% .$fraction_positive_predicted %>% mean(),
							  avg.target_prediction_performances_discrete_fisher=
								  target_prediction_performances_discrete_fisher %>% unlist() %>% mean()
							  )

	}
	
	return(list("p_ligand_target_network"=p_ligand_target_network,"p_ligand_pearson"=p_ligand_pearson,
				"p_ligand_receptor_network"=p_ligand_receptor_network,
				"gene.oi"=gene.oi,
				"gene.bg"=gene.bg,
				"eval.list"=eval.out.list,
				"gene.unexplained"=gene.unexplained))

}

run.gene.analysis.volcano <- function(gene.desc.top,out.prefix,plot.type="png",vector.friendly=F)
{
	dir.create(sprintf("%s/all",dirname(out.prefix)),F,T)
	dir.create(sprintf("%s/TF",dirname(out.prefix)),F,T)
	dir.create(sprintf("%s/kinase",dirname(out.prefix)),F,T)
	dir.create(sprintf("%s/regulation.epigenetic",dirname(out.prefix)),F,T)
	dir.create(sprintf("%s/drug.target",dirname(out.prefix)),F,T)
	dir.create(sprintf("%s/secreted",dirname(out.prefix)),F,T)

	l_ply(unique(gene.desc.top$meta.cluster),function(mcls){
		###mcls <- "CD8.c08c09.Tex.CXCL13"
		o.tb <- gene.desc.top[order(-comb.ES),][meta.cluster==mcls,]
		write.table(o.tb,file=sprintf("%s.DEG.%s.all.txt",out.prefix,mcls),row.names=F,sep="\t",quote=F)

		############ all genes
		out.prefix.plot <- sprintf("%s/all/%s",dirname(out.prefix),basename(out.prefix))
		set.seed(123456)
		o.plot.tb <- o.tb
		o.plot.tb[,geneSet.drug.target:=geneSet.drug.target.FDA|geneSet.drug.target.potential]
		o.plot.tb[comb.ES>2,comb.ES:=2]
		#o.plot.tb[comb.padj<1e-15,comb.padj:=1e-15]
		o.plot.tb[,sig:=factor(sig,levels=c("TRUE","FALSE"))]
		p <- ggplot(o.plot.tb,aes(y=-log10(comb.padj),x=comb.ES)) + 
				geom_point(aes(color=sig,shape=sig.cate)) +
				scale_shape_manual(values=c("notSig"=16,"Tier1"=17,"Tier2"=15,"Tier3"=3)) +
				scale_color_manual(values=c("TRUE"="#F8766D","FALSE"="#00BFC4")) +
				xlab("Combined Effect Size") + 
				ylab("-log10(Adj.CombinedPValue)") +
		#		coord_cartesian(ylim=c(0,20)) +
				theme_bw()
		if(!vector.friendly){
			p <- p + geom_text_repel(data=o.plot.tb[sig==T,head(.SD,n=30)],aes(label=geneSymbol),
								size=if(vector.friendly) 3.5 else 2)
		}else{
		  	  ####p <- Seurat::AugmentPlot(p,width=width,height=height)
			  tmpfilename <- sprintf("%s.DEG.%s.all.tmp.png",out.prefix.plot,mcls)
			  #print(tmpfilename)
			  ggsave(filename=tmpfilename,
					 plot = p + theme_void() + theme(legend.position = "none",
									  axis.line.x = element_blank(), axis.line.y = element_blank(),
									  axis.title.x = element_blank(), axis.title.y = element_blank(),
									  axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
									  plot.title = element_blank()),
					 width=5.0,height=5.0)
			  pbuild.params <- ggplot_build(plot = p)$layout$panel_params[[1]]
			  range.values <- c( pbuild.params$x.range, pbuild.params$y.range)
			  img <- png::readPNG(source = tmpfilename)
			  blank <- ggplot(data = p$data,mapping = aes(y=-log10(comb.padj),x=comb.ES)) +
						  geom_blank() +
						  xlab("Combined Effect Size") +
						  ylab("-log10(Adj.CombinedPValue)")
			  blank <- blank + p$theme + coord_cartesian(xlim = range.values[1:2], ylim = range.values[3:4], expand = F)
			  blank <- blank + annotation_raster(raster = img,
												 xmin = range.values[1], xmax = range.values[2],
												 ymin = range.values[3], ymax = range.values[4])
			  blank <- blank + geom_hline(yintercept=seq(0,100,25),linetype=2,alpha=0.2) +
				  geom_vline(xintercept=seq(-2,2,0.5),linetype=2,alpha=0.2)
			  blank <- blank + geom_text_repel(data=o.plot.tb[sig==T,head(.SD,n=30)],aes(label=geneSymbol),
								size=if(vector.friendly) 3.5 else 2,segment.size=0.25)
			  legend.blank <- cowplot::get_legend(p)
			  p <- cowplot::plot_grid(blank, legend.blank, rel_widths = c(3.0,1))
			  #p <- blank
			  file.remove(tmpfilename)
		}
		ggsave(sprintf("%s.DEG.%s.all.%s",out.prefix.plot,mcls,plot.type),width=5.5,height=4,useDingbats=F)

		############ 
		for(a.cate in c("TF","kinase","regulation.epigenetic","secreted","drug.target"))
		{
			out.prefix.plot <- sprintf("%s/%s/%s",dirname(out.prefix),a.cate,basename(out.prefix))
			set.seed(123456)
			f.show.all <- o.plot.tb[[sprintf("geneSet.%s",a.cate)]]==T
			f.show.hl <- o.plot.tb[[sprintf("geneSet.%s",a.cate)]]==T & o.plot.tb[["sig"]]==T

			p <- ggplot(o.plot.tb[f.show.all,],aes(y=-log10(comb.padj),x=comb.ES)) + 
					geom_point(aes(color=sig,shape=sig.cate)) +
					scale_shape_manual(values=c("notSig"=16,"Tier1"=17,"Tier2"=15,"Tier3"=3)) +
					scale_color_manual(values=c("TRUE"="#F8766D","FALSE"="#00BFC4")) +
					xlab("Combined Effect Size") + 
					ylab("-log10(Adj.CombinedPValue)") +
			#		coord_cartesian(ylim=c(0,20)) +
					theme_bw()
			if(!vector.friendly){
				p <- p + geom_text_repel(data=o.plot.tb[f.show.hl,][,head(.SD,n=30)],
									aes(label=geneSymbol),size=2)
			}else{
		  	  ####p <- Seurat::AugmentPlot(p,width=width,height=height)
			  tmpfilename <- sprintf("%s.DEG.%s.%s.tmp.png",out.prefix.plot,mcls,a.cate)
			  #print(tmpfilename)
			  ggsave(filename=tmpfilename,
					 plot = p + theme_void() + theme(legend.position = "none",
									  axis.line.x = element_blank(), axis.line.y = element_blank(),
									  axis.title.x = element_blank(), axis.title.y = element_blank(),
									  axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
									  plot.title = element_blank()),
					 width=5.0,height=5.0)
			  pbuild.params <- ggplot_build(plot = p)$layout$panel_params[[1]]
			  range.values <- c( pbuild.params$x.range, pbuild.params$y.range)
			  img <- png::readPNG(source = tmpfilename)
			  blank <- ggplot(data = p$data,mapping = aes(y=-log10(comb.padj),x=comb.ES)) +
						  geom_blank() +
						  xlab("Combined Effect Size") +
						  ylab("-log10(Adj.CombinedPValue)")
			  blank <- blank + p$theme + coord_cartesian(xlim = range.values[1:2], ylim = range.values[3:4], expand = F)
			  blank <- blank + annotation_raster(raster = img,
												 xmin = range.values[1], xmax = range.values[2],
												 ymin = range.values[3], ymax = range.values[4])
			  blank <- blank + geom_hline(yintercept=seq(0,100,25),linetype=2,alpha=0.2) +
				  geom_vline(xintercept=seq(-2,2,0.5),linetype=2,alpha=0.2)
			  blank <- blank + geom_text_repel(data=o.plot.tb[f.show.hl,][,head(.SD,n=30)],
									aes(label=geneSymbol),size=3,segment.size=0.25)
			  legend.blank <- cowplot::get_legend(p)
			  p <- cowplot::plot_grid(blank, legend.blank, rel_widths = c(3.0,1))
			  #p <- blank
			  file.remove(tmpfilename)
			}
			ggsave(sprintf("%s.DEG.%s.%s.%s",out.prefix.plot,mcls,a.cate,plot.type),width=5.5,height=4,useDingbats=F)

		}

	},.parallel=T)

}

run.gene.analysis.signatureAmongDatasetPlot <- function(obj,gene.desc.top,out.prefix)
{

	dir.create(sprintf("%s/signatureAmongDatasetPlot.Tier1",dirname(out.prefix)),F,T)
	dir.create(sprintf("%s/signatureAmongDatasetPlot.TF",dirname(out.prefix)),F,T)
	l_ply(unique(gene.desc.top$meta.cluster),function(mcls){
	    out.prefix.plot <- sprintf("%s/signatureAmongDatasetPlot.Tier1/%s",
								   dirname(out.prefix),basename(out.prefix))
		signatureAmongDatasetPlot(obj,
								  gene.desc.top=gene.desc.top[sig.cate=="Tier1" & meta.cluster==mcls,],
								  mcls=mcls,th.universal=0.5,
								  out.prefix=sprintf("%s.DEG.Tier1",out.prefix.plot),assay.name="sig",k.column=2)
	    out.prefix.plot <- sprintf("%s/signatureAmongDatasetPlot.TF/%s",
								   dirname(out.prefix),basename(out.prefix))
		signatureAmongDatasetPlot(obj,
								  gene.desc.top=gene.desc.top[geneSet.TF==TRUE & sig==T & meta.cluster==mcls,],
								  mcls=mcls,th.universal=0.5,
								  out.prefix=sprintf("%s.DEG.TF",out.prefix.plot),assay.name="sig",k.column=2)
	},.parallel=T)

}

sigGeneBublePlot <- function(dat.plot.tb,out.prefix,ntop=30,nplot=1000,my.seed=123456,
								 pdf.width=7,pdf.height=4,text.hl=NULL)
{

	dat.plot.tb[,Z:=median.meanScale]
	dat.plot.tb[median.meanScale >  1,Z:=1]
	dat.plot.tb[median.meanScale < -1,Z:=-1]
	dat.plot.text.tb <- head(dat.plot.tb[order(median.rank),][freq.sig>=0.1,],n=ntop)
	dat.plot.text.tb[,gene:=sprintf("italic('%s')",geneSymbol)]
	### highlight text
	if(!is.null(text.hl)){
		dat.plot.hl.tb <- dat.plot.tb[geneSymbol %in% text.hl,]
		dat.plot.hl.tb[,gene:=sprintf("italic('%s')",geneSymbol)]
		dat.plot.tb <- dat.plot.tb[!(geneSymbol %in% dat.plot.hl.tb$geneSymbol),]
		dat.plot.text.tb <- dat.plot.text.tb[!(geneSymbol %in% dat.plot.hl.tb$geneSymbol),]
	}else{
		dat.plot.hl.tb <- NULL
	}

	p <- ggplot(head(dat.plot.tb,n=nplot),aes(median.rank,freq.sig)) +
			geom_point(aes(color=Z,size=freq.sig*20),alpha=0.8)+
			scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))+
			geom_text_repel(data=dat.plot.text.tb,aes(label=gene),parse=T,
							segment.alpha=0.3,seed=my.seed,
							size=3,force=28)
	if(!is.null(dat.plot.hl.tb)){
		p <- p + geom_point(aes(color=Z,size=freq.sig*20),alpha=0.8,data=dat.plot.hl.tb)+
			#scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))+
			geom_text_repel(data=dat.plot.hl.tb,aes(label=gene),parse=T,colour="darkgreen",segment.alpha=0.3,seed=my.seed, size=3.5,force=28)
	}
	p <- p + guides(size = FALSE)+
			theme_bw()
	ggsave(sprintf("%s.pdf",out.prefix),width=pdf.width,height=pdf.height)
}

sigGenePPPlot <- function(dat.plot.x.tb,dat.plot.y.tb,idx.col="comb.ES",out.prefix=NULL,
						  ntop.x=10,ntop.y=10,ntop.both=10,
						  ntop.x.TF=10,ntop.y.TF=10,ntop.both.TF=10,
						  nplot=1000,my.seed=123456,
						  clamp.x=NULL,clamp.y=NULL,text.force=1,
						  min.freq.sig=0.1,sig.tb=NULL,lab.x="X",lab.y="Y",
						  plot.type="png",vector.friendly=F,
						  pdf.width=5,pdf.height=5,text.hl.tb=NULL,do.plot=T)
{

	dat.plot.tb <- merge(dat.plot.x.tb,dat.plot.y.tb,by=c("geneID","geneSymbol"))
	dat.plot.tb$Z.x <- dat.plot.tb[[sprintf("%s.x",idx.col)]]
	dat.plot.tb$Z.y <- dat.plot.tb[[sprintf("%s.y",idx.col)]]
	dat.plot.tb[,Z:=ifelse(Z.x>Z.y,Z.x,Z.y)]
	dat.plot.tb[,freq.sig:=ifelse(Z.x>Z.y,freq.sig.x,freq.sig.y)]
	#dat.plot.tb[median.meanScale.x >  1,Z.x:=1]
	#dat.plot.tb[median.meanScale.y < -1,Z.y:=-1]
	#dat.plot.tb[,Z:=ifelse(median.rank.x<median.rank.y,Z.x,Z.y)]
	#dat.plot.tb[,freq.sig:=ifelse(median.rank.x<median.rank.y,freq.sig.x,freq.sig.y)]
	dat.plot.tb[,gene:=sprintf("italic('%s')",geneSymbol)]
	dat.plot.tb[geneSet.TF.y==T,gene:=sprintf("bold(italic('%s'))",geneSymbol)]

	dat.plot.tb[,sig.cate2:="notSig"]
	dat.plot.tb[sig.x==T & sig.y==T,sig.cate2:="Both"]
	dat.plot.tb[sig.x==T & sig.y==F,sig.cate2:=lab.x]
	dat.plot.tb[sig.x==F & sig.y==T,sig.cate2:=lab.y]
	if(!is.null(clamp.x)){
		dat.plot.tb[Z.x < clamp.x[1], Z.x:=clamp.x[1]]
		dat.plot.tb[Z.x > clamp.x[2], Z.x:=clamp.x[2]]
	}
	if(!is.null(clamp.y)){
		dat.plot.tb[Z.y < clamp.y[1], Z.y:=clamp.y[1]]
		dat.plot.tb[Z.y > clamp.y[2], Z.y:=clamp.y[2]]
	}
	dat.plot.tb[,sig.cate2:=factor(sig.cate2,levels=c("notSig","Both",lab.x,lab.y))]

	dat.plot.text.x.TF.tb <- head(dat.plot.tb[sig.x==T & sig.y==F & geneSet.TF.x==T,][order(-Z.x),],n=ntop.x.TF)
	dat.plot.text.y.TF.tb <- head(dat.plot.tb[sig.y==T & sig.x==F & geneSet.TF.y==T,][order(-Z.y),],n=ntop.y.TF)
	dat.plot.text.both.TF.tb <- head(dat.plot.tb[sig.y==T & sig.x==T & geneSet.TF.x==T & geneSet.TF.y==T,][order(-Z.y),],n=ntop.both.TF)
	dat.plot.text.x.tb <- head(dat.plot.tb[sig.x==T & sig.y==F,][!(geneSymbol %in% dat.plot.text.x.TF.tb$geneSymbol),][order(-Z.x),],n=ntop.x)
	dat.plot.text.y.tb <- head(dat.plot.tb[sig.y==T & sig.x==F,][!(geneSymbol %in% dat.plot.text.y.TF.tb$geneSymbol),][order(-Z.y),],n=ntop.y)
	dat.plot.text.both.tb <- head(dat.plot.tb[sig.y==T & sig.x==T,][!(geneSymbol %in% dat.plot.text.both.TF.tb$geneSymbol),][order(-Z.y),],n=ntop.both)
	text.hl.tb <- rbind(dat.plot.text.x.TF.tb,dat.plot.text.y.TF.tb,dat.plot.text.both.TF.tb,
						dat.plot.text.x.tb,dat.plot.text.y.tb,dat.plot.text.both.tb)
#	dat.plot.text.x.tb <- head(dat.plot.tb[freq.sig.x>=min.freq.sig & freq.sig.y < min.freq.sig,
#							   ][order(median.rank.x),],n=ntop.x)
#	dat.plot.text.y.tb <- head(dat.plot.tb[freq.sig.y>=min.freq.sig & freq.sig.x < min.freq.sig,
#							   ][order(median.rank.y),],n=ntop.y)

	set.seed(my.seed)
	dat.plot.tb <- dat.plot.tb[sample(nrow(dat.plot.tb)),]

#	if(!is.null(sig.tb)){
#		p <- ggplot(dat.plot.tb[!(geneID %in% sig.tb$geneID),],aes(median.rank.x,median.rank.y)) +
#				geom_point(aes(size=freq.sig*20),color="lightgray",alpha=0.8)+
#				geom_point(aes(color=Z,size=freq.sig*20),alpha=0.8,
#						   data=dat.plot.tb[(geneID %in% sig.tb$geneID),])+
#				scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))
#	}else{
#		p <- ggplot(dat.plot.tb,aes(median.rank.x,median.rank.y)) +
#				geom_point(aes(color=Z,size=freq.sig*20),alpha=0.8)+
#				scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))
#	}
	v.pretty <- pretty(c(dat.plot.tb$Z.x,dat.plot.tb$Z.y))
	print(v.pretty)
	p <- ggplot(dat.plot.tb,aes(Z.x,Z.y)) +
			#geom_point(aes(color=Z,size=freq.sig*20),alpha=0.8)+
			#scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))
			geom_point(aes(color=sig.cate2),alpha=0.5)+
			##scale_colour_manual(values=c("notSig"="lightgray","Both"="#4DAF4A","X"="#377EB8","Y"="#E41A1C"))
			scale_colour_manual(labels=c("notSig","Both",lab.x,lab.y),
								values=c("lightgray","#4DAF4A","#377EB8","#E41A1C")) +
			coord_cartesian(xlim=c(v.pretty[1],v.pretty[length(v.pretty)]), ylim=c(v.pretty[1],v.pretty[length(v.pretty)])) +
			xlab(lab.x)+
			ylab(lab.y)
	if(vector.friendly){
		tmpfilename <- sprintf("%s.tmp.png",if(is.null(out.prefix)) "" else out.prefix)
		ggsave(filename=tmpfilename,
			 plot = p + theme_void() + theme(legend.position = "none",
							  axis.line.x = element_blank(), axis.line.y = element_blank(),
							  axis.title.x = element_blank(), axis.title.y = element_blank(),
							  axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
							  plot.title = element_blank()),
			 width=5.0,height=5.0)
		pbuild.params <- ggplot_build(plot = p)$layout$panel_params[[1]]
		range.values <- c( pbuild.params$x.range, pbuild.params$y.range)
		img <- png::readPNG(source = tmpfilename)
		blank <- ggplot(data = p$data,mapping = aes(Z.x,Z.y)) +
				  geom_blank() +
				  xlab(sprintf("Combined Effect Size (%s)",lab.x)) +
				  ylab(sprintf("Combined Effect Size (%s)",lab.y))
		blank <- blank + p$theme + coord_cartesian(xlim = range.values[1:2], ylim = range.values[3:4], expand = F)
		blank <- blank + annotation_raster(raster = img,
										 xmin = range.values[1], xmax = range.values[2],
										 ymin = range.values[3], ymax = range.values[4])
		#blank <- blank + geom_hline(yintercept=seq(0,100,25),linetype=2,alpha=0.2) +
		#	geom_vline(xintercept=seq(-2,2,0.5),linetype=2,alpha=0.2)
		text.hl.tb.debug <<- text.hl.tb
		dat.plot.tb.debug <<- dat.plot.tb
		col.tex <- structure(c("lightgray","#4DAF4A","#377EB8","#E41A1C"),
							 names=c("notSig","Both",lab.x,lab.y))
		blank <- blank + geom_text_repel(data=text.hl.tb,
							aes(label=gene,color=sig.cate2,size=geneSet.TF.y),parse=T,
							#aes(label=gene,size=geneSet.TF.y),parse=T,color="green",
							segment.alpha=0.8,seed=my.seed,segment.size=0.25,
							show.legend=F,
							force=text.force) +
			scale_size_manual(values=c("TRUE"=3,"FALSE"=2)) +
			scale_colour_manual(values=col.tex) +
			guides(size = FALSE)+
			geom_abline(slope=1,intercept=0,linetype=2,size=0.8,alpha=0.2) +
			theme_bw()
		legend.blank <- cowplot::get_legend(p)
		p <- cowplot::plot_grid(blank, legend.blank, rel_widths = c(2.2,1))
		#p <- blank
		file.remove(tmpfilename)
	}else{
		p <- p + 
			geom_text_repel(data=text.hl.tb,
							aes(label=gene,color=sig.cate2,size=geneSet.TF.y),parse=T,
							segment.alpha=0.8,seed=my.seed,
							show.legend=F,
							force=text.force) +
			scale_size_manual(values=c("TRUE"=3,"FALSE"=2)) +
			guides(size = FALSE)+
			geom_abline(slope=1,intercept=0,linetype=2,size=0.8,alpha=0.2) +
			theme_bw()
	}
	if(do.plot){
		ggsave(sprintf("%s.%s",out.prefix,plot.type),width=pdf.width,height=pdf.height)
	}else{
		return(p)
	}
}


sigGeneVennPlot <- function(v.list,background.list,out.prefix,col.venn=NULL,fill.venn=NULL,venn.alpha=c(1,0.7,0.7))
{

    for(i in seq_along(v.list)){
	v.list[[i]] <- intersect(v.list[[i]],background.list)
    }
    n.sample <- length(unique(background.list))
    cat(sprintf("n.sample: %d\n",n.sample))

    library("VennDiagram")
    library("gridBase")
    pdf(sprintf("%s.venn.pdf",out.prefix),width=7,height=7)
    if(length(v.list)==2){
        opar <- par(mar=c(7,7,7,7),cex.lab=1.5,cex.main=1.5,xpd=T)
    }else{
        opar <- par(mar=c(2,2,2,2),cex.lab=1.5,cex.main=1.5,xpd=T)
    }
    plot.new()
    #title(main="venn",sub=sprintf("p=%s",p.value),cex=1.2)
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    ##venn.plot <- venn.diagram(v.list,filename = NULL, cat.cex=1.5,cex=1.5,
    venn.plot <- venn.diagram(v.list,filename = NULL, cat.cex=1.5,cex=1.5,
                              hyper.test=T,total.population=n.sample,lower.tail=F,sub.pos=c(0.5,0.25),sub.cex=1.5,
							  col=if(is.null(col.venn)) "black" else col.venn,fill=col.venn,
							  alpha=venn.alpha[seq_len(length(v.list))],
                              margin=0.2,cat.dist=0.15,na = "remove")
    grid.draw(venn.plot)
    par(opar)
    dev.off()
}

sigGeneVennTable <- function(gene.tb,cmp,only.sig=T)
{
	col.out <- c("comb.ES","comb.ES.sd","comb.Z","comb.p","comb.padj","sig","sig.cate")
	mapping.TF.tb <- unique(gene.tb[,c("geneID","geneSet.TF"),with=F])
	mapping.TF.vec <- structure(mapping.TF.tb$geneSet.TF,names=mapping.TF.tb$geneID)
	out.tb <- dcast(gene.tb[meta.cluster %in% cmp,c("geneID","meta.cluster",col.out),with=F],
					geneID~meta.cluster,value.var=col.out)
	out.tb[,geneSet.TF:=mapping.TF.vec[geneID]]
	f.gene <- apply(out.tb[,grepl("^sig_",colnames(out.tb),perl=T),with=F],1,
					function(x){ any(x==T,na.rm=T) })
	if(only.sig==T){
		out.tb <- out.tb[f.gene,]
	}
	return(out.tb)
}

sigGeneHeatmap <- function(out.prefix,gene.desc.top,sce.pb,gene.to.show.tb,value.var="comb.ES",colSet=list(),...)
{
    ### show combined ES
    g.plot.tb <- dcast(gene.desc.top[geneID %in% gene.to.show.tb$geneID,c("geneID","meta.cluster","comb.ES"),with=F],
					       geneID~meta.cluster, value.var=value.var)

    g.plot.mtx <- as.matrix(g.plot.tb[,-1])
    rownames(g.plot.mtx) <- g.plot.tb[[1]]
    g.plot.mtx <- g.plot.mtx[gene.to.show.tb$geneID,]
    g.plot.mtx[g.plot.mtx < -1] <- -1
    g.plot.mtx[g.plot.mtx >  1] <- 1
    cuts <- c(-1,-0.5,-0.15,0.15,0.5,1)
    bin.values <- c("< -0.5","[0.5,-0.15)","[-0.15,0.15]","(0.15,0.5]",">0.5")
    bin.values <- factor(bin.values,levels=bin.values)
    ##g.plot.mtx.bin <- matrix(bin.values[findInterval(g.plot.mtx, cuts,rightmost.closed=T)],
    g.plot.mtx.bin <- matrix(findInterval(g.plot.mtx, cuts,rightmost.closed=T),
						       ncol=ncol(g.plot.mtx))
    colnames(g.plot.mtx.bin) <- colnames(g.plot.mtx)
    rownames(g.plot.mtx.bin) <- rownames(g.plot.mtx)
    head(g.plot.mtx,n=3)
    head(g.plot.mtx.bin,n=3)
    sscVis::plotMatrix.simple(g.plot.mtx.bin,out.prefix=sprintf("%s.slim.bin",out.prefix),
			      col.ht=structure(rev(brewer.pal(5,name="RdBu")), names=1:5),
			      par.legend=list(labels=rev(bin.values), at=5:1),
			      row.split=gene.to.show.tb$Group,
			      pdf.width=if(ncol(g.plot.mtx)<20) 6.5 else 8,
			      par.heatmap=list(cex.column=0.8,border=T,
					       row_gap = unit(0, "mm"),
					       row_title_rot=0),
			      exp.name="comb.ES",...)
    sscVis::plotMatrix.simple(g.plot.mtx,out.prefix=sprintf("%s.slim",out.prefix),
			      row.split=gene.to.show.tb$Group,
			      palatte=rev(brewer.pal(n = 7,name = "RdBu")),
			      z.lo=-0.6,z.hi=0.6,
			      par.legend=list(at = seq(-0.6,0.6,0.15)),
			      pdf.width=if(ncol(g.plot.mtx)<20) 6.5 else 8,
			      par.heatmap=list(cex.column=0.8,border=T,
					       row_gap = unit(0, "mm"),
					       row_title_rot=0),
			      exp.name="comb.ES",...)

    ### show every study
    ssc.plot.heatmap(sce.pb,out.prefix=out.prefix,
		     columns="meta.cluster",columns.order="meta.cluster",
		     colSet=colSet,
		     gene.desc=gene.to.show.tb,
		     row.split=gene.to.show.tb$Group,
		     column.split=sce.pb$meta.cluster,
		     par.heatmap=list(row_gap = unit(0, "mm"),
				      column_gap = unit(0, "mm"),
				      row_title_rot = 0,
				      column_title_gp = gpar(fontsize = 0),
				      border = TRUE),
		     pdf.width=20,pdf.height=10,do.scale=F,
		     palette.name="RdBu",
		     z.lo=-0.5,z.hi=0.5,z.step=0.25,
		     do.clustering.row=F,
		     do.clustering.col=T,...)
}

changeSomeNames <- function(obj,col.mcls="meta.cluster",col.ctype="cancerType",col.dataset="dataset")
{
    #### change cluster name here
    if(!is.null(obj[[col.mcls]])){
#		obj[[sprintf("%s.old",col.mcls)]] <- obj[[sprintf("%s",col.mcls)]]
#		obj[[col.mcls]][ obj[[sprintf("%s.old",col.mcls)]]=="CD8.c05.Tem.TNF" ] <- "CD8.c05.Tact.TNF"
#		obj[[col.mcls]][ obj[[sprintf("%s.old",col.mcls)]]=="CD4.c09.Trm.CAPG" ] <- "CD4.c09.Tm.CAPG"
#		obj[[col.mcls]][ obj[[sprintf("%s.old",col.mcls)]]=="CD4.c04.Tm.TNF" ] <- "CD4.c04.Tact.TNF"
    }

    #### change cancerType here
    if(!is.null(obj[[col.ctype]])){
	    obj[[sprintf("%s.old",col.ctype)]] <- obj[[sprintf("%s",col.ctype)]]
	    obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="BC" ] <- "BRCA"
	    ##obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="LUNG" ] <- "NSCLC"
	    obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="LUNG" ] <- "LC"
	    obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="Melanoma" ] <- "MELA"
	    obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="BM" & obj[[col.dataset]]=="AML.PeterVanGalen2019" ] <- "AML"
	    obj[[col.ctype]][ obj[[col.dataset]]=="MELA.MichaelADurante2020" ] <- "UVM"
    }

    #### change dataset here
    if(!is.null(obj[[col.dataset]])){
	    obj[[sprintf("%s.old",col.dataset)]] <- obj[[sprintf("%s",col.dataset)]]
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CRC.ZiyiLi10X" ] <- "CRC.LeiZhang2020.10X"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="HCC.YaoHeSS2" ] <- "HCC.QimingZhang2019.SS2"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="HCC.YaoHe10X" ] <- "HCC.QimingZhang2019.10X"
	    ######## ???
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CHOL.YaoHe10X" ] <- "HCC.QimingZhang2019.10X"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CHOL.LichunMa2019" ] <- "LIHC.LichunMa2019"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="HCC.LichunMa2019" ] <- "LIHC.LichunMa2019"
	    #########
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.Diether2018" ] <- "LC.DietherLambrechts2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.QianqianSong2019" ] <- "LC.QianqianSong2019"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.RapolasZilionis2019" ] <- "LC.RapolasZilionis2019"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.AshleyMaynard2020" ] <- "LC.AshleyMaynard2020"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="Melanoma.LivnatJerby-Arnon2018" ] <- "MELA.LivnatJerby-Arnon2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="Melanoma.MosheSade-Feldman2018" ] <- "MELA.MosheSade-Feldman2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="Melanoma.HanjieLi2018" ] <- "MELA.HanjieLi2019"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="MELA.MichaelADurante2020" ] <- "UVM.MichaelADurante2020"
	    ##obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="STAD.BoxiKang2019" ] <- "STAD.BoxiKang2020"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="STAD.BoxiKang2019" ] <- "STAD.BoxiKang2021"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="NPC.YangLiu2020" ] <- "NPC.YangLiu2021"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="HCC.zhangLabSS2" ] <- "HCC.ChunhongZheng2017"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.zhangLabSS2" ] <- "LC.XinyiGuo2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CRC.zhangLabSS2" ] <- "CRC.LeiZhang2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.Elham2018.10X" ] <- "BRCA.ElhamAzizi2018.10X"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.Elham2018.Indrop" ] <- "BRCA.ElhamAzizi2018.InDrop"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.Peter2018" ] <- "BRCA.PeterSavas2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.zhangLab5P" ] <- "BRCA.zhangLab5P"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CHOL.zhangLabSS2" ] <- "CHOL.thisStudy"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="NPC.XiliangWang2019" ] <- "NPC.YangLiu2021"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.AshleyMLaughney2020" ] <- "LC.AshleyMLaughney2020"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.NayoungKim2020" ] <- "LC.NayoungKim2020"
	    obj[[col.dataset]] <- gsub("zhangLab5P","thisStudy",obj[[col.dataset]])
    }
    
    return(obj)
}


correctCellInfo <- function(cellInfo.tb)
{
	#cellInfo.tb[,sampleID:=""]
	#cellInfo.tb[,treatment:=""]
	cellInfo.tb$sampleID <- ""
	cellInfo.tb$treatment <- ""
	cellInfo.tb[dataset=="AML.PeterVanGalen2019",sampleID:=batchV]
	cellInfo.tb[dataset=="AML.PeterVanGalen2019" & grepl("-D0",batchV),treatment:="baseline"]
	cellInfo.tb[dataset=="AML.PeterVanGalen2019" & !grepl("-D0",batchV),treatment:="post.treatment"]
	cellInfo.tb[dataset=="AML.PeterVanGalen2019" & grepl("^BM",batchV),treatment:="normal"]

	cellInfo.tb[grepl("^BCC.KathrynEYost2019",dataset),sampleID:=libraryID]
	cellInfo.tb[grepl("^BCC.KathrynEYost2019",dataset) & grepl("\\.pre\\b",sampleID),
				treatment:="baseline"]
	cellInfo.tb[grepl("^BCC.KathrynEYost2019",dataset) & grepl("\\.post\\b",sampleID),
				treatment:="post.treatment"]

	cellInfo.tb[grepl("^SCC.KathrynEYost2019",dataset),sampleID:=libraryID]
	cellInfo.tb[grepl("^SCC.KathrynEYost2019",dataset) & grepl("\\.pre\\b",sampleID),
				treatment:="baseline"]
	cellInfo.tb[grepl("^SCC.KathrynEYost2019",dataset) & grepl("\\.post\\b",sampleID),
				treatment:="post.treatment"]

	cellInfo.tb[grepl("^LC.RapolasZilionis2019",dataset),sampleID:=sprintf("%s.%s",patient,loc)]
	cellInfo.tb[grepl("^LC.RapolasZilionis2019",dataset) & patient=="p2",treatment:="post.treatment"]
	cellInfo.tb[grepl("^LC.RapolasZilionis2019",dataset) & patient!="p2",treatment:="baseline"]

	cellInfo.tb[grepl("^LIHC.LichunMa2019",dataset),sampleID:=patient]
	cellInfo.tb[grepl("^LIHC.LichunMa2019",dataset) & 
				patient %in% c("C25","C39","C56","C60","C66",
							   "H21","H23","H28","H30","H38","H34","H65"),
				treatment:="baseline"]
	cellInfo.tb[grepl("^LIHC.LichunMa2019",dataset) &
				patient %in% c("C26","C29","C35","C42","C46","H18","H37"),
				treatment:="post.treatment"]

	cellInfo.tb[grepl("^BLCA.DavidY.Oh2020",dataset) &
		    grepl("^(Anti|Chemo)",patient,perl=T),
		    treatment:="post.treatment"]
	cellInfo.tb[grepl("^BLCA.DavidY.Oh2020",dataset) &
		    grepl("^Untreated",patient,perl=T),
		    treatment:="baseline"]
	cellInfo.tb[grepl("^BLCA.DavidY.Oh2020",dataset),
		    sampleID:=sprintf("%s.%s",patient,loc)]

	tmp.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/LUNG.AshleyMaynard2020/LUNG.AshleyMaynard2020.sinfo.tb")
	tmp.mapping.tb <- structure(tmp.tb$treatment,names=tmp.tb$sampleID)
	cellInfo.tb[dataset=="LC.AshleyMaynard2020",sampleID:=libraryID]
	cellInfo.tb[dataset=="LC.AshleyMaynard2020",treatment:=tmp.mapping.tb[sampleID] ]

	tmp.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/LUNG.AshleyMLaughney2020/LUNG.AshleyMLaughney2020.sampleInfo.txt")
	tmp.mapping.tb <- structure(tmp.tb$treatment,names=tmp.tb$sampleID)
	print(cellInfo.tb[dataset=="LC.AshleyMLaughney2020",][1:2,])
	cellInfo.tb[dataset=="LC.AshleyMLaughney2020",sampleID:=libraryID]
	cellInfo.tb[dataset=="LC.AshleyMLaughney2020",treatment:=tmp.mapping.tb[sampleID] ]

	########### ???
	#f.dataset <- grepl("^LIHC.LichunMa2019",cellInfo.tb$dataset) & cellInfo.tb$cancerType=="CHOL"
	#cellInfo.tb$dataset[f.dataset] <- "CHOL.LichunMa2019"
	#f.dataset <- grepl("^LIHC.LichunMa2019",cellInfo.tb$dataset) & cellInfo.tb$cancerType=="HCC"
	#cellInfo.tb$dataset[f.dataset] <- "HCC.LichunMa2019"
	###########

	tmp.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Melanoma.HanjieLi2018/Melanoma.HanjieLi2018.cellInfo.txt")
	tmp.id.mapping <- structure(tmp.tb$sampleID,names=tmp.tb$cellID)
	cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset),sampleID:=tmp.id.mapping[cellID] ]
	cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset) & 
				patient %in% c("p1","p11","p12","p13","p15","p16","p17",
							   "p18","p19","p21","p23","p24","p25","p26","p3"),
			treatment:="baseline"]
	cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset) & 
				patient %in% c("p10","p2","p20","p27","p28","p4","p5","p6","p8","p9"),
			treatment:="post.treatment"]
	### primary or metastasis
	cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset),loc:="M"]
	cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset) & patient %in% c("p25","p26"),loc:="T"]

	tmp.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Melanoma.LivnatJerby-Arnon2018/Melanoma.LivnatJerby-Arnon.cellInfo.txt")
	tmp.tb[treatment.group=="treatment.naive",treatment.group:="baseline"]
	tmp.treat.mapping <- structure(tmp.tb$treatment.group,names=tmp.tb$cellID)
	cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset),sampleID:=patient]
	cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset),treatment:=tmp.treat.mapping[cellID]]
	cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset) & patient %in% c("Mel129pa","Mel129pb"),patient:="Mel129"]
	cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset),loc:="M"]
	cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset) &
		    patient %in% c("Mel84","Mel129","Mel105"),loc:="T"]

	tmp.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Melanoma.MosheSade-Feldman2018/sampleInfo.txt")
	tmp.id.mapping <- structure(tmp.tb$oID,names=tmp.tb$title)
	tmp.tb[,treatment.group:=""]
	tmp.tb[Timepoint=="Pre",treatment.group:="baseline"]
	tmp.tb[Timepoint=="Post",treatment.group:="post.treatment"]
	tmp.treat.mapping <- structure(tmp.tb$treatment.group,names=tmp.tb$title)
	cellInfo.tb[grepl("^MELA.MosheSade-Feldman2018",dataset),sampleID:=tmp.id.mapping[cellID] ]
	cellInfo.tb[grepl("^MELA.MosheSade-Feldman2018",dataset),treatment:=tmp.treat.mapping[cellID]]
	#### all metastasis ? (immune cells isolated from 48 tumor biopsies taken from 32 metastatic melanoma patients)
	cellInfo.tb[grepl("^MELA.MosheSade-Feldman2018",dataset),loc:="M" ]

	######## treatment naive only samples
	f.cell <- cellInfo.tb[,(grepl("^(HNSCC.SidharthVPuram2017|CHOL.thisStudy|CRC.LeiZhang2018|HCC.ChunhongZheng2017|LC.XinyiGuo2018)",dataset))]
	##cellInfo.tb[f.cell ,libraryID:=sprintf("%s.%s",patient,loc)]
	cellInfo.tb[f.cell, sampleID:=sprintf("%s%s",patient,loc)]
	cellInfo.tb[f.cell, treatment:="baseline"]

	f.cell <- cellInfo.tb[, grepl("^BRCA.ElhamAzizi2018",dataset) & patient %in% c("BC9","BC10","BC11")]
	cellInfo.tb[f.cell,libraryID:=sprintf("%s%s",patient,loc)]
							
	f.cell <- cellInfo.tb[,grepl("^(BRCA.ElhamAzizi2018|BRCA.PeterSavas2018|CRC.LeiZhang2020.10X|HCC.QimingZhang2019.10X|HCC.QimingZhang2019|NPC.YangLiu2021|LC.DietherLambrechts2018|LC.QianqianSong2019|PACA.JunyaPeng2019|RC.MatthewDYoung2018|STAD.BoxiKang2021|UVM.MichaelADurante2020|CRC.Hae-OckLee2020.KUL3|CRC.Hae-OckLee2020.SMC|LC.NayoungKim2020|MM.OksanaZavidij2020|MULT.MichalSlyper2020|MULT.ThomasD.Wu2020)",dataset) | 
			( grepl("thisStudy",dataset) & dataset!="CHOL.thisStudy" ) ]
	cellInfo.tb[f.cell, sampleID:=libraryID]
	cellInfo.tb[f.cell, treatment:="baseline"]


	#### correct cancerType
	cellInfo.tb[dataset=="HCC.QimingZhang2019.10X" & patient=="D20171215",cancerType:="CHOL"]
	####### ???
	##cellInfo.tb[dataset=="HCC.QimingZhang2019.10X" & patient=="D20171215",dataset:="CHOL.QimingZhang2019.10X"]
	#######

	#### add stype
	cellInfo.tb[,stype:=sapply(strsplit(as.character(meta.cluster),"\\."),"[",1)]
	#### add patient.uid
	cellInfo.tb[,patient.uid:=sprintf("%s.%s",dataset,patient)]
	#### samples for frequency analysis
	cellInfo.tb[,usedForFreq:="Y"]
	cellInfo.tb[dataset=="CRC.LeiZhang2020.10X" & patient %in% c("P0408","P0613","P1025","P1026"),usedForFreq:="Y"]
	cellInfo.tb[dataset=="CRC.LeiZhang2020.10X" & patient %in% c("P0410","P0104"),usedForFreq:="SortingMyeloid"]
	cellInfo.tb[dataset=="HCC.QimingZhang2019.SS2" & patient=="D20171109", usedForFreq:="Dup.CrossPlatformBenchmark"]
	cellInfo.tb[grepl("^(CHOL.thisStudy|HCC.ChunhongZheng2017|LC.XinyiGuo2018|CRC.LeiZhang2018)",dataset) &
				stype=="CD4",usedForFreq:="SortingTreg"]
	cellInfo.tb[dataset=="MULT.MichalSlyper2020" & sampleID %in% c("HTAPP-735-SMP-3841_fresh_channel2","HTAPP-285-SMP-751_fresh_channel1"),
		    usedForFreq:="Biopsy"]
	#cellInfo.tb[dataset=="SCC.KathrynEYost2019" & patient=="su010", patient:="" ]
	cellInfo.tb[usedForFreq!="Y",.N,by=c("dataset","patient","stype","usedForFreq","cancerType")][order(dataset,stype),]
## CRC.LeiZhang2020.10X, samples below contain all immune cells (CD45+): "CRC.P0408", "CRC.P0613", "CRC.P1025", "CRC.P1026"; other samples were sorted by other enrichment criteria

## D20171109 present in both HCC.QimingZhang2019.10X and HCC.QimingZhang2019.SS2. It was used to platform benchmark in the original study. Remove D20171109 in HCC.QimingZhang2019.SS2.
## cellInfo.tb[patient=="D20171109",table(dataset)]

## su010 present in both BCC.KathrynEYost2019 and SCC.KathrynEYost2019
## cellInfo.tb[patient=="su010",table(dataset)]

	{
		cellInfo.tb[grepl("AML.PeterVanGalen2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment),]
		cellInfo.tb[grepl("^BCC.KathrynEYost2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment),]
		cellInfo.tb[grepl("^SCC.KathrynEYost2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment),]
		cellInfo.tb[grepl("^LC.RapolasZilionis2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment),]
		cellInfo.tb[grepl("^LIHC.LichunMa2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment,cancerType,patient),]
		cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment,cancerType,patient),]
		cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment,cancerType,patient),]
		cellInfo.tb[grepl("^MELA.MosheSade-Feldman2018",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment,cancerType,patient),]

	}

	return(cellInfo.tb)
}


resetSig <- function(obj)
{
	##### redefine sig
	zs <- assay(obj,"dprime")
	adjP <- assay(obj,"adj.P.Val")
	sig <- assay(obj,"sig")
	zs[1:4,1:3]
	sig[1:4,1:3]
	adjP[1:4,1:3]
	sig[!(adjP < 0.01 & zs > 0.15)] <- 0
	sig[adjP < 0.01 & zs > 0.15] <- 1
	assay(obj,"sig") <- sig
	return(obj)
}


getCancerTypeSpecificGene <- function(obj,cType,mcls)
{
	#cType <- "OV"
	#mcls <- "CD8.c08c09.Tex.CXCL13"
	obj.case <- obj[,obj$cancerType==cType & obj$meta.cluster==mcls]
	obj.control <- obj[,obj$cancerType!=cType & obj$meta.cluster==mcls]
	if(!(ncol(obj.case)>0 && ncol(obj.control)>0)){
		return(NULL)
	}
	f.pos <- rowSums(assay(obj.case,"adj.P.Val") < 0.01)==ncol(obj.case) &
				rowSums(assay(obj.case,"dprime") > 0.15)==ncol(obj.case)
	f.neg <- rowSums( assay(obj.control,"adj.P.Val") > 0.05 | assay(obj.control,"dprime") < 0.1 )==ncol(obj.control)
	f.sel <- f.pos & f.neg
	f.sel.name <- rownames(obj.case)[f.sel]
	rmeans <- rowMeans(assay(obj.case,"dprime"))[f.sel]
	###rmeans.control <- apply(assay(obj.control,"dprime"),1,median)[f.sel]
	es.combi <- sscVis::directEScombi(assay(obj.control,"dprime"), assay(obj.control,"vardprime"))
	f.sel.name.sort <- f.sel.name[order(-rmeans)]
	ret.tb <- data.table(cancerType=cType,meta.cluster=mcls,geneID=f.sel.name.sort,
						 ES.case=rmeans[f.sel.name.sort],
						 ES.control=es.combi[f.sel.name.sort,"comb.ES"])
	return(ret.tb)
}

signatureAmongDatasetPlot <- function(obj,gene.desc.top,mcls,out.prefix,assay.name="sig",k.column=1,
									  col.freq.sig="freq.sig",th.universal=0.75,
									  pmode="",f.gene=NULL,row.gp=NULL,ES.case=NULL,ES.control=NULL)
{
	#mcls <- "CD8.c08c09.Tex.CXCL13"
	o.tb <- gene.desc.top[order(-comb.ES),][meta.cluster==mcls,]
	o.tb$freq.sig <- o.tb[[col.freq.sig]]
	if(nrow(o.tb)<2){ return(NULL) }
	#### heatmap (Tier1)
	##f.gene <- head(o.tb[sig.cate=="Tier1" ,][["geneID"]],n=30)
	#f.gene <- o.tb[sig.cate=="Tier1" ,][["geneID"]]
	if(is.null(f.gene)){
		f.gene <- o.tb[["geneID"]]
	}
	f.sample <- obj$meta.cluster==mcls
	if(assay.name=="dprime")
	{
		###zs <- assay(obj,"dprime")[f.gene,f.sample]
		zs <- assay(obj,assay.name)[f.gene,f.sample,drop=F]
		print(zs[1:4,1:3])
		colnames(zs) <- obj$dataset[f.sample]
		sscVis::plotMatrix.simple(zs,exp.name="ES",
									  mytitle=mcls,
									  z.lo=-1.0,z.hi=1.0,
									  palatte=rev(brewer.pal(n = 7,name = "RdBu")),
									  #palatte=(brewer.pal(n = 7,name = "Blues")),
									  #out.prefix=sprintf("%s.DEG.%s.Tier1.ht.dprime",out.prefix,mcls),
									  out.prefix=sprintf("%s.ht.dprime.%s",out.prefix,mcls),
									  show.number=F,
									  pdf.width = 8, pdf.height = 10,
									  #clust.column="cutreeDynamic",
									  #clust.row="cutreeDynamic",
									  #clust.column="cutree",
									  #clust.row="cutree",
									  waterfall.col=T,
									  row_dend_width = unit(3, "cm"),
									  column_dend_height = unit(3, "cm"),
									  #par.warterfall=list(method.distance="",minClusterSize=5),
									  #par.warterfall=list(method.distance="",k=2),
									  par.warterfall=list(score.alpha = 1.1,do.norm=F),
									  show.dendrogram=T)
	}else if(assay.name=="sig")
	{
		gene.universal <- o.tb[sig==T & freq.sig > th.universal,]
		##gene.universal <- o.plot.tb[sig.cate=="Tier1",]
		print(str(f.gene))
		print(summary(f.sample))
		zs <- assay(obj,assay.name)[f.gene,f.sample,drop=F]
		print(dim(zs))
		#print(zs[1:4,1:3])
		row.ann.tb <- data.frame(geneID=rownames(zs),
								 TF.cate=ifelse(rownames(zs) %in% gene.universal$geneID,"universal","diverse"),
								 sig.freq=o.tb$freq.sig[match(rownames(zs),o.tb$geneID)],
								 comb.ES=o.tb$comb.ES[match(rownames(zs),o.tb$geneID)],
								 comb.ES.sd=o.tb$comb.ES.sd[match(rownames(zs),o.tb$geneID)])
		row.ann.tb$TF.cate <- factor(row.ann.tb$TF.cate,levels=c("universal","diverse"))
		rownames(row.ann.tb) <- rownames(zs)
		colnames(zs) <- obj$dataset[f.sample]
		if(pmode=="cancerSpecific")
		{
			#f.sort <- order(rowMax(zs))
			#row.ann.tb <- row.ann.tb[f.sort,]
			#zs <- zs[f.sort,]
			dd <- assay(obj,"dprime")[rownames(zs),,drop=F]
			zs.test <<- zs
			sscVis::plotMatrix.simple(zs,exp.name="sig",
										  mytitle=mcls,
										  #z.lo=-1.0,z.hi=1.0,
										  #palatte=rev(brewer.pal(n = 7,name = "RdBu")),
										  palatte=(brewer.pal(n = 7,name = "Blues")),
										  out.prefix=sprintf("%s.ht.sig.%s",out.prefix,mcls),
										  show.number=F,par.legend=list(color_bar="discrete",at=c(0,1)),
										  pdf.width = 10, pdf.height = 10,
										  #row.split=row.ann.tb$TF.cate,
										  row.split=row.gp,row_title_rot = 0,
										  right_annotation=rowAnnotation(ES.case=anno_barplot(ES.case ,bar_width=1,ylim=c(0,1),
																							  width = unit(2, "cm"),
																							  gp = gpar(fill = "#74C476",col=NA)),
																		 ES.control=anno_barplot(ES.control ,bar_width=1,ylim=c(-0.5,0.5),
																							  width = unit(2, "cm"),
																							  gp = gpar(fill = "#00BFC4",col=NA))),
										  column_gap = unit(0, "mm"),row_gap = unit(0, "mm"),border=T,
										  column_split=if(k.column>=2) k.column else NULL,
										  #clust.column="cutreeDynamic",
										  #clust.row="cutreeDynamic",
										  #clust.column="cutree",
										  #clust.row="cutree",
										  waterfall.col=T,
										  row_dend_width = unit(3, "cm"),
										  column_dend_height = unit(3, "cm"),
										  #par.warterfall=list(method.distance="",minClusterSize=5),
										  #par.warterfall=list(method.distance="",k=k.column),
										  par.warterfall=list(score.alpha = 15,do.norm=F),
										  show.dendrogram=F)
		}else{
			f.sort <- order(-row.ann.tb$sig.freq)
			row.ann.tb <- row.ann.tb[f.sort,]
			zs <- zs[f.sort,]
			
			sscVis:::plotMatrix.simple(zs,exp.name="sig",
										  mytitle=mcls,
										  #z.lo=-1.0,z.hi=1.0,
										  #palatte=rev(brewer.pal(n = 7,name = "RdBu")),
										  palatte=(brewer.pal(n = 7,name = "Blues")),
										  out.prefix=sprintf("%s.ht.sig.%s",out.prefix,mcls),
										  show.number=F,par.legend=list(color_bar="discrete",at=c(0,1)),
										  pdf.width = 8, pdf.height = 10,
										  row.split=row.ann.tb$TF.cate,
										  right_annotation=rowAnnotation(sig.freq=anno_barplot(rowSums(zs)/ncol(zs),
																							   bar_width=1,
																							  width = unit(2, "cm"),
																							  gp = gpar(fill = "#74C476",
																										col=NA))),
										  column_gap = unit(0, "mm"),border=T,
										  column_split=if(k.column>=2) k.column else NULL,
										  #clust.column="cutreeDynamic",
										  #clust.row="cutreeDynamic",
										  clust.column="cutree",
										  #clust.row="cutree",
										  #waterfall.col=T,
										  row_dend_width = unit(3, "cm"),
										  column_dend_height = unit(3, "cm"),
										  #par.warterfall=list(method.distance="",minClusterSize=5),
										  par.warterfall=list(method.distance="",k=k.column),
										  #par.warterfall=list(score.alpha = 1.1,do.norm=F),
										  show.dendrogram=T)
		}
	}
	
}


do.plot.freq.heatmap <- function(dat.plot.a,colSet,mapping.vec,group.var="meta.cluster",
				 dat.plot.b=NULL,out.prefix,k.column=6,k.row=1,
				 fltBySD=F,TH.cor=0.35,ann.bar.r=1,k.corMat=NULL,
				 pdf.width=8,pdf.height=8,...)
{
    library("ComplexHeatmap")
    library("dendextend")

    ht.a.tb <- dcast(dat.plot.a,group.var~donor.var,value.var="freq",fill=0)
    ht.a <- as.matrix(ht.a.tb[,-1])
    rownames(ht.a) <- ht.a.tb[[1]]
    print(head(ht.a[,1:3]))
    f.rm.a <- colSums(ht.a) < 1
    cat(sprintf("samples with meta.clusters' frequency sum < 1 (data A):\n"))
    print(colnames(ht.a)[f.rm.a])
    ht.a <- ht.a[,!f.rm.a]
    dat.in.a <- ht.a
    dat.in <- dat.in.a

    if(!is.null(dat.plot.b))
    {
	ht.b.tb <- dcast(dat.plot.b,group.var~donor.var,value.var="freq",fill=0)
	ht.b <- as.matrix(ht.b.tb[,-1])
	rownames(ht.b) <- ht.b.tb[[1]]
	print(head(ht.b[,1:3]))
	f.rm.b <- colSums(ht.b) < 1
	cat(sprintf("samples with meta.clusters' frequency sum < 1 (data B):\n"))
	print(colnames(ht.b)[f.rm.b])
	ht.b <- ht.b[,!f.rm.b]
	f.samples <- intersect(colnames(dat.in.a),colnames(ht.b))
	dat.in.a <- dat.in.a[,f.samples,drop=F]
	dat.in.b <- ht.b[,f.samples,drop=F]
	dat.in <- rbind(dat.in.a,dat.in.b)
	mcls.color.b <- colSet[[group.var]][rownames(ht.b)]
    }


    if(fltBySD)
    {
	rsd <- rowSds(dat.in)
	th.sd <- quantile(rsd,0.5)
	rn.flt <- rownames(dat.in)[ rsd >= th.sd ]
	dat.in <- dat.in[rn.flt,,drop=F]
	dat.in.a <- dat.in.a[rownames(dat.in.a) %in% rn.flt,,drop=F]
	if(!is.null(dat.plot.b)){
		dat.in.b <- dat.in.b[rownames(dat.in.b) %in% rn.flt,,drop=F]
	}
    }

    mcls.color <- colSet[[group.var]][rownames(dat.in)]
    mcls.color.a <- colSet[[group.var]][rownames(dat.in.a)]
    if(!is.null(dat.plot.b)){
	    mcls.color.b <- colSet[[group.var]][rownames(dat.in.b)]
    }
    
    mcls.color.debug <<- mcls.color
    dat.in.debug <<- dat.in

    #### try PCA
    {
	sce.f <- ssc.build(dat.in)
	rowData(sce.f)$rsd <- rowSds(dat.in)
	rowData(sce.f)$gene.all <- T
	rowData(sce.f)
	sce.f$cancerType <- mapping.vec[colnames(sce.f)]
	sce.f <- ssc.reduceDim(sce.f,method="pca",method.vgene="gene.all",pca.npc=5,
						   tSNE.perplexity=30,
						   seed=12345)
	p <- ssc.plot.tsne(sce.f,columns="cancerType",colSet=colSet,reduced.name="pca")
	ggsave(file=sprintf("%s.pca.cancerType.png",out.prefix),width=6,height=4)
	#p <- ssc.plot.tsne(sce.f,gene=c("CD8.c12.Tex.CXCL13","CD8.c10.Trm.ZNF683","CD8.c15.ISG.IFIT1"),
	#				   par.geneOnTSNE = list(scales = "free",pt.alpha=1),
	#				   colSet=colSet,reduced.name="pca")
	#ggsave(file=sprintf("%s.pca.gene.png",out.prefix),width=10,height=3)
	#p <- ssc.plot.tsne(sce.f,gene=c("CD8.c12.Tex.CXCL13","CD8.c10.Trm.ZNF683","CD8.c15.ISG.IFIT1"),
	#				   par.geneOnTSNE = list(scales = "free",pt.alpha=1),
	#				   colSet=colSet,reduced.name="pca.tsne")
	#ggsave(file=sprintf("%s.pca.gene.tsne.png",out.prefix),width=10,height=3)
    }

    ######## correlation and select variables
    {
	cor.dat <- cor(t(dat.in),method="spearman")
	cor.hclust <- run.cutreeDynamic(dat.in,deepSplit=2, minClusterSize=2)
	cor.hclust$branch <- dendextend::set(cor.hclust$branch,"branches_lwd", 1.5)
	sscVis:::plotMatrix.simple(cor.dat,show.dendrogram=T,
				   out.prefix=sprintf("%s.cor",out.prefix),
				   show.number=F,
				   z.lo = -0.4, z.hi=0.4,
				   clust.column=cor.hclust$branch,
				   clust.row=cor.hclust$branch,
				   exp.name=expression(Corr),
				   row_dend_width = unit(2.8, "cm"),
				   par.heatmap=list(
						    column_dend_height = unit(2.8, "cm")
						    ),
				   palatte=rev(brewer.pal(n = 7, name = "RdBu")),
				   pdf.width=8,pdf.height=8)

	#### filter by TH.cor
	mcls.cor.max <- apply(cor.dat,1,function(x){ max(abs(x[x!=1])) })
	f.cor <- mcls.cor.max >= TH.cor
	cat(sprintf("mcls.cor.max:"))
	print(sort(mcls.cor.max))
	cor.dat.flt <- cor.dat[f.cor,f.cor]
	if(is.null(k.corMat)){
	    cor.hclust.flt <- run.cutreeDynamic(dat.in[f.cor,],deepSplit=2, minClusterSize=2)
	}else{
	    #cor.hclust.flt <- run.cutree(dat.in[f.cor,],method.distance="spearman",k=k.corMat)
	    #cor.hclust.flt <- run.cutree(dat.in[f.cor,],method.distance="cosine",k=k.corMat)
	    #cor.hclust.flt <- run.cutree(dat.in[f.cor,],method.distance="cosine",k=k.corMat,method.hclust="complete")
	    cor.hclust.flt <- run.cutree(dat.in[f.cor,],method.distance="cosine",k=k.corMat)
	}
	cor.hclust.flt$branch <- dendextend::set(cor.hclust.flt$branch,"branches_lwd", 1.5)
	sscVis::plotMatrix.simple(cor.dat.flt,show.dendrogram=T,
				  out.prefix=sprintf("%s.cor.flt",out.prefix),
				  show.number=F,
				  clust.column=cor.hclust.flt$branch,
				  clust.row=cor.hclust.flt$branch,
				  exp.name=expression(Corr),
				  #z.lo = -0.5, z.hi=0.5,
				  #par.legend=list(at=seq(-0.5,0.5,0.25)),
				  z.lo = -0.4, z.hi=0.4,
				  par.heatmap=list(row_dend_width = unit(1.5, "cm"),
						   row_names_gp=gpar(fontsize=10),
						   column_names_gp=gpar(fontsize=10),
						   column_dend_height = unit(1.5, "cm")),
				  palatte=rev(brewer.pal(n = 7, name = "RdBu")),
				  pdf.width=8,pdf.height=7.7)

    }

    if(!is.null(dat.plot.b)){
	axis_param.column <- default_axis_param("column")
	axis_param.column$gp <- gpar(fontsize=12)
	ha.col <- HeatmapAnnotation(cancerType=mapping.vec[colnames(dat.in)],
				    CD4=anno_barplot(t(dat.in.b),gp = gpar(fill = mcls.color.b,col=NA),
						     axis_param=axis_param.column,border=F,
						     bar_width=1, height = unit(3.0*ann.bar.r, "cm")),
				    CD8=anno_barplot(t(dat.in.a), gp = gpar(fill = mcls.color.a,col=NA),
						     axis_param=axis_param.column,border=F,
						     bar_width=1, height = unit(3.0*ann.bar.r, "cm")),
				    col = colSet["cancerType"],
				    #annotation_legend_param = list(),
				    show_legend = T,
				    ###simple_anno_size = unit(1, "cm"))
				    simple_anno_size = unit(0.5, "cm"))
    }else{
	ha.col <- HeatmapAnnotation(cancerType=mapping.vec[colnames(dat.in)],
				    CD8=anno_barplot(t(dat.in.a), gp = gpar(fill = mcls.color.a,col=NA),
						     bar_width=1, height = unit(2.5, "cm")),
				    col = colSet["cancerType"],
				    #annotation_legend_param = list(),
				    show_legend = T,
				    simple_anno_size = unit(1, "cm"))
    }

    makeHTPlot <- function(dat.in,out.prefix)
    {

	###obj.hclust.col <- run.cutree(t(dat.in),k=k.column)
	##obj.hclust.col <- run.cutree(t(dat.in),k=k.column,method.distance="spearman")
	###obj.hclust.col <- run.cutreeDynamic(t(dat.in),deepSplit=1, minClusterSize=2)
	obj.hclust.col <- run.cutree(t(dat.in),k=k.column,method.distance="")
	obj.hclust.row <- run.cutree(dat.in,k=k.row)
	#obj.hclust.col <- run.cutree(t(dat.in),k=k.column,method.distance="cosine")
	#obj.hclust.row <- run.cutree(dat.in,k=k.row,method.distance="cosine")

	obj.hclust.col$branch <- dendextend::set(obj.hclust.col$branch,"branches_lwd", 2)
	obj.hclust.row$branch <- dendextend::set(obj.hclust.row$branch,"branches_lwd", 2.5)
	##sscClust:::plot.matrix.simple(t(scale(t(dat.in))), show.dendrogram=T,
	dat.in.bin <- floor(floor(dat.in * 100)/5)
	dat.in.bin[ dat.in.bin > 9 ] <- 9
	dat.in.bin[1:3,1:4]
	#print(str(dat.in))
	dat.in.z <- t(scale(t(dat.in)))
	dat.in.z[dat.in.z > 3] <- 3
	dat.in.z[dat.in.z < -3] <- -3

	obj.hclust.row.debug <<- obj.hclust.row
	ha.row = rowAnnotation(overall=anno_boxplot(dat.in,outline = FALSE,width = unit(4, "cm"),
						    gp = gpar(fill = mcls.color[rownames(dat.in)])))

	sscVis::plotMatrix.simple(
				  #dat.in.z,
				  #palatte=rev(brewer.pal(n = 7,name = "RdBu")),
				  #par.legend=list(ncol = 1,at=seq(-3,3,1.5)),
				  dat.in.bin,
				  col.ht=rev(structure(viridis::magma(10),names=0:9 )),
				  par.legend=list(ncol = 1,labels=rev(sprintf("%s%%~%s%%",5*(0:9),c(5*(1:9),100) ))),
				  ####col.ht=rev(structure(colorRampPalette(brewer.pal(9,name="Blues"))(10),names=0:9 )),
				  show.dendrogram=T,
				  out.prefix=sprintf("%s",out.prefix),
				  show.number=F,
				  clust.row=obj.hclust.row$branch,
				  clust.column=obj.hclust.col$branch,
				  exp.name=expression(Freq),
				  par.heatmap=list(
						   column_split=if(k.column>=2) k.column else NULL,
						   row_names_gp=gpar(fontsize=12),
						   right_annotation = ha.row,
						   top_annotation = ha.col,
						   column_gap = unit(0.8, "mm"),border=FALSE,
						   row_dend_width = unit(1.5, "cm"),
						   column_dend_height = unit(1.5, "cm")
						   ),
				  #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
				  #palatte=viridis::magma(7),
				  #palatte=rev(brewer.pal(n = 7,name = "RdBu")),
				  pdf.width = pdf.width, pdf.height = pdf.height,...)
	return(list("obj.hclust.col"=obj.hclust.col,"obj.hclust.row"=obj.hclust.row,"cor.dat.flt"=cor.dat.flt))
    }

    ret.list.allRow <- makeHTPlot(dat.in,out.prefix)
    ret.list.fltRow <- makeHTPlot(dat.in[f.cor,,drop=F],sprintf("%s.flt",out.prefix))

    return(list("all"=ret.list.allRow,"flt"=ret.list.fltRow))
}


calProliferationScore <- function(obj,assay.name,gene.prol,out.prefix=NULL,method="mean")
{
    f.gene <- rowData(obj)[,"display.name"] %in% gene.prol
    exp.sub <- as.matrix(assay(obj[f.gene,],assay.name))
    f.zero <- apply(exp.sub,1,function(x){ all(x==0) })
    if(sum(f.zero) > 0){
	    cat(sprintf("Number of gene with value zero in all cells: %d\n",sum(f.zero)))
    }
    exp.sub <- exp.sub[!f.zero,]

    if(method=="mean")
    {
	score.prol <- colMeans(exp.sub)
	dat.score <- classify.outlier(score.prol,out.prefix=out.prefix)
	out.tb <- as.data.table(dat.score$score.cls.tb)
	colnames(out.tb)[1] <- "cellID"
    }else if(method=="AUCell")
    {
	require("AUCell")
	#####
	#f.gene <- rowData(obj)[,"display.name"] %in% gene.prol
	#exp.sub <- as.matrix(assay(obj,assay.name))
	exp.sub <- (assay(obj,assay.name))
	rownames(exp.sub) <- rowData(obj)[,"display.name"]
	f.zero <- apply(exp.sub,1,function(x){ all(x==0) })
	if(sum(f.zero) > 0){
	    cat(sprintf("Number of gene with value zero in all cells: %d\n",sum(f.zero)))
	}
	exp.sub <- exp.sub[!f.zero,]

	#####
	pdf(sprintf("%s.buildRankings.1.pdf",out.prefix),width=7,height=7)
	cells_rankings <- AUCell_buildRankings(exp.sub)
	dev.off()

	geneSets <- list("prol"=intersect(rownames(exp.sub),gene.prol))
	#### geneSets <- GSEABase::GeneSet(genes, setName="geneSet1") # alternative
	cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)

	pdf(sprintf("%s.exploreThresholds.1.pdf",out.prefix),width=7,height=7)
	#par(mfrow=c(3,3))
	set.seed(123)
	cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
	dev.off()

	if("L_k2" %in% rownames(cells_assignment$prol$aucThr$thresholds)){
	    th.prol <- cells_assignment$prol$aucThr$thresholds["L_k2","threshold"]
	}else{
	    th.prol <- 0.09
	}

	##geneSetName <- rownames(cells_AUC)[grep("prol", rownames(cells_AUC))]
	pdf(sprintf("%s.exploreThresholds.2.pdf",out.prefix),width=7,height=7)
	AUCell_plotHist(cells_AUC["prol",], aucThr=th.prol)
	abline(v=th.prol)
	dev.off()
	
	out.tb <- data.table(cellID=colnames(cells_AUC),
			     proliferationScore.bin=as.integer(getAUC(cells_AUC)["prol",]>th.prol),
			     proliferationScore=getAUC(cells_AUC)["prol",])
	out.tb$classification <- out.tb$proliferationScore.bin
	dat.score <- cells_AUC

	my.dat.score <- classify.outlier(getAUC(cells_AUC)["prol",],out.prefix=out.prefix)
	out.tb[,myCls:=my.dat.score$score.cls.tb[cellID,"classification"]]

    }
    return(list("out.tb"=out.tb,"detail"=dat.score))
}


dataOnRDPlot <- function(seu,sce,out.prefix,rd="umap",graph.name="RNA_snn",resolution.vec=seq(0.1,5,0.1))
{

	dir.create(dirname(out.prefix),F,T)
	#resolution.vec <- seq(0.5,3.2,0.1)
	#resolution.vec <- seq(0.1,5,0.1)

	p <- DimPlot(seu, reduction = rd, pt.size = .5, group.by = "dataset")
	ggsave(sprintf("%s.%s.dataset.png",out.prefix,rd),width=11,height=5)
	p <- DimPlot(seu, reduction = rd, pt.size = .5, group.by = "dataset.tech")
	ggsave(sprintf("%s.%s.dataset.tech.png",out.prefix,rd),width=9.5,height=5)
	p <- DimPlot(seu, reduction = rd, pt.size = .5,
				 group.by = "dataset",split.by="dataset",ncol=6)
	ggsave(sprintf("%s.%s.dataset.split.png",out.prefix,rd),width=18,height=14)
	p <- DimPlot(seu, reduction = rd, pt.size = .5,
				 group.by = "dataset.tech",split.by="dataset.tech",ncol=5)
	ggsave(sprintf("%s.%s.dataset.tech.split.png",out.prefix,rd),width=14,height=11)

	################ 
	plot.resolution.list <- list()
	for(t.res in resolution.vec){
		cate.res <- sprintf("%s_res.%s",graph.name,t.res)
		plot.resolution.list[[cate.res]] <- DimPlot(seu, reduction = rd,
													pt.size=0.1,label.size=2,
													label = TRUE,
													group.by=cate.res) +
			NoLegend()
	}
	for(i in seq_len(length(plot.resolution.list)/4))
	{
		pp <- plot_grid(plotlist=plot.resolution.list[((i-1)*4+1):(i*4)],
						ncol = 2,align = "hv")
		save_plot(sprintf("%s.%s.res.%d.png",out.prefix,rd,i),pp, ncol = 2, base_aspect_ratio=0.55)
	}

	plot.resolution.list <- list()
	for(t.res in resolution.vec){
		cate.res <- sprintf("%s_res.%s",graph.name,t.res)
		plot.resolution.list[[cate.res]] <- ssc.plot.tsne(sce,columns = cate.res,
														  reduced.name = sprintf("%s",rd),
														  colSet=list(),size=0.1,label=2,
														  base_aspect_ratio = 1.2)
	}
	for(i in seq_len(length(plot.resolution.list)/4))
	{
		pp <- plot_grid(plotlist=plot.resolution.list[((i-1)*4+1):(i*4)],
						ncol = 2,align = "hv")
		save_plot(sprintf("%s.%s.res.sceStyle.%d.png",out.prefix,rd,i),pp, ncol = 2,
				  base_aspect_ratio=0.9,base_height=5.5)
	}

	makeGeneOnTSNEPlot(sce,rd,out.prefix,geneOnUmap.list=g.GeneOnUmap.list)

	#### density
	ssc.plot.tsne(sce,plotDensity=T,reduced.name=sprintf("%s",rd),
			  out.prefix=sprintf("%s.%s",out.prefix,rd))

}

makeGeneOnTSNEPlot <- function(sce,rd,out.prefix,geneOnUmap.list=g.GeneOnUmap.list,
							   plot.ncol=NULL,plot.nrow=NULL,plot.type="png",
							   plot.width=NULL,plot.height=NULL,...)
{
	dir.create(dirname(out.prefix),F,T)
	## gene on umap
	l_ply(seq_along(geneOnUmap.list),function(i){
		gene.tmp <- intersect(geneOnUmap.list[[i]],rowData(sce)$display.name)
		if(is.null(plot.ncol)){
			plot.ncol <- if(length(gene.tmp)>3) floor(sqrt(length(gene.tmp))+0.5) else 3
		}
		if(is.null(plot.nrow)){
			plot.nrow <- ceiling(length(gene.tmp)/plot.ncol)
		}
		if(is.null(plot.width)){
			plot.width <- if(plot.ncol > 3) 14 else if(plot.ncol>2) 10 else if(plot.ncol>1) 7  else 3.5
		}
		if(is.null(plot.height)){
			plot.height <- if(plot.nrow>3) 11 else if(plot.nrow>2) 8 else if(plot.nrow>1) 5.4 else 2.7
		}
		if(length(gene.tmp)>0){
			p <- ssc.plot.tsne(sce,assay.name="exprs",adjB=NULL,
						  gene=gene.tmp,clamp=c(-0.5,1.5),
						  ##gene=gene.tmp,clamp=c(-0.5,0.5),par.legend=list(breaks=c(-0.5,-0.25,0,0.25,0.5)),
						  p.ncol=plot.ncol,
						  ##par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha = 0.5),
						  par.geneOnTSNE=list(scales="fixed",pt.order="random",pt.alpha = 0.5),
						  reduced.name=sprintf("%s",rd),...)
			ggsave(sprintf("%s.%s.marker.%s.%s",
						   out.prefix,rd,names(geneOnUmap.list)[i],plot.type),
				   width=plot.width,
				   height=plot.height)
				   ####height=if(length(gene.tmp)>9) 11 else if(length(gene.tmp)>6) 8 else if(length(gene.tmp)>3) 5.4 else 2.7)
		}
	},.parallel=T)
}


makeViolinPlot <- function(sce,out.prefix,assay.name="exprs",col.group="ClusterID",
					   geneOnUmap.list=g.GeneOnUmap.list,plot.type="png",clamp=c(-2.5,5),
					   plot.width=NULL,plot.height=NULL,...)
{
	dir.create(dirname(out.prefix),F,T)
	nCls <- length(table(sce[[col.group]]))
	l_ply(seq_along(geneOnUmap.list),function(i){
		gene.tmp <- intersect(geneOnUmap.list[[i]],rowData(sce)$display.name)
		nCol <- if(length(gene.tmp)>14) 2 else 1
		if(is.null(plot.width)){ plot.width <- if(nCls<=12) 6*nCol else 8*nCol }
		if(is.null(plot.height)){ plot.height <- if(length(gene.tmp)>10) 20 else 8 }
		if(length(gene.tmp)>0){
			p <- ssc.plot.violin(sce,assay.name=assay.name,adjB=NULL,
								 p.ncol=nCol,
								 clamp = clamp, gene=gene.tmp, group.var=col.group,...)
			ggsave(sprintf("%s.violin.marker.%s.%s.%s",
						   out.prefix,col.group,names(geneOnUmap.list)[i],plot.type),
				   width=plot.width,
				   height=plot.height)
		}
	},.parallel=T)

}


LISIDistPlot <- function(lisi.pca,lisi.harmony,group=NULL)
{
    lisi.pca$rd <- "PCA"
    colnames(lisi.pca)[3] <- "ClusterID"
    lisi.harmony$rd <- "Harmony"
    colnames(lisi.harmony)[3] <- "ClusterID"
    dat.plot <- rbind(lisi.pca,lisi.harmony)

    dat.plot$cellID <- rownames(dat.plot)
    dat.plot <- as.data.table(dat.plot)
    dat.plot <- melt(dat.plot,id.vars=c("cellID","rd"))
    if(!is.null(group)){
	    dat.plot <- dat.plot[variable %in% group,]
    }
    p <- ggplot(dat.plot,aes(value)) +
		    geom_density(aes(color=rd,group=rd,fill=rd),alpha=0.5,size=0.5) +
		    xlab("LISI") +
		    theme_bw()
    if(length(unique(dat.plot$variable))>1){
	    p <- p + facet_wrap(~variable,ncol=1,scales="free")
    }
    return(p)
}

saveSCEPerDataSet <- function(exp.list.table,meta.tb,out.prefix)
{
	sce.list <- llply(seq_len(nrow(exp.list.table)),function(i){
	  sce <- readRDS(exp.list.table$scefile[i])
	  sce$dataset <- exp.list.table$data.id[i]
	  if(!("norm_exprs" %in% assayNames(sce)))
	  {
		  assay(sce,"norm_exprs") <- assay(sce,"log2TPM")
	  }
	  if(exp.list.table$data.id[i] %in% c("HCC.YaoHe10X","HCC.YaoHeSS2") ){
		  f.gene.tmp <- which(rowData(sce)[,"display.name"]=="MALAT1-ENSG00000251562")
		  rowData(sce)[f.gene.tmp,"display.name"] <- "MALAT1"
	  }
	  sce
	  ##assay.name <- "norm_exprs"
	},.parallel=T)
	names(sce.list) <- exp.list.table$data.id

	#### add meta cluster info to sce.list
	dir.sce <- sprintf("%s/sce",dirname(out.prefix))
	dir.create(dir.sce,F,T)
	#### core
	cn.tb <- as.data.table(ldply(seq_along(sce.list),function(i){
				obj <- sce.list[[i]]
				dataset.id <- names(sce.list)[i]
				loginfo(sprintf("saving sce (%s):",dataset.id))
				meta.search.tb <- as.data.frame(meta.tb[dataset==dataset.id,])
				rownames(meta.search.tb) <- meta.search.tb$cellID
				f.cell <- intersect(colnames(obj),rownames(meta.search.tb))
				obj <- obj[,f.cell]
				##rownames(meta.search.df) <- meta.search.df$cellID
				colData(obj) <- DataFrame(meta.search.tb[f.cell,])
				saveRDS(obj,file=sprintf("%s/%s.sce.rds",dir.sce,dataset.id))

				##### special case
				if(dataset.id=="LIHC.LichunMa2019"){
					for(ctype in c("HCC","CHOL")){
						obj.t <- obj[,obj$cancerType==ctype]
						saveRDS(obj.t,file=sprintf("%s/%s.LichunMa2019.sce.rds",dir.sce,ctype))
					}
				}
				if(dataset.id=="HCC.YaoHe10X"){
					obj$cancerType[obj$patient=="D20171215"] <- "CHOL"
					for(ctype in c("HCC","CHOL")){
						obj.t <- obj[,obj$cancerType==ctype]
						saveRDS(obj.t,file=sprintf("%s/%s.YaoHe10X.mod.sce.rds",dir.sce,ctype))
					}
				}
				if(dataset.id=="MULT.MichalSlyper2020"){
				    for(ctype in c("LC","BC")){
					obj.t <- obj[,obj$cancerType==ctype]
					saveRDS(obj.t,file=sprintf("%s/%s.MichalSlyper2020.sce.rds",dir.sce,ctype))
				    }
				}
				if(dataset.id=="MULT.ThomasD.Wu2020"){
				    for(ctype in c("CRC","UCEC","LC","RC")){
					obj.t <- obj[,obj$cancerType==ctype]
					saveRDS(obj.t,file=sprintf("%s/%s.ThomasD.Wu2020.sce.rds",dir.sce,ctype))
				    }
				}
				#####cellInfo.tb[dataset=="HCC.QimingZhang2019.10X" & patient=="D20171215",cancerType:="CHOL"]

				#obj <- readRDS(file=sprintf("%s/%s.sce.rds",dir.sce,dataset.id))
				data.table(dataset=dataset.id,ncell=ncol(obj))
							},.parallel=T))
	cn.tb
}


zDistPlot <- function(exp.list.table,out.prefix,gene.top.tb,ncores=12)
{
	RhpcBLASctl::omp_set_num_threads(1)
	doParallel::registerDoParallel(cores = ncores)

	sce.list <- llply(seq_len(nrow(exp.list.table)),function(i){
	  sce <- readRDS(exp.list.table$efile[i])
	  sce$dataset <- exp.list.table$data.id[i]
	  if(!("norm_exprs" %in% assayNames(sce)))
	  {
		  assay(sce,"norm_exprs") <- assay(sce,"log2TPM")
	  }
	  sce
	  ##assay.name <- "norm_exprs"
	},.parallel=T)
	names(sce.list) <- exp.list.table$data.id

	gene.list <- gene.top.tb[,c("geneSymbol","meta.cluster")][!duplicated(geneSymbol),]
	gene.list.vec <- structure(gene.list$meta.cluster,names=gene.list$geneSymbol)

	loginfo("merge data ...")
	dat.merged <- as.data.table(ldply(seq_len(nrow(exp.list.table)),function(i){
							data.id <- exp.list.table$data.id[i]
							#sce <- sce.list[[i]]
							#sce <- sce[rowData(sce)$display.name %in% gene.list$geneSymbol,]
							sce <- ssc.scale(sce.list[[i]],gene.symbol=gene.list$geneSymbol,
											 adjB="batchV",do.scale=T)
							dat.x <- assay(sce,"norm_exprs.scale")
							rownames(dat.x) <- rowData(sce)$display.name
							dat.x.tb <- data.table(geneSymbol=rownames(dat.x),
												   meta.cluster=gene.list.vec[rownames(dat.x)],
												   dataset=data.id)
							dat.x.tb <- cbind(dat.x.tb,dat.x)
							dat.x.tb[1:4,1:5]
							dat.x.tb <- melt(dat.x.tb,id.vars=c("geneSymbol","meta.cluster","dataset"))
							dat.x.tb
	},.parallel=ncores))
	loginfo("merge data done.")

	l_ply(seq_len(nrow(gene.list)),function(i){
			  gene.test <- gene.list$geneSymbol[i]
			  mcls <- gene.list$meta.cluster[i]
			  dat.plot <- dat.merged[geneSymbol==gene.test,]
			  p <- ggplot(dat.plot,aes(value)) +
					geom_density(aes(color=dataset,group=dataset)) +
					coord_cartesian(xlim=c(-4,8)) +
					xlab("z score") +
					theme_bw()
			  ggsave(sprintf("%s.%s.%s.png",out.prefix,mcls,gene.test),width=12,height=4)
	},.parallel=T)


}

fetchMetaClusterID2CusterFullName <- function(col.use="cluster.name.full")
{
    ##ret.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/zhangLab.10X/inte.metaClust.20200111/list/name.conversion.txt")
    ret.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/zhangLab.10X/inte.metaClust.20201223/list/name.conversion.txt")
    ret.vec <- structure(ret.tb[[col.use]],names=ret.tb$meta.cluster)
    return(ret.vec)
}

fetechGeneSetList <- function(gmode=0)
{
    hpa.tb <- fread("/lustre1/zeminz_pkuhpc/00.database/transcriptionalRegulatory/humanProteinAtlas/version.20200107/proteinatlas.gid.tsv")
    
    colnames(hpa.tb) <- make.names(colnames(hpa.tb))
    gset.list <- list("TF"=fread("/lustre1/zeminz_pkuhpc/zhenglt/00.database/marker/transcription.factor.activity.GO0003700.list")$geneSymbol,
		      "kinase"=fread("/lustre1/zeminz_pkuhpc/zhenglt/00.database/marker/kinase.activity.GO0016301.list")$geneSymbol,
		      "regulation.epigenetic"=c(fread("/lustre1/zeminz_pkuhpc/zhenglt/00.database/marker/regulation.of.gene.expression.epigenetic.GO0040029.list")$geneSymbol,
						fread("/lustre1/zeminz_pkuhpc/zhenglt/00.database/marker/regulation.of.molecular.function.epigenetic.GO0040030.list")$geneSymbol),
		      "secreted"=hpa.tb[grepl("secreted proteins",Protein.class,perl=T),][["Gene"]],
		      "membrane"=hpa.tb[grepl("membrane proteins",Protein.class,perl=T),][["Gene"]],
		      "enzyme"=hpa.tb[grepl("Enzymes",Protein.class,perl=T),][["Gene"]],
		      "drug.target.FDA"=hpa.tb[grepl("FDA approved drug targets",Protein.class,perl=T),][["Gene"]],
		      "drug.target.potential"=hpa.tb[grepl("Potential drug targets",Protein.class,perl=T),][["Gene"]]
		      )
    if(gmode==1){
	    gset.list[["cytokine"]] <- fread("/lustre1/zeminz_pkuhpc/zhenglt/00.database/marker/cytokine.activity.GO0005125.list")$geneSymbol
	    gset.list[["cytokineReceptor"]] <- fread("/lustre1/zeminz_pkuhpc/zhenglt/00.database/marker/cytokine.receptor.activity.GO0004896.list")$geneSymbol
	    gset.list[["chemokine"]] <- fread("/lustre1/zeminz_pkuhpc/zhenglt/00.database/marker/chemokine.activity.GO0008009.list")$geneSymbol
	    gset.list[["chemokineReceptor"]] <- fread("/lustre1/zeminz_pkuhpc/zhenglt/00.database/marker/chemokine.receptor.activity.GO0004950.list")$geneSymbol
	    gset.list[["CD.molecular"]] <- fread("/lustre1/zeminz_pkuhpc/zhenglt/00.database/marker/HGNC.CD.molecules.txt")$Approved.symbol
    }
    return(gset.list)
}

gen.gsea.script <- function(gene.desc.top,sh.dir,out.prefix,db.file,
							group.var="meta.cluster",rank.var="comb.ES",
							gsea.scoring.scheme="classic",
							gsea.max=500)
{
	for(mcls in unique(sort(gene.desc.top[[group.var]])) ){
		filename.geneSymbol <- sprintf("%s.geneSymbol.%s.rnk",out.prefix,mcls)
		dir.create(dirname(filename.geneSymbol),F,T)
		f.mcls <- gene.desc.top[[group.var]]==mcls
		write.table(gene.desc.top[f.mcls,c("geneSymbol",rank.var),with=F],
					file = filename.geneSymbol,
					quote = F,sep = "\t",row.names = F,col.names = F)

		odir <- sprintf("%s.%s.%s.GSEA",out.prefix,mcls,gsea.scoring.scheme)
		dir.create(odir,showWarnings = F,recursive = T)
		dir.create(sh.dir,showWarnings = F,recursive = T)
		sink(file = sprintf("%s/external.GSEA.%s.job",sh.dir,mcls))
		cat("#!/bin/bash\n")
		cat("#SBATCH -p cn-long\n")
		cat("#SBATCH -N 1\n")
		cat("#SBATCH --ntasks-per-node=1\n")
		cat(sprintf("#SBATCH -o gsea.%s.%%j.out\n",mcls))
		cat(sprintf("#SBATCH -e gsea.%s.%%j.err\n",mcls))
		cat("#SBATCH --no-requeue\n")
		cat("#SBATCH -A zeminz_g1\n")
		cat("#SBATCH --qos=zeminzcnl\n")
		cat("source /lustre1/zeminz_pkuhpc/zhenglt/.bashrc_zhenglt\n")
		cmd.str <- sprintf("/lustre1/zeminz_pkuhpc/01.bin/GSEA/GSEA_Linux_4.0.3/gsea-cli.sh GSEAPreranked -gmx %s -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk %s -scoring_scheme %s -rpt_label %s -create_svgs true -include_only_symbols true -make_sets true -plot_top_x 100 -rnd_seed timestamp -set_max %s -set_min 5 -zip_report false -out %s \n",
						   db.file,
						   filename.geneSymbol,
						   gsea.scoring.scheme,
						   sprintf("%s",mcls),
						   gsea.max,odir)
		cat(cmd.str)
		#ret <- system(cmd.str)
		sink()
	}
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
    count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
    sum.col <- colSums(count.dist)
    sum.row <- rowSums(count.dist)
    count.dist.tb <- as.data.frame(count.dist)
    setDT(count.dist.tb,keep.rownames=T)
    count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
    colnames(count.dist.melt.tb) <- c("rid","cid","count")
    count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
						   this.row <- count.dist.melt.tb$rid[i]
						   this.col <- count.dist.melt.tb$cid[i]
						   this.c <- count.dist.melt.tb$count[i]
						   other.col.c <- sum.col[this.col]-this.c
						   this.m <- matrix(c(this.c,
								      sum.row[this.row]-this.c,
								      other.col.c,
								      sum(sum.col)-sum.row[this.row]-other.col.c),
								    ncol=2)
						   res.test <- fisher.test(this.m)
						   data.frame(rid=this.row,
							      cid=this.col,
							      p.value=res.test$p.value,
							      OR=res.test$estimate)
					       }))
    count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
								    by=c("rid","cid"))
    count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
    #count.dist.melt.ext.tb[adj.p.value < 0.05,]
    #count.dist.melt.ext.tb[p.value < 0.05 & OR > 0,]
    
    return(count.dist.melt.ext.tb)

}

do.tissueDist <- function(cellInfo.tb,out.prefix,pdf.width=3,pdf.height=5,verbose=0)
{
    library("Startrac")
    dir.create(dirname(out.prefix),F,T)

    cellInfo.tb[,meta.cluster:=as.character(meta.cluster)]
    loc.avai.vec <- unique(cellInfo.tb[["loc"]])
    loc.avai.vec <- intersect(c("P","N","T"),loc.avai.vec)
    count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
    freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
    freq.dist.bin <- floor(freq.dist * 100 / 10)
    print(freq.dist.bin)

    {
	count.dist.melt.ext.tb <- test.dist.table(count.dist)
	p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
	OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
	OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
	rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
    }

    startrac.dist <- unclass(calTissueDist(cellInfo.tb,colname.cluster="meta.cluster"))
    startrac.dist <- startrac.dist[,loc.avai.vec]
    
    cuts <- c(0, 0.8, 1.2,Inf)
    startrac.dist.bin.values <- factor(c("-", "+/-", "+"),levels=c("-", "+/-", "+"))
    startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
							       ncol=ncol(startrac.dist))
    colnames(startrac.dist.bin) <- colnames(startrac.dist)
    rownames(startrac.dist.bin) <- rownames(startrac.dist)

#	sscClust:::plot.matrix.simple(freq.dist.bin,
#								  col.ht=rev(structure(colorRampPalette(brewer.pal(9,name="Blues"))(10),
#												   names=0:9 )),
#								  par.legend=list(labels=rev(sprintf("%s%%~%s%%",10*(0:9),c(10*(1:9),100) )),
#												  at=0:9),
#								  out.prefix=sprintf("%s.freq.dist",out.prefix),
#								  show.number=F,clust.row=T,exp.name=expression(italic(Freq)),
#								  #palatte=(brewer.pal(n = 7,name = "Blues")),
#								  #palatte=viridis::viridis(7),
#								  pdf.width = 4.5, pdf.height = pdf.height)

#	sscClust:::plot.matrix.simple(startrac.dist.bin,
#								  col.ht=rev(structure(viridis::viridis(3),
#													   names=levels(startrac.dist.bin.values))),
#								  out.prefix=sprintf("%s.startrac.dist.bin",out.prefix),
#								  show.number=F,clust.row=T,exp.name=expression(italic(R)[o/e]),
#								  pdf.width = pdf.width, pdf.height = pdf.height)

    sscVis::plotMatrix.simple(startrac.dist,
							      out.prefix=sprintf("%s.startrac.dist",out.prefix),
							      show.number=F,
							      clust.row=T,
							      #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
							      exp.name=expression(italic(R)[o/e]),
							      z.hi=2,
							      #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
							      palatte=viridis::viridis(7),
							      pdf.width = 4, pdf.height = pdf.height)

    sscVis::plotMatrix.simple(OR.dist.mtx,
							      out.prefix=sprintf("%s.OR.dist",out.prefix),
							      show.number=F,
							      #clust.row=T,
							      #par.legend=list(color_bar = "discrete",at=seq(0,4,0.5)),
							      waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
							      exp.name=expression(italic(OR)),
							      z.hi=4,
							      #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
							      palatte=viridis::viridis(7),
							      pdf.width = 4, pdf.height = pdf.height)
    if(verbose==1){
	return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
		    "p.dist.tb"=p.dist.tb,
		    "OR.dist.tb"=OR.dist.tb,
		    "OR.dist.mtx"=OR.dist.mtx))
    }else{
	return(OR.dist.mtx)
    }

}

run.harmony <- function(exp.list.table,out.prefix,gene.exclude.file,ncores=12,npc=15,TH.gene.occ=1,
							   contamination.vec=NULL,cor.MALAT1=F,cor.DIG=T)
{
	RhpcBLASctl::omp_set_num_threads(1)
	doParallel::registerDoParallel(cores = ncores)

	env.misc <- loadToEnv(gene.exclude.file)
	cat(sprintf("black gene list:\n"))
	env.misc$all.gene.ignore.df %>% head %>% print

	cat(sprintf("contamination.vec:\n"))
	print(str(contamination.vec))

	loginfo("read data ...")
	seu.list <- llply(seq_len(nrow(exp.list.table)),function(i){
						  seu <- readRDS(exp.list.table$seufile[i])
						  data.id <- exp.list.table$data.id[i]
						  seu$dataset <- data.id
						  seu 
	},.parallel=T)
	names(seu.list) <- exp.list.table$data.id

	sce.list <- llply(seq_len(nrow(exp.list.table)),function(i){
	  sce <- readRDS(exp.list.table$efile[i])
	  sce$dataset <- exp.list.table$data.id[i]
	  if(!("norm_exprs" %in% assayNames(sce)))
	  {
		  assay(sce,"norm_exprs") <- assay(sce,"log2TPM")
	  }
	  if(exp.list.table$data.id[i] %in% c("HCC.YaoHe10X","HCC.YaoHeSS2") ){
		  f.gene.tmp <- which(rowData(sce)[,"display.name"]=="MALAT1-ENSG00000251562")
		  rowData(sce)[f.gene.tmp,"display.name"] <- "MALAT1"
	  }
	  sce
	  ##assay.name <- "norm_exprs"
	},.parallel=T)
	names(sce.list) <- exp.list.table$data.id

	#### get genes used for PCA. input from  sce.list and gene.de.list
	{

		gene.occ <- table(unlist(sapply(sce.list,function(x){ unname(rowData(x)[,"display.name"]) })))
		gene.occ <- gene.occ/length(sce.list)
		gene.occ <- sort(gene.occ,decreasing=T)
		gene.common <- names(gene.occ)[gene.occ >= TH.gene.occ]
		gene.common.all <- names(gene.occ)[gene.occ >= 0.85]
		loginfo(sprintf("total %d common genes obtain ",length(gene.common)))

		l_ply(names(sce.list),function(x){
				  ##cat(sprintf("gene %s in dataset %s: %s\n","CX3CR1",x,"CX3CR1" %in% rowData(sce.list[[x]])[,"display.name"]))
				  #cat(sprintf("gene %s in dataset %s: %s\n","MALAT1",x,"MALAT1" %in% rowData(sce.list[[x]])[,"display.name"]))
				  #cat(sprintf("gene %s in dataset %s: %s\n","IL17A",x,"IL17A" %in% rowData(sce.list[[x]])[,"display.name"]))
				  cat(sprintf("gene %s in dataset %s: %s\n","CCR8",x,"CCR8" %in% rowData(sce.list[[x]])[,"display.name"]))
		},.parallel=T)

		gene.de.list <- list()
		for(i in seq_len(nrow(exp.list.table))){
			id.d <- exp.list.table$data.id[i]
			ifile <- exp.list.table$ifile[i]
			de.out <- readRDS(ifile)
			gene.de.list[[id.d]] <- de.out$all
			gene.de.list[[id.d]]$geneID <- gene.de.list[[id.d]]$geneSymbol
		}
		names(gene.de.list) <- exp.list.table$data.id

		gene.rank.tb <- as.data.table(ldply(names(gene.de.list),function(x){
												ret.tb <- unique(gene.de.list[[x]][,c("geneID","F.rank")])
												ret.tb$dataset.id <- x
												ret.tb <- ret.tb[geneID %in% gene.common,]
												return(ret.tb)
														 }))
		gene.rank.tb <- dcast(gene.rank.tb,geneID~dataset.id,value.var="F.rank",fill=1)
		gene.rank.tb$median.F.rank <- rowMedians(as.matrix(gene.rank.tb[,-c("geneID"),with=F]))
		gene.rank.tb <- gene.rank.tb[order(median.F.rank),]
		##gene.rank.tb <- gene.rank.tb[geneID %in% gene.common,]
		#rowData(sce.pb)$median.F.rank <- gene.rank.tb[["median.F.rank"]][match(rownames(sce.pb),gene.rank.tb$geneID)]

		#### ASH1L MALAT1  CLDN15  ZNF117 NF1 ZNF546 HLA-A HLA-B HLA-C HLA-E PDIA3 IL2RG PTPRC B2M 
		f.gene.blackList <- (gene.rank.tb$geneID %in% env.misc$all.gene.ignore.df$geneSymbol) |
								grepl("^RP[LS]",gene.rank.tb$geneID,perl=T) |
								gene.rank.tb$geneID=="MALAT1"
		gene.de.common <- head(gene.rank.tb[!f.gene.blackList,][["geneID"]],n=1500)
		###gene.de.common <- head(gene.rank.tb$geneID,n=1500)

		cat(sprintf("RP gene in gene.de.common:\n"))
		print(gene.de.common[grepl("^RP[LS]",gene.de.common,perl=T)])
		
	}

	#### cellInfo
	{
		seu.res <- 2
		meta.tb <- ldply(seq_len(nrow(exp.list.table)),function(i){
								m.tb <- fread(exp.list.table$cfile[i])
								m.tb$dataset <- exp.list.table$data.id[i]
								a.res <- ifelse(exp.list.table$platform[i]=="SmartSeq2",
												sprintf("RNA_snn_res.%s",seu.res),
												sprintf("SCT_snn_res.%s",seu.res))
								m.tb$ClusterID <- sprintf("%s.%s",m.tb$dataset,m.tb[[a.res]])
								m.tb$dataset <- exp.list.table$data.id[i]
								o.tb <- m.tb[,c("patient","cellID","libraryID","cancerType","loc","batchV","TCR","dataset","ClusterID")]
								##o.tb <- m.tb[,.(NCell=.N),by=c("cancerType","dataset","ClusterID")]
								return(o.tb)
							   })
		rownames(meta.tb) <- sprintf("%s.%s",meta.tb$dataset,meta.tb$cellID)
		meta.tb$dataset.tech <- gsub("^.+?\\.","",meta.tb$dataset)
	}

	loginfo("merge data ...")
	dat.merged <- do.call(cbind,llply(seq_len(nrow(exp.list.table)),function(i){
							data.id <- exp.list.table$data.id[i]
							seu <- seu.list[[i]]
							sce <- sce.list[[i]]
							cat(sprintf("all(colnames(seu)==colnames(sce)) in %s: %s\n",data.id,
										all(colnames(seu)==colnames(sce))))
							dat.x <- assay(sce,"norm_exprs")
							if(!"seu.id" %in% colnames(rowData(sce))){
								rowData(sce)[,"seu.id"] <- gsub("_","-",rowData(sce)[,"display.name"])
							}
							rownames(dat.x) <- rowData(sce)[,"seu.id"]
							seu.x <- CreateSeuratObject(dat.x,project="panC", meta.data=seu[[]])
							adj.cov <- c("S.Score","G2M.Score","DIG.Score1")
							if("percent.mito" %in% colnames(seu.x[[]])){
								adj.cov <- c("percent.mito",adj.cov)
							}
						    nBatch <- length(table(seu.x$batchV))
							if(nBatch>1){
								adj.cov <- c("batchV",adj.cov)
							}
							seu.x <- ScaleData(seu.x,do.scale=T,features=gene.de.common,vars.to.regress = adj.cov,verbose=F)
							dat.x <- GetAssayData(seu.x,"scale.data")
							colnames(dat.x) <- sprintf("%s.%s",data.id,colnames(dat.x))
							#dat.x[1:4,1:5]
							return((dat.x[gene.de.common,]))
	},.parallel=ncores))

	rownames(dat.merged) <- gsub("_","-",rownames(dat.merged))
	dim(dat.merged)
	dat.merged[1:4,1:5]

	seu.list <- NULL
	sce.list <- NULL
	gc()

	loginfo("begin Seurat pipeline ...")
	##### seurat object
	seu.merged <- CreateSeuratObject((dat.merged),project="panC",meta.data=meta.tb[colnames(dat.merged),])
	##seu.merged$dataset.tech <- gsub("^.+?\\.","",seu.merged$dataset)
	seu.merged <- SetAssayData(seu.merged,"scale.data",dat.merged)

	seu.merged <- RunPCA(seu.merged,features = rownames(seu.merged), npcs = npc, verbose = FALSE)
	seu.merged <- RunUMAP(seu.merged,reduction="pca",dims=1:npc)
	seu.merged <- RunHarmony(seu.merged, c("dataset"))
	seu.merged <- RunUMAP(seu.merged,reduction = "harmony",reduction.name = "harmony.umap", dims = 1:npc)
	#seu.merged <- RunTSNE(seu.merged,reduction ="pca",tsne.method="FIt-SNE",dims=1:npc,reduction.name="tsne")
	#seu.merged <- RunTSNE(seu.merged,reduction ="harmony",tsne.method="FIt-SNE",dims=1:npc,reduction.name="harmony.tsne")

	loginfo("begin Rtsne ...")
	seu.merged <- RunTSNE(seu.merged,reduction ="pca",tsne.method="Rtsne",dims=1:npc,reduction.name="tsne.Rtsne")
	seu.merged <- RunTSNE(seu.merged,reduction ="harmony",tsne.method="Rtsne",dims=1:npc,reduction.name="harmony.tsne.Rtsne")
	loginfo("end Rtsne ...")

	resolution.vec <- seq(0.5,2.4,0.1)
	seu.merged <- FindNeighbors(seu.merged,reduction = "harmony", dims = 1:npc)
	#seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,algorithm=3)
	seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,algorithm=1)

	seu.merged <- FindNeighbors(seu.merged,reduction = "pca", dims = 1:npc,graph.name="RNA_pca_snn")
	#seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,algorithm=3, graph.name="RNA_pca_snn")
	seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,algorithm=1, graph.name="RNA_pca_snn")

	sce.merged <- ssc.build((dat.merged))
	colData(sce.merged) <- DataFrame(meta.tb[colnames(sce.merged),])
	rowData(sce.merged)$gene.de.common <- T

	all(colnames(seu.merged)==colnames(sce.merged))
	reducedDim(sce.merged,"pca") <- Embeddings(seu.merged,"pca")
	reducedDim(sce.merged,"umap") <- Embeddings(seu.merged,"umap")
	reducedDim(sce.merged,"harmony") <- Embeddings(seu.merged,"harmony")
	reducedDim(sce.merged,"harmony.umap") <- Embeddings(seu.merged,"harmony.umap")
	#reducedDim(sce.merged,"tsne") <- Embeddings(seu.merged,"tsne")
	#reducedDim(sce.merged,"harmony.tsne") <- Embeddings(seu.merged,"harmony.tsne")
	reducedDim(sce.merged,"tsne.Rtsne") <- Embeddings(seu.merged,"tsne.Rtsne")
	reducedDim(sce.merged,"harmony.tsne.Rtsne") <- Embeddings(seu.merged,"harmony.tsne.Rtsne")

	for(a.res in resolution.vec){
		sce.merged[[sprintf("RNA_pca_snn_res.%s",a.res)]] <- seu.merged@meta.data[,sprintf("RNA_pca_snn_res.%s",a.res)]
		sce.merged[[sprintf("RNA_snn_res.%s",a.res)]] <- seu.merged@meta.data[,sprintf("RNA_snn_res.%s",a.res)]
	}

	loginfo("save result ...")
	saveRDS(seu.merged,file=sprintf("%s.seu.merged.rds",out.prefix))
	saveRDS(sce.merged,file=sprintf("%s.sce.merged.rds",out.prefix))
	
	return(list("seu.merged"=seu.merged,"sce.merged"=sce.merged))
}

run.raw <- function(exp.list.table,out.prefix,gene.exclude.file,ncores=12,npc=15,TH.gene.occ=1,
							   contamination.vec=NULL,cor.MALAT1=F,cor.DIG=T)
{
	RhpcBLASctl::omp_set_num_threads(1)
	doParallel::registerDoParallel(cores = ncores)

	env.misc <- loadToEnv(gene.exclude.file)
	cat(sprintf("black gene list:\n"))
	env.misc$all.gene.ignore.df %>% head %>% print

	cat(sprintf("contamination.vec:\n"))
	print(str(contamination.vec))

	loginfo("read data ...")
	seu.list <- llply(seq_len(nrow(exp.list.table)),function(i){
						  seu <- readRDS(exp.list.table$seufile[i])
						  data.id <- exp.list.table$data.id[i]
						  seu$dataset <- data.id
						  seu 
	},.parallel=T)
	names(seu.list) <- exp.list.table$data.id

	sce.list <- llply(seq_len(nrow(exp.list.table)),function(i){
	  sce <- readRDS(exp.list.table$efile[i])
	  sce$dataset <- exp.list.table$data.id[i]
	  if(!("norm_exprs" %in% assayNames(sce)))
	  {
		  assay(sce,"norm_exprs") <- assay(sce,"log2TPM")
	  }
	  if(exp.list.table$data.id[i] %in% c("HCC.YaoHe10X","HCC.YaoHeSS2") ){
		  f.gene.tmp <- which(rowData(sce)[,"display.name"]=="MALAT1-ENSG00000251562")
		  rowData(sce)[f.gene.tmp,"display.name"] <- "MALAT1"
	  }
	  sce
	  ##assay.name <- "norm_exprs"
	},.parallel=T)
	names(sce.list) <- exp.list.table$data.id

	#### get genes used for PCA. input from  sce.list and gene.de.list
	{

		gene.occ <- table(unlist(sapply(sce.list,function(x){ unname(rowData(x)[,"display.name"]) })))
		gene.occ <- gene.occ/length(sce.list)
		gene.occ <- sort(gene.occ,decreasing=T)
		gene.common <- names(gene.occ)[gene.occ >= TH.gene.occ]
		gene.common.all <- names(gene.occ)[gene.occ >= 0.85]
		loginfo(sprintf("total %d common genes obtain ",length(gene.common)))

		l_ply(names(sce.list),function(x){
				  ##cat(sprintf("gene %s in dataset %s: %s\n","CX3CR1",x,"CX3CR1" %in% rowData(sce.list[[x]])[,"display.name"]))
				  #cat(sprintf("gene %s in dataset %s: %s\n","MALAT1",x,"MALAT1" %in% rowData(sce.list[[x]])[,"display.name"]))
				  #cat(sprintf("gene %s in dataset %s: %s\n","IL17A",x,"IL17A" %in% rowData(sce.list[[x]])[,"display.name"]))
				  cat(sprintf("gene %s in dataset %s: %s\n","CCR8",x,"CCR8" %in% rowData(sce.list[[x]])[,"display.name"]))
		},.parallel=T)

		gene.de.list <- list()
		for(i in seq_len(nrow(exp.list.table))){
			id.d <- exp.list.table$data.id[i]
			ifile <- exp.list.table$ifile[i]
			de.out <- readRDS(ifile)
			gene.de.list[[id.d]] <- de.out$all
			gene.de.list[[id.d]]$geneID <- gene.de.list[[id.d]]$geneSymbol
		}
		names(gene.de.list) <- exp.list.table$data.id

		gene.rank.tb <- as.data.table(ldply(names(gene.de.list),function(x){
												ret.tb <- unique(gene.de.list[[x]][,c("geneID","F.rank")])
												ret.tb$dataset.id <- x
												ret.tb <- ret.tb[geneID %in% gene.common,]
												return(ret.tb)
														 }))
		gene.rank.tb <- dcast(gene.rank.tb,geneID~dataset.id,value.var="F.rank",fill=1)
		gene.rank.tb$median.F.rank <- rowMedians(as.matrix(gene.rank.tb[,-c("geneID"),with=F]))
		gene.rank.tb <- gene.rank.tb[order(median.F.rank),]
		##gene.rank.tb <- gene.rank.tb[geneID %in% gene.common,]
		#rowData(sce.pb)$median.F.rank <- gene.rank.tb[["median.F.rank"]][match(rownames(sce.pb),gene.rank.tb$geneID)]

		#### ASH1L MALAT1  CLDN15  ZNF117 NF1 ZNF546 HLA-A HLA-B HLA-C HLA-E PDIA3 IL2RG PTPRC B2M 
		f.gene.blackList <- (gene.rank.tb$geneID %in% env.misc$all.gene.ignore.df$geneSymbol) |
								grepl("^RP[LS]",gene.rank.tb$geneID,perl=T) |
								gene.rank.tb$geneID=="MALAT1"
		gene.de.common <- head(gene.rank.tb[!f.gene.blackList,][["geneID"]],n=1500)
		###gene.de.common <- head(gene.rank.tb$geneID,n=1500)

		cat(sprintf("RP gene in gene.de.common:\n"))
		print(gene.de.common[grepl("^RP[LS]",gene.de.common,perl=T)])
		
	}

	#### cellInfo
	{
		seu.res <- 2
		meta.tb <- ldply(seq_len(nrow(exp.list.table)),function(i){
								m.tb <- fread(exp.list.table$cfile[i])
								m.tb$dataset <- exp.list.table$data.id[i]
								a.res <- ifelse(exp.list.table$platform[i]=="SmartSeq2",
												sprintf("RNA_snn_res.%s",seu.res),
												sprintf("SCT_snn_res.%s",seu.res))
								m.tb$ClusterID <- sprintf("%s.%s",m.tb$dataset,m.tb[[a.res]])
								m.tb$dataset <- exp.list.table$data.id[i]
								o.tb <- m.tb[,c("patient","cellID","libraryID","cancerType","loc","batchV","TCR","dataset","ClusterID")]
								##o.tb <- m.tb[,.(NCell=.N),by=c("cancerType","dataset","ClusterID")]
								return(o.tb)
							   })
		rownames(meta.tb) <- sprintf("%s.%s",meta.tb$dataset,meta.tb$cellID)
		meta.tb$dataset.tech <- gsub("^.+?\\.","",meta.tb$dataset)
	}

	loginfo("merge data ...")
	dat.merged <- do.call(cbind,llply(seq_len(nrow(exp.list.table)),function(i){
							data.id <- exp.list.table$data.id[i]
							seu <- seu.list[[i]]
							sce <- sce.list[[i]]
							cat(sprintf("all(colnames(seu)==colnames(sce)) in %s: %s\n",data.id,
										all(colnames(seu)==colnames(sce))))
							dat.x <- assay(sce,"norm_exprs")
							if(!"seu.id" %in% colnames(rowData(sce))){
								rowData(sce)[,"seu.id"] <- gsub("_","-",rowData(sce)[,"display.name"])
							}
							rownames(dat.x) <- rowData(sce)[,"seu.id"]
							seu.x <- CreateSeuratObject(dat.x,project="panC", meta.data=seu[[]])
							adj.cov <- c("S.Score","G2M.Score","DIG.Score1")
							if("percent.mito" %in% colnames(seu.x[[]])){
								adj.cov <- c("percent.mito",adj.cov)
							}
							nBatch <- length(table(seu.x$batchV))
							if(nBatch>1){
								adj.cov <- c("batchV",adj.cov)
							}
							#seu.x <- ScaleData(seu.x,do.scale=T,features=gene.de.common,vars.to.regress = adj.cov,verbose=F)
							#dat.x <- GetAssayData(seu.x,"scale.data")
							dat.x <- GetAssayData(seu.x,"data")
							#colnames(dat.x) <- sprintf("%s.%s",data.id,colnames(dat.x))
							colnames(dat.x) <- unname(sprintf("%s.%s",data.id,colnames(dat.x)))
							#dat.x[1:4,1:5]
							return((dat.x[gene.de.common,]))
	},.parallel=ncores))

	#rownames(dat.merged) <- gsub("_","-",rownames(dat.merged))
	rownames(dat.merged) <- unname(gsub("_","-",rownames(dat.merged)))
	print(dim(dat.merged))
	print(dat.merged[1:4,1:5])
	print(class(dat.merged))
	print(str(dat.merged))
	print(str(colnames(dat.merged)))
	print(str(rownames(dat.merged)))

	seu.list <- NULL
	sce.list <- NULL
	gc()

	loginfo("begin Seurat pipeline ...")
	##### seurat object
	seu.merged <- CreateSeuratObject((dat.merged),project="panC",meta.data=meta.tb[colnames(dat.merged),])
	####seu.merged$dataset.tech <- gsub("^.+?\\.","",seu.merged$dataset)
	#seu.merged <- SetAssayData(seu.merged,"scale.data",dat.merged)
	seu.merged <- ScaleData(seu.merged,do.scale=T,verbose=F)

	seu.merged <- RunPCA(seu.merged,features = rownames(seu.merged), npcs = npc, verbose = FALSE)
	seu.merged <- RunUMAP(seu.merged,reduction="pca",dims=1:npc)
	seu.merged <- RunHarmony(seu.merged, c("dataset"))
	seu.merged <- RunUMAP(seu.merged,reduction = "harmony",reduction.name = "harmony.umap", dims = 1:npc)
	#seu.merged <- RunTSNE(seu.merged,reduction ="pca",tsne.method="FIt-SNE",dims=1:npc,reduction.name="tsne")
	#seu.merged <- RunTSNE(seu.merged,reduction ="harmony",tsne.method="FIt-SNE",dims=1:npc,reduction.name="harmony.tsne")

	loginfo("begin Rtsne ...")
	seu.merged <- RunTSNE(seu.merged,reduction ="pca",tsne.method="Rtsne",dims=1:npc,reduction.name="tsne.Rtsne")
	seu.merged <- RunTSNE(seu.merged,reduction ="harmony",tsne.method="Rtsne",dims=1:npc,reduction.name="harmony.tsne.Rtsne")
	loginfo("end Rtsne ...")

	resolution.vec <- seq(0.5,2.4,0.1)
	seu.merged <- FindNeighbors(seu.merged,reduction = "harmony", dims = 1:npc)
	#seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,algorithm=3)
	seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,algorithm=1)

	seu.merged <- FindNeighbors(seu.merged,reduction = "pca", dims = 1:npc,graph.name="RNA_pca_snn")
	#seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,algorithm=3, graph.name="RNA_pca_snn")
	seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,algorithm=1, graph.name="RNA_pca_snn")

	saveRDS(seu.merged,file=sprintf("%s.seu.merged.rds",out.prefix))

	sce.merged <- ssc.build((dat.merged))
	colData(sce.merged) <- DataFrame(meta.tb[colnames(sce.merged),])
	rowData(sce.merged)$gene.de.common <- T

	all(colnames(seu.merged)==colnames(sce.merged))
	reducedDim(sce.merged,"pca") <- Embeddings(seu.merged,"pca")
	reducedDim(sce.merged,"umap") <- Embeddings(seu.merged,"umap")
	reducedDim(sce.merged,"harmony") <- Embeddings(seu.merged,"harmony")
	reducedDim(sce.merged,"harmony.umap") <- Embeddings(seu.merged,"harmony.umap")
	#reducedDim(sce.merged,"tsne") <- Embeddings(seu.merged,"tsne")
	#reducedDim(sce.merged,"harmony.tsne") <- Embeddings(seu.merged,"harmony.tsne")
	reducedDim(sce.merged,"tsne.Rtsne") <- Embeddings(seu.merged,"tsne.Rtsne")
	reducedDim(sce.merged,"harmony.tsne.Rtsne") <- Embeddings(seu.merged,"harmony.tsne.Rtsne")

	for(a.res in resolution.vec){
		sce.merged[[sprintf("RNA_pca_snn_res.%s",a.res)]] <- seu.merged@meta.data[,sprintf("RNA_pca_snn_res.%s",a.res)]
		sce.merged[[sprintf("RNA_snn_res.%s",a.res)]] <- seu.merged@meta.data[,sprintf("RNA_snn_res.%s",a.res)]
	}

	loginfo("save result ...")
	saveRDS(seu.merged,file=sprintf("%s.seu.merged.rds",out.prefix))
	saveRDS(sce.merged,file=sprintf("%s.sce.merged.rds",out.prefix))
	
	return(list("seu.merged"=seu.merged,"sce.merged"=sce.merged))
}


plotNightingaleRose <- function(dat.plot.NightingaleRose,empty_bar=2,
					    y.pretty=c(0,20000,40000,60000),
					    y.lim.min=-20000,
					    exp.tick=0.5,
					    my.title=NULL,
					    y.colum="numOfCell",
					    x.colum="cancerType",
					    group.colum="pub",
					    pallette.group=structure(scales::viridis_pal(option = "viridis")(5)[c(2,3)],
								     names=c("published","thisStudy")))
{
    dat.plot.NightingaleRose$y.value <- dat.plot.NightingaleRose[[y.colum]]
    dat.plot.NightingaleRose$x.value <- dat.plot.NightingaleRose[[x.colum]]
    #empty_bar <- 1
    to_add <- matrix(NA, empty_bar, ncol(dat.plot.NightingaleRose))
    colnames(to_add) <- colnames(dat.plot.NightingaleRose)
    dat.plot.NightingaleRose <- rbind(dat.plot.NightingaleRose, to_add)
    dat.plot.NightingaleRose[,id:=as.integer(x.value)]
    dat.plot.NightingaleRose[is.na(x.value),
			     id:=seq(max(dat.plot.NightingaleRose$id,na.rm=T)+1,
				     max(dat.plot.NightingaleRose$id,na.rm=T)+empty_bar)]
    dat.plot.NightingaleRose[,id.factor:=as.factor(id)]

    ##label_data <- dat.plot.NightingaleRose[,.(numOfCell=sum(.SD$numOfCell)),by=c("cancerType","id")]
    label_data <- dat.plot.NightingaleRose[,.(y.value=sum(.SD$y.value)),by=c(x.colum,"id")]
    label_data$x.value <- label_data[[x.colum]]
    number_of_bar <- nrow(label_data)
    #label_data[,id:=number_of_bar-as.integer(cancerType)+1]
    #label_data[,id:=as.integer(cancerType)]
    angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar
    #angle <-  90 - 360 * (label_data$id+0.5) /number_of_bar
    label_data$hjust<-ifelse( angle < -90, 1, 0)
    label_data$angle<-ifelse(angle < -90, angle+180, angle)
    if(is.null(y.pretty)){
	y.pretty <- label_data[,pretty(y.value)]
    }
    if(is.null(y.lim.min)){
	y.lim.min <- -y.pretty[2]
    }
    nticks <- length(y.pretty)

    p <- ggbarplot(dat.plot.NightingaleRose,x="id.factor",y="y.value",group=group.colum,
			   fill=group.colum,
			   color=NA,width=0.8,
			   title=my.title,
			   #ylab="number of cells",xlab="cancer types",
			   legend="right") + 
	#facet_wrap(~tech.cate,scales="free_x") +
	#coord_flip() +
	scale_fill_manual(values=pallette.group)+
	### axis
	geom_segment(data=data.table(x=rep(max(dat.plot.NightingaleRose$id)+0.5,1),
				     y.start=0,
				     y.end=y.pretty[nticks]),
		     aes(x = x, y = y.start, xend = x, yend = y.end),
		     colour = "grey", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
	### ticks
	geom_segment(data=data.table(x.start=rep(max(dat.plot.NightingaleRose$id)+0.5,nticks),
				     x.end=rep(max(dat.plot.NightingaleRose$id)+0.5,
					       nticks)-seq(0.3,0.2,-0.1/(nticks-1))*exp.tick,
				     y.start=y.pretty,
				     y.end=y.pretty),
		     aes(x = x.start, y = y.start, xend = x.end, yend = y.end),
		     colour = "grey", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
	### axis.y
	annotate("text", x = rep(max(dat.plot.NightingaleRose$id)+0.5,nticks),
		 y = y.pretty,
		 label = y.pretty,
		 color="grey", size=5 , angle=0, hjust=1.2,vjust=0.5) +
	ylim(y.lim.min,y.pretty[nticks]) +
	coord_polar(start=0,direction=1) +
	### label
	geom_text(data=label_data, aes(x=id, y=y.value+y.pretty[2]/10,
				       label=x.value,
				       angle=angle,
				       hjust=hjust),
		  color="black", alpha=0.6, size=4.2,
		  #angle= label_data$angle,
		  inherit.aes = FALSE ) +
	theme(
	      plot.title = element_text(hjust = 0.5),
	      axis.ticks=element_blank(),axis.line=element_blank(),
	      panel.grid = element_blank(),plot.margin = unit(rep(0,4), "cm"),
	      axis.text = element_blank(), axis.title = element_blank())
    return(p)
}

#############################

cal.CSI <- function(incMat,lowerBoundary=0.05,do.parallel=T,ret.type="matrix")
{
	###reg.cor <- sscVis:::cor.BLAS(reg.inc.mtx)
	reg.cor <- sscVis:::cor.BLAS(incMat)
	cor.pairs <- as.data.table(t(combn(colnames(reg.cor),2)))
	cor.pairs$pcc <- reg.cor[lower.tri(reg.cor)]
	cor.pairs <- rbind(cor.pairs,data.table(V1=colnames(reg.cor),
											V2=colnames(reg.cor),
											pcc=1))
	csi.tb <- as.data.table(ldply(seq_len(nrow(cor.pairs)),function(i){
									  g1 <- cor.pairs$V1[i]
									  g2 <- cor.pairs$V2[i]
									  th.pcc <- cor.pairs$pcc[i] - lowerBoundary
									  go <- setdiff(colnames(reg.cor),c(g1,g2))
									  f.go <- rowSums(reg.cor[go,c(g1,g2)] < th.pcc)==2
									  index <- sum(f.go)/ncol(reg.cor)
									  if(g1!=g2){
										  return(data.table(V1=c(g1,g2),
													 V2=c(g2,g1),
													 index=rep(index,2)))
									  }else{
										  return(data.table(V1=g1,V2=g2,index=index))
									  }
							},.parallel=do.parallel))
	csi.tb <- csi.tb[order(V1,V2),]
	if(ret.type=="data.table"){
		return(csi.tb)
	}else if(ret.type=="matrix"){
		csi.tmp <- dcast(csi.tb,V2~V1,value.var="index")
		csi.mtx <- as.matrix(csi.tmp[,-1])
		rownames(csi.mtx) <- csi.tmp[[1]]
		return(csi.mtx)
	}else{
		return(NULL)
	}
}

#### custom theme
{
	library("grid")
	element_custom <- function(...) {
	  structure(list(...), class = c("element_custom", "element_blank"))
	}
	element_grob.element_custom <- function(element, label, x, y, ...)  {
	  tg <- textGrob(label, x=x,
					 y = unit(0.0,"npc"),
					 gp=gpar(col=element$colour,fontsize=element$size),
					 rot = 90,just="left")
	  padding <- unit(0.0,"line")
	  rg <- rectGrob(x=x,
					 y = unit(0.5,"npc"),
					 height=grobHeight(tg)+unit(0.5,"char"), width=unit(0.68,"line")+padding,
					 gp=gpar(fill = element$fill, col=NA, alpha=0.75))
	  gTree(children=gList(rg, tg), height=grobHeight(tg) + padding, cl="custom_axis")
	}
	heightDetails.custom_axis <- function(y) y$height + unit(2,"mm") # fudge

	#qplot(1:10,1:10) + theme(axis.text.x = element_custom(fill=1:3,size=200))
	#ggsave(sprintf("%s.test.png",out.prefix))
	
}


