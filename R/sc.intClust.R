####suppressPackageStartupMessages(library("sscClust"))
####suppressPackageStartupMessages(library("Seurat"))
####suppressPackageStartupMessages(library("harmony"))
####suppressPackageStartupMessages(library("tictoc"))
####suppressPackageStartupMessages(library("plyr"))
####suppressPackageStartupMessages(library("doParallel"))
####suppressPackageStartupMessages(library("Matrix"))
####suppressPackageStartupMessages(library("data.table"))
####suppressPackageStartupMessages(library("R.utils"))
####suppressPackageStartupMessages(library("gplots"))
####suppressPackageStartupMessages(library("ggpubr"))
####suppressPackageStartupMessages(library("ggplot2"))
####suppressPackageStartupMessages(library("ggrepel"))
####suppressPackageStartupMessages(library("cowplot"))
####suppressPackageStartupMessages(library("limma"))
####suppressPackageStartupMessages(library("lisi"))
####suppressPackageStartupMessages(library("RColorBrewer"))
####suppressPackageStartupMessages(library("ComplexHeatmap"))
####library("metacell")

####source("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/PanC.T/lib/plot.func.R")

#' highly variable genes finding (method using percetile rank values from F test)
#' @importFrom data.table data.table as.data.table `:=`
#' @importFrom plyr ldply
#' @importFrom stats median
#' @param gene.rank.tb data.table object
#' @param n.common integer; number of common top genes. (default: 1000)
#' @param n.specific integer; number of top genes with high dataset-specificities. (default: 1000)
#' @param th.rank double; threshold for rank. (default: 0.1)
#' @return a data.table object
#' @details For specific genes, the top n.specific genes with percetile rank values < th.rank are selected.
HVG.From.GeneRankTb <- function(gene.rank.tb,n.common=1000,n.specific=1000,th.rank=0.1)
{
    gene.rank.tb <- gene.rank.tb[order(median.F.rank),]
    gene.rank.tb[,rank:=seq_len(nrow(gene.rank.tb))]
    gene.rank.mat <- as.matrix(gene.rank.tb[,-c("geneID","median.F.rank","rank")])
    rownames(gene.rank.mat) <- gene.rank.tb$geneID
    #print(gene.rank.mat[1:4,1:min(5,ncol(gene.rank.mat))])
    gene.rank.comp.mat <- 1 - gene.rank.mat
    
    ###
    f.common <- gene.rank.tb[rank <= n.common,][["geneID"]]
    gene.rank.tail.mat <- gene.rank.mat[!(rownames(gene.rank.mat) %in% f.common),]
    f.specific <- NULL
    if(nrow(gene.rank.tail.mat) > 0){
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
        #print(f.specific.tb[geneID %in% c("IL17A","IL17F","IL23R","RORC","IL26","IL22","IL21","CD4","CD8A","CD8B"),])
        f.specific <- f.specific.tb[rank <= n.specific,][["geneID"]]
    }

    ####
#    sscVis:::plotMatrix.simple(gene.rank.tail.mat[f.specific,],
#			       out.prefix=sprintf("%s.gene.rank.datasetSpecific",out.prefix),
#			       show.dendrogram=T,clust.row=T,clust.column=T,
#			       pdf.width = 8,pdf.height = 10,
#			       exp.name="1-percRank")
    ret.tb <- data.table(geneID=f.common,category="common")
    if(length(f.specific) > 0){
        ret.tb <- rbind(ret.tb,data.table(geneID=f.specific,category="specific"))
    }
    return(ret.tb)
    ##return(c(f.common,f.specific))
}


#' wrapper for running leiden clustering
#' @importFrom scran buildSNNGraph
#' @importFrom leiden leiden
#' @param dat.pca data.frame; row for sample, column for variable
#' @param SNN.k integer; parameter k of scran::buildSNNGraph. (default: 20)
#' @param myseed integer; seed for random number generation. (default: 123456)
#' @param ... ; other parameters passed to leiden::leiden
#' @return a object returned from leiden::leiden
#' @details implementation in Seurat (FindClusters) can only detect singletons :-(
#' @export
run.Leiden <- function(dat.pca,SNN.k=20,myseed=123456,...)
{
	####x <- Embeddings(seu.merged,"harmony")[,1:15]
	snn.gr <- scran::buildSNNGraph(t(dat.pca), k=SNN.k,d=NA)
	clust.res <- leiden::leiden(snn.gr,seed=myseed,...)
	return(clust.res)
}

#' merge data file file table
#' @importFrom SingleCellExperiment rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom Seurat CreateSeuratObject ScaleData RunPCA FindNeighbors FindClusters
#' @importFrom harmony RunHarmony
#' @importFrom sscVis loginfo ssc.build ssc.average.cell
#' @importFrom plyr llply
#' @param exp.list.table data.table; one line for a dataset
#' @param gene.de.common character vector; common genes across multiple datasets
#' @param seu.list list; list containing seu objects
#' @param sce.list list; list containing sce objects
#' @param res.hi integer; high resolution used for mini-clusters identification
#' @param method.clustering character; clustering method for mini-clusters identification. (default: "louvain")
#' @param ncores integer; number of CPU cores to use. (default: 12)
#' @param cor.var character vector; subset of c("S.Score","G2M.Score","DIG.Score1","ISG.Score1","score.MALAT1")
#' @param contamination.vec character vector; cells to be excluded. (default: NULL)
#' @param use.harmony logical; use harmony to correct for batch effect (batches are defined in column batchV). (default: FALSE)
#' @return a list containing 3 components: dat.avg, meta.extra.tb and clust
#' @details For each dataset, the function first identify mini-clusters, then calculate the average expressions of mini-clusters.
mergeDataFromFileTable <- function(exp.list.table,gene.de.common,seu.list,sce.list,
				   res.hi,method.clustering="louvain",
				   #cor.cellCycle=T,cor.MALAT1=F,cor.DIG=T,cor.ISG=F,
                   ncores=12,
				   cor.var=c("S.Score","G2M.Score","DIG.Score1"),
				   contamination.vec=NULL,use.harmony=F)
{
    ret.list <- llply(seq_len(nrow(exp.list.table)),function(i){
		    data.id <- exp.list.table$data.id[i]
		    ##### seu and sce --> seu.x
		    seu <- seu.list[[i]]
		    sce <- sce.list[[i]]
		    loginfo(sprintf("all(colnames(seu)==colnames(sce)) in %s: %s",data.id,
					    all(colnames(seu)==colnames(sce))))
		    dat.x <- assay(sce,"norm_exprs")
		    if(!"seu.id" %in% colnames(rowData(sce))){
			    rowData(sce)[,"seu.id"] <- gsub("_","-",rowData(sce)[,"display.name"])
		    }
		    rownames(dat.x) <- rowData(sce)[,"seu.id"]
		    #### pad zero
		    {
			gene.pad <- setdiff(gene.de.common,rownames(dat.x))
			loginfo(sprintf("For dataset %s, those genes will be imputed by 0: %s",
						data.id, 
						paste(gene.pad,collapse=", ")))
			if(length(gene.pad)>0){
				value.pad <-  0
				mtx.pad <- matrix(rep(value.pad,length(gene.pad)*ncol(dat.x)),nrow=length(gene.pad))
				rownames(mtx.pad) <- gene.pad
				dat.x <- rbind(dat.x,mtx.pad)
			}
		    }

            score.MALAT1 <- NULL
            if("MALAT1" %in% rowData(sce)$display.name){
                score.MALAT1 <- assay(sce,"norm_exprs")[rowData(sce)[,"display.name"]=="MALAT1",]
            }else{
                cor.var <- setdiff(cor.var,"score.MALAT1")
            }

		    #if(data.id %in% c("HCC.YaoHe10X","HCC.YaoHeSS2") ){
            #    gene.MALAT1.vec <- c("MALAT1-ENSG00000251562")
		    #}else{
            #    gene.MALAT1.vec <- c("MALAT1")
		    #}
		    ##seu <- AddModuleScore(seu, features=list("score.MALAT1"=gene.MALAT1.vec), name="score.MALAT1",
		    ##					  pool = NULL, nbin = 24, ctrl = 100)

		    sce <- NULL
		    if(!is.null(contamination.vec)){
                f.cell.cont <- sprintf("%s.%s",data.id,colnames(dat.x)) %in% contamination.vec
                loginfo(sprintf("Number of contaminated cells (%s): %d.",data.id,sum(f.cell.cont)))
                dat.x <- dat.x[,!f.cell.cont]
                score.MALAT1 <- score.MALAT1[colnames(dat.x)]
		    }
		    seu.x <- CreateSeuratObject(dat.x,project="panC", meta.data=seu[[]][colnames(dat.x),])
		    ### seu$score.MALAT11, two "1"
		    #score.MALAT1 <- seu$score.MALAT11
            if(!is.null(score.MALAT1)){
		        seu.x$score.MALAT1 <- score.MALAT1
            }
		    seu <- NULL
		    gc()
		    ##### regression
		    meta.extra.tb <- data.table(cellID.uniq=sprintf("%s.%s",data.id,colnames(seu.x)),
									    S.Score=seu.x$S.Score,
									    G2M.Score=seu.x$G2M.Score,
									    Phase=seu.x$Phase,
									    DIG.Score1=seu.x$DIG.Score1,
									    score.MALAT1=score.MALAT1,
									    ISG.Score1=seu.x$ISG.Score1,
									    percent.mito=NA,
									    stringsAsFactors=F)
		    adj.cov <- cor.var
	        ###### patch ######
	        if(!"percent.mito" %in% colnames(seu.x[[]])) {
                idx.pmt <- grep("^(percent.mito|percent.mt)$",colnames(seu.x[[]]),value=T,perl=T)
                if(length(idx.pmt) > 0){
                    seu.x$percent.mito <- seu.x[[]][,idx.pmt[1]]
                }
	        }
	        if(!"patient" %in% colnames(seu.x[[]])) {
                idx.patient <- grep("^(patient|Patient)$",colnames(seu.x[[]]),value=T,perl=T)
                if(length(idx.patient) > 0){
                    seu.x$patient <- seu.x[[]][,idx.patient[1]]
                }else{
                    seu.x$patient <- "PXX"
                }
	        }
	        if(!"batchV" %in% colnames(seu.x[[]])) { seu.x$batchV <- seu.x$patient  }
	        ###################

	        if("percent.mito" %in% colnames(seu.x[[]])){
                adj.cov <- c("percent.mito",adj.cov)
                meta.extra.tb$percent.mito <- seu.x$percent.mito
		    }

		    nBatch <- length(table(seu.x$batchV))
		    ###adj.cov <- c()
		    if(nBatch>1 && !use.harmony){
		        adj.cov <- c("batchV",adj.cov)
		    }
		    loginfo(sprintf("ScaleData on dataset %s, with adj.cov : %s ",
				    data.id,
				    paste(adj.cov,collapse=", ")))
		    seu.x <- ScaleData(seu.x,do.scale=T,features=gene.de.common,
						       vars.to.regress = adj.cov,verbose=F)
		    
		    seu.x <- RunPCA(seu.x,features=rownames(seu.x),npcs= 15,verbose = FALSE)

            loginfo(sprintf("use harmony to correct for batch effect in each dataset: %s",use.harmony))
            rd.use <- "pca"
            if(nBatch>1 && use.harmony){
                seu.x <- RunHarmony(seu.x, c("batchV"),verbose=F)
                rd.use <- "harmony"
            }

		    loginfo(sprintf("find mini-clusters on dataset %s ",data.id))
		    ###### Seurat high resolution ####
		    if(method.clustering=="leiden"){
                clust.x <- run.Leiden(Embeddings(seu.x,rd.use),SNN.k=20,myseed=123456,resolution_parameter=res.hi)
                seu.x@meta.data[[sprintf("RNA_snn_res.%d",res.hi)]] <- clust.x
		    }else if(method.clustering=="louvain"){
                seu.x <- FindNeighbors(seu.x, reduction = rd.use,
                                       #k=if(ncol(seu.x)<500) 5 else 10,
                                       #k=if(ncol(seu.x)<500) 5 else 5,
                                       k=if(ncol(seu.x)<500) 10 else 10,
                                       #k=if(ncol(seu.x)<500) 10 else 20,
                                       dims = 1:15,nn.eps=0,force.recalc=T)
                ####res.hi <- if(ncol(seu.x)<500){ 25 } else res.hi
                res.hi <- if(ncol(seu.x) < 100) { 5 }else if(ncol(seu.x)<500){ 25 } else res.hi
                seu.x <- FindClusters(seu.x,resolution =res.hi, algorithm=1,verbose=F)
		    }

		    ###### single cell #####
		    loginfo(sprintf("calculate mini-cluster level expression (dataset %s)",data.id))
		    sce.x <- ssc.build(GetAssayData(seu.x,"scale.data"))
		    loginfo(sprintf("all(colnames(sce.x)==colnames(seu.x))? %s",
				    all(colnames(sce.x)==colnames(seu.x))))
		    sce.x$ClusterID <- sprintf("%s.C%04d",data.id,
					       as.integer(as.character(seu.x@meta.data[[sprintf("RNA_snn_res.%d",res.hi)]])))
		    dat.avg <- ssc.average.cell(sce.x,column="ClusterID",ret.type="data.mtx")
		    colnames(dat.avg) <- sprintf("%s", (colnames(dat.avg)))
		    meta.extra.tb$miniCluster <- sce.x$ClusterID
		    ### todo: save the new PCA result
		    return(list("dat.avg"=dat.avg[gene.de.common,],
				"meta.extra.tb"=meta.extra.tb,
				"clust"=structure(sce.x$ClusterID, names=sprintf("%s.%s",data.id,colnames(sce.x)))))
	    },.parallel=ncores)
    return(ret.list)
}

#' collapse gene by cell expression data to gene by mini-cluster expression data, with input stored in SingleCellExperiment object
#' @importFrom SingleCellExperiment colData rowData
#' @importFrom sscVis loginfo ssc.scale ssc.average.cell
#' @importFrom plyr llply
#' @param exp.list.table data.table; one line for a dataset
#' @param gene.common character vector; common genes across multiple datasets
#' @param sce.list list; list containing sce objects
#' @param group.vec character; cell id to mini-cluster mapping vector
#' @param ncores integer; number of CPU cores to use
#' @param contamination.vec character vector; cells to be excluded. (default: NULL)
#' @return a list, each component of which is expression matrix of gene.common by min-clusters
#' @details For each dataset, the function first identify mini-clusters, then calculate the average expressions of mini-clusters.
mergeSCEDataFromFileTable <- function(exp.list.table,gene.common,sce.list,group.vec,ncores=6,contamination.vec=NULL)
{
    ret.list <- llply(seq_len(nrow(exp.list.table)),function(i){
                    data.id <- exp.list.table$data.id[i]
                    ##### seu and sce --> seu.x
                    sce <- sce.list[[i]]
                    colnames(sce) <- sprintf("%s.%s",data.id,colnames(sce))
                    if(!is.null(contamination.vec)){
                        f.cell.cont <- colnames(sce) %in% contamination.vec
                        loginfo(sprintf("Number of contaminated cells (%s, mergeSCEDataFromFileTable): %d",
                                    data.id,sum(f.cell.cont)))
                        sce <- sce[,!f.cell.cont]
                    }
                    sce$miniCluster <- group.vec[colnames(sce)]
                    ########### patch
                    if(!"patient" %in% colnames(colData(sce))) {
                        idx.patient <- grep("^(patient|Patient)$",colnames(colData(sce)),value=T,perl=T)
                        if(length(idx.patient) > 0){
                            sce$patient <- colData(sce)[,idx.patient[1]]
                        }else{
                            sce$patient <- "PXX"
                        }
                    }
                    if(!"batchV" %in% colnames(colData(sce))) { sce$batchV <- sce$patient  }
                    #########################
                    sce <- ssc.scale(sce,gene.symbol=gene.common,assay.name="norm_exprs",
                                     adjB="batchV",do.scale=T)
                    dat.avg <- ssc.average.cell(sce,"norm_exprs.scale",
                                    column="miniCluster",ret.type="data.mtx")
                    rownames(dat.avg) <- rowData(sce)$display.name
                    #### pad zero
                    {
                        gene.pad <- setdiff(gene.common,rownames(dat.avg))
                        loginfo(sprintf("For dataset %s, %d genes will be imputed by 0",
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

#' run the integration pipeline
#' @importFrom SingleCellExperiment reducedDim colData rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom Seurat CreateSeuratObject SetAssayData RunPCA RunUMAP RunTSNE FindNeighbors FindClusters  Embeddings
#' @importFrom sscVis loginfo ssc.build
#' @importFrom harmony RunHarmony
#' @importFrom data.table as.data.table dcast
#' @importFrom matrixStats rowMedians
#' @importFrom plyr llply ldply
#' @importFrom R.utils loadToEnv
#' @importFrom utils head str write.table
#' @importFrom ggplot2 ggplot ggsave theme element_text
#' @importFrom ggpubr ggboxplot
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @param exp.list.table data.table; one line for a dataset
#' @param out.prefix character; output prefix
#' @param gene.exclude.file character; file contains the genes to be excluded
#' @param nGene.common integer; number of common genes. (default: 1500)
#' @param nGene.specific integer; number of dataset specific genes. (default: 0)
#' @param ncores integer; number of CPU cores to use. (default: 12)
#' @param npc integer; number of principal components to use. (default: 15)
#' @param TH.gene.occ double; range from 0 to 1. genes present in >= TH.gene.occ datasets are used. (default: 1)
#' @param res.hi integer; high resolution used for mini-clusters identification. (default: 50)
#' @param method.clustering character; clustering method for mini-clusters identification. (default: "louvain")
#' @param cor.var character vector; subset of c("S.Score","G2M.Score","DIG.Score1","ISG.Score1","score.MALAT1")
#' @param use.harmony logical; use harmony to correct for batch effect (batches are defined in column batchV). (default: FALSE)
#' @param contamination.vec character vector; cells to be excluded. (default: NULL)
#' @return a list containing 3 components: sce.merged, seu.merged and meta.tb
#' @details For each dataset, the function first identify mini-clusters, then calculate the average expressions of mini-clusters. The gene by mini-cluster expression data will pass the pipeline: PCA, harmony, UMAP/Clustering.
#' @export
run.inte.metaClust <- function(exp.list.table,
			       out.prefix,
			       gene.exclude.file,
			       nGene.common=1500,nGene.specific=0,
			       ncores=12,npc=15,TH.gene.occ=1,
			       res.hi=50,method.clustering="louvain",
			       #cor.cellCycle=T,cor.MALAT1=F,cor.DIG=T,cor.ISG=F,
			       cor.var=c("S.Score","G2M.Score","DIG.Score1"),
                   use.harmony=F,
			       contamination.vec=NULL)
{
    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = ncores)

    env.misc <- loadToEnv(gene.exclude.file)
    loginfo(sprintf("A total of %d genes in the black list, the top ones:",nrow(env.misc$all.gene.ignore.df)))
    print(head(env.misc$all.gene.ignore.df))

    loginfo(sprintf("contamination.vec:"))
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
	gene.common.all <- gene.common
	###gene.common.all <- names(gene.occ)[gene.occ >= 0.85]
	loginfo(sprintf("total %d common genes obtained.",length(gene.common)))

	l_ply(names(sce.list),function(x){
	  ##cat(sprintf("gene %s in dataset %s: %s\n","CX3CR1",x,"CX3CR1" %in% rowData(sce.list[[x]])[,"display.name"]))
	  #cat(sprintf("gene %s in dataset %s: %s\n","MALAT1",x,"MALAT1" %in% rowData(sce.list[[x]])[,"display.name"]))
	  #cat(sprintf("gene %s in dataset %s: %s\n","IL17A",x,"IL17A" %in% rowData(sce.list[[x]])[,"display.name"]))
	  loginfo(sprintf("gene %s in dataset %s? %s.","CCR8",x,"CCR8" %in% rowData(sce.list[[x]])[,"display.name"]))
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
					      n.common=nGene.common,n.specific=nGene.specific,th.rank=0.1)
	write.table(gene.de.common.tmp.tb,
		    file=sprintf("%s.gene.de.common%d.specific%s.tb",out.prefix,nGene.common,nGene.specific),
		    row.names=F,sep="\t",quote=F)
	gene.de.common <- gene.de.common.tmp.tb$geneID

#	    gene.de.common.tmp.tb <- HVG.From.GeneRankTb(gene.rank.tb[!f.gene.blackList,],
#						  ##n.common=1500,n.specific=0,th.rank=0.1)
#						  ##n.common=750,n.specific=750,th.rank=0.1)
#						  n.common=1000,n.specific=1000,th.rank=0.1)
#	    write.table(gene.de.common.tmp.tb,
#			file=sprintf("%s.gene.de.common.1000.1000.tb",out.prefix),
#			row.names=F,sep="\t",quote=F)
#	    gene.de.common <- gene.de.common.tmp.tb$geneID

	cat(sprintf("RP gene in gene.de.common:\n"))
	print(gene.de.common[grepl("^RP[LS]",gene.de.common,perl=T)])
	    
    }

    #### cellInfo
    {
	seu.res <- 2
	meta.tb <- ldply(seq_len(nrow(exp.list.table)),function(i){
                m.tb <- as.data.table(colData(sce.list[[i]]))
                m.tb$dataset <- exp.list.table$data.id[i]
                if(sprintf("SCT_snn_res.%s",seu.res) %in% colnames(m.tb))
                {
                    a.res <- sprintf("SCT_snn_res.%s",seu.res)
                }else{
                    a.res <- sprintf("RNA_snn_res.%s",seu.res)
                }
                m.tb$ClusterID <- sprintf("%s.%s",m.tb$dataset,m.tb[[a.res]])
                m.tb$dataset <- exp.list.table$data.id[i]
                if(!"patient" %in% colnames(m.tb))
                {
                    idx.patient <- grep("^(patient|Patient)$",colnames(m.tb),value=T,perl=T)
                    if(length(idx.patient) > 0){
                        m.tb$patient <- m.tb[[idx.patient[1]]]
                    }else{
                        m.tb$patient <- "PXX"
                    }
                }
                if(!"cellID" %in% colnames(m.tb)) { m.tb$cellID <- colnames(sce.list[[i]]) }
                if(!"libraryID" %in% colnames(m.tb)) { m.tb$libraryID <- "LXX"  }
                if(!"cancerType" %in% colnames(m.tb)) {
                    if("DiseaseType" %in% colnames(m.tb)){
                        m.tb$cancerType <- m.tb$DiseaseType
                    }else{
                        m.tb$cancerType <- "UNKNOWN"
                    }
                }
                if(!"loc" %in% colnames(m.tb)) {
                    if("Tissue" %in% colnames(m.tb)){
                        m.tb$loc <- m.tb$Tissue
                    }else{
                        m.tb$loc <- "T"
                    }
                }
                if(!"batchV" %in% colnames(m.tb)) { m.tb$batchV <- m.tb$patient  }
                if(!"TCR" %in% colnames(m.tb)) { m.tb$TCR <- ""  }
                o.tb <- m.tb[,c("patient","cellID","libraryID","cancerType","loc",
                                "batchV","TCR","dataset","ClusterID")]
                        ##o.tb <- m.tb[,.(NCell=.N),by=c("cancerType","dataset","ClusterID")]
                return(o.tb)
			})
	rownames(meta.tb) <- sprintf("%s.%s",meta.tb$dataset,meta.tb$cellID)
	meta.tb$dataset.tech <- gsub("^.+?\\.","",meta.tb$dataset)
    #print(as.data.table(meta.tb)[,.N,by=c("dataset","patient")])
    #return()
    }

    loginfo("begin merge data ...")
    dat.merged.list <- mergeDataFromFileTable(exp.list.table,gene.de.common,seu.list,sce.list,res.hi,
					      cor.var=cor.var,
                          ncores=ncores,
					      contamination.vec=contamination.vec)
    loginfo("end merge data.")
    loginfo(sprintf("A total of %d datasets are combined while the number of input datasets is %d.",
                    length(dat.merged.list),
                    length(sce.list)))
    
    dat.merged.mtx <- do.call(cbind,llply(dat.merged.list,function(x){ x[["dat.avg"]] }))
    rownames(dat.merged.mtx) <- gsub("_","-",rownames(dat.merged.mtx))
    loginfo("dim(dat.merged.mtx):")
    print(dim(dat.merged.mtx))
    print(dat.merged.mtx[1:4,1:3])

    meta.extra.tb <- as.data.frame(ldply(dat.merged.list,function(x){ x[["meta.extra.tb"]]}))
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
    ggsave(sprintf("%s.miniCluster.size.png",out.prefix),width=10,height=5)

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
    seu.merged <- RunUMAP(seu.merged,reduction="pca",dims=1:npc,verbose=F)
    seu.merged <- RunHarmony(seu.merged, c("dataset"),verbose=F)
    seu.merged <- RunUMAP(seu.merged,reduction = "harmony",reduction.name = "harmony.umap", dims = 1:npc,verbose=F)
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
	    seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,verbose=F)
	    seu.merged <- FindNeighbors(seu.merged,reduction = "pca", dims = 1:npc,graph.name="RNA_pca_snn")
	    seu.merged <- FindClusters(seu.merged,resolution = resolution.vec,graph.name="RNA_pca_snn",verbose=F)
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
    loginfo(sprintf("prepare sce.merged containing z-score expression of genes in gene.common.all"))
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

#' convert the limma result to gene by meta-cluster data stored in an SingleCellExperiment object
#' @importFrom SingleCellExperiment rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors `metadata<-` metadata
#' @importFrom sscVis loginfo ssc.build ssc.plot.heatmap
#' @importFrom sscClust integrate.by.avg effectsize rank.de.gene
#' @importFrom data.table as.data.table data.table `:=`
#' @importFrom plyr llply l_ply
#' @importFrom stats pt qnorm
#' @importFrom utils head write.table
#' @param de.limma.tb data.table; one line for a dataset. columns "data.id", "platform", and "dfile" are required
#' @param out.prefix character; output prefix
#' @param ncores integer; number of CPU cores to use. (default: 8)
#' @param min.ncells integer; only meta-clusters with number of cells > min.ncells are used. (default: 30)
#' @param min.ncellsStudy integer; only datasets with number of cells > min.ncellsStudy are used. (default: 200)
#' @param gset.list list; list containing gene sets. (default: NULL)
#' @param direct.sig logical; if TRUE, genes exibit significance in all datasets will be assigned "sig" directly, irrespective of combined ES and combined adjusted p vlaue. (default: TRUE)
#' @param TH.gene.occ double; range from 0 to 1. genes present in >= TH.gene.occ datasets are used. (default: 1)
#' @param de.mode character; mode of differential expression analysis. (default: "multiAsTwo")
#' @param column.exp character; convert the column of limma result to assay data. (default: "meanScale")
#' @param gene.used character; only keep genes in gene.used. (default: NULL)
#' @param colSet list; mapping iterms in the names to colors in the values. (default: list())
#' @return a SingleCellExperiment object
#' @details convert the limma result to gene by meta-cluster data stored in an SingleCellExperiment object
#' @export
convertLimmaToSCE <- function(de.limma.tb,out.prefix,ncores=8,
				    min.ncells=30,min.ncellsStudy=200,
				    gset.list=NULL,
                    direct.sig=T,TH.gene.occ=1,
				    de.mode="multiAsTwo",column.exp="meanScale",
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
            gene.de.list[[id.d]] <- de.out
            gene.de.list[[id.d]]$all$geneID <- gene.de.list[[id.d]]$all$geneSymbol
        }

        if(is.null(gene.used)){
            gene.occ <- table(unlist(sapply(gene.de.list,function(x){ unique(x$all$geneSymbol) })))
            gene.occ <- gene.occ/length(gene.de.list)
            gene.occ <- sort(gene.occ,decreasing=T)
            gene.common <- names(gene.occ)[gene.occ >= TH.gene.occ]
            loginfo(sprintf("total %d common genes obtained, with TH.gene.occ==%s.",length(gene.common),TH.gene.occ))
            gene.used <- gene.common
        }

        for(i in seq_len(nrow(de.limma.tb)))
        {
            id.d <- de.limma.tb$data.id[i]
            #dfile <- de.limma.tb$dfile[i]
            #de.out <- readRDS(dfile)
            de.out <- gene.de.list[[id.d]]
            #### number of cells in each ClusterID
            #de.out$all[1:2,]
            ncells.control.vec <- NULL
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
                ncell.control.df <- unique(as.data.table(de.out$all[,c("cluster","length._control"),with=F]))
                ncell.df <- merge(ncell.df,ncell.control.df,by="cluster")
                colnames(ncell.df) <- c("ClusterID","ncells","ncells.control")
                ncells.vec <- structure(ncell.df$ncells,names=ncell.df$ClusterID)
                ncells.control.vec <- structure(ncell.df$ncells.control,names=ncell.df$ClusterID)
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

                                 ##print(all(rownames(ES)==names(zp)))
                                 out.tb <- data.table(geneSymbol=geneID.mapping.vec[rownames(ES)],cluster=group.id)
                                 out.tb <- cbind(out.tb,ES,zp)
                                 return(out.tb)
                            },.parallel=T))

            #### information from de.out$all
            ###gene.de.list[[id.d]] <- de.out$all
            ###gene.de.list[[id.d]]$geneID <- gene.de.list[[id.d]]$geneSymbol
            mnames <- c("logFC","t","P.Value","adj.P.Val","meanExp","meanScale","freq._case","freq._control",
                        "OR","OR.p.value","OR.adj.pvalue")

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
                                        #### padding values for these ars also 1: OR, OR.p.value, OR.adj.pvalue
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
            if(!is.null(ncells.control.vec)) { sce.obj$nCellsOther <- ncells.control.vec[sce.obj$ClusterID] }
            ###sce.obj <- sce.obj[,sce.obj$nCellsCluster >= min.ncells]
            sce.obj <- sce.obj[,sce.obj$nCellsCluster >= min.ncells & sce.obj$nCellsStudy >= min.ncellsStudy]

            ################
            sce.list[[id.d]] <- sce.obj
            cat(sprintf("%d, %s\n",i,de.limma.tb$data.id[i]))
            ####if(length(gene.common)==0)
            ####    ##gene.common <- rownames(sce.list[[id.d]])
            ####    gene.common <- rowData(sce.list[[id.d]])$display.name
            ####else
            ####    gene.common <- intersect(gene.common,rowData(sce.list[[id.d]])$display.name)
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
	    cat(sprintf("%s, %d genes\n",aid,length(unique(gene.de.list[[aid]]$all$geneID))))
    })

    .tmp.gene.de.list <- llply(names(gene.de.list),function(x){ gene.de.list[[x]]$all })
    names(.tmp.gene.de.list) <- names(gene.de.list)

    #####################################
    sce.pb <- integrate.by.avg(sce.list,sprintf("%s",out.prefix),
				     assay.name=column.exp, is.avg=T,n.downsample=NULL,
				     ncores=ncores,gene.de.list=.tmp.gene.de.list, de.thres=1000,
				     do.clustering=F,
				     avg.by="ClusterID")
    #sce.debug <<- sce.pb

    #### some meta info
    #mm <- regexec(".+\\.(CD[48]\\.c\\d+.+)$",colnames(sce.pb),perl=T)
    #sce.pb$meta.cluster <- sapply(regmatches(colnames(sce.pb),mm),"[",2)
    m <- regexec("^(.+?)\\.(.+?)\\.(.+)$",colnames(sce.pb),perl=T)
    mm <- regmatches(colnames(sce.pb),m)
    sce.pb$dataset.id <- NULL
    sce.pb$meta.cluster <- sapply(mm,"[",4)
    sce.pb$cancerType <- sapply(mm,"[",2)
    sce.pb$dataset <- sprintf("%s.%s",sapply(mm,"[",2),sapply(mm,"[",3))
    #sce.pb$cancerType <- sce.pb$dataset.id
    #sce.pb$dataset.old <- gsub("\\.CD[48]\\..+$","",sce.pb$ClusterID,perl=T)
    #sce.pb$dataset <- sce.pb$dataset.old
    dataset.mapping <- structure(de.limma.tb$platform,names=de.limma.tb$data.id)
    sce.pb$tech <- ifelse(dataset.mapping[sce.pb$dataset]=="SmartSeq2","SmartSeq2","Droplet")
    sce.pb$weight.tech <- ifelse(sce.pb$tech=="SmartSeq2",1,0.5)
    #### HCC and LIHC
    #sce.pb <- changeSomeNames(sce.pb)
    metadata(sce.pb)$ssc$colSet <- colSet

    ### in limma output, the significance is set by logFC and p-value
    ### here, reset significance by effect size and p-value
    sce.pb <- resetSig(sce.pb)
    

    ##### core signature genes
    {
        ################# sig genes table
        gene.desc.top <- rank.de.gene(sce.pb,group="meta.cluster",
                             weight.adj="weight.tech", group.2nd="cancerType",ncores=ncores)

        if(direct.sig){
            gene.desc.top[,sig:=(comb.padj < 0.01 & comb.ES>0.15) | freq.sig==1]
        }else{
            gene.desc.top[,sig:=comb.padj < 0.01 & comb.ES>0.15]
        }
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
#####        ntop <- 50
#####        sig.prevelance <- 0.4
#####        gene.top.mrank.tb <- gene.desc.top[median.rank<0.01 & freq.sig >= 0.15,][order(meta.cluster),]
#####        gene.top.mrank.tb <- gene.top.mrank.tb[,head(.SD,n=ntop),by=c("meta.cluster")]
#####        if(!is.null(out.prefix)){
#####            write.table(gene.top.mrank.tb,sprintf("%s.gene.top.mrank.n%d.prev%4.2f.mrank%4.2f.txt",
#####                                                  out.prefix,ntop,0.15,0.01),
#####                        row.names=F,sep="\t",quote=F)
#####        }
#####        gene.top.mrank.plot.tb <- gene.top.mrank.tb[freq.sig >=sig.prevelance, ][!duplicated(geneID),]
#####
#####        ##### heatmap
#####        ssc.plot.heatmap(sce.pb,out.prefix=sprintf("%s.gene.top%d.mrank.prev%4.2f.mrank%4.2f",
#####                               out.prefix,ntop,sig.prevelance,0.01),
#####                 columns="meta.cluster",columns.order="meta.cluster",colSet=colSet,
#####                 gene.desc=gene.top.mrank.plot.tb,
#####                 row.split=gene.top.mrank.plot.tb$Group,
#####                 column.split=sce.pb$meta.cluster,
#####                 #border = TRUE,
#####                 #row_title_rot = 0,
#####                 par.heatmap=list(row_title_rot=0,border=T,
#####                          row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
#####                          column_title_gp = gpar(fontsize = 0)),
#####                 mytitle=sprintf("Core Genes (top%d.prev%4.2f)",ntop,sig.prevelance),
#####                 pdf.width=20,pdf.height=10,do.scale=F,
#####                 palette.name="RdBu",
#####                 z.lo=-0.5,z.hi=0.5,z.step=0.25,
#####                 do.clustering.row=F,
#####                 do.clustering.col=T)
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
        if("geneSet.TF" %in% colnames(gene.desc.top)){
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

#' make heatmap for top signature genes
#' @importFrom sscVis plotMatrix.simple ssc.plot.heatmap
#' @importFrom data.table dcast
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid gpar unit
#' @param out.prefix character; output prefix
#' @param gene.desc.top data.table; gene info table
#' @param sce.pb object of SingleCellExperiment; gene by meta-cluster values stored in this object
#' @param gene.to.show.tb data.table; infomation of genes to be shown
#' @param value.var character; stored value to use. (default: "comb.ES")
#' @param colSet list; color set to use. (default: list())
#' @param ... ; parameters passed to sscVis::plotMatrix.simple
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
    #head(g.plot.mtx,n=3)
    #head(g.plot.mtx.bin,n=3)
    if(F){
        plotMatrix.simple(g.plot.mtx.bin,out.prefix=sprintf("%s.slim.bin",out.prefix),
                      col.ht=structure(rev(brewer.pal(5,name="RdBu")), names=1:5),
                      par.legend=list(labels=rev(bin.values), at=5:1),
                      row.split=gene.to.show.tb$Group,
                      pdf.width=if(ncol(g.plot.mtx)<20) 6.5 else 8,
                      par.heatmap=list(cex.column=0.8,border=T,
                               row_gap = unit(0, "mm"),
                               row_title_rot=0),
                      exp.name="comb.ES",...)
    }
    plotMatrix.simple(g.plot.mtx,out.prefix=sprintf("%s.slim",out.prefix),
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
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.Diether2018" ] <- "NSCLC.DietherLambrechts2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.QianqianSong2019" ] <- "NSCLC.QianqianSong2019"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.RapolasZilionis2019" ] <- "NSCLC.RapolasZilionis2019"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="Melanoma.LivnatJerby-Arnon2018" ] <- "MELA.LivnatJerby-Arnon2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="Melanoma.MosheSade-Feldman2018" ] <- "MELA.MosheSade-Feldman2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="Melanoma.HanjieLi2018" ] <- "MELA.HanjieLi2019"
	    ##obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="STAD.BoxiKang2019" ] <- "STAD.BoxiKang2020"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="STAD.BoxiKang2019" ] <- "STAD.BoxiKang2021"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="HCC.zhangLabSS2" ] <- "HCC.ChunhongZheng2017"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.zhangLabSS2" ] <- "NSCLC.XinyiGuo2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CRC.zhangLabSS2" ] <- "CRC.LeiZhang2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.Elham2018.10X" ] <- "BRCA.ElhamAzizi2018.10X"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.Elham2018.Indrop" ] <- "BRCA.ElhamAzizi2018.InDrop"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.Peter2018" ] <- "BRCA.PeterSavas2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.zhangLab5P" ] <- "BRCA.zhangLab5P"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CHOL.zhangLabSS2" ] <- "CHOL.thisStudy"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="NPC.XiliangWang2019" ] <- "NPC.YangLiu2021"
	    obj[[col.dataset]] <- gsub("zhangLab5P","thisStudy",obj[[col.dataset]])
    }
    
    return(obj)
}

#' reset the significance
#' @importFrom SummarizedExperiment assay
#' @param obj object of \code{SingleCellExperiment} class
#' @return a SingleCellExperiment object
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

#' make a long type gene table, given a SingleCellExperiment object
#' @importFrom sscVis ssc.toLongTable collapseEffectSizeLong
#' @importFrom plyr llply
#' @param obj object of \code{SingleCellExperiment} class
#' @param out.prefix output prefix
#' @param th.adj.P double; threshold for adjusted p-value. (default: 0.01)
#' @param th.dprime double; threshold for dprime. (default: 0.15)
#' @export
make.geneTableLong <- function(obj,out.prefix,th.adj.P=0.01,th.dprime=0.15)
{
    out.list <- (llply(unique(sort(obj[["meta.cluster"]])),function(x){
        dat.long <- ssc.toLongTable(obj[,obj$meta.cluster==x],gene.id=NULL,
                        assay.name=c("meanExp","meanScale","logFC",
                             "zp","dprime","vardprime",
                             "P.Value","adj.P.Val","sig","freq._case"),
                        col.idx=colnames(colData(obj)))
        dat.collapse <- collapseEffectSizeLong(dat.long,mode.collapse="comb",
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


#' render a knitr report
#' @importFrom rmarkdown render
#' @param template.file character; template file
#' @param out.file character; output file
#' @param par.list list; parameter list for rmarkdown::render
#' @return a SingleCellExperiment object
#' @export
render_KnitReport <- function(template.file,out.file, par.list=list() ) {
  render(template.file, params = par.list, output_file = out.file)
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
						    label = TRUE, group.by=cate.res) +
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
							  colSet=list(),size=0.1,label=2, base_aspect_ratio = 1.2)
    }
    for(i in seq_len(length(plot.resolution.list)/4))
    {
	pp <- plot_grid(plotlist=plot.resolution.list[((i-1)*4+1):(i*4)],
					ncol = 2,align = "hv")
	save_plot(sprintf("%s.%s.res.sceStyle.%d.png",out.prefix,rd,i),pp, ncol = 2,
			  base_aspect_ratio=0.9,base_height=5.5)
    }

    makeGeneOnTSNEPlot(sce,rd,out.prefix,geneOnUmap.list=g.geneOnUmap.list)

    #### density
    ssc.plot.tsne(sce,plotDensity=T,reduced.name=sprintf("%s",rd),
		      out.prefix=sprintf("%s.%s",out.prefix,rd))

}

makeGeneOnTSNEPlot <- function(sce,rd,out.prefix,
			       geneOnUmap.list=g.geneOnUmap.list,
			       plot.ncol=NULL,plot.nrow=NULL,plot.type="png",
			       plot.width=NULL,plot.height=NULL,do.parallel=T,...)
{
    if(!is.null(out.prefix)){ dir.create(dirname(out.prefix),F,T) }
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
            if(!is.null(out.prefix)){
            ggsave(sprintf("%s.%s.marker.%s.%s", out.prefix,rd,
                       names(geneOnUmap.list)[i],plot.type),
                   width=plot.width, height=plot.height)
            }else{
            print(p)
            }
        }
    },.parallel=do.parallel)
}

makeViolinPlot <- function(sce,out.prefix,assay.name="exprs",col.group="ClusterID",
					   geneOnUmap.list=g.geneOnUmap.list,plot.type="png",clamp=c(-2.5,5),
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

#############################







