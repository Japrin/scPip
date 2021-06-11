
#require("Seurat")
#require("tictoc")
#require("tictoc")
#require("dplyr")
#require("cowplot")
#require("plyr")
#require("ggplot2")
#suppressPackageStartupMessages(library("harmony"))

#' Wraper for highly variable genes finding
#' @importFrom Seurat FindVariableFeatures `VariableFeatures<-`
#' @importFrom tibble rownames_to_column
#' @importFrom utils head str
#' @importFrom magrittr `%>%`
#' @importFrom dplyr arrange
#' @param seu object of \code{Seurat}
#' @param gene.exclude.df data.frame; gene blak list. Required column: seu.id.
#' @param n.top integer; number of top genes. (default: 1500)
#' @param measurement character; "counts", "TPM" or "cpm". (default: "counts")
#' @return a Seurat object
#' @export
run.HVG <- function(seu,gene.exclude.df,n.top=1500,measurement="counts")
{
    seu <- FindVariableFeatures(object = seu)
    if(measurement=="TPM"){
	#### TPM
	hvg.gene.info <- seu@assays$RNA@meta.features %>% 
		tibble::rownames_to_column(var="geneSymbol") %>%
		arrange(-vst.variance)
    }else{
	#### counts, cpm
	hvg.gene.info <- seu@assays$RNA@meta.features %>% 
		tibble::rownames_to_column(var="geneSymbol") %>%
		arrange((-vst.variance.standardized))
    }
    f.hvg <- !(hvg.gene.info$geneSymbol %in% gene.exclude.df[["seu.id"]]) &
	!(grepl("^RP[LS]",hvg.gene.info$geneSymbol,perl=T)) &
	!(hvg.gene.info$geneSymbol %in% c("MALAT1","MALAT1-ENSG00000251562",
					  "MALAT1-ENSG00000279576","MALAT1-ENSG00000278217"))
    hvg.gene.info.flt <- hvg.gene.info[f.hvg,]
    if(measurement=="TPM"){
        hvg.gene.info.flt$rank.perc <- rank(-hvg.gene.info.flt$vst.variance)/nrow(hvg.gene.info.flt)
    }else {
        hvg.gene.info.flt$rank.perc <- rank(-hvg.gene.info.flt$vst.variance.standardized)/nrow(hvg.gene.info.flt)
    }
    VariableFeatures(seu) <- head(hvg.gene.info.flt,n=n.top)$geneSymbol
    print(str(seu@assays$RNA@var.features))
    print(head(hvg.gene.info.flt))
    seu@assays$RNA@misc$meta.features.flt <- hvg.gene.info.flt
    return(seu)
}

#' Wraper for running Seurat3 pipeline
#' @importFrom Seurat CreateSeuratObject SetAssayData GetAssayData CellCycleScoring AddModuleScore ScaleData SCTransform RunPCA ProjectDim RunUMAP RunTSNE FindNeighbors FindClusters Embeddings DimPlot NoLegend
#' @importFrom SingleCellExperiment `reducedDim<-`
#' @importFrom SummarizedExperiment assayNames assay `assay<-` rowData `rowData<-` colData `colData<-`
#' @importFrom harmony RunHarmony
#' @importFrom data.table data.table fread
#' @importFrom S4Vectors DataFrame
#' @importFrom sscVis ssc.build loginfo ssc.plot.tsne ssc.plot.violin
#' @importFrom sscClust ssc.DEGene.limma
#' @importFrom ggplot2 ggsave
#' @importFrom cowplot plot_grid save_plot
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom utils head str data
#' @importFrom tictoc tic toc
#' @param seu object of \code{Seurat}
#' @param sce object of \code{SingleCellExperiment}
#' @param out.prefix character; output prefix
#' @param gene.exclude.df data.frame; gene blak list. Required column: seu.id.
#' @param n.top integer; number of top genes. (default: 1500)
#' @param measurement character; "counts", "TPM" or "cpm". (default: "counts")
#' @param platform character; "10X", "SmartSeq2", "InDrop" etc.. (default: "10X")
#' @param opt.res character; optimal resolution (default: "1")
#' @param use.sctransform logical; whether use scTransform method (default: FALSE)
#' @param aid character; an ID (default: "PRJ")
#' @param plot.rd character vector; reducedDimNames used for plots (default: c("umap"))
#' @param opt.npc integer; optimal number of principal componets to use (default: 15)
#' @param ncores integer; number of CPU cores to use (default: 16)
#' @param cor.var character vector; Subset of c("S.Score","G2M.Score","DIG.Score1","ISG.Score1") or NULL. (default: c("S.Score","G2M.Score","DIG.Score1")). If certain variables are corrected (!is.null(cor.var) || cor.var!="NULL"), the pipeline will also correct for batchV and percent.mito.
#' @param ncell.deg integer; number of cell to downsample. used in the differentially expressed gene analysis. (default: 1500)
#' @param do.deg logical; whether perform the differentially expressed gene analysis. (default: FALSE)
#' @param do.scale logical; whether scale the expression data. (default: FALSE)
#' @param use.harmony logical; whether use the harmony method. (default: FALSE)
#' @param gene.mapping.table data.table; used for gene ID conversion. (default: NULL)
#' @param res.addition character vector; additional resolution parameters. Internally, resolutions from 0.1 to 2.4 will be used. (default: NULL)
#' @param run.stage integer; running stage. (default: 100)
#' @return a list contain a Seurat object and a SingleCellExperiment object
#' @details run the Seurat3 pipeline
#' @export
run.Seurat3 <- function(seu,sce,out.prefix,gene.exclude.df,n.top=1500,
			measurement="counts",platform="10X",
			opt.res="1",use.sctransform=F,aid="PRJ",plot.rd=c("umap"),
			opt.npc=15,ncores=16,
            #do.adj=T,
            cor.var=c("S.Score","G2M.Score","DIG.Score1"),
            ncell.deg=1500,do.deg=F,do.scale=F,
            use.harmony=F,
			gene.mapping.table=NULL,res.addition=NULL,
			run.stage=100)
{
    if(use.sctransform && platform!="SmartSeq2"){
	    opt.res <- sprintf("SCT_snn_res.%s",opt.res)
    }else{
	    opt.res <- sprintf("RNA_snn_res.%s",opt.res)
    }
    cat(sprintf("opt.res: %s\n",opt.res))
    if(measurement=="TPM"){
        #### TPM
        assay.name <- "log2TPM"
    }else{
        #### counts, cpm
        assay.name <- "norm_exprs"
    }

    ###
    if(is.null(seu) && !is.null(sce)){
        if("counts" %in% assayNames(sce)){
        exp.mat <- assay(sce,"counts")
        rownames(exp.mat) <- rowData(sce)[["seu.id"]]
        seu <- CreateSeuratObject(exp.mat, min.cells = 0, min.features = 0,
                  meta.data = as.data.frame(colData(sce)), project = aid)
        exp.mat <- assay(sce,assay.name)
        rownames(exp.mat) <- rowData(sce)[["seu.id"]]
        seu <- SetAssayData(seu, assay = "RNA", slot = 'data',new.data = exp.mat)
        rm(exp.mat)
        }else{
        ### norm_exprs, i.e. log2(x+1), expected
        exp.mat <- assay(sce,assay.name)
        rownames(exp.mat) <- rowData(sce)[["seu.id"]]
        seu <- CreateSeuratObject(exp.mat, min.cells = 0, min.features = 0,
                  meta.data = as.data.frame(colData(sce)), project = aid)
        }
    }

    if(is.null(gene.mapping.table)){
        gene.mapping.table <- data.table(geneID=as.character(rownames(seu)),
                         seu.id=as.character(rownames(seu)),
                         display.name=as.character(rownames(seu)))
    }

    if(is.null(sce) && !is.null(seu)){
	exp.mat <- GetAssayData(seu, slot = "counts")
	exp.mat[1:4,1:5]
	rownames(exp.mat) <- gene.mapping.table$geneID[match(rownames(exp.mat),
							     gene.mapping.table$seu.id)]
	sce <- ssc.build(exp.mat,assay.name="counts")
	rowData(sce)$geneID <- rownames(sce)
	rowData(sce)$display.name <- gene.mapping.table[match(rownames(sce),
							      gene.mapping.table$geneID)
							][["display.name"]]
	names(rowData(sce)$display.name) <- rownames(sce)
	rowData(sce)$seu.id <- gene.mapping.table[match(rownames(sce),
							gene.mapping.table$geneID)][["seu.id"]]
	print(head(rowData(sce)))

	print(all(colnames(seu)==rownames(colData(sce))))
	print(all(rownames(seu)==rowData(sce)$seu.id))
	colData(sce) <- DataFrame(seu[[]])
	exp.mat <- GetAssayData(seu, slot = "data")
	rownames(exp.mat) <- gene.mapping.table[match(rownames(exp.mat),
						      gene.mapping.table$seu.id)][["geneID"]]
	assay(sce,"norm_exprs") <- exp.mat
	#metadata(sce)$ssc[["variable.gene"]][["var.seurat"]] <- VariableFeatures(seu)
	#reducedDim(sce,"seurat.pca") <- Embeddings(seu, reduction = "pca")
	#reducedDim(sce,"seurat.umap") <- Embeddings(seu, reduction = "umap")
    }

    if(!all(colnames(sce)==colnames(seu))){
	    warning("seu and sce are not consistent!")
    }

    ###### patch ######
    if(!"percent.mito" %in% colnames(seu[[]]))
    {
        idx.pmt <- grep("^(percent.mito|percent.mt)$",colnames(seu[[]]),value=T,perl=T)
        if(length(idx.pmt) > 0){
        seu$percent.mito <- seu[[]][,idx.pmt[1]]
        sce$percent.mito <- seu$percent.mito
        }
    }
    if(!"patient" %in% colnames(seu[[]]))
    {
        idx.patient <- grep("^(patient|Patient|PatientID)$",colnames(seu[[]]),value=T,perl=T)
        if(length(idx.patient) > 0){
        seu$patient <- seu[[]][,idx.patient[1]]
        }else{
        seu$patient <- "PXX"
        }
    }
    if(!"batchV" %in% colnames(seu[[]])) { seu$batchV <- seu$patient  }
    sce$patient <- seu$patient
    sce$batchV <- seu$batchV
    ###################

    nBatch <- length(table(sce$batchV))
    loginfo(sprintf("Total batchs: %d",nBatch))
    print(str(table(sce$batchV)))

    #####
    
    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = ncores)
    
    plot.all <- function(rd="umap",resolution.vec=seq(0.1,2.4,0.1))
    {
        if(!is.null(seu@meta.data[["loc"]])){
        p <- DimPlot(seu,reduction=rd, group.by = "loc")
        ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"loc"),width=5,height=4)
        }

        if(!is.null(seu@meta.data[["stype"]])){
            p <- DimPlot(seu,reduction=rd, group.by = "stype")
            ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"stype"),width=5,height=4)
        }
        
        if(!is.null(seu@meta.data[["cancerType"]])){
        p <- DimPlot(seu,reduction=rd, group.by = "cancerType")
        ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"cancerType"),width=5,height=4)
        }

        if(!is.null(seu@meta.data[["libraryID"]])){
            if(length(unique(seu$libraryID))>20){
                p <- DimPlot(seu,reduction=rd,label=F,label.size=2,
                 group.by = "libraryID") + NoLegend()
                ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"libraryID"),width=4.5,height=4)
            }else{
                p <- DimPlot(seu,reduction=rd,label=T,label.size=2, group.by = "libraryID")
                ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"libraryID"),width=6.5,height=4)
            }
        }

        if(!is.null(seu@meta.data[["patient"]])){
            if(length(unique(seu$patient))>20){
                p <- DimPlot(seu,reduction=rd, group.by = "patient") + NoLegend()
                ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"patient"),width=4.5,height=4)
            }else{
                p <- DimPlot(seu,reduction=rd, group.by = "patient")
                ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"patient"),width=6.5,height=4)
            }
        }

        plot.resolution.list <- list()
        for(t.res in resolution.vec){
            if(use.sctransform && platform!="SmartSeq2"){
                cate.res <- sprintf("SCT_snn_res.%s",t.res)
            }else{
                cate.res <- sprintf("RNA_snn_res.%s",t.res)
            }
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
            save_plot(sprintf("%s.%s.res.%d.png",out.prefix,rd,i),pp,
                  ncol = 2, base_aspect_ratio=0.55)
        }

        plot.resolution.list <- list()
        for(t.res in resolution.vec){
            if(use.sctransform && platform!="SmartSeq2"){
                cate.res <- sprintf("SCT_snn_res.%s",t.res)
            }else{
                cate.res <- sprintf("RNA_snn_res.%s",t.res)
            }
            plot.resolution.list[[cate.res]] <- ssc.plot.tsne(sce,columns = cate.res,
                        reduced.name = sprintf("seurat.%s",rd),
                        colSet=list(),size=0.1,label=3,
                        par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                        base_aspect_ratio = 1.2)
        }

        #plot.resolution.list.debug <<- plot.resolution.list
        #out.prefix.debug <<- out.prefix
        for(i in seq_len(length(plot.resolution.list)/4))
        {
            pp <- plot_grid(plotlist=plot.resolution.list[((i-1)*4+1):(i*4)],
                    ncol = 2,align = "hv")
            save_plot(sprintf("%s.%s.res.sceStyle.%d.png",out.prefix,rd,i),pp,
                  ncol = 2, base_aspect_ratio=0.9,base_height=5.5)
        }

        ## gene on umap
        loginfo(sprintf("bein plotting geneOnUmap ..."))
        l_ply(seq_along(g.geneOnUmap.list),function(i){
            gene.tmp <- intersect(g.geneOnUmap.list[[i]],rowData(sce)$display.name)
            loginfo(sprintf("(begin) geneSet %s",names(g.geneOnUmap.list)[i]))
            if(length(gene.tmp)>0){
                p <- ssc.plot.tsne(sce,assay.name=assay.name,adjB=if(nBatch>1) "batchV" else NULL,
                       gene=gene.tmp,
                       par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                      reduced.name=sprintf("seurat.%s",rd))
                ggsave(sprintf("%s.seurat.%s.marker.%s.png",
                   out.prefix,rd,names(g.geneOnUmap.list)[i]),
                   width=10,
                   height=if(length(gene.tmp)>9) 11 else if(length(gene.tmp)>6) 8 else if(length(gene.tmp)>3) 5.4 else 2.7)
            }
            loginfo(sprintf("(end) geneSet %s",names(g.geneOnUmap.list)[i]))
        },.parallel=T)
        loginfo(sprintf("end plotting geneOnUmap."))

        #### density
        ssc.plot.tsne(sce,plotDensity=T,reduced.name=sprintf("seurat.%s",rd),
                  out.prefix=sprintf("%s.seurat.%s",out.prefix,rd))

        colSet <- list()

        ssc.plot.tsne(sce, columns = "percent.mito", 
                      reduced.name = sprintf("seurat.%s",rd),
                      colSet=colSet,size=0.03,
                      par.geneOnTSNE = list(scales = "free",pt.order = "random"),
                      vector.friendly=T,
                      out.prefix = sprintf("%s.seurat.%s.groupBy.%s",out.prefix,rd,"percent.mito"),
                      base_aspect_ratio = 1.30)

        ssc.plot.tsne(sce, columns = "nCount_RNA", 
                      reduced.name = sprintf("seurat.%s",rd),
                      colSet=colSet,size=0.03,
                      par.geneOnTSNE = list(scales = "free",pt.order = "random"),
                      vector.friendly=T,
                      out.prefix = sprintf("%s.seurat.%s.groupBy.%s",out.prefix,rd,"nCount_RNA"),
                      base_aspect_ratio = 1.30)

        ssc.plot.tsne(sce, columns = "nFeature_RNA", 
                      reduced.name = sprintf("seurat.%s",rd),
                      colSet=colSet,size=0.03,
                      par.geneOnTSNE = list(scales = "free",pt.order = "random"),
                      vector.friendly=T,
                      out.prefix = sprintf("%s.seurat.%s.groupBy.%s",out.prefix,rd,"nFeature_RNA"),
                      base_aspect_ratio = 1.30)
    }

    old.par <- NULL
    if(file.exists(sprintf("%s.run.par.txt",out.prefix))){
	    old.par <- fread(sprintf("%s.run.par.txt",out.prefix))
    }
    write.table(data.frame(opt.res=opt.res),
		file=sprintf("%s.run.par.txt",out.prefix),row.names=F,quote=F,sep="\t")

    if(file.exists(sprintf("%s.sce.rds",out.prefix)))
    {
        loginfo(sprintf("loading pre-calculated result ..."))
        sce <- readRDS(sprintf("%s.sce.rds",out.prefix))
        seu <- readRDS(sprintf("%s.seu.rds",out.prefix))
        if(old.par$opt.res==opt.res){
            return(list("seu"=seu,"sce"=sce))
        }
    }else
    {
        loginfo(sprintf("running Seurat pipeline ..."))

        seu <- run.HVG(seu,gene.exclude.df,n.top=n.top,measurement=measurement)

        adj.cov <- cor.var
        if(adj.cov[1]=="NULL") { adj.cov <- NULL }
        if(!is.null(adj.cov)){
        loginfo(sprintf("CellCycleScoring ..."))
        data("cc.genes",package="Seurat")
        seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes,
                                g2m.features = cc.genes$g2m.genes,
                                set.ident = FALSE)
        
        loginfo(sprintf("AddModuleScore ..."))
        dat.ext.dir <- system.file("extdata",package="scPip")
        glist.HSP <- fread(sprintf("%s/byZhangLab.stress.glist.gz",dat.ext.dir))$geneSymbol
        #glist.HSP <- intersect(glist.HSP,rownames(seu))
        seu <- AddModuleScore(seu, features=list("DIG.Score"=glist.HSP), name="DIG.Score",
                      pool = NULL, nbin = 24, ctrl = 100)

        glist.ISG <- fread(sprintf("%s/ISG.MSigDB.BROWNE_INTERFERON_RESPONSIVE_GENES.detected.glist.gz",dat.ext.dir))$geneSymbol
        seu <- AddModuleScore(seu, features=list("ISG.Score"=glist.ISG), name="ISG.Score",
                      pool = NULL, nbin = 24, ctrl = 100)

        ### if correct something, always correct for batchV and percent.mito
        if("percent.mito" %in% colnames(seu[[]])){
            adj.cov <- c(adj.cov,"percent.mito")
        }
        if(nBatch>1){
            adj.cov <- c("batchV",adj.cov)
        }
        }
        loginfo(sprintf("adj: %s\n",paste(adj.cov,collapse=",")))
        print(head(seu[[]]))
        loginfo(sprintf("do.scale: %s\n",do.scale))

        loginfo(sprintf("Scale ..."))
        if(use.sctransform && platform!="SmartSeq2"){
        cat(sprintf("SCTransform:\n"))
        seu <- SCTransform(seu, variable.features.n=n.top, do.scale=do.scale, vars.to.regress = adj.cov)
        }else{
        cat(sprintf("ScaleData:\n"))
        seu <- ScaleData(object = seu,do.scale=do.scale,vars.to.regress = adj.cov)
        }

        if(run.stage==0){
        return(list("seu"=seu,"sce"=sce))
        }

        #### PCA
        seu <- RunPCA(object = seu)
        #print(x = seu[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
        #p <- VizDimLoadings(seu)
        #ggsave(sprintf("%s.pca.01.pdf",out.prefix),width=7,height=14)

        seu <- ProjectDim(object = seu)

        #pdf(sprintf("%s.pca.02.pdf",out.prefix),width=14,height=18)
        #DimHeatmap(object = seu, dims = 1:15, cells = 500, balanced = TRUE,fast = TRUE)
        #dev.off()

        #p <- ElbowPlot(object = seu,ndims=50)
        #ggsave(sprintf("%s.pca.03.pdf",out.prefix),width=5,height=4)

        loginfo(sprintf("use.harmony: %s",use.harmony))
        use.rd <- "pca"
        if(use.harmony){
            seu <- RunHarmony(seu, c("batchV"),verbose=F)
            use.rd <- "harmony"
        }

        ######### UMAP
        tic("RunUMAP...")
        seu <- RunUMAP(object = seu, reduction = use.rd,dims = 1:opt.npc,verbose=F)
        toc()
        
        ######### tSNE 
        if("tsne" %in% plot.rd){
        tic("RunTSNE...")
        seu <- RunTSNE(object = seu, reduction = use.rd,dims = 1:opt.npc,verbose=F)
        toc()
        }

        #######################################
        
        #### clustring
        seu <- FindNeighbors(object = seu, reduction = use.rd, dims = 1:opt.npc)

        resolution.vec <- seq(0.1,2.4,0.1)
        seu <- FindClusters(object = seu,resolution = c(resolution.vec,res.addition),verbose=F)

        for(t.res in resolution.vec){
        if(use.sctransform && platform!="SmartSeq2"){
            print(table(seu[[sprintf("SCT_snn_res.%s",t.res)]]))
            aa.res <- sprintf("SCT_snn_res.%s",t.res)
        }else{
            print(table(seu[[sprintf("RNA_snn_res.%s",t.res)]]))
            aa.res <- sprintf("RNA_snn_res.%s",t.res)
        }
        #sce[[aa.res]] <- seu.merged@meta.data[,aa.res]
        }

        #### for sscClust
        reducedDim(sce,"seurat.pca") <- Embeddings(seu, reduction = "pca")
        reducedDim(sce,"seurat.umap") <- Embeddings(seu, reduction = "umap")
        if("tsne" %in% names(seu@reductions)){
        reducedDim(sce,"seurat.tsne") <- Embeddings(seu, reduction = "tsne")
        }

        if(use.sctransform && platform!="SmartSeq2"){
        idx.colRes <- grep("^SCT_snn_res",colnames(seu[[]]),value=T)
        }else{
        idx.colRes <- grep("^RNA_snn_res",colnames(seu[[]]),value=T)
        }
        for(idx in idx.colRes){
        colData(sce)[[idx]] <- sprintf("C%02d",as.integer(as.character(seu[[]][,idx])))
        }

        ### patch
        colData(sce)[is.na(sce$majorCluster) | sce$majorCluster=="unknown","majorCluster"] <- ""
        ### 

        sce$ClusterID <- colData(sce)[[opt.res]]
        print("all(colnames(sce)==colnames(seu))?")
        print(all(colnames(sce)==colnames(seu)))
        seu$ClusterID <- sce$ClusterID

        for(i.rd in plot.rd){
        plot.all(rd=i.rd)
        }
    }

    sce$ClusterID <- colData(sce)[[opt.res]]
    print("all(colnames(sce)==colnames(seu))?")
    print(all(colnames(sce)==colnames(seu)))
    seu$ClusterID <- sce$ClusterID

    for(i.rd in plot.rd){
	ssc.plot.tsne(sce, columns = "ClusterID",
		      reduced.name = sprintf("seurat.%s",i.rd),
		      colSet=list(),size=0.03,label=3,
		      par.geneOnTSNE = list(scales = "free",pt.order = "random"),
		      vector.friendly=T,
		      out.prefix = sprintf("%s.seurat.%s.groupBy.%s",out.prefix,i.rd,"ClusterID"),
		      base_aspect_ratio = 1.30)
    }



    ##### save result
    #write.table(colData(sce), sprintf("%s.clusterInfo.txt",out.prefix),sep = "\t",
    #			row.names = F, quote = F)
    saveRDS(seu, file = sprintf("%s.seu.rds",out.prefix))
    saveRDS(sce, file = sprintf("%s.sce.rds",out.prefix))
    #seu <- readRDS(file = sprintf("%s.seu.rds",out.prefix))
    #sce <- readRDS(file = sprintf("%s.sce.rds",out.prefix))
    cat(sprintf("data saved!\n"))

    #############################
    ## violin
    if(!all(colnames(seu)==colnames(sce))){
	warning("seu and sce are not consistent !")
    }
    
    nCls <- length(table(sce$ClusterID))
    l_ply(seq_along(g.geneOnUmap.list),function(i){
        gene.tmp <- intersect(g.geneOnUmap.list[[i]],rowData(sce)$display.name)
        if(length(gene.tmp)>0){
            p <- ssc.plot.violin(sce,assay.name=assay.name,
                     adjB=if(nBatch>1) "batchV" else NULL,
                     clamp = c(-4, 8), gene=gene.tmp, group.var="ClusterID")
            ggsave(sprintf("%s.seurat.violin.marker.%s.png",
                   out.prefix,names(g.geneOnUmap.list)[i]),
                   width=if(nCls<=12) 6 else 8,height=8)
        }
    },.parallel=T)

### DE Genes
####	tic("DE genes")
####	de.out <- ssc.clusterMarkerGene(sce, assay.name=assay.name,
####					ncell.downsample=if(ncol(sce)>30000) 1000 else NULL,
####					group.var="ClusterID",
####					batch=if(nBatch>1) "batchV" else NULL,
####					out.prefix=sprintf("%s.de.opt.res",out.prefix),
####					n.cores=ncores, do.plot=T,
####					F.FDR.THRESHOLD=0.01,
####					pairwise.P.THRESHOLD=0.01,
####					pairwise.FC.THRESHOLD=0.25, verbose=F)
####	toc()

    if(do.deg){
        dir.create(sprintf("%s/limma",dirname(out.prefix)),F,T)

        tic("limma")
        de.out <- ssc.DEGene.limma(sce,assay.name=assay.name,ncell.downsample=ncell.deg,
                       group.var="ClusterID",batch=if(nBatch>1) "batchV" else NULL,
                       out.prefix=sprintf("%s/limma/%s",
                                  dirname(out.prefix),basename(out.prefix)),
                       n.cores=ncores,verbose=3, group.mode="multiAsTwo",
                       T.logFC=if(platform=="SmartSeq2") 1 else 0.25)
        toc()

        saveRDS(de.out,file=sprintf("%s.de.out.limma.rda",
                        sprintf("%s/limma/%s",
                            dirname(out.prefix),
                            basename(out.prefix))))
    }

    return(list("seu"=seu,"sce"=sce))
}

#' Identification of gamma delta T cells (Fred's method)
#' @importFrom SingleCellExperiment rowData
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom sscClust classify.outlier
#' @param obj object of \code{SingleCellExperiment}; assay "counts" is required.
#' @param GSx character vector; gdT marker gene list. (default: c("CD3D","CD3E","TRDC","TRGC1","TRGC2")).
#' @param GSy character vector; CD8 T cell marker gene list. (default: c("CD8A","CD8B")).
#' @param col.name character; prefix of column names, to be added to obj. (default: "Score.gammaDeltaT")
#' @param th.score double; threshold of the score. (default: 0.35)
#' @param out.prefix character; output prefix. (default: NULL)
#' @return a SingleCellExperiment object
#' @export
cal.signatureScore.gdT.Fred <- function(obj,GSx=c("CD3D","CD3E","TRDC","TRGC1","TRGC2"),
					GSy=c("CD8A","CD8B"),
					col.name="Score.gammaDeltaT",th.score=0.35,out.prefix=NULL)
{
    if(!all(GSx %in% rowData(obj)$display.name) || !all(GSy %in% rowData(obj)$display.name)){
	    warning(sprintf("not all genes in GSx and GSy in the obj !!\n"))
	    return(obj)
    }

    if(!("counts" %in% assayNames(obj))){
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

#' Identification of gamma delta T cells (simple average-threshold method)
#' @importFrom SingleCellExperiment rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom data.table setDT
#' @importFrom ggplot2 ggplot geom_density theme_bw geom_vline facet_wrap ggsave aes_string
#' @param obj object of \code{SingleCellExperiment}
#' @param out.prefix character; output prefix. (default: NULL)
#' @param assay.name character vector; which assay to use. (default: "norm_exprs").
#' @param vis.v double vector; for vertical lines in visulization. (default: c(0.25,0.5)).
#' @param Th.CD3 double; threshold for T cell signature. (default: 0.25)
#' @param Th.DC double; threshold for delta receptor constant chain. (default: 0.25)
#' @param Th.GC1 double; threshold for gamma receptor constant chain 1. (default: 0.25)
#' @param Th.GC2 double; threshold for gamma receptor constant chain 2. (default: 0.25)
#' @return a SingleCellExperiment object
#' @export
inSilico.TGammaDelta <- function(obj,out.prefix=NULL,assay.name="norm_exprs",vis.v=c(0.25,0.5),
                           Th.CD3=0.25,Th.DC=0.25,Th.GC1=0.25,Th.GC2=0.25)
{
    #library("data.table")
    #library("ggplot2")
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

    if(!is.null(out.prefix)){
	dir.create(dirname(out.prefix),F,T)
	p <- ggplot(dat.plot.melt, aes_string(x="norm_exprs", fill = "gene", colour = "gene")) +
	    geom_density(alpha = 0.1) +
		    theme_bw() +
	    geom_vline(xintercept = vis.v,linetype=2) +
	    facet_wrap(~gene,ncol=3,scales="free_y")
	ggsave(sprintf("%s.inSiliso.marker.density.pdf",out.prefix),width=7,height=3.5)
    }

    obj$TGammaDelta <- "unknown"
    obj$TGammaDelta[dat.plot[["CD3"]] > Th.CD3 & dat.plot[["TRDC"]] > Th.DC ] <- "undetermined"
    obj$TGammaDelta[dat.plot[["CD3"]] > Th.CD3 & dat.plot[["TRDC"]] > Th.DC & ( dat.plot[["TRGC1"]] > Th.GC1 & dat.plot[["TRGC2"]] < Th.GC2 ) ] <- "delta2"
    obj$TGammaDelta[dat.plot[["CD3"]] > Th.CD3 & dat.plot[["TRDC"]] > Th.DC & ( dat.plot[["TRGC1"]] < Th.GC1 & dat.plot[["TRGC2"]] > Th.GC2 ) ] <- "delta1"
    table(obj$TGammaDelta)

    #bin.tb <- sscClust::classify.outlier(dat.plot[["TRDC"]],out.prefix=out.prefix,e.name="TRDC",th.score=0.25)
    ##write.table(colData(obj),sprintf("%s.cellInfo.txt",out.prefix),row.names=F,sep="\t",quote=F)
    return(obj)
}

#' Identification of T cells (simple average-threshold method)
#' @importFrom SingleCellExperiment rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom sscVis ssc.scale
#' @importFrom data.table setDT
#' @importFrom ggpubr ggdensity
#' @importFrom ggplot2 ggplot geom_density theme_bw theme geom_vline facet_wrap ggsave
#' @param sce object of \code{SingleCellExperiment}
#' @param out.prefix character; output prefix. (default: NULL)
#' @param assay.name character vector; which assay to use. (default: "norm_exprs").
#' @param vis.v double vector; for vertical lines in visulization. (default: c(0.25,0.5)).
#' @param Th.CD3 double; threshold for T cell signature (mean of CD3D, CD3G). (default: 0.25)
#' @param Th.CD8 double; threshold for CD8 (mean of CD8A, CD8B). (default: 0.25)
#' @param Th.CD4 double; threshold for CD4 (expression of CD4). (default: 0.25)
#' @param Th.TH double; threshold for Thelper signature (mean of CD4, CD40LG). (default: 0.25)
#' @param Th.TR double; threshold for Treg signature (mean of CD4, FOXP3). (default: 0.25)
#' @param Th.DC double; threshold for delta receptor constant chain. (default: 0.25)
#' @param Th.GC1 double; threshold for gamma receptor constant chain 1. (default: 0.25)
#' @param Th.GC2 double; threshold for gamma receptor constant chain 2. (default: 0.25)
#' @param do.zscore logical; whether use zscore for calculation. (default: FALSE)
#' @param do.rescue logical; whether use "rescue" mode. (default: FALSE)
#' @return a SingleCellExperiment object
#' @details columns stype and gdType will be added to the sce
#' @export
inSilico.TCell <- function(sce, out.prefix, assay.name="norm_exprs",vis.v=c(0.25,0.5),
                           Th.CD3=0.25,Th.CD8=0.5,Th.CD4=0.5,Th.TH=0.25,Th.TR=0.25,do.zscore=F,
			   Th.DC=0.25,Th.GC1=0.25,Th.GC2=0.25,
			   do.rescue=T)
{
    #### in silico classification
    #library("data.table")
    #library("ggplot2")
    #library("ggpubr")
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

    if(!file.exists(dirname(out.prefix))){
        dir.create(dirname(out.prefix),F,T)
    }
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

#' calculate the signature score of one potential contaminated cell types (simple average-threshold method)
#' @importFrom SummarizedExperiment assay `colData<-`
#' @importFrom Seurat GetAssayData
#' @importFrom data.table data.table melt
#' @importFrom ggplot2 ggplot geom_density geom_vline facet_wrap ggsave aes_string
#' @param obj object of \code{SingleCellExperiment} or \code{Seurat}
#' @param out.prefix character; output prefix. (default: NULL)
#' @param assay.name character vector; which assay to use. (default: "norm_exprs").
#' @param g.name character; signature name of the cell type. (default: "plasmaB").
#' @param g.test character vector; signature genes of the cell type. (default: c("CD79A", "JCHAIN", "SDC1")).
#' @param score.t double; threshold for signature score (mean of g.test). (default: 1)
#' @param vis.v double vector; for vertical lines in visulization. (default: c(0.25,0.5,1)).
#' @return an object of \code{SingleCellExperiment} or \code{Seurat}
#' @details Two columns will be added to the obj, "%s.score" and "%s.class, where %s is the g.name.
#' @export
fill.contamination <- function(obj,out.prefix,assay.name="norm_exprs",
			       g.name="plasmaB",g.test=c("CD79A", "JCHAIN", "SDC1"),
			       score.t=1,vis.v=c(0.25,0.5,1))
{
    #require("ggplot2")
    #require("data.table")
    ##### plasma contamination
    if(class(obj)[1]=="SingleCellExperiment"){
        exp.data <- assay(obj,assay.name)[g.test,,drop=F]
        B.score.tb <- data.table(cell=colnames(obj),B.score= colMeans(exp.data))
        B.score.tb <- cbind(B.score.tb,as.matrix(t(exp.data)))
        colData(obj)[[sprintf("%s.score",g.name)]] <- B.score.tb$B.score
        colData(obj)[[sprintf("%s.class",g.name)]] <- obj[[sprintf("%s.score",g.name)]] > score.t
    }else if(class(obj)[1]=="Seurat"){
        exp.data <- GetAssayData(obj,"data")[g.test,,drop=F]
        B.score.tb <- data.table(cell=colnames(obj),B.score= colMeans(exp.data))
        B.score.tb <- cbind(B.score.tb,as.matrix(t(exp.data)))
        obj[[sprintf("%s.score",g.name)]] <- B.score.tb$B.score
        obj[[sprintf("%s.class",g.name)]] <- obj[[sprintf("%s.score",g.name)]] > score.t
    }else{
        return(obj)
    }

    B.score.tb.melt <- melt(B.score.tb,id="cell")
    colnames(B.score.tb.melt) <- c("cell","gene","norm_exprs")

    if(!file.exists(dirname(out.prefix))){
        dir.create(dirname(out.prefix),F,T)
    }
    p <- ggplot(B.score.tb.melt, aes_string(x="norm_exprs", fill = "gene", colour = "gene")) +
		    geom_density(alpha = 0.1) +
		    geom_vline(xintercept = vis.v,linetype=2) +
		    facet_wrap(~gene,ncol=2,scales="free_y")
    ggsave(sprintf("%s.marker.%s.density.pdf",out.prefix,g.name),width=8,height=8)
    return(obj)
}

#' calculate the proliferation score
#' @importFrom SingleCellExperiment rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom sscClust classify.outlier
#' @importFrom data.table as.data.table data.table
#' @param obj object of \code{SingleCellExperiment}
#' @param gene.prol character vector; genes to use.
#' @param assay.name character vector; which assay to use. (default: "norm_exprs").
#' @param out.prefix character; output prefix. (default: NULL)
#' @param method character; method to use. (default: "mean").
#' @return a list
#' @details calculate the proliferation score.
#' @export
calProliferationScore <- function(obj,gene.prol,assay.name="norm_exprs",out.prefix=NULL,method="mean")
{
    f.gene <- rowData(obj)[,"display.name"] %in% gene.prol
    exp.sub <- as.matrix(assay(obj[f.gene,],assay.name))
    f.zero <- apply(exp.sub,1,function(x){ all(x==0) })
    if(sum(f.zero) > 0){
	    cat(sprintf("Number of gene with value zero in all cells: %d\n",sum(f.zero)))
    }
    exp.sub <- exp.sub[!f.zero,]

    dat.score <- NULL
    out.tb <- NULL
    if(method=="mean")
    {
        score.prol <- colMeans(exp.sub)
        dat.score <- classify.outlier(score.prol,out.prefix=out.prefix)
        out.tb <- as.data.table(dat.score$score.cls.tb)
        colnames(out.tb)[1] <- "cellID"
    }
####    else if(method=="AUCell") {
####        require("AUCell")
####        #####
####        #f.gene <- rowData(obj)[,"display.name"] %in% gene.prol
####        #exp.sub <- as.matrix(assay(obj,assay.name))
####        exp.sub <- (assay(obj,assay.name))
####        rownames(exp.sub) <- rowData(obj)[,"display.name"]
####        f.zero <- apply(exp.sub,1,function(x){ all(x==0) })
####        if(sum(f.zero) > 0){
####            cat(sprintf("Number of gene with value zero in all cells: %d\n",sum(f.zero)))
####        }
####        exp.sub <- exp.sub[!f.zero,]
####
####        #####
####        pdf(sprintf("%s.buildRankings.1.pdf",out.prefix),width=7,height=7)
####        cells_rankings <- AUCell_buildRankings(exp.sub)
####        dev.off()
####
####        geneSets <- list("prol"=intersect(rownames(exp.sub),gene.prol))
####        #### geneSets <- GSEABase::GeneSet(genes, setName="geneSet1") # alternative
####        cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
####
####        pdf(sprintf("%s.exploreThresholds.1.pdf",out.prefix),width=7,height=7)
####        #par(mfrow=c(3,3))
####        set.seed(123)
####        cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
####        dev.off()
####
####        if("L_k2" %in% rownames(cells_assignment$prol$aucThr$thresholds)){
####            th.prol <- cells_assignment$prol$aucThr$thresholds["L_k2","threshold"]
####        }else{
####            th.prol <- 0.09
####        }
####
####        ##geneSetName <- rownames(cells_AUC)[grep("prol", rownames(cells_AUC))]
####        pdf(sprintf("%s.exploreThresholds.2.pdf",out.prefix),width=7,height=7)
####        AUCell_plotHist(cells_AUC["prol",], aucThr=th.prol)
####        abline(v=th.prol)
####        dev.off()
####        
####        out.tb <- data.table(cellID=colnames(cells_AUC),
####                     proliferationScore.bin=as.integer(getAUC(cells_AUC)["prol",]>th.prol),
####                     proliferationScore=getAUC(cells_AUC)["prol",])
####        out.tb$classification <- out.tb$proliferationScore.bin
####        dat.score <- cells_AUC
####
####        my.dat.score <- classify.outlier(getAUC(cells_AUC)["prol",],out.prefix=out.prefix)
####        out.tb[,myCls:=my.dat.score$score.cls.tb[cellID,"classification"]]
####
####    }
    return(list("out.tb"=out.tb,"detail"=dat.score))
}


