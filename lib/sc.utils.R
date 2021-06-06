
require("Seurat")
require("tictoc")
require("tictoc")
require("dplyr")
require("cowplot")
require("plyr")
require("ggplot2")

###################### variable genes
run.HVG <- function(seu,gene.exclude.df,n.top=1500,measurement="counts")
{
    seu <- FindVariableFeatures(object = seu)
    if(measurement=="TPM"){
	#### TPM
	hvg.gene.info <- seu@assays$RNA@meta.features %>% 
		tibble::rownames_to_column(var="geneSymbol") %>%
		arrange(desc(vst.variance))
    }else{
	#### counts, cpm
	hvg.gene.info <- seu@assays$RNA@meta.features %>% 
		tibble::rownames_to_column(var="geneSymbol") %>%
		arrange(desc(vst.variance.standardized))
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
    hvg.gene.info.flt %>% head %>% print
    seu@assays$RNA@misc$meta.features.flt <- hvg.gene.info.flt
    return(seu)
}

### special genes
{

    g.geneOnUmap.list <- list(
			      "epithelial.01"=c("EPCAM", "SFN", "KRT8", "KRT10","KRT18",
						"KRT19","KRTCAP3","KRT15",""),
			      "endothelial"=c("PECAM1","VWF","CDH5","PLVAP","CLDN5",
					      "SLCO2A1","RAMP2","EGFL7","EMCN"),
			      "CAF"=c("DCN", "COL1A2", "COL1A1", "COL3A1", "COL6A2", "COL6A1",
				      "PCOLCE", "C1S", "LUM"),
			      "SMC"=c("MYH11", "ACTA2", "TAGLN", "CNN1", "MYL9",
				      "DES", "TPM2", "MYLK", "PLN"),
			      "B.01"=c("BLK", "CD19", "CD70", "CD79A","CD79B",
				       "MS4A1","FCER2","SDC1","JCHAIN"),
			      "Mast"=c("ENPP3", "KIT", "PTPRC", "TPSB2","SLC18A2",
				       "FCER1A", "TPSAB1", "CPA3", "HPGDS"),
			      "Neutrophils"=c("FUT4","CD3E","TLR2","CD14","CD33","FCGR1A","CSF3R","",""),
			      "Monocyte"=c("CD14","VCAN", "S100A9","FCGR3A","FCEG1G",
					   "LST1","ITGAM", "HLA-DQA1","PTPRC"),
			      "DC.01"=c("LILRA4","CD1C", "CLEC9A","CLEC10A","FCER1A",
					"RSG2","LAMP3", "THBD", "IDO1"),
			      "DC.02"=c("ITGAX","IRF4", "IRF8", "CD274", "LTB",
				  "CCL19","CCL21","IL3RA","NRP1","CLE4C"),
			      "Macrophase.01"=c("CD68","THBS1","CD163", "C1QA","APOE",
						"FCGR3A","MARCO", "S100A4","S100A9"),
			      "Macrophase.02"=c("HLA-DQA1", "CST3", "LYZ", "ISGF6","FUT4",
						"SEPP1","IFI30","FCGR1A","GPNMB"),
			      "Macrophase.03"=c("CTSK", "MMP9", "SIGLEC15", "LYVE1","SPP1",
						"FOLR2","C1QC","TREM2","LPL"),
			      "Macrophase.04"=c("MT1G", "CCL18", "CD36", "CD52","FABP4",
						"C3","CXCL8","PDK4","HK2"),
			      "Macrophase.05"=c("FCN1", "IL1B", "NLRP3", "PLTP","CCR2",
						"CCL2","CCL18","MMP9","SLC2A1"),
			      "Macrophase.06"=c("GCHFR","CTSD","CTSL","CTSB", "APOE","APOC1",
						"CCL2","CCL18","SLC2A1"),
			      "Macrophase.07"=c("S100A8","S100A9", "CXCL10", "CXCL11", "OLR1",
					      "IDO1","CALHM6","FN1","SERPINA1"),
			      "Macrophase.08"=c("CD14","FCGR3A", "CX3CR1", "MT1G", "CCL2",
						"SLC2A1","FABP4","FOLR2","LYVE1"),
			      "Macrophase.09"=c("FCN1","IL1B","GPNMB","C1QC","SPP1","TREM2","LYZ","CSF1R","CD163"),
			      ##"Granulocyte"=c("CEACAM1","CEACAM6","CEACAM3","ENPP3","CD3E","","","",""),
			      "NK"=c("FCGR3A", "NCR1","NCAM1","CD3E", "NKG7",
				     "KLRD1","KIR2DL3","KIR2DL4","KLRB1"),
			      #"ILC"=c("IL4I1","KIT","KRT86","KRT81","ZFP36L1", "IL23R","CD83","NR4A1","EEF1A1"),
			      "T.01"=c("CD4", "CD8A", "CD8B", "CD2", "CD3D",
				       "CD3G", "NKG7", "GNLY", "FGFBP2"),
			      "T.02"=c("CX3CR1", "GZMK", "ZNF683","CD40LG","FOXP3",
				       "TNFRSF9","IL7R","IL2RA","TYROBP"),
			      "T.03"=c("TRGC1","TRGC2","TRBC1","TRBC2","TRDC",
				       "TRAC","TRGV9","TRDV2","TRAV1-2"),
			      "T.04"=c("CCR7","GZMK","LAYN","HAVCR2","CXCL13",
				       "CX3CR1","ZNF683","SLC4A10","CD160"),
			      "T.05"=c("IL17A","RORC","IL23R","CCR6","BCL6",
				       "CXCR5", "CXCR6","XCL1","XCL2"),
			      "T.06"=c("TCF7","MAL","CCL5","CCL4","CCL4L2",
				       "GZMA","GPR183","EOMES","KLRG1"),
			      "T.07"=c("PDCD1","CTLA4","IFNG","TNF","GZMB",
				       "RGS1","S1PR1","ITGAE","MALAT1"),
			      "cellCycle"=c("TYMSOS","PCLAF","TYMS","ZWINT","MCM2",
					    "MIR3917","UHRF1","TK1","MKI67"),
			      "ISG"=c("ISG15","ISG20","IFIT1","IFIT2",
				      "IRF7", "OAS1","OAS2","STAT1"),
			      "DIG"=c("DNAJB1","DUSP1","FOS","JUN","HSP90AA1",
				      "HSP90AB1","HSPA1A","HSPA6")
						    )
}

run.Seurat3 <- function(seu,sce,out.prefix,gene.exclude.df,n.top=1500,
			measurement="counts",platform="10X",
			opt.res="1",use.sctransform=F,aid="PRJ",plot.rd=c("umap","tsne"),
			opt.npc=15,ncores=16,do.adj=T,ncell.deg=1500,do.deg=F,do.scale=F,
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
        gene.mapping.table <- data.table(geneID=rownames(seu),seu.id=rownames(seu),display.name=rownames(seu))
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
	    rowData(sce) %>% head %>% print

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
        idx.patient <- grep("^(patient|Patient)$",colnames(seu[[]]),value=T,perl=T)
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

	    plot.resolution.list.debug <<- plot.resolution.list
	    out.prefix.debug <<- out.prefix
	    for(i in seq_len(length(plot.resolution.list)/4))
	    {
            pp <- plot_grid(plotlist=plot.resolution.list[((i-1)*4+1):(i*4)],
                            ncol = 2,align = "hv")
            save_plot(sprintf("%s.%s.res.sceStyle.%d.png",out.prefix,rd,i),pp,
                      ncol = 2, base_aspect_ratio=0.9,base_height=5.5)
	    }

	    ## gene on umap
	    l_ply(seq_along(g.geneOnUmap.list),function(i){
            gene.tmp <- intersect(g.geneOnUmap.list[[i]],rowData(sce)$display.name)
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
	    },.parallel=T)

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
	    if(is.character(do.adj)){
            adj.cov <- c("S.Score","G2M.Score","DIG.Score1",do.adj)
            ####adj.cov <- c("PROL.Score","DIG.Score1",do.adj)
            do.adj <- T
	    }else{
            adj.cov <- c("S.Score","G2M.Score","DIG.Score1")
	    }

	    if(do.adj){
            loginfo(sprintf("CellCycleScoring ..."))
            seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
            
            loginfo(sprintf("AddModuleScore ..."))
            ##glist.HSP <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/geneSet/exclude/byHanjieLi.stress.glist")$geneSymbol
            src.dir <- sprintf("%s/../",this.dir())
            glist.HSP <- fread(sprintf("%s/data/geneSet/byZhangLab.stress.glist.gz",src.dir))$geneSymbol
            #glist.HSP <- intersect(glist.HSP,rownames(seu))
            seu <- AddModuleScore(seu, features=list("DIG.Score"=glist.HSP), name="DIG.Score",
                                  pool = NULL, nbin = 24, ctrl = 100)

            glist.ISG <- fread(sprintf("%s/data/geneSet/ISG.MSigDB.BROWNE_INTERFERON_RESPONSIVE_GENES.detected.glist.gz",src.dir))$geneSymbol
            seu <- AddModuleScore(seu, features=list("ISG.Score"=glist.ISG), name="ISG.Score",
                                  pool = NULL, nbin = 24, ctrl = 100)

            #glist.prol <- fread(sprintf("%s/data/geneSet/byDE.cellCycle.hc.glist.gz",src.dir))$geneSymbol
            #glist.prol <- head(glist.prol,n=30)
            #seu <- AddModuleScore(seu, features=list("PROL.Score"=glist.prol), name="PROL.Score",
            #				      pool = NULL, nbin = 24, ctrl = 100)

            if("percent.mito" %in% colnames(seu[[]])){
                adj.cov <- c(adj.cov,"percent.mito")
            }
            if(nBatch>1){
                adj.cov <- c("batchV",adj.cov)
            }
            cat(sprintf("adj: %s\n",paste(adj.cov,collapse=",")))
            print(head(seu[[]]))
	    }else{
            adj.cov <- NULL
	    }
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

	    ######### UMAP
	    tic("RunUMAP...")
	    seu <- RunUMAP(object = seu, reduction = "pca",dims = 1:opt.npc,verbose=F)
	    toc()
	    
	    ######### tSNE 
	    if("tsne" %in% plot.rd){
            tic("RunTSNE...")
            seu <- RunTSNE(object = seu, reduction = "pca",dims = 1:opt.npc,verbose=F)
            toc()
	    }

	    #######################################
	    
	    #### clustring
	    seu <- FindNeighbors(object = seu, reduction = "pca", dims = 1:opt.npc)

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
	    colData(sce)[sce$majorCluster=="unknown","majorCluster"] <- ""
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
	    saveRDS(de.out,file=sprintf("%s.de.out.limma.rda",
					sprintf("%s/limma/%s",
						dirname(out.prefix),
						basename(out.prefix))))
	    toc()
	}

	return(list("seu"=seu,"sce"=sce))
}

cal.signatureScore.gdT.Fred <- function(obj,GSx=c("CD3D","CD3E","TRDC","TRGC1","TRGC2"),
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



