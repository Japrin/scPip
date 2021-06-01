
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
			 "B.01"=c("BLK", "CD19", "CD70", "CD79A","CD79B",
				     "MS4A1","FCER2","SDC1","JCHAIN"),
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
			"Monocyte"=c("CD14","VCAN", "S100A9","FCGR3A","FCEG1G",
				     "LST1","ITGAM", "HLA-DQA1","PTPRC"),
			"DC.01"=c("LILRA4","CD1C", "CLEC9A","CLEC10A","FCER1A",
				  "RSG2","LAMP3", "THBD", "IDO1"),
			"DC.02"=c("ITGAX","IRF4", "IRF8", "CD274", "LTB",
				  "CCL19","CCL21","IL3RA","NRP1","CLE4C"),
			"Mast"=c("ENPP3", "KIT", "PTPRC", "TPSB2","SLC18A2",
				 "FCER1A", "TPSAB1", "CPA3", "HPGDS"),
			"Granulocyte"=c("CEACAM1","CEACAM6","CEACAM3","ENPP3","CD3E","","","",""),
			"Neutrophils"=c("FUT4","CD3E","TLR2","CD14","CD33","FCGR1A","CSF3R","",""),
			"NK"=c("FCGR3A", "NCR1","NCAM1","CD3E", "NKG7",
			       "KLRD1","KIR2DL3","KIR3DL2","KLRB1"),
			"ILC"=c("IL4I1","KIT","KRT86","KRT81","ZFP36L1",
				"IL23R","CD83","NR4A1","EEF1A1"),
			"T.01"=c("CD4", "CD8A", "CD8B", "CD2", "CD3D",
				 "CD3G", "NKG7", "GNLY", "FGFBP2"),
			"T.02"=c("CX3CR1", "GZMK", "ZNF683", "IL7R","GZMB",
				 "CD40LG","FOXP3","CCR7","MKI67"),
			"T.03"=c("TRGC1","TRGC2","TRBC1","TRBC2","TRDC",
				 "TRAC","FOXP3","IL2RA","CD40LG"),
			"T.04"=c("CCR7","GZMK","LAYN","HAVCR2","CXCL13",
				 "CX3CR1","ZNF683","S1PR1","RGS1"),
			"T.05"=c("CXCR6","SLC4A10","CCL5","CD160","CD6",
				 "TMIGD2","MAL","SELL","KIR2DL4"),
			"T.06"=c("IL17A","GZMA","CXCR6","HAVCR2","CXCR5",
				 "BCL6","RORC","IL23R","IFNG"),
			"T.07"=c("SELL","MAL","TCF7","CTLA4","IL10",
				 "IL2RA","ANXA1","ANXA2","TNFRSF9"),
			"T.08"=c("CCL4L2","CCL4","XCL1","XCL2","CRTAM",
				 "RGCC","LMNA","FCER1G","KRT81"),
			"T.09"=c("ENTPD1","LAYN","CTLA4","TNFRSF9", "MYO1E",
				 "MYO7A", "KRT86", "RBPJ", "TIGIT"),
			"T.10"=c("SIRPG", "GEM", "PHLDA1", "TNFRSF18", "AFAP1L2",
				 "DUSP4", "MIR155HG","PDCD1","TOX"),
			##"CAF"=c("COL1A1","COL1A2","COL6A1","COL6A2","COL6A3",
			##	    "DCN","FAP","THY1","LUM"),
			"CAF"=c("DCN", "COL1A2", "COL1A1", "COL3A1", "COL6A2", "COL6A1",
				   "PCOLCE", "C1S", "LUM"),
			"endothelial"=c("PECAM1","VWF","CDH5","PLVAP","CLDN5",
					"SLCO2A1","RAMP2","EGFL7","EMCN"),
			"SMC"=c("MYH11", "ACTA2", "TAGLN", "CNN1", "MYL9",
				"DES", "TPM2", "MYLK", "PLN"),
			"cellCycle"=c("TYMSOS","PCLAF","TYMS","ZWINT","MCM2",
				      "MIR3917","UHRF1","TK1","MKI67"),
			"IFN"=c("ISG15","ISG20","IFIT1","IFIT2","IRF7",
				"OAS1","OAS2","STAT1"),
			"stress"=c("DNAJB1","DUSP1","FOS","JUN","HSP90AA1",
				   "HSP90AB1","HSPA1A","HSPA6"),
			"other.01"=c("PTPRC", "RGS1", "S1PR1", "ITGAE","HAVCR2",
				     "CXCL13","TNFRSF9","CTLA4","MALAT1")
						    )
}

### see that defined in cancer/lib/scRNAToolKit.R
####inSilico.TCell <- function(sce, out.prefix, assay.name="norm_exprs",vis.v=c(0.25,0.5),
####                           Th.CD3=0.25,Th.CD8=0.5,Th.CD4=0.5,Th.TH=0.25,Th.TR=0.25,rescue.CD4=T)
####{
####    #### in silico classification
####    library("data.table")
####    library("ggplot2")
####    ##gene.to.test <- c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G","CD40LG","FOXP3","IL2RA")
####    gene.to.test <- c("CD4","CD8A","CD8B","CD3D","CD3G","CD40LG","FOXP3","IL2RA")
####    f.gene <- which(rowData(sce)$display.name %in% gene.to.test)
####    gene.to.test <- structure(rowData(sce)$display.name[f.gene],names=rownames(sce)[f.gene])
####    gene.to.test <- gene.to.test[order(gene.to.test)]
####
####    dat.plot <- as.data.frame(t(as.matrix(assay(sce,assay.name)[names(gene.to.test),])))
####    colnames(dat.plot) <- gene.to.test
####    ##dat.plot[,"CD3"] <- apply(dat.plot[,c("CD3D","CD3E","CD3G")],1,mean)
####    dat.plot[,"CD3"] <- apply(dat.plot[,c("CD3D","CD3G")],1,mean)
####    dat.plot[,"CD8"] <- apply(dat.plot[,c("CD8A","CD8B")],1,mean)
####    dat.plot[,"TH"] <- apply(dat.plot[,c("CD4","CD40LG")],1,mean)
####    dat.plot[,"TR"] <- apply(dat.plot[,c("CD4","FOXP3")],1,mean)
####    dat.plot.melt <- melt(as.matrix(dat.plot))
####    colnames(dat.plot.melt) <- c("cell","gene","norm_exprs")
####
####    p <- ggplot(dat.plot.melt, aes(norm_exprs, fill = gene, colour = gene)) +
####        geom_density(alpha = 0.1) +
####        geom_vline(xintercept = vis.v,linetype=2) +
####        facet_wrap(~gene,ncol=3,scales="free_y")
####    ggsave(sprintf("%s.inSiliso.marker.density.pdf",out.prefix),width=7,height=8)
####
####    sce$stype <- "unknown"
####    ####sce$stype[dat.plot[,"CD3"] < 0.25] <- "noCD3"
####    sce$stype[dat.plot[,"CD3"] > Th.CD3 & dat.plot[,"CD8"] > Th.CD8 & dat.plot[,"CD4"] < Th.CD4] <- "CD8"
####    sce$stype[dat.plot[,"CD3"] > Th.CD3 & dat.plot[,"CD8"] > Th.CD8 & dat.plot[,"CD4"] > Th.CD4] <- "DP"
####    sce$stype[dat.plot[,"CD3"] > Th.CD3 & dat.plot[,"CD8"] < Th.CD8 & dat.plot[,"CD4"] > Th.CD4] <- "CD4"
####    sce$stype[dat.plot[,"CD3"] > Th.CD3 & dat.plot[,"CD8"] < Th.CD8 & dat.plot[,"CD4"] < Th.CD4] <- "DN"
####    if(rescue.CD4==T){
####	### todo: add requirement: DN and not gamma delta T
####	sce$stype[sce$stype=="DN" & (dat.plot[,"TH"] > Th.TH | dat.plot[,"TR"] > Th.TR)] <- "CD4"
####    }
####    ##sce$stype[dat.plot[,"CD3"] > 0.25 & dat.plot[,"CD8"] < 0.5 & dat.plot[,"CD4"] < 0.5 & dat.plot[,"TH"] > 0.25 ] <- "TH"
####    ##sce$stype[dat.plot[,"CD3"] > 0.25 & dat.plot[,"CD8"] < 0.5 & dat.plot[,"CD4"] < 0.5 & dat.plot[,"TR"] > 0.25 ] <- "TR"
####    table(sce$stype)
####
####    write.table(colData(sce),sprintf("%s.cellInfo.txt",out.prefix),row.names=F,sep="\t",quote=F)
####    return(sce)
####}
####


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
	if(is.null(sce) && !is.null(seu) && !is.null(gene.mapping.table)){
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
	    rowData(sce) %>% head

	    all(colnames(seu)==rownames(colData(sce)))
	    all(rownames(seu)==rowData(sce)$seu.id)
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
	nBatch <- length(table(sce$batchV))

	#####
	
	RhpcBLASctl::omp_set_num_threads(1)
	doParallel::registerDoParallel(cores = ncores)
	
	plot.all <- function(rd="umap",resolution.vec=seq(0.1,2.4,0.1))
	{
	    p <- DimPlot(seu,reduction=rd, group.by = "loc")
	    ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"loc"),width=5,height=4)

	    if(!is.null(seu@meta.data[["stype"]])){
		    p <- DimPlot(seu,reduction=rd, group.by = "stype")
		    ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"stype"),width=5,height=4)
	    }
	    
	    p <- DimPlot(seu,reduction=rd, group.by = "cancerType")
	    ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"cancerType"),width=5,height=4)

	    if(length(unique(seu$libraryID))>20){
		    p <- DimPlot(seu,reduction=rd,label=F,label.size=2,
				 group.by = "libraryID") + NoLegend()
		    ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"libraryID"),width=4.5,height=4)
	    }else{
		    p <- DimPlot(seu,reduction=rd,label=T,label.size=2, group.by = "libraryID")
		    ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"libraryID"),width=6.5,height=4)
	    }

	    if(length(unique(seu$patient))>20){
		    p <- DimPlot(seu,reduction=rd, group.by = "patient") + NoLegend()
		    ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"patient"),width=4.5,height=4)
	    }else{
		    p <- DimPlot(seu,reduction=rd, group.by = "patient")
		    ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"patient"),width=6.5,height=4)
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
		glist.HSP <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/geneSet/exclude/byZhangLab.stress.glist")$geneSymbol
		#glist.HSP <- intersect(glist.HSP,rownames(seu))
		seu <- AddModuleScore(seu, features=list("DIG.Score"=glist.HSP), name="DIG.Score",
						      pool = NULL, nbin = 24, ctrl = 100)

		glist.ISG <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/geneSet/exclude/ISG.MSigDB.BROWNE_INTERFERON_RESPONSIVE_GENES.detected.glist")$geneSymbol
		seu <- AddModuleScore(seu, features=list("ISG.Score"=glist.ISG), name="ISG.Score",
						      pool = NULL, nbin = 24, ctrl = 100)

		#glist.prol <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/geneSet/exclude/byDE.cellCycle.hc.glist")$geneSymbol
		#glist.prol <- head(glist.prol,n=30)
		#seu <- AddModuleScore(seu, features=list("PROL.Score"=glist.prol), name="PROL.Score",
		#				      pool = NULL, nbin = 24, ctrl = 100)

		###adj.cov <- c("S.Score","G2M.Score","DIG.Score1")
		####adj.cov <- c()
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
	    print(x = seu[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
	    p <- VizDimLoadings(seu)
	    ggsave(sprintf("%s.pca.01.pdf",out.prefix),width=7,height=14)

	    seu <- ProjectDim(object = seu)

	    pdf(sprintf("%s.pca.02.pdf",out.prefix),width=14,height=18)
	    DimHeatmap(object = seu, dims = 1:15, cells = 500, balanced = TRUE,fast = TRUE)
	    dev.off()

	    p <- ElbowPlot(object = seu,ndims=50)
	    ggsave(sprintf("%s.pca.03.pdf",out.prefix),width=5,height=4)

	    ######### UMAP
	    tic("RunUMAP...")
	    seu <- RunUMAP(object = seu, reduction = "pca",dims = 1:opt.npc)
	    toc()
	    
	    ######### tSNE 
	    if("tsne" %in% plot.rd){
		tic("RunTSNE...")
		seu <- RunTSNE(object = seu, reduction = "pca",dims = 1:opt.npc)
		toc()
	    }

	    ################# in silico FACS #####
####	    sce <- inSilico.TCell(sce, out.prefix, assay.name=assay.name,vis.v=c(0.25,0.5),
####				  Th.CD3=0.25,Th.CD8=0.5,Th.CD4=0.5,Th.TH=0.25,Th.TR=0.25)
####	    table(sce$stype)
####	    seu$stype <- sce$stype
####	    table(seu$stype)
	    #######################################
	    
	    #### clustring
	    seu <- FindNeighbors(object = seu, reduction = "pca", dims = 1:opt.npc)

	    ##resolution.vec <- c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6)
	    #resolution.vec <- c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,2.1,2.2,2.3,2.4)
	    ##resolution.vec <- seq(0.5,2.4,0.1)
	    #resolution.vec <- seq(0.1,0.4,0.1)
	    resolution.vec <- seq(0.1,2.4,0.1)
	    seu <- FindClusters(object = seu,resolution = c(resolution.vec,res.addition))
	    #seu <- FindClusters(object = seu,resolution = c(resolution.vec,5,10,25))

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
				       n.cores=ncores,verbose=3, group.mode="multi",
				       T.logFC=if(platform=="SmartSeq2") 1 else 0.25)
	    saveRDS(de.out,file=sprintf("%s.de.out.limma.rda",
					sprintf("%s/limma/%s",
						dirname(out.prefix),
						basename(out.prefix))))
	    toc()
	}

	return(list("seu"=seu,"sce"=sce))
}

divid.CD4CD8 <- function(sce,out.prefix,gene.use,
			measurement="counts",platform="10X",
			opt.res="0.6",use.sctransform=F,aid="PRJ",plot.rd=c("umap","tsne"),
			opt.npc=5,ncores=16,do.adj=F,
			gene.mapping.table=NULL,res.addition=NULL)
{

    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = ncores)

    {
	gene.CD4CD8.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/zhangLab.10X/divid.CD4.CD8/OUT.sigGene/T.CD4CD8.by.Tn.Temra.txt.gz")
	gene.CD4CD8.tb[,geneID:=as.character(geneID)]
	gene.DNDP.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/human.thymic.development.JongEunPark2020/OUT.JongEunPark2020/JongEunPark2020.deg.T.set4.top.txt")
	gene.DNDP.tb <- gene.DNDP.tb[cluster %in% c("DP","DN"),]
	gene.DNDP.tb[,geneID:=as.character(geneID)]
	gene.DNDP.tb <- gene.DNDP.tb[,c("geneID","cluster","geneSymbol")]
	colnames(gene.DNDP.tb)[2] <- "category"
	gene.div.tb <- rbind(gene.CD4CD8.tb,gene.DNDP.tb)
	gene.div.tb$ENSG <- rownames(sce)[match(gene.div.tb$geneSymbol,rowData(sce)$display.name)]
	gene.div.tb <- gene.div.tb[!is.na(ENSG),]
	f.gene <- intersect(gene.div.tb$geneSymbol,rowData(sce)$display.name)
	gene.div.tb <- gene.div.tb[geneSymbol %in% f.gene,]

	geneOnMap <- list("T01"=c("CD4", "CD8A", "CD8B","CCR7","SELL","CX3CR1","CD40LG","KLRC4","CTSW"),
			       "T02"=c("CD4", "CD8A", "CD8B","TRAT1","RNASET2","CTSB","CD40LG","APBA2","KLRC3"),
			       "T.DN"=head(gene.div.tb[category=="DN",][["geneSymbol"]],n=9),
			       "T.DP"=head(gene.div.tb[category=="DP",][["geneSymbol"]],n=9)
			       )
	gene.use <- unique(c(gene.div.tb$geneSymbol,sapply(geneOnMap,function(x){ x  })))
	sce.tmp <- sce[rowData(sce)$display.name %in% gene.use,]
	rowData(sce.tmp)$gene.div <- rowData(sce.tmp)$display.name %in% gene.div.tb[category %in% c("CD4.hi","CD8.hi"),][["geneSymbol"]]

	sce.tmp <- ssc.reduceDim(sce.tmp,pca.npc=5,
			     assay.name = "norm_exprs",method="pca",method.vgene="gene.div",
			     autoTSNE=F,dim.name="div.pca")

	
	vgene <- rowData(sce.tmp)[["gene.div"]]
        proj_data <- uwot::umap(as.matrix(BiocGenerics::t(assay(sce.tmp[vgene,], "norm_exprs"))),
				init = "spca", pca = 5, n_threads = ncores)
	reducedDim(sce.tmp,"div.umap") <- proj_data

        proj_data <- uwot::umap(reducedDim(sce.tmp,"div.pca")[,1:5],
				metric="cosine",
				init = "spca", pca = NULL, n_threads = ncores)
	reducedDim(sce.tmp,"div.umap") <- proj_data


	l_ply(c("div.pca","div.umap"),function(rd) {
	    ssc.plot.tsne(sce.tmp, columns = "stype",
				  reduced.name = rd,
				  par.geneOnTSNE = list(scales = "free",pt.order = "random"),
				  colSet=list(),size=0.01,label=NULL,
				  base_aspect_ratio = 1.50)
	    ggsave(file = sprintf("%s.%s.groupBy.%s.tttt.png",out.prefix,rd,"stype"),width=5.0,height=4)
	})


	#sce.tmp <- ssc.reduceDim(sce.tmp,pca.npc=5,
	#		     assay.name = "norm_exprs",method="umap",method.vgene="gene.div",
	#		     autoTSNE=F,dim.name="div.umap",ncore=ncores)

	#### using seurat
	rownames(sce.tmp) <- rowData(sce.tmp)$seu.id
	seu <- as.Seurat(sce.tmp,counts="counts",data = "norm_exprs")

    }

    ###
    plot.all <- function(rd="umap",resolution.vec=seq(0.5,2.4,0.1))
    {
	p <- DimPlot(seu,reduction=rd, group.by = "loc")
	ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"loc"),width=5,height=4)

	if(!is.null(seu@meta.data[["stype"]])){
		p <- DimPlot(seu,reduction=rd, group.by = "stype")
		ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"stype"),width=5,height=4)
	}
	
	p <- DimPlot(seu,reduction=rd, group.by = "cancerType")
	ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"cancerType"),width=5,height=4)

	if(length(unique(seu$libraryID))>20){
		p <- DimPlot(seu,reduction=rd,label=F,label.size=2,
			     group.by = "libraryID") + NoLegend()
		ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"libraryID"),width=4.5,height=4)
	}else{
		p <- DimPlot(seu,reduction=rd,label=T,label.size=2, group.by = "libraryID")
		ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"libraryID"),width=6.5,height=4)
	}

	if(length(unique(seu$patient))>20){
		p <- DimPlot(seu,reduction=rd, group.by = "patient") + NoLegend()
		ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"patient"),width=4.5,height=4)
	}else{
		p <- DimPlot(seu,reduction=rd, group.by = "patient")
		ggsave(sprintf("%s.%s.%s.png",out.prefix,rd,"patient"),width=6.5,height=4)
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

	ssc.plot.tsne(sce, columns = "nFeature_RNA", 
				  reduced.name = sprintf("seurat.%s",rd),
				  colSet=colSet,size=0.03,
				  par.geneOnTSNE = list(scales = "free",pt.order = "random"),
				  vector.friendly=T,
				  out.prefix = sprintf("%s.seurat.%s.groupBy.%s",out.prefix,rd,"nFeature_RNA"),
				  base_aspect_ratio = 1.30)
    }

    
    {
	loginfo(sprintf("running Seurat pipeline ..."))

	VariableFeatures(seu) <- gene.div.tb[category %in% c("CD4.hi","CD8.hi"),][["geneSymbol"]]
	print((seu@assays$RNA@var.features))

	adj.cov <- NULL
	loginfo(sprintf("Scale ..."))
	
	seu <- ScaleData(object = seu,do.scale=F,vars.to.regress = adj.cov)

	#### PCA
	seu <- RunPCA(object = seu)
	print(x = seu[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
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
	seu <- RunUMAP(object = seu, reduction = "pca",dims = 1:opt.npc)
	toc()

	p <- DimPlot(seu,reduction="umap", group.by = "stype")
	ggsave(sprintf("%s.%s.%s.ttt.png",out.prefix,"umap","stype"),width=5,height=4)

	######### tSNE 
	if("tsne" %in% plot.rd){
	    tic("RunTSNE...")
	    seu <- RunTSNE(object = seu, reduction = "pca",dims = 1:opt.npc)
	    toc()
	}
	
	#### clustring
	seu <- FindNeighbors(object = seu, reduction = "pca", dims = 1:opt.npc)
	resolution.vec <- seq(0.1,2.4,0.1)
	seu <- FindClusters(object = seu,resolution = c(resolution.vec,res.addition))
	#seu <- FindClusters(object = seu,resolution = c(resolution.vec,5,10,25))

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
       

    return(list("seu"=seu,"sce"=sce))
}


