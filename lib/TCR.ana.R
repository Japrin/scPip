
#### visualization of clonotype
vis.clonotype <- function(out.prefix,
			  clone.x,sce.list,cinfo.clone.tb,
			  gene.sig.list,
			  sig.pretty=seq(-3,3,1),
			  gene.use=c("CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4",
				     "IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
			  mcls.plot=c("CD8.c15.ISG.IFIT1",
				      "CD8.c12.Tex.CXCL13"),
			  mcls.sig=NULL,
			  pdf.width=16.5,pdf.height=6,
			  prob.mat=NULL,sortByProb=F)
{
    #####x <- "UCEC-P20190312_C001978:59"
    ##clone.perCell.tb <- cinfo.clone.tb[clone.id==clone.x,]
    if(is.null(mcls.plot)){
	clone.perCell.tb <- cinfo.clone.tb[clone.id==clone.x,]
	mcls.plot <- as.character(sort(unique(clone.perCell.tb$meta.cluster)))
    }else{
	clone.perCell.tb <- cinfo.clone.tb[clone.id==clone.x & meta.cluster %in% mcls.plot,]
    }
    data.id <- clone.perCell.tb$dataset.old[1]
    cat(sprintf("Processing clonotype %s (%s)\n",clone.x,data.id))

    if(is.null(sce.list[[data.id]])){
	sce.list[[data.id]] <- readRDS(sprintf("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/zhangLab.10X/inte.metaClust.20200111/OUT.int.S35.CD8.k10.10.v1.occ1/sce/%s.sce.rds",data.id))
    }

    #### mean   
    ##f.gene.zz <- match(unique(gene.sig.tb$geneID),rowData(sce.list[[data.id]])$display.name)
    f.gene.zz <- match(unique(c(unname(do.call(c,gene.sig.list)),gene.use)),
                       rowData(sce.list[[data.id]])$display.name)
    f.gene.zz <- f.gene.zz[!is.na(f.gene.zz)]    

    sce.z <- ssc.assay.zscore(obj=sce.list[[data.id]][f.gene.zz,],
						      assay.name="norm_exprs")
    sce.plot.avg <- ssc.average.gene(sce.z,assay.name="norm_exprs.z",gene=gene.sig.list)
    sce.plot.avg <- sce.plot.avg[,clone.perCell.tb$Cell_Name]
    print("all(clone.perCell.tb$Cell_Name==colnames(sce.plot.avg))")
    print(all(clone.perCell.tb$Cell_Name==colnames(sce.plot.avg)))
    colData(sce.plot.avg) <- DataFrame(clone.perCell.tb[,c("Cell_Name","meta.cluster","clone.id"),
				       with=F])
    sce.plot.avg$meta.cluster <- as.character(sce.plot.avg$meta.cluster)
    colnames(sce.plot.avg) <- sce.plot.avg$Cell_Name
    sce.plot.avg$cluster.name <- fetchMetaClusterID2CusterFullName()[sce.plot.avg$meta.cluster]

    #sig.pretty <- seq(-1,1,0.5)
    #####
    if(is.null(mcls.sig)){
	mcls.sig <- mcls.plot
    }
    sig.plot.vec <- fetchMetaClusterID2CusterFullName()[mcls.sig]
    #print(str(sig.plot.vec))
    #print(sce.plot.avg)
    
    sce.plot <- sce.z
    sce.plot <- sce.plot[,clone.perCell.tb$Cell_Name]
    print("all(clone.perCell.tb$Cell_Name==colnames(sce.plot))")
    print(all(clone.perCell.tb$Cell_Name==colnames(sce.plot)))
    colData(sce.plot) <- cbind(colData(sce.plot),clone.perCell.tb[,c("meta.cluster","clone.id"),with=F])
    sce.plot$meta.cluster <- as.character(sce.plot$meta.cluster)
    sce.plot$cluster.name <- fetchMetaClusterID2CusterFullName()[sce.plot$meta.cluster]

    print("all(colnames(sce.plot)==colnames(sce.plot.avg))")
    print(all(colnames(sce.plot)==colnames(sce.plot.avg)))
    for(a.sig in names(sig.plot.vec)){
	#print(sig.plot.vec[a.sig])
	a.dat <- assay(sce.plot.avg,"exprs")[sig.plot.vec[a.sig],]
	#a.dat <- assay(sce.plot.avg,"exprs")[a.sig,]
	a.dat[a.dat < sig.pretty[1] ] <- sig.pretty[1]
	a.dat[a.dat > sig.pretty[length(sig.pretty)] ] <- sig.pretty[length(sig.pretty)]
	colData(sce.plot)[[sprintf("sig.%s",a.sig)]] <- a.dat
	### add cluster assignment probability
	if(!is.null(prob.mat)){
	    colData(sce.plot)[[sprintf("prob.%s",a.sig)]] <- prob.mat[colnames(sce.plot),a.sig]
	}
    }

    {
	dat.for.sort <- cbind(cellID=colnames(sce.plot.avg),
			      meta.cluster=sce.plot.avg$meta.cluster,
			      as.data.table(t(assay(sce.plot.avg[sig.plot.vec,]))))
	if(!is.null(prob.mat)){
	    .tmp.colData <- colData(sce.plot)
	    dat.for.sort <- as.data.table(cbind(dat.for.sort,.tmp.colData[,grepl("prob.",colnames(.tmp.colData))]))
	}
	#print(head(dat.for.sort))
	colnames(dat.for.sort)[3:(2+length(sig.plot.vec))] <- names(sig.plot.vec)

	dat.tmp <- dat.for.sort[!meta.cluster %in% names(sig.plot.vec),]
	for(mcls.xx in mcls.sig)
	{
	    dat.tmp <- rbind(dat.tmp,setorderv(dat.for.sort[meta.cluster==mcls.xx,],mcls.xx))
	}
	if(!is.null(prob.mat) & sortByProb){
	    dat.tmp <- dat.for.sort[!meta.cluster %in% names(sig.plot.vec),]
	    for(mcls.xx in mcls.sig)
	    {
		dat.tmp <- rbind(dat.tmp,setorderv(dat.for.sort[meta.cluster==mcls.xx,],sprintf("prob.%s",mcls.xx)))
	    }
	}
	dat.for.sort <- dat.tmp
    }


    #gene.use <- intersect(c("HAVCR2","CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "CTLA4","PDCD1","IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
#			  rowData(sce.plot)$display.name)

    if(is.null(gene.use)){
	gene.use <- unique(unname(do.call(c,gene.sig.list)))
    }
    print(str(gene.use))

    gene.use <- intersect(gene.use,
			  rowData(sce.plot)$display.name)
    ###print(str(gene.use))
    f.gene.use <- match(gene.use,rowData(sce.plot)$display.name)
    sce.plot.comb <- sce.plot[f.gene.use,
			      colnames(ssc.assay.hclust(sce.plot.avg[(sig.plot.vec),],
							"exprs",order.col=T,
							clustering.distance="cosine",
							clustering.method="ward.D2"))]
    annotation_legend_param <- llply(mcls.sig,function(x){ list(at=sig.pretty) })
    names(annotation_legend_param) <- sprintf("sig.%s",mcls.sig)
    annotation_legend_param_prob <- llply(mcls.sig,function(x){ list(at=seq(0,1,0.2),color_bar = "continuous", legend_direction = "horizontal",legend_width = unit(4, "cm"), legend_height = unit(2,"cm")) })
    names(annotation_legend_param_prob) <- sprintf("prob.%s",mcls.sig)
    colSet.plot <- g.colSet
    if(!is.null(prob.mat)){
	annotation_legend_param <- c(annotation_legend_param,annotation_legend_param_prob)
	colSet.prob <- llply(mcls.sig,function(x){ colorRamp2(seq(0,1,length=9),sscVis:::getColorPaletteFromNameContinuous("Greys"),space="LAB") })
	names(colSet.prob) <- sprintf("prob.%s",mcls.sig)
	colSet.plot <- c(colSet.plot,colSet.prob)
    }
    #print(sprintf("%s.example.%s.comb",out.prefix,gsub(":","_",clone.x,perl=T))) 
    #print((sce.plot.comb[,dat.for.sort$cellID]))
    ssc.plot.heatmap(sce.plot.comb[,dat.for.sort$cellID],assay.name="norm_exprs.z",
		     out.prefix=sprintf("%s.example.%s.comb",out.prefix,gsub(":","_",clone.x,perl=T)),
		     columns=if(is.null(prob.mat)) c("cluster.name",sprintf("sig.%s",names(sig.plot.vec))) else c("cluster.name",sprintf("sig.%s",names(sig.plot.vec)),sprintf("prob.%s",names(sig.plot.vec))),
		     columns.order=c("cluster.name"),
		     clustering.distance="cosine",
		     clustering.method="ward.D2",
		     colSet=colSet.plot,do.scale=F,
		     do.clustering.row=F,
		     do.clustering.col=F,
		     par.heatmap=list(cex.row=1.2,cex.column=0),
		     #palette.name="RdBu",
		     palette.name="cividis",
		     palette.ann.numeric="RdBu",
		     Y.level.ann.numeric=sig.pretty,
		     annotation_legend_param=annotation_legend_param,
		     z.lo = -3.0, z.hi = 3.0, z.step = 1,
		     ###pdf.width = 16.5, pdf.height = 6,
		     pdf.width = pdf.width, pdf.height = pdf.height,
		     mytitle = clone.x)

}

ana.clonotypeAcrossMcls.moreThanTwo <- function(object,
			    in.dat,
			    out.prefix,
			    aid="PathTrm",
			    lim.LLR=10,
			    par.barplot=list("pdf.width.byPatientF"=7,"pdf.height.byPatientF"=7.5,
					     "pdf.width.byPatientT"=7,"pdf.height.byPatientT"=8),
			    par.ggbarplot=NULL,
			    mcls.moi=c("CD8.c02.Tm.IL7R","CD8.c10.Trm.ZNF683","CD8.c12.Tex.CXCL13"))
{

    dir.create(dirname(out.prefix),F,T)

    dat.block <- object@clonotype.dist.cluster[,mcls.moi]
    ### clone.info.flt.tb
    dat.block.clone.index <- Startrac:::mrow.entropy(dat.block)
    dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
    clone.info.index.tb <- data.table(cloneID=names(dat.block.clone.index),
				      pIndex.tran=dat.block.clone.index)
    f.crossAll <- rowSums(dat.block > 0)==ncol(dat.block)
    clone.info.index.tb <- clone.info.index.tb[pIndex.tran > 0 & f.crossAll,]
    clone.info.patient.tb <- in.dat[,.N,by=c("cancerType","dataset","dataset.old","patient","cloneID")]
    clone.info.flt.tb <- merge(clone.info.index.tb,clone.info.patient.tb,by="cloneID")

    clone.LLR.tb <- ldply(seq_len(nrow(clone.info.flt.tb)),function(i){
			calCloneLLR(in.dat[majorCluster %in% mcls.moi,],
			    prob.mat=prob.svm.mat.list[[ clone.info.flt.tb$dataset.old[i] ]],
			    cloneID.x=clone.info.flt.tb$cloneID[i],verbose=F)
			})
    clone.info.flt.tb <- merge(clone.info.flt.tb,clone.LLR.tb,by="cloneID")
    clone.info.flt.tb <- clone.info.flt.tb[order(-LLR,pIndex.tran),]
    clone.info.flt.tb <- clone.info.flt.tb[LLR > 1 & G.best==G.obs,]

    if(nrow(clone.info.flt.tb) < 0){
	ret.list <- list("dat.block.flt"=NULL,
			 "dat.block.flt.freq"=NULL,
			 "clone.info.flt.tb"=NULL)
	warning(sprintf("no clonotype found !"))
	return(ret.list)
    }

    ## dat.block
    dat.block.flt <- dat.block[clone.info.flt.tb$cloneID,,drop=F]
    dat.block.flt.freq <- sweep(dat.block.flt,1,rowSums(dat.block.flt),"/")


    do.plot.bar <- function(dat.block.flt.freq,clone.info.flt.tb,
			    out.prefix,
			    pdf.width=4,
			    pdf.height=9,
			    splitByPatient=F)
    {

	dat.plot.tb <- cbind(data.table(cloneID=rownames(dat.block.flt.freq)),
			     dat.block.flt.freq)
	dat.plot.tb <- merge(clone.info.flt.tb,dat.plot.tb)
	dat.plot.tb <- dat.plot.tb[order(LLR,-pIndex.tran),]
	dat.plot.melt <- melt(dat.plot.tb[,c("cloneID",mcls.moi,"patient","cancerType"),with=F],
			      variable.name = "meta.cluster", value.name = "freq",
			      id.vars=c("cloneID","patient","cancerType"))
	dat.plot.melt[,cloneID.fac:=factor(cloneID,levels=dat.plot.tb$cloneID)]
	#dat.plot.melt[,meta.cluster:=factor(meta.cluster,levels=sort(unique(as.character(meta.cluster)),
	#							     decreasing=T))]
	dat.plot.melt[,meta.cluster:=factor(meta.cluster,levels=rev(mcls.moi))]
	dat.plot.tb[,cloneID.fac:=factor(cloneID,levels=cloneID)]
	if(!is.null(lim.LLR)){
	    dat.plot.tb[LLR>lim.LLR,LLR:=lim.LLR]
	}

	### freq
	p1 <- do.call(ggbarplot,c(list(data=dat.plot.melt,x="cloneID.fac",y="freq",
				       color=NA, fill="meta.cluster"),
				  par.ggbarplot))
	p1 <- p1 +
		labs(x="",y="Frequency")+
		scale_fill_manual(values=g.colSet$meta.cluster) +
		theme(axis.text.y=element_blank(),legend.position="right",
		      axis.text.x=element_text(size=10),
		      strip.placement = "outside",
		      strip.background=element_blank(),
		      strip.text.y.left=element_text(angle=0)) +
		coord_flip()
	if(splitByPatient){
	    p1 <- p1 + facet_grid(patient~.,scales = "free", space = "free",switch="y")
	}

	### pIndex.tran
	p2 <- do.call(ggbarplot,c(list(data=dat.plot.tb,x="cloneID.fac",y="pIndex.tran",
				       color=NA, fill="steelblue"),
				  par.ggbarplot))
	p2 <- p2 + 
		labs(x="",y="Entropy")+
		geom_hline(yintercept=0.1,linetype="dashed") +
		theme(
		      axis.text.y=element_blank(),
		      axis.text.x=element_text(size=10),
		      legend.position="right",
		      #strip.placement = "outside",
		      strip.background=element_blank(),
		      strip.text.y=element_blank()
		      #strip.text.y.left=element_text(angle=0)
		      ) +
		coord_flip()
	if(splitByPatient){
	    #p2 <- p2 + theme(axis.text.y=element_blank())
	    p2 <- p2 + facet_grid(patient~.,scales = "free", space = "free",switch="y")
	}

	### LLR
	p3 <- do.call(ggbarplot,c(list(data=dat.plot.tb,x="cloneID.fac",y="LLR",
				       color=NA, fill="#82B446"),
				  par.ggbarplot))
	p3 <- p3 +
		labs(x="")+
		geom_hline(yintercept=1,linetype="dashed") +
		theme(
		      axis.text.y=element_blank(),
		      axis.text.x=element_text(size=10),
		      legend.position="right",
		      #strip.placement = "outside",
		      strip.background=element_blank(),
		      strip.text.y=element_blank()
		      #strip.text.y.left=element_text(angle=0)
		      ) +
		coord_flip()
	if(splitByPatient){
	    ##p3 <- p3 + theme(axis.text.y=element_blank())
	    p3 <- p3 + facet_grid(patient~.,scales = "free", space = "free",switch="y")
	}

	p00 <- cowplot::plot_grid(p1+theme(legend.position="none"),
				p2,
				p3,
				align="h",nrow=1,
				rel_widths=if(splitByPatient) c(1,0.4,0.4) else c(1,0.4,0.4)
				)
	pout <- cowplot::plot_grid(p00,get_legend(p1),ncol=1,rel_heights=c(0.7,0.3))
	ggsave(sprintf("%s.test.pdf",out.prefix),plot=pout,width=pdf.width,height=pdf.height)

	p00 <- cowplot::plot_grid(p1+theme(legend.position="none"),
				p3,
				align="h",nrow=1,
				rel_widths=if(splitByPatient) c(4.5,1.78) else c(4.5,1.78)
				)
	pout <- cowplot::plot_grid(p00,get_legend(p1),ncol=1,rel_heights=c(0.7,0.3))
	ggsave(sprintf("%s.pdf",out.prefix),plot=pout,width=6.28,height=pdf.height)
	cat(sprintf("----- width: %4.2f, height: %4.2f\n",pdf.width,pdf.height))
    }

    do.plot.bar(dat.block.flt.freq,
		clone.info.flt.tb,
		out.prefix=sprintf("%s.bar.freq.byPatientF",out.prefix),
		pdf.width=par.barplot$pdf.width.byPatientF,pdf.height=par.barplot$pdf.height.byPatientF,splitByPatient=F)
    do.plot.bar(dat.block.flt.freq,
		clone.info.flt.tb,
		out.prefix=sprintf("%s.bar.freq.byPatientT",out.prefix),
		pdf.width=par.barplot$pdf.width.byPatientT,pdf.height=par.barplot$pdf.height.byPatientT,splitByPatient=T)

    #### specific clonotype
    {
	clone.example.vec <- head(clone.info.flt.tb$cloneID,n=10)

	gene.sig.slim.list <- list("CD8.c02.Tm.IL7R"=c("IL7R","ZFP36L2","CXCR4","ZFP36","ANXA1","GPR183"),
				   "CD8.c10.Trm.ZNF683"=c("ZNF683","HOPX","CAPG"),
				   "CD8.c05.Tem.CXCR5"=c("GZMK", "CD74", "CST7"),
				   "CD8.c06.Tem.GZMK"=c("GZMK", "CD74", "CST7"),
				   "CD8.c11.Tex.PDCD1"=c("GZMK", "CD74", "CST7"),
				   "CD8.c15.ISG.IFIT1"=c("IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
				   "CD8.c12.Tex.CXCL13"=c("CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4"),
				   "CD8.c14.Tex.TCF7"=c("CD200","GNG4","CXCL13","TNFRSF4","TCF7"),
				   "CD8.c17.Tm.NME1"=c("BIRC5","RRM2","CDC20","NME1","CENPW")
				   )
	gene.sig.slim.list <- gene.sig.slim.list[mcls.moi]
	names(gene.sig.slim.list) <- fetchMetaClusterID2CusterFullName()[names(gene.sig.slim.list)]
	print(str(gene.sig.slim.list))

	l_ply(seq_along(clone.example.vec),function(i){
	    dataset.old.id <- clone.info.flt.tb[cloneID==clone.example.vec[i],][["dataset.old"]][1]
	    vis.clonotype(out.prefix=sprintf("%s.%s", out.prefix,"exampleClone"),
			  clone.x=clone.example.vec[i],
			  sce.list=sce.list,
			  cinfo.clone.tb=in.dat[stype=="CD8",],
			  ##gene.sig.list=gene.sig.slim.list, ### manually picked
			  gene.sig.list=gene.sig.list, ### top 50
			  sig.pretty=seq(-1,1,0.5),
			  gene.use=unique(do.call(c,unname(gene.sig.slim.list))),
			  prob.mat=prob.svm.mat.list[[dataset.old.id]],
			  mcls.plot=mcls.moi)
				 })
    }
    ret.list <- list("dat.block.flt"=dat.block.flt,
		"dat.block.flt.freq"=dat.block.flt.freq,
		"clone.info.flt.tb"=clone.info.flt.tb)
    saveRDS(ret.list,sprintf("%s.ret.list.rds",out.prefix))
    return(ret.list)
	
}





