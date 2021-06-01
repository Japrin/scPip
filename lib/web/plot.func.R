
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



makeFig.ExampleGeneBarplot <- function(gene.to.plot,
                                       mcls.plot,
                                       gene.long.tb,
                                       gene.desc.top,
                                       gene.long.collapsed.tb=NULL,
                                       mod.sort=2,th.dprime=0.15,
                                       colSet.cancerType=NULL,ncol=1)
{
    prepare.data.for.plot <- function(dat.long,gene.desc.top,mcls,a.gene,mod.sort=3)
    {
	if(mod.sort==1 || mod.sort==2){
	    mapping.dataset.tb <- unique(data.table(dataset=gene.long.tb$dataset,dataset.old=gene.long.tb$dataset.old))
	    mapping.dataset.tb[dataset.old=="CHOL.YaoHe10X",dataset:="CHOL.QimingZhang2019.10X"]
	    mapping.dataset.tb[dataset.old=="CHOL.LichunMa2019",dataset:="CHOL.LichunMa2019"]
	    mapping.dataset.vec <- structure(mapping.dataset.tb$dataset,names=mapping.dataset.tb$dataset.old)
	}
	
	dat.plot <- dat.long[meta.cluster==mcls & geneID==a.gene,]
	dat.meta <- gene.desc.top[meta.cluster==mcls & geneSymbol==a.gene,]

	if(mod.sort==1 || mod.sort==2){
	    dat.plot[,dataset:=gsub("\\.CD[48].+$","",aid)]
	    dat.plot[,dataset:=mapping.dataset.vec[dataset]]
	    dat.plot <- dat.plot[,c("aid","geneID","P.Value","adj.P.Val","sig",
							    "dprime","vardprime","dataset",
							    "meta.cluster","cancerType"),with=F]
	    if(mod.sort==1){
		    tmp.order <- dat.plot[order(-meta.cluster,dprime),]
	    }else if(mod.sort==2){
		    tmp.order <- dat.plot[order(-meta.cluster,-dprime),]
		    dat.plot[,cancerType:=factor(cancerType,levels=rev(unique(tmp.order$cancerType)))]
		    tmp.order <- dat.plot[order(-meta.cluster,cancerType,dprime),]
	    }
	    dat.plot[,dataset:=factor(dataset,levels=c(unique(tmp.order$dataset),"combined"))]
	    dat.plot[,dataset.cate:="perStudies"]
	    dat.meta.plot <- data.table(aid="combined",geneID=a.gene,
							    P.Value=dat.meta$comb.p,
							    adj.P.Val=dat.meta$comb.padj,
							    sig=dat.meta$sig,
							    dprime=dat.meta$comb.ES,
							    vardprime=dat.meta$comb.ES.sd^2,
							    dataset="combined",
							    meta.cluster=mcls,
							    cancerType="panC",
							    dataset.cate="combined")
	}else if(mod.sort==3){
		dat.plot <- dat.plot[,c("geneID","cancerType","dprime","vardprime",
								"P.Value","adj.P.Val","sig",
								"meta.cluster"),with=F]
		tmp.order <- dat.plot[order(-meta.cluster,sig,dprime),]
		dat.plot[,cancerType:=factor(cancerType,levels=c(unique(tmp.order$cancerType),"panC"))]
		dat.meta.plot <- data.table(geneID=a.gene,
									cancerType="panC",
									dprime=dat.meta$comb.ES,
									vardprime=dat.meta$comb.ES.sd^2,
									P.Value=dat.meta$comb.p,
									adj.P.Val=dat.meta$comb.padj,
									sig=dat.meta$sig,
									meta.cluster=mcls)
	}
	dat.plot <- rbind(dat.plot,dat.meta.plot)
	dat.plot[,lower:=dprime-sqrt(vardprime)]
	dat.plot[,upper:=dprime+sqrt(vardprime)]
	dat.plot[,sig:=as.logical(sig)]
	#### for display purpose
	dat.plot[sig==F,adj.P.Val:=1]
	dat.plot[sig==T & adj.P.Val > 0.01,adj.P.Val:=0.01]
	return(dat.plot)
    }

    

    ### gene by gene
    dat.fig.list <- llply(gene.to.plot,function(a.gene){
	    if(mod.sort==3){
		    dat.plot <- prepare.data.for.plot(gene.long.collapsed.tb,gene.desc.top,mcls.plot,a.gene,mod.sort=mod.sort)
	    }else{
		    dat.plot <- prepare.data.for.plot(gene.long.tb,gene.desc.top,mcls.plot,a.gene,mod.sort=mod.sort)
	    }
	    ###dat.plot.debug <<- dat.plot
	    vline.x <- NULL
	    if(mod.sort==3){
		    vline.x <- dat.plot[sig==T,min(as.integer(cancerType))]
	    }
    
	    
	    ###print(vline.x)
	    p <- ggplot(dat.plot, aes_string(if(mod.sort==3) "cancerType" else "dataset", "dprime")) + 
		    geom_col(aes(fill=adj.P.Val),color=NA,position = "dodge2",width=0.8) +
		    scale_fill_distiller(palette = "Blues",breaks=c(0,0.002,0.004,0.006,0.008,0.01),
							     limits=c(0,0.01),na.value="lightgray",
							     guide = guide_colorbar(label.theme = element_text(angle = 45, hjust = 1),
													    direction = "horizontal")) +
		    ####geom_col(aes(fill=sig),color=NA,position = "dodge2") +
		    ####scale_fill_manual(values=c("TRUE"="#F8766D","FALSE"="#00BFC4")) +
		    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
		    geom_hline(yintercept=th.dprime,linetype="dashed",alpha=0.8,color="lightgray") +
		    theme_classic() +
		    ##labs(x="dataset",y="Effect Size",title=a.gene) +
		    labs(x=if(mod.sort==3) "cancerType" else "dataset",y="Effect Size",title=a.gene) +
		    theme(plot.title = element_text(hjust = 0.5))
		    ##theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5),
		    #theme(axis.text.x=element_custom(fill=colSet$cancerType[as.character(levels(dat.plot$cancerType))]),
	    if(mod.sort==3){
		    p <- p +theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5))
	    }else{
		    mapping.d2c.tb <- unique(dat.plot[,c("dataset","cancerType"),with=F])
		    mapping.d2c.tb[,dataset:=as.character(dataset)]
		    mapping.d2c.tb[,cancerType:=as.character(cancerType)]
		    mapping.d2c <- structure(mapping.d2c.tb$cancerType,names=mapping.d2c.tb$dataset)
		    par.fill <- colSet.cancerType[mapping.d2c[as.character(intersect(levels(dat.plot$dataset),unique(dat.plot$dataset)))]]
		    #print(par.fill)
		    #print(str(par.fill))
		    p <- p + theme(axis.text.x=element_custom(size=10, fill=par.fill))
		    #p <- p + theme(axis.text.x=element_text(size=10, angle = 90))
		    
	    }
	    if(!is.null(vline.x) && is.finite(vline.x)){
		    p <- p + geom_vline(xintercept=vline.x-0.5,linetype="dashed",color="red",alpha=0.8)
	    }
	    # ggsave(sprintf("%s.cmp.gene.example.%s.pdf",out.prefix,a.gene),
	    # 	   width=if(mod.sort==3) 6.5 else 7,height=if(mod.sort==3) 2.5 else 4.0)
	    return(p)
			      },.parallel=F)
    #return(dat.fig.list)
    cowplot::plot_grid(cowplot::plot_grid(plotlist=llply(dat.fig.list,
							 function(x){ x + theme(legend.position = "none") }),
					  ncol = ncol),
		       cowplot::get_legend(dat.fig.list[[1]]),
		       ncol = 1,rel_heights = c(length(dat.fig.list),0.5))
		       
    #saveRDS(dat.fig.list,file=sprintf("%s.cmp.gene.example.fig.rds",out.prefix))

}

