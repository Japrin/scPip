#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-a", "--aFile", type="character", required=TRUE, help="input h5ad file")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
parser$add_argument("-q", "--specie", type="character", default="human", help="one of human, mouse [default %(default)s]")
parser$add_argument("-c", "--stype", type="character", help="only analyze stype specified (default all)")
parser$add_argument("-u", "--geneIDFile", type="character", help="gene id mapping file")
parser$add_argument("-d", "--npc", type="integer",default=15L, help="[default %(default)s]")
parser$add_argument("-n", "--ncores", type="integer",default=16L, help="[default %(default)s]")
parser$add_argument("-m", "--ntop", type="integer",default=1500L, help="top HVG genes[default %(default)s]")
parser$add_argument("-r", "--resolution",type="character",default="1",help="best resolution [default %(default)s]")
parser$add_argument("-s", "--resTest", type="character",default="seq(0.1,2.4,0.1)", help="resolutions to test [default %(default)s]")
parser$add_argument("-y", "--harmony",action="store_true",default=FALSE,help="[default %(default)s]")
parser$add_argument("-i", "--integration",type="character",help="[default %(default)s]")
parser$add_argument("-g", "--deg",action="store_true",default=FALSE,help="[default %(default)s]")
parser$add_argument("-w", "--ncellDEG",type="integer",default=1500,
                    help="number of cells to downsample to for each group. used in DEG analysis. [default %(default)s]")
#parser$add_argument("-s", "--scale",action="store_true",default=FALSE,help="[default %(default)s]")
parser$add_argument("-j", "--corVar", type="character", default="S_score,G2M_score,DIG.Score,percent.mito",
		    help="subset of S_score,G2M_score,DIG.Score,percent.mito, or NULL. [default %(default)s]")
parser$add_argument("-f", "--filterout",type="character",help="Format is COLUMN_ID:COLUMN_VAL_1,COLUMN_VAL_2,COLUMN_VAL_3. Filter out cells with COLUMN_ID in one of COLUMN_VAL_1, COLUMN_VAL_2, and COLUMN_VAL_3.")
parser$add_argument("-k", "--keep",type="character",help="Format is COLUMN_ID:COLUMN_VAL_1,COLUMN_VAL_2,COLUMN_VAL_3. Keep only cells with COLUMN_ID in one of COLUMN_VAL_1, COLUMN_VAL_2, and COLUMN_VAL_3.")
parser$add_argument("-z", "--markerFile",type="character",help="file define cell type marker. used for --removeContamination. [default %(default)s]")
parser$add_argument("-x", "--removeContamination",type="character",help="comma separated string indicates subset of predefined signature (plasmaB, caf, epi, T, cd8, and mac) will be calculated and cells with high signature scores will be removed. For example, use plasmaB:0.75,caf:0.75,epi:0.75,T:0.25 for myeloid cell analysis. The numer after colon is the threshold")
args <- parser$parse_args()
print(args)

############## tune parametrs  ########

adata.file <- args$aFile
out.prefix <- args$outPrefix
###opt.npc <- args$npc
opt.npc <- as.integer(eval(parse(text = args$npc)))
opt.ncores <- as.integer(args$ncores)
opt.ntop <- as.integer(args$ntop)
opt.resTest <-  eval(parse(text=args$resTest))
opt.stype <- args$stype
opt.resolution <- args$resolution
opt.harmony <- args$harmony
opt.integration <- args$integration
opt.doDEG <- args$deg
opt.ncell.deg <- args$ncellDEG
#opt.scale <- args$scale
opt.cor.var <- if(args$corVar=="") c("") else unlist(strsplit(args$corVar,",",perl=T))
opt.filterout <- args$filterout
opt.keep <- args$keep
opt.geneIDFile <- args$geneIDFile
opt.markerFile <- args$markerFile
opt.removeContamination <- args$removeContamination
opt.specie <- args$specie

dir.create(dirname(out.prefix),F,T)

saveRDS(args,file=sprintf("%s.args.rds",out.prefix))
###args <- readRDS(file=sprintf("%s.args.rds",out.prefix))

############## tune parametrs  ########
suppressMessages(library("sscVis"))
suppressMessages(library("anndata"))
#suppressMessages(library("scanpy"))
#suppressMessages(library("sscClust"))
#suppressMessages(library("Seurat"))
suppressMessages(library("tictoc"))
suppressMessages(library("plyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("tibble"))
suppressMessages(library("doParallel"))
suppressMessages(library("Matrix"))
suppressMessages(library("data.table"))
suppressMessages(library("R.utils"))
suppressMessages(library("gplots"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("cowplot"))
suppressMessages(library("limma"))
suppressMessages(library("reticulate"))
suppressMessages(library("scPip"))
suppressMessages(library("anndata"))
options(stringsAsFactors = FALSE)
sc <- import("scanpy")

dat.ext.dir <- system.file("extdata",package="scPip")
if(opt.specie=="human"){
    ###gene.exclude.file <- sprintf("%s/exclude.gene.misc.human.v3.RData",dat.ext.dir)
    gene.exclude.file <- sprintf("%s/exclude.gene.misc.human.v4.RData",dat.ext.dir)
}else if(opt.specie=="mouse"){
    ###gene.exclude.file <- sprintf("%s/exclude.gene.misc.mouse.v3.RData",dat.ext.dir)
    gene.exclude.file <- sprintf("%s/exclude.gene.misc.mouse.v4.RData",dat.ext.dir)
}else{
    gene.exclude.file <- NULL
}

if(!is.null(gene.exclude.file)){
    env.misc <- loadToEnv(gene.exclude.file)
    g.all.gene.ignore.df <- env.misc$all.gene.ignore.df
    g.all.gene.ignore.df %>% head
}else{
    g.all.gene.ignore.df <- NULL
}

######################

gene.mapping.table <- NULL
if(!is.null(opt.geneIDFile) && file.exists(opt.geneIDFile) && grepl("\\.rds$",opt.geneIDFile,perl=T)){
    gene.mapping.table <- readRDS(opt.geneIDFile)
}

adata <- anndata::read_h5ad(adata.file)

##### filter out doublets predictd by Scrublet
if(!is.null(adata) && "pred.th" %in% colnames(adata$obs) ){
    adata <- adata[adata$obs$pred.th==FALSE]
}

##### check & clean sce #####
if(!is.null(adata) && "percent.mito" %in% colnames(adata$obs)){
    if(all(is.na(adata$obs$percent.mito))){
        adata$obs$percent.mito <- NULL
    }
}

if(!is.null(adata) && ("libraryID" %in% colnames(adata$obs)) && all(is.na(adata$obs$libraryID))){
    adata$obs$libraryID <- "Unk"
}

##### filter out cells with hgih percent.mito
if(!is.null(adata) && "percent.mito" %in% colnames(adata$obs) ){
	if(max(adata$obs$percent.mito) <1){
        adata <- adata[adata$obs$percent.mito < 0.1]
	}else{
        adata <- adata[adata$obs$percent.mito < 10]
	}
}

if(!is.null(adata) && !is.null(opt.filterout)){
	if(!file.exists(opt.filterout)){
		col.filter <- unlist(strsplit(opt.filterout,":"))[1]
		col.value <- unlist(strsplit(unlist(strsplit(opt.filterout,":"))[2],","))
        if(col.filter %in% colnames(adata$obs)){
            cat(sprintf("filter cells with %s in c(%s)\n",col.filter,paste(col.value,collapse=",")))
            if(grepl("^(>|<|=)",col.value)){
                f.cell <- eval(parse(text=sprintf("adata$obs[,\'%s\'] %s",col.filter,col.value)))
            }else{
                f.cell <- adata$obs[,col.filter] %in% col.value
            }

            print(summary(f.cell))
            adata <- adata[!f.cell]
        }else{
            warning(sprintf("The meta-data doesnot contain %s\n",col.filter))
        }
	}
}

if(!is.null(adata) && !is.null(opt.keep)){
	if(!file.exists(opt.keep)){
		col.keep <- unlist(strsplit(opt.keep,":"))[1]
		col.value <- unlist(strsplit(unlist(strsplit(opt.keep,":"))[2],","))
        if(col.keep %in% colnames(adata$obs)){
            cat(sprintf("keep only cells with %s in c(%s)\n",col.keep,paste(col.value,collapse=",")))
            if(grepl("^(>|<|=)",col.value)){
                f.cell <- eval(parse(text=sprintf("adata$obs[,\'%s\'] %s",col.keep,col.value)))
            }else{
                f.cell <- adata$obs[,col.keep] %in% col.value
            }
            print(summary(f.cell))
            adata <- adata[f.cell]
        }else{
            warning(sprintf("The meta-data doesnot contain %s\n",col.keep))
        }
    }
}

if(!is.null(opt.removeContamination)) {
    ### opt.removeContamination <- "plasmaB:0.75,caf:0.75,epi:0.75,T:0.25"
    if(is.null(opt.markerFile)){
        g.cont.gene.list <- list("plasmaB"=c("JCHAIN"),
                                 "B"=c("CD79A","MS4A1"),
                                 "caf"=c("COL1A2", "COL1A1", "COL3A1","LUM"),
                                 "epi"=c("KRT18","KRT19","EPCAM"),
                                 "T"=c("CD3D","CD3G"),
                                 "cd8"=c("CD8A","CD8B"),
                                 "mac"=c("LYZ","C1QA","C1QB","CD68"))
    }else{
        marker.rmCont.tb <- fread(opt.markerFile,head=F)
        g.cont.gene.list <- llply(split(marker.rmCont.tb,by=c("V1")),function(x){ x$V2 })
    }

    sig.vec <- unlist(strsplit(opt.removeContamination,","))
    sig.name <- unname(sapply(sig.vec,function(x){ unlist(strsplit(x,":"))[1] }))
    sig.thre <- unname(sapply(sig.vec,function(x){ as.numeric(unlist(strsplit(x,":"))[2]) }))

    for(i in seq_along(sig.name)){
        loginfo(sprintf("calculate signature score of %s cells ...",sig.name[i]))
        adata <- fill.contamination(adata,out.prefix,
                                  g.name=sig.name[i],
                                  g.test=g.cont.gene.list[[sig.name[i]]],
                                  score.t=sig.thre[i],
                                  vis.v=c(0.25,0.5,0.75,1))
    }
   
    idx.sig.class <- intersect(sprintf("%s.class",sig.name),colnames(adata$obs))
    f.cont.mtx <- adata$obs[,idx.sig.class,drop=F]
    f.cont <- rowSums(f.cont.mtx) > 0
    loginfo(sprintf("A total number of potential contamination: %d\n",sum(f.cont)))
    print(colSums(f.cont.mtx==T))
    adata <- adata[!f.cont]

}



tic("run.scanpy")
obj.list <- run.scanpy(adata,out.prefix,
                        gene.exclude.df=g.all.gene.ignore.df,
                        n.top=opt.ntop,
						opt.res=opt.resolution,
                        plot.rd=c("umap"),
                        opt.npc=opt.npc,ncores=opt.ncores,
                        res.test=opt.resTest,
			            cor.var=opt.cor.var,
                        use.harmony=opt.harmony,
                        method.integration=opt.integration,
                        specie=opt.specie,
                        do.deg=opt.doDEG,
                        ncell.deg=opt.ncell.deg,
                        ###do.adj=T,do.scale=F,
                        gene.mapping.table=gene.mapping.table)
toc()


