#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-a", "--aFile", type="character", required=TRUE, help="input seu file list")
parser$add_argument("-b", "--bFile", type="character", required=TRUE, help="input sce file list")
parser$add_argument("-c", "--stype", type="character", help="only analyze stype specified (default all)")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
#parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id")
parser$add_argument("-d", "--npc", type="integer",default=15L, help="[default %(default)s]")
parser$add_argument("-n", "--ncores", type="integer",default=16L, help="[default %(default)s]")
parser$add_argument("-m", "--measurement",type="character",default="counts",help="[default %(default)s]")
parser$add_argument("-r", "--resolution",type="character",default="2",help="[default %(default)s]")
parser$add_argument("-t", "--scTransform",action="store_true",default=FALSE,help="[default %(default)s]")
parser$add_argument("-g", "--deg",action="store_true",default=FALSE,help="[default %(default)s]")
parser$add_argument("-f", "--filterout",type="character",help="filterout cells")
parser$add_argument("-p", "--platform",type="character",required=TRUE,help="platform such as 10X, SmartSeq2")
args <- parser$parse_args()
print(args)

############## tune parametrs  ########

seu.file <- args$aFile
sce.file <- args$bFile
out.prefix <- args$outPrefix
opt.npc <- args$npc
opt.ncores <- args$ncores
opt.measurement <- args$measurement
opt.platform <- args$platform
opt.stype <- args$stype
opt.resolution <- args$resolution
opt.scTransform <- args$scTransform
opt.doDEG <- args$deg

dir.create(dirname(out.prefix),F,T)

saveRDS(args,file=sprintf("%s.args.rds",out.prefix))
###args <- readRDS(file=sprintf("%s.args.rds",out.prefix))

############## tune parametrs  ########
library("sscClust")
library("Seurat")
library("tictoc")
library("plyr")
library("dplyr")
library("tibble")
library("doParallel")
library("sscClust")
library("Matrix")
library("data.table")
library("R.utils")
library("gplots")
library("ggplot2")
library("ggpubr")
library("cowplot")
library("limma")
library("reticulate")
library("this.path")
options(stringsAsFactors = FALSE)

g.src.dir <- sprintf("%s/../",this.dir())
source(sprintf("%s/lib/sc.utils.R",g.src.dir))
gene.exclude.file <- sprintf("%s/data/geneSet/exclude.gene.misc.misc.v3.RData",g.src.dir)
#print(gene.exclude.file)

######################

env.misc <- loadToEnv(gene.exclude.file)
env.misc$all.gene.ignore.df %>% head

seu <- NULL
sce <- NULL

if(seu.file!="-" && file.exists(seu.file)){
    if(grepl("\\.rds$",seu.file)){
        seu <- readRDS(seu.file)
    }else{
        env.a <- loadToEnv(seu.file)
        obj.name.a <- names(env.a)[1]
        seu <- env.a[[obj.name.a]]
        rm(env.a)
    }
}
if(sce.file!="-" && file.exists(sce.file)){
    if(grepl("\\.rds$",sce.file)){
        sce <- readRDS(sce.file)
    }else{
        env.b <- loadToEnv(sce.file)
        obj.name.b <- names(env.b)[1]
        sce <- env.b[[obj.name.b]]
        rm(env.b)
    }
}

if(!is.null(seu) && !is.null(opt.stype)){
	seu <- seu[,seu$stype==opt.stype]
}

if(!is.null(seu) && "percent.mito" %in% colnames(seu[[]])){
	if(max(seu$percent.mito) <1){
		seu  <- subset(seu, subset = percent.mito<0.1)
	}else{
		seu  <- subset(seu, subset = percent.mito<10)
	}
}

if(!is.null(seu) && !is.null(args$filterout)){
	if(!file.exists(args$filterout)){
		col.filter <- unlist(strsplit(args$filterout,":"))[1]
		col.value <- unlist(strsplit(unlist(strsplit(args$filterout,":"))[2],","))
		cat(sprintf("filter cells with %s in c(%s)\n",col.filter,paste(col.value,collapse=",")))
		f.cell <- seu[[]][,col.filter] %in% col.value
		print(summary(f.cell))
		seu <- seu[,!f.cell]
	}
}

if(!is.null(seu) && !is.null(sce)){
    f.cell <- intersect(colnames(seu),colnames(sce))
    seu <- seu[,f.cell]
    sce <- sce[,colnames(seu)]
}

tic("run.Seurat3")
obj.list <- run.Seurat3(seu,sce,out.prefix,
                        gene.exclude.df=env.misc$all.gene.ignore.df,
                        n.top=1500,
						measurement=opt.measurement,platform=opt.platform,
                        use.sctransform=opt.scTransform,
                        do.deg=opt.doDEG,
                        ###do.adj=T,do.scale=F,
                        plot.rd=c("umap"),
						opt.res=opt.resolution,
                        opt.npc=opt.npc,ncores=opt.ncores)
toc()


