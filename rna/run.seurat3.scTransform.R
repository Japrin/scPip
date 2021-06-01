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
parser$add_argument("-r", "--resolution",type="character",default="RNA_snn_res.1",help="[default %(default)s]")
parser$add_argument("-f", "--filterout",type="character",help="filterout cells")
parser$add_argument("-p", "--platform",type="character",required=TRUE,help="platform such as 10X, SmartSeq2")
args <- parser$parse_args()
print(args)

############## tune parametrs  ########
#seu.file <- "/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/merge.phase02/seurat/rmIFNStress/OUT.seurat.CD8.wMAIT/patientT.scaleF/panC_CD8.wMAIT.seu.forInte.RData"
#sce.file <- "/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/merge.phase02/seurat/rmIFNStress/OUT.seurat.CD8.wMAIT/patientT.scaleF/panC_CD8.wMAIT.sce.forInte.RData"
#out.prefix <- "/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/zhangLab.10X/A20190225/inte.bbknn.panC/OUT.byDataset/test/byDataset.CD8"
#opt.npc <- 15L
#opt.ncores <- 16L
#opt.measurement <- "counts"

seu.file <- args$aFile
sce.file <- args$bFile
out.prefix <- args$outPrefix
opt.npc <- args$npc
opt.ncores <- args$ncores
opt.measurement <- args$measurement
opt.platform <- args$platform
opt.stype <- args$stype
opt.resolution <- args$resolution

##gene.exclude.file <- "/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/geneSet/exclude/exclude.gene.misc.misc.RData"
gene.exclude.file <- "/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/geneSet/exclude/exclude.gene.misc.misc.v3.RData"


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

options(stringsAsFactors = FALSE)

#####source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/pipeline/lib/sc.func.seurat3.R")

######################

env.misc <- loadToEnv(gene.exclude.file)
env.misc$all.gene.ignore.df %>% head

if(grepl("\\.rds$",seu.file)){
	seu <- readRDS(seu.file)
	sce <- readRDS(sce.file)
}else{
	env.a <- loadToEnv(seu.file)
	env.b <- loadToEnv(sce.file)
	obj.name.a <- names(env.a)[1]
	obj.name.b <- names(env.b)[1]
	seu <- env.a[[obj.name.a]]
	sce <- env.b[[obj.name.b]]
	rm(env.a)
	rm(env.b)
}

if(!is.null(opt.stype)){
	seu <- seu[,seu$stype==opt.stype]
}

if("percent.mito" %in% colnames(seu[[]])){
	if(max(seu$percent.mito) <1){
		seu  <- subset(seu, subset = percent.mito<0.1)
	}else{
		seu  <- subset(seu, subset = percent.mito<10)
	}
}

if(!is.null(args$filterout)){
	if(!file.exists(args$filterout)){
		col.filter <- unlist(strsplit(args$filterout,":"))[1]
		col.value <- unlist(strsplit(unlist(strsplit(args$filterout,":"))[2],","))
		cat(sprintf("filter cells with %s in c(%s)\n",col.filter,paste(col.value,collapse=",")))
		f.cell <- seu[[]][,col.filter] %in% col.value
		print(summary(f.cell))
		seu <- seu[,!f.cell]
	}
}


f.cell <- intersect(colnames(seu),colnames(sce))
seu <- seu[,f.cell]
sce <- sce[,colnames(seu)]

tic("run.Seurat3")
obj.list <- run.Seurat3(seu,sce,out.prefix,gene.exclude.df=env.misc$all.gene.ignore.df,n.top=1500,
						measurement=opt.measurement,platform=opt.platform,use.sctransform=T,
						opt.res=opt.resolution, opt.npc=opt.npc,ncores=opt.ncores)
toc()


