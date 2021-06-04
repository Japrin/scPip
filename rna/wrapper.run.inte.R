#!/usr/bin/env Rscript

####suppressPackageStartupMessages(library("argparse"))
####
####parser <- ArgumentParser()
####parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="input seu file list")
####parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
######parser$add_argument("-d", "--npc", type="integer",default=15L, help="[default %(default)s]")
####args <- parser$parse_args()
####print(args)
####
####exp.list.file <- args$inFile
####out.prefix <- args$outPrefix
#####gene.file <- args$geneFile
####
####dir.create(dirname(out.prefix),F,T)
####saveRDS(args,file=sprintf("%s.args.rds",out.prefix))
######args <- readRDS(sprintf("%s.args.rds",out.prefix))

####
exp.list.file <- "/workspace/zhengliangtao/work/proj_BRCA/ana/inte.metaClust/list/obj.T.list.r1.list"
out.prefix <- "/workspace/zhengliangtao/work/proj_BRCA/ana/inte.metaClust/OUT.int.S2.T/int.S2.T"
dir.create(dirname(out.prefix),F,T)

g.src.dir <- sprintf("%s/../",this.dir())
source(sprintf("%s/lib/inte.comb.miniClust.lib.R",g.src.dir))
#source("/workspace/zhengliangtao/02.pipeline/pipeline/lib/inte.comb.miniClust.lib.R")

ncores <- 12
gene.exclude.file <- sprintf("%s/data/geneSet/exclude.gene.misc.misc.v3.RData",g.src.dir)
#gene.exclude.file <- "/workspace/zhengliangtao/02.pipeline/pipeline/data/geneSet/exclude.gene.misc.misc.v3.RData"

options(stringsAsFactors = FALSE)

#env.misc <- loadToEnv(gene.exclude.file)

exp.list.table <- fread(cmd=sprintf("awk '!/^#/' %s",exp.list.file))

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = ncores)

tic("run.inte.metaClust")
ret.list <- run.inte.metaClust(exp.list.table,out.prefix,gene.exclude.file,ncores=12,npc=15,res.hi=50,TH.gene.occ=0.85)
toc()

seu.merged <- ret.list[["seu.merged"]]
sce.merged <- ret.list[["sce.merged"]]
meta.tb <- ret.list[["meta.tb"]]

#seu.merged <- readRDS(file=sprintf("%s.seu.merged.rds",out.prefix))
#sce.merged <- readRDS(file=sprintf("%s.sce.merged.rds",out.prefix))
#meta.tb <- readRDS(file=sprintf("%s.meta.tb.rds",out.prefix))

#seu.merged$dataset %>% table
#seu.merged$dataset.tech %>% table


dataOnRDPlot(seu.merged,sce.merged,
			 sprintf("%s/%s/%s",dirname(out.prefix),"umap.algorithm1",basename(out.prefix)),
			 rd="umap",graph.name="RNA_pca_snn")
dataOnRDPlot(seu.merged,sce.merged,
			 sprintf("%s/%s/%s",dirname(out.prefix),"harmony.umap.algorithm1",basename(out.prefix)),
			 rd="harmony.umap")
#dataOnRDPlot(seu.merged,sce.merged,
#			 sprintf("%s/%s/%s",dirname(out.prefix),"tsne.Rtsne.algorithm1",basename(out.prefix)),
#			 rd="tsne.Rtsne",graph.name="RNA_pca_snn")
#dataOnRDPlot(seu.merged,sce.merged,
#			 sprintf("%s/%s/%s",dirname(out.prefix),"harmony.tsne.Rtsne.algorithm1",basename(out.prefix)),
#			 rd="harmony.tsne.Rtsne")



