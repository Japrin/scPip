#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("this.path"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="input files list")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
parser$add_argument("-j", "--corVar", type="character", default="S.Score,G2M.Score,DIG.Score1",
		    help="subset of S.Score,G2M.Score,DIG.Score1,ISG.Score1,score.MALAT1 [default %(default)s]")
parser$add_argument("-f", "--excludeCells", type="character",
		    help="vector containing cells to be excluded, stored in a .rds file ")
##parser$add_argument("-d", "--npc", type="integer",default=15L, help="[default %(default)s]")
args <- parser$parse_args()
cat(sprintf("############# arguments #############\n"))
print(args)
cat(sprintf("#####################################\n"))

exp.list.file <- args$inFile
out.prefix <- args$outPrefix
#gene.file <- args$geneFile
opt.cor.var <- unlist(strsplit(args$corVar,",",perl=T))
opt.excludeCells.file <- args$excludeCells

dir.create(dirname(out.prefix),F,T)
saveRDS(args,file=sprintf("%s.args.rds",out.prefix))
##args <- readRDS(sprintf("%s.args.rds",out.prefix))

####
#exp.list.file <- "list/obj.T.list.r1.list"
#out.prefix <- "./OUT.int.S3.T/int.S3.T"
dir.create(dirname(out.prefix),F,T)

###opt.cor.var <- c("S.Score","G2M.Score", "DIG.Score1")

g.src.dir <- sprintf("%s/../",this.dir())
source(sprintf("%s/lib/inte.comb.miniClust.lib.R",g.src.dir))
#source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/scPip/lib/inte.comb.miniClust.lib.R")

ncores <- 12
gene.exclude.file <- sprintf("%s/data/geneSet/exclude.gene.misc.misc.v3.RData",g.src.dir)
#gene.exclude.file <- "/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/scPip/data/geneSet/exclude.gene.misc.misc.v3.RData"

options(stringsAsFactors = FALSE)

#env.misc <- loadToEnv(gene.exclude.file)

#### data.id measurement     platform        defile  scefile seufile
exp.list.table <- fread(cmd=sprintf("awk '!/^#/' %s",exp.list.file))

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = ncores)

contamination.vec <- NULL
if(!is.null(opt.excludeCells.file) && file.exists(opt.excludeCells.file)){
    contamination.vec <- readRDS(opt.excludeCells.file)
}

tic("run.inte.metaClust")
ret.list <- run.inte.metaClust(exp.list.table, out.prefix, gene.exclude.file,
			       #cor.cellCycle=T,cor.MALAT1=F,cor.DIG=T,cor.ISG=T,
			       cor.var=opt.cor.var,
			       contamination.vec=contamination.vec,
			       ncores=12,npc=15,res.hi=50,TH.gene.occ=0.85)
toc()

seu.merged <- ret.list[["seu.merged"]]
sce.merged <- ret.list[["sce.merged"]]
meta.tb <- ret.list[["meta.tb"]]

#seu.merged <- readRDS(file=sprintf("%s.seu.merged.rds",out.prefix))
#sce.merged <- readRDS(file=sprintf("%s.sce.merged.rds",out.prefix))
#meta.tb <- readRDS(file=sprintf("%s.meta.tb.rds",out.prefix))

#dataOnRDPlot(seu.merged,sce.merged,
#			 sprintf("%s/%s/%s",dirname(out.prefix),"umap.algorithm1",basename(out.prefix)),
#			 rd="umap",graph.name="RNA_pca_snn")
dataOnRDPlot(seu.merged,sce.merged,
			 sprintf("%s/%s/%s",dirname(out.prefix),"harmony.umap.algorithm1",basename(out.prefix)),
			 rd="harmony.umap")
#dataOnRDPlot(seu.merged,sce.merged,
#			 sprintf("%s/%s/%s",dirname(out.prefix),"tsne.Rtsne.algorithm1",basename(out.prefix)),
#			 rd="tsne.Rtsne",graph.name="RNA_pca_snn")
#dataOnRDPlot(seu.merged,sce.merged,
#			 sprintf("%s/%s/%s",dirname(out.prefix),"harmony.tsne.Rtsne.algorithm1",basename(out.prefix)),
#			 rd="harmony.tsne.Rtsne")



