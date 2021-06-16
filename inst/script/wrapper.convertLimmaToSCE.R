#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--limmaFile", type="character", required=TRUE, help="input limma file list")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
parser$add_argument("-n", "--ncores", type="integer",default=16L, help="[default %(default)s]")
#parser$add_argument("-m", "--measurement",type="character",default="counts",help="[default %(default)s]")
args <- parser$parse_args()
print(args)


in.file <- args$limmaFile
out.prefix <- args$outPrefix
opt.ncores <- args$ncores
#opt.measurement <- args$measurement

dir.create(dirname(out.prefix),F,T)

############## tune parametrs  ########
library("scPip")
library("data.table")
library("tictoc")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = opt.ncores)
options(stringsAsFactors = FALSE)

######################
de.limma.tb <- fread(in.file)

tic("convertLimmaToSCE")
sce.pb <- convertLimmaToSCE(de.limma.tb,out.prefix,ncores=8,
				    min.ncells=30,min.ncellsStudy=200,
				    gset.list=NULL,
				    de.mode="multiAsTwo",column.exp="meanScale",
				    gene.used=NULL,colSet=list())
toc()


#saveRDS(de.out,file=sprintf("%s.de.out.rda",out.prefix))


