#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--limmaFile", type="character", required=TRUE, help="input limma file list")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
parser$add_argument("-c", "--occ", type="double", default=0.85, help="genes detected in >= OCC datasets will be kept [default %(default)s]")
parser$add_argument("-n", "--ncores", type="integer",default=16L, help="[default %(default)s]")
#parser$add_argument("-m", "--measurement",type="character",default="counts",help="[default %(default)s]")
args <- parser$parse_args()
print(args)


in.file <- args$limmaFile
out.prefix <- args$outPrefix
opt.gene.occ <- args$occ
opt.ncores <- args$ncores
#opt.measurement <- args$measurement

dir.create(dirname(out.prefix),F,T)

############## tune parametrs  ########
library("scPip")
library("data.table")
library("tictoc")
library("sscVis")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = opt.ncores)
options(stringsAsFactors = FALSE)

######################
de.limma.tb <- fread(in.file)

sce.pb.file <- sprintf("%s.sce.pb.rds",out.prefix)
if(file.exists(sce.pb.file)){
    loginfo(sprintf("load sce.pb file..."))
    sce.pb <- readRDS(sce.pb.file)
}else{
    tic("convertLimmaToSCE")
    sce.pb <- convertLimmaToSCE(de.limma.tb,out.prefix,ncores=8,
                        min.ncells=30,min.ncellsStudy=200,
                        gset.list=NULL,
                        de.mode="multiAsTwo",column.exp="meanScale",
                        TH.gene.occ=opt.gene.occ,
                        gene.used=NULL,colSet=list())
    toc()
}

tic("make.geneTableLong()")
make.geneTableLong(sce.pb,out.prefix,th.adj.P=0.01,th.dprime=0.15)
toc()

