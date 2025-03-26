#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressMessages(library("anndata"))
suppressMessages(library("scibet"))
#suppressMessages(library(tidyverse))
suppressMessages(library("data.table"))
suppressMessages(library("plyr"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="input h5ad file. 'cellID' is required in the obs")
parser$add_argument("-d", "--dataset", type="character", default="DataSet01", help="dataset name. [default %(default)s]")
parser$add_argument("-m", "--modelPath", type="character", help="model path. If not specified, use a ESCA model stored under the directory of the package ")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="out prefix")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
out.prefix <- args$outPrefix
opt.dataset <- args$dataset
model_path <- args$modelPath
#in.file = c('OUT.ann.all/ESCA.HuyQDinh2021.scanpy.h5ad')
#out.prefix = './OUT.scibet/PanC'
#opt.dataset = 'HuyQDinh2021'
dir.create(dirname(out.prefix),F,T)


##############Set larger connection buffer for reading###################
#Sys.setenv ("VROOM_CONNECTION_SIZE" = 131072 * 2) 
if(is.null(model_path)){
    model_path <- system.file("extdata/scibet/model_ESCA",package="scPip")
}

############## Prediction ################################################
{
    
    model_all = readRDS(file.path(model_path,'model.major.all.rds'))
    cellSubtype.vec <- c("T8","Th","Treg","ILC","B","Plasma",
                         "Neutro","pDC","DC","M","Mast",
                         "Endo","Fibro","SMC","Epi","Glia")
    model_sub_list <- sapply(cellSubtype.vec, function(x){
                                readRDS(file.path(model_path, paste0("model.sub.", x, ".rds")))
                             })
    names(model_sub_list) <- cellSubtype.vec
    
    #### read and process expression matrix
    dat_test <- read_h5ad(in.file)
    dat_test_mtx = expm1(as.matrix(dat_test$X))
    
    #### predict major cell types
    prd.vec = model_all(dat_test_mtx)
    prd.tb <- data.table("cellID"=rownames(dat_test_mtx),dataset=opt.dataset,"scibetMajor"=prd.vec)
    
    ##### within each major group, predict sub types
    prd.full.tb <- as.data.table(ldply(cellSubtype.vec,function(ct_sel) {
        prd.x.tb <- prd.tb[scibetMajor==ct_sel,]
        dat_test_mtx_x = dat_test_mtx[prd.x.tb$cellID,]
        prd.x.vec <- model_sub_list[[ct_sel]](dat_test_mtx_x)
        prd.x.tb$scibetSub <- prd.x.vec
        print(sprintf("cellSubtype done (%s)!",ct_sel))
        rm(dat_test_mtx_x)
        return(prd.x.tb)
    }))
    rm(dat_test_mtx,dat_test)
    gc()
}

saveRDS(prd.full.tb,file=sprintf("%s.ann.tb.rds",out_path))

