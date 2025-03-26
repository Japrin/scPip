#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="input h5ad file. 'cellID' is required in the obs")
#parser$add_argument("-d", "--dataset", type="character", default="DataSet01", help="dataset name. [default %(default)s]")
parser$add_argument("-a", "--annFile", type="character", required=TRUE, help="annotation file by run.scibet.R")
#parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="out prefix")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
#out.prefix <- args$outPrefix
ann.file <- args$annFile
out.prefix = gsub(".h5ad$","",in.file)
#dir.create(dirname(out.prefix),F,T)

############ settings ##########
{
    library("R.utils")
    library("data.table")
    library("sscVis")
    library("tictoc")
    library("plyr")
    library("ggplot2")
    library("ggpubr")
    library("anndata")

    require("reticulate")
    sc <- import("scanpy")
    plt <- import("matplotlib.pyplot")

    sc$set_figure_params(dpi=300,dpi_save=300,fontsize=12)

    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = 8)

}

############ function ##########
{
    split_umap <- function(adata,split_by,plot.ncol=4,plot.nrow=NULL,...)
    {

        my.xlim <- pretty(adata$obsm$X_umap[,1])
        my.ylim <- pretty(adata$obsm$X_umap[,2])
        my.xlim <- my.xlim[c(1,length(my.xlim))]
        my.ylim <- my.ylim[c(1,length(my.ylim))]

        categories <- levels(adata$obs[[split_by]])
        if(is.null(plot.nrow)){
            plot.nrow <- ceiling(length(categories) / plot.ncol)
        }
        fig_axs <- plt$subplots(as.integer(plot.nrow),
                                as.integer(plot.ncol),
                                figsize=c(5*plot.ncol, 4*plot.nrow))
        names(fig_axs) <- c("fig","axs")
        axs <- matrix(fig_axs[["axs"]],nrow=1)
        for(i in seq_along(categories)){
            ax <- axs[[i]]
            sc$pl$umap(adata[adata$obs[[split_by]] == categories[i] ], ax=ax, show=FALSE, title=categories[i],...)
            ax$set_xlim(my.xlim[1],my.xlim[2])
            ax$set_ylim(my.ylim[1],my.ylim[2])
        }
        plt$tight_layout()

    }

}

############ annotation ##########
{

    adata <- anndata::read_h5ad(in.file)
    ann.tb <- readRDS(ann.file)

    setkey(ann.tb,"cellID")
    ann.tb <- ann.tb[adata$obs_names,]
    print("all(adata$obs_names == ann.tb$cellID) ?")
    print(all(adata$obs_names == ann.tb$cellID))
    adata$obs$scibetMajor <- ann.tb$scibetMajor
    adata$obs$scibetSub <- ann.tb$scibetSub

    ###
    {

        my.xlim <- pretty(adata$obsm$X_umap[,1])
        my.ylim <- pretty(adata$obsm$X_umap[,2])
        my.xlim <- my.xlim[c(1,length(my.xlim))]
        my.ylim <- my.ylim[c(1,length(my.ylim))]

        sp <- plt$subplots(figsize=c(4,4))
        ax <- sp[[2]]
        sc$pl$umap(adata, color="scibetMajor",legend_loc='on data',ax=sp[[2]],legend_fontsize=10)
        ax$set_xlim(my.xlim[1],my.xlim[2]);
        ax$set_ylim(my.ylim[1],my.ylim[2]);
        plt$tight_layout()
        plt$savefig(sprintf("%s.%s.%s.png",out.prefix,"umap","scibetMajor"))

        split_umap(adata,"scibetMajor",plot.ncol=4,color="scibetMajor",legend_loc='on data',legend_fontsize=10)
        plt$savefig(sprintf("%s.%s.%s.split.00.png",out.prefix,"umap","scibetMajor"))
    
    }

    tic("write_h5ad ...")
    adata$write_h5ad(sprintf("%s.scibet.h5ad",out.prefix))
    toc()
    #adata <- anndata::read_h5ad(sprintf("%s.scibet.h5ad",out.prefix))

    m.tb <- as.data.table(adata$obs)
    saveRDS(m.tb,file=sprintf("%s.m.tb.rds",out.prefix))
    
}

