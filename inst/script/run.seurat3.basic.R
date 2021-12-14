#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-a", "--aFile", type="character", required=TRUE, help="input seu file list")
parser$add_argument("-b", "--bFile", type="character", required=TRUE, help="input sce file list")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
parser$add_argument("-c", "--stype", type="character", help="only analyze stype specified (default all)")
parser$add_argument("-u", "--geneIDFile", type="character", help="gene id mapping file")
parser$add_argument("-d", "--npc", type="integer",default=15L, help="[default %(default)s]")
parser$add_argument("-n", "--ncores", type="integer",default=16L, help="[default %(default)s]")
parser$add_argument("-m", "--measurement",type="character",default="counts",help="[default %(default)s]")
parser$add_argument("-r", "--resolution",type="character",default="2",help="[default %(default)s]")
parser$add_argument("-t", "--scTransform",action="store_true",default=FALSE,help="[default %(default)s]")
parser$add_argument("-y", "--harmony",action="store_true",default=FALSE,help="[default %(default)s]")
parser$add_argument("-g", "--deg",action="store_true",default=FALSE,help="[default %(default)s]")
parser$add_argument("-w", "--ncellDEG",type="integer",default=1500,
                    help="number of cells to downsample to for each group. used in DEG analysis. [default %(default)s]")
parser$add_argument("-s", "--scale",action="store_true",default=FALSE,help="[default %(default)s]")
parser$add_argument("-j", "--corVar", type="character", default="S.Score,G2M.Score,DIG.Score1",
		    help="subset of S.Score,G2M.Score,DIG.Score1,ISG.Score1, or NULL. If correct something, always correct for batchV and percent.mito. [default %(default)s]")
parser$add_argument("-f", "--filterout",type="character",help="Format is COLUMN_ID:COLUMN_VAL_1,COLUMN_VAL_2,COLUMN_VAL_3. Filter out cells with COLUMN_ID in one of COLUMN_VAL_1, COLUMN_VAL_2, and COLUMN_VAL_3.")
parser$add_argument("-k", "--keep",type="character",help="Format is COLUMN_ID:COLUMN_VAL_1,COLUMN_VAL_2,COLUMN_VAL_3. Keep only cells with COLUMN_ID in one of COLUMN_VAL_1, COLUMN_VAL_2, and COLUMN_VAL_3.")
parser$add_argument("-x", "--removeContamination",type="character",help="comma separated string indicates subset of predefined signature (plasmaB, caf, epi, T, cd8, and mac) will be calculated and cells with high signature scores will be removed. For example, use plasmaB:0.75,caf:0.75,epi:0.75,T:0.25 for myeloid cell analysis. The numer after colon is the threshold")
parser$add_argument("-p", "--platform",type="character",default="10X",
                    help="platform such as 10X, SmartSeq2 [default %(default)s)]")
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
opt.harmony <- args$harmony
opt.doDEG <- args$deg
opt.ncell.deg <- args$ncellDEG
opt.scale <- args$scale
opt.cor.var <- if(args$corVar=="") c("") else unlist(strsplit(args$corVar,",",perl=T))
opt.filterout <- args$filterout
opt.keep <- args$keep
opt.geneIDFile <- args$geneIDFile
opt.removeContamination <- args$removeContamination

dir.create(dirname(out.prefix),F,T)

saveRDS(args,file=sprintf("%s.args.rds",out.prefix))
###args <- readRDS(file=sprintf("%s.args.rds",out.prefix))

############## tune parametrs  ########
suppressMessages(library("sscVis"))
suppressMessages(library("sscClust"))
suppressMessages(library("Seurat"))
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
options(stringsAsFactors = FALSE)

dat.ext.dir <- system.file("extdata",package="scPip")
gene.exclude.file <- sprintf("%s/exclude.gene.misc.misc.v3.RData",dat.ext.dir)

######################

gene.mapping.table <- NULL
if(!is.null(opt.geneIDFile) && file.exists(opt.geneIDFile) && grepl("\\.rds$",opt.geneIDFile,perl=T)){
    gene.mapping.table <- readRDS(opt.geneIDFile)
}

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

##### check & clean sce #####
if(!is.null(sce) && "percent.mito" %in% colnames(colData(sce))){
    if(all(is.na(sce$percent.mito))){
        sce$percent.mito <- NULL
    }
}

if(!is.null(sce) && !("seu.id" %in% colnames(rowData(sce)))){
    rowData(sce)[["seu.id"]] <- gsub("_","-",rowData(sce)[["display.name"]])
}

if(!is.null(sce) && ("libraryID" %in% colnames(colData(sce))) && all(is.na(sce$libraryID))){
    sce$libraryID <- "Unk"
}



if(!is.null(seu) && !is.null(opt.filterout)){
	if(!file.exists(opt.filterout)){
		col.filter <- unlist(strsplit(opt.filterout,":"))[1]
		col.value <- unlist(strsplit(unlist(strsplit(opt.filterout,":"))[2],","))
        if(col.filter %in% colnames(seu[[]])){
            cat(sprintf("filter cells with %s in c(%s)\n",col.filter,paste(col.value,collapse=",")))
            f.cell <- seu[[]][,col.filter] %in% col.value
            print(summary(f.cell))
            seu <- seu[,!f.cell]
        }else{
            warning(sprintf("The meta-data doesnot contain %s\n",col.filter))
        }
	}
}

if(!is.null(seu) && !is.null(opt.keep)){
	if(!file.exists(opt.keep)){
		col.keep <- unlist(strsplit(opt.keep,":"))[1]
		col.value <- unlist(strsplit(unlist(strsplit(opt.keep,":"))[2],","))
        if(col.keep %in% colnames(seu[[]])){
            cat(sprintf("keep only cells with %s in c(%s)\n",col.keep,paste(col.value,collapse=",")))
            f.cell <- seu[[]][,col.keep] %in% col.value
            print(summary(f.cell))
            seu <- seu[,f.cell]
        }else{
            warning(sprintf("The meta-data doesnot contain %s\n",col.keep))
        }
    }
}

if(!is.null(opt.removeContamination)){
    ### opt.removeContamination <- "plasmaB:0.75,caf:0.75,epi:0.75,T:0.25"

    g.cont.gene.list <- list("plasmaB"=c("JCHAIN"),
                             "caf"=c("COL1A2", "COL1A1", "COL3A1","LUM"),
                             "epi"=c("KRT18","KRT19","EPCAM"),
                             "T"=c("CD3D","CD3G"),
                             "cd8"=c("CD8A","CD8B"),
                             "mac"=c("LYZ","C1QA","C1QB","CD68"))

    sig.vec <- unlist(strsplit(opt.removeContamination,","))
    sig.name <- unname(sapply(sig.vec,function(x){ unlist(strsplit(x,":"))[1] }))
    sig.thre <- unname(sapply(sig.vec,function(x){ as.numeric(unlist(strsplit(x,":"))[2]) }))

    for(i in seq_along(sig.name)){
        loginfo(sprintf("calculate signature score of %s cells ...",sig.name[i]))
        seu <- fill.contamination(seu,out.prefix,
                                  g.name=sig.name[i],
                                  g.test=g.cont.gene.list[[sig.name[i]]],
                                  score.t=sig.thre[i],
                                  vis.v=c(0.25,0.5,0.75,1))
    }
   
    idx.sig.class <- intersect(sprintf("%s.class",sig.name),colnames(seu[[]]))
    f.cont.mtx <- seu[[]][,idx.sig.class,drop=F]
    f.cont <- rowSums(f.cont.mtx) > 0
    loginfo(sprintf("A total number of potential contamination: %d\n",sum(f.cont)))
    print(colSums(f.cont.mtx==T))
    seu <- seu[,!f.cont]

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
                        use.harmony=opt.harmony,
                        do.deg=opt.doDEG,
                        ncell.deg=opt.ncell.deg,
                        gene.mapping.table=gene.mapping.table,
                        ###do.adj=T,do.scale=F,
			            cor.var=opt.cor.var,
                        do.scale=opt.scale,
                        plot.rd=c("umap"),
						opt.res=opt.resolution,
                        opt.npc=opt.npc,ncores=opt.ncores)
toc()


