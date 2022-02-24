#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
parser$add_argument("-e", "--expFile", type="character", required=TRUE, help="input seu file list")
parser$add_argument("-a", "--aid", type="character", default="SAMPLE_X", help="aid")
parser$add_argument("-n", "--ncores", type="integer",default=8L, help="[default %(default)s]")
parser$add_argument("-k", "--keep",type="character",help="Format is COLUMN_ID:COLUMN_VAL_1,COLUMN_VAL_2,COLUMN_VAL_3. Keep only cells with COLUMN_ID in one of COLUMN_VAL_1, COLUMN_VAL_2, and COLUMN_VAL_3.")
parser$add_argument("-d", "--distance", type="character", default="euclidean", help="euclidean, pearson, or spearman [default %(default)s]")
args <- parser$parse_args()
print(args)

########

out.prefix <- args$outPrefix
exp.file <- args$expFile
opt.keep <- args$keep
opt.ncores <- args$ncores
opt.aid <- args$aid
opt.distance <- args$distance

#out.prefix <- "./OUT.test/test.P032"
#exp.file <- "OUT.exp/ESCA.A20220210.seu.doubletFlt.all.rds"
#opt.keep <- "patient:P032"
#opt.ncores <- 4
#opt.aid <- "P032"
#opt.distance <- "euclidean"

dir.create(dirname(out.prefix),F,T)

library("Seurat")
library("copykat")
library("tictoc")


tic("readRDS(exp.file)")
seu <- readRDS(exp.file)
toc()

#### 
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

{

    exp.rawdata <- as.matrix(seu@assays$RNA@counts)

    cDir <- getwd()
    setwd(dirname(out.prefix))

    tic("copykat")
    copykat.test <- copykat(rawmat=exp.rawdata, id.type="S",
                            sam.name=opt.aid,
                            distance=opt.distance,
                            n.cores=opt.ncores,
                            output.seg="FLASE")
    toc()

    setwd(cDir)

    saveRDS(copykat.test,file=sprintf("%s.copykat.rds",out.prefix))

}

#############
#copykat.test <- readRDS(sprintf("%s.copykat.rds",out.prefix))

pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

##### heatmap plot 01
{

####    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
####
####    chr <- as.numeric(CNA.test$chrom) %% 2+1
####    rbPal1 <- colorRampPalette(c('black','grey'))
####    CHR <- rbPal1(2)[as.numeric(chr)]
####    chr1 <- cbind(CHR,CHR)
####
####    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
####    com.preN <- pred.test$copykat.pred
####    pred <- rbPal5(2)[as.numeric(factor(com.preN))]
####
####    cells <- rbind(pred,pred)
####    col_breaks = c(seq(-1,-0.4,length=50),
####                   seq(-0.4,-0.2,length=150),
####                   seq(-0.2,0.2,length=600),
####                   seq(0.2,0.4,length=150),
####                   seq(0.4, 1,length=50))
####
####    png(sprintf("%s.heatmap.01.png",out.prefix),width=800,height=700)
####    heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),
####              dendrogram="r",
####              distfun = function(x) parallelDist::parDist(x,threads = opt.ncores, method = "euclidean"),
####              hclustfun = function(x) hclust(x, method="ward.D2"),
####              ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
####              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
####              keysize=1, density.info="none", trace="none",
####              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
####              symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
####
####    legend("topright",
####           paste("pred.",names(table(com.preN)),sep=""),
####           pch=15,
####           col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1],
####           cex=0.6, bty="n")
####    dev.off()

}

##### heatmap plot 02
{

    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

    chr <- as.numeric(CNA.test$chrom) %% 2+1
    rbPal1 <- colorRampPalette(c('black','grey'))
    CHR <- rbPal1(2)[as.numeric(chr)]
    chr1 <- cbind(CHR,CHR)

    tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
    tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% make.names(tumor.cells))]
    hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads = opt.ncores, method = "euclidean"), method = "ward.D2")
    hc.umap <- cutree(hcc,2)

    rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
    subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
    cells <- rbind(subpop,subpop)

    col_breaks = c(seq(-1,-0.4,length=50),
                   seq(-0.4,-0.2,length=150),
                   seq(-0.2,0.2,length=600),
                   seq(0.2,0.4,length=150),
                   seq(0.4, 1,length=50))

    png(sprintf("%s.heatmap.02.png",out.prefix),width=800,height=700)
    heatmap.3(t(tumor.mat),
              dendrogram="r",
              distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"),
              hclustfun = function(x) hclust(x, method="ward.D2"),
              ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
              notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              keysize=1, density.info="none", trace="none",
              cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

    legend("topright",
           c("c1","c2"),
           pch=15,
           col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4],
           cex=0.9, bty='n')
    dev.off()

}



