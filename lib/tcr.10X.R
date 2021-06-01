#!/usr/bin/env Rscript

library("magrittr")
library("Biostrings")
library("rjson")
library("plyr")
library("data.table")
library("Matrix")
library("ggplot2")
library("ggpubr")



##### extend consensus 
CRegionSeq.list <- list()
CRegionSeq.list[["mouse"]] <- list("TRA"=DNAString("ACATCCAGAACCCAGAACCTGCTGTGTACCAGTTAAAAGATCCTCGGTCTCAGGACAGCACCCTCTGCCTGTTCACCGACTTTGACTCCCAAATCAATGTGCCGAAAACCATGGAATCTGGAACGTTCATCACTGACAAAACTGTGCTGGACATGAAAGCTATGGATTCCAAGAGCAATGGGGCCATTGCCTGGAGCAACCAGACAAGCTTCACCTGCCAAGATATCTTCAAAGAGACCAACGCCACCTACCCCAGTTCAGACGTTCCCTGTGATGCCACGTTGACTGAGAAAAGCTTTGAAACAGATATGAACCTAAACTTTCAAAACCTGTCAGTTATGGGACTCCGAATCCTCCTGCTGAAAGTAGCCGGATTTAACCTGCTCATGACGCTGAGGCTGTGGTCCAGTTGA"),
                        "TRB"=DNAString("AGGATCTGAGAAATGTGACTCCACCCAAGGTCTCCTTGTTTGAGCCATCAAAAGCAGAGATTGCAAACAAACAAAAGGCTACCCTCGTGTGCTTGGCCAGGGGCTTCTTCCCTGACCACGTGGAGCTGAGCTGGTGGGTGAATGGCAAGGAGGTCCACAGTGGGGTCAGCACGGACCCTCAGGCCTACAAGGAGAGCAATTATAGCTACTGCCTGAGCAGCCGCCTGAGGGTCTCTGCTACCTTCTGGCACAATCCTCGAAACCACTTCCGCTGCCAAGTGCAGTTCCATGGGCTTTCAGAGGAGGACAAGTGGCCAGAGGGCTCACCCAAACCTGTCACACAGAACATCAGTGCAGAGGCCTGGGGCCGAGCAGACTGTGGAATCACTTCAGCATCCTATCATCAGGGGGTTCTGTCTGCAACCATCCTCTATGAGATCCTACTGGGGAAGGCCACCCTATATGCTGTGCTGGTCAGTGGCCTGGTGCTGATGGCCATGGTCAAGAAAAAAAATTCCTGAGACAAACTTTTATGCATCCTGAGCCGTTCTTCACCCTGGCCATAGATTTTCCTGCACCTTCTCTAATTCCTGTTCCTAAGAACTTGTCTCTTCTTCCTCCATGGATATCCATCCTTCCTCGTTGACACCTTGACTCTGAAACAGACTAAATCAATAAAAACATGGAGTTAA"))

CRegionSeq.list[["human"]] <- list("TRA"=DNAString("ATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTCTGCCTATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAAGTAAGGATTCTGATGTGTATATCACAGACAAAACTGTGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGGCCTGGAGCAACAAATCTGACTTTGCATGTGCAAACGCCTTCAACAACAGCATTATTCCAGAAGACACCTTCTTCCCCAGCCCAGAAAGTTCCTGTGATGTCAAGCTGGTCGAGAAAAGCTTTGAAACAGATACGAACCTAAACTTTCAAAACCTGTCAGTGATTGGGTTCCGAATCCTCCTCCTGAAAGTGGCCGGGTTTAATCTGCTCATGACGCTGCGGCTGTGGTCC"),
                                   "TRB"=DNAString("AGGACCTGAAAAACGTGTTCCCACCCAAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTACCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACAGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCTGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTCACCTCCGAGTCTTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCTTGCTAGGGAAGGCCACCTTGTATGCCGTGCTGGTCAGTGCCCTCGTGCTGATGGCCATGGTCAAGAGAAAGGATTCCAGAGGCTAG"))

####### 

get.ext.consensus <- function(x,species="mouse"){
  ann.regionType <- sapply(seq_along(x$annotations),function(i){
    x$annotations[[i]]$feature$region_type
  })
  names(x$annotations) <- ann.regionType
  
  dna.seq <- DNAString(x$sequence)
  
  if(!"C-REGION" %in% ann.regionType)
  {
    return(list("AA"="NO_C-REGION","neucletide"="NO_C-REGION"))
  }else{
    chain <- x$annotations$`C-REGION`$feature$chain
    ### start: 0-based coordinate in the json file
    contig_match_start <- x$annotations$`C-REGION`$contig_match_start+1
    contig_match_end <- x$annotations$`C-REGION`$contig_match_end
    ## C-region sequence in consensus sequence:
    ##dna.seq[contig_match_start:contig_match_end]
    dna.ext.seq <- c(dna.seq[1:(contig_match_start-1)],CRegionSeq.list[[species]][[chain]])
    start_codon_pos <- x$start_codon_pos
    aa.seq <- translate(dna.ext.seq[(start_codon_pos+1):length(dna.ext.seq)])
    ### original aa sequence
    ###x$aa_sequence
    return(list("AA"=as.character(aa.seq),
                "neucletide"=as.character(dna.ext.seq[(start_codon_pos+1):length(dna.ext.seq)])))
  }
}


get.annotation.match.start.leftmost <- function(x){
  ann.regionType <- sapply(seq_along(x$annotations),function(i){
    x$annotations[[i]]$feature$region_type
  })
  names(x$annotations) <- ann.regionType
  
  if("5'UTR" %in% ann.regionType){
    return(x$annotations[["5'UTR"]][["annotation_match_start"]])
  }else if("L-REGION+V-REGION" %in% ann.regionType){
    return(x$annotations[["L-REGION+V-REGION"]][["annotation_match_start"]])
  }else{
    return("Unexpected region_type")  
  }
}

get.annotation.feature <- function(x,region_type){
  ann.regionType <- sapply(seq_along(x$annotations),function(i){
    x$annotations[[i]]$feature$region_type
  })
  names(x$annotations) <- ann.regionType
  
  if(region_type %in% ann.regionType){
    return(x$annotations[[region_type]][["feature"]][["gene_name"]])
  }else{
    return("")  
  }
}



##### re-define clonotypes and output tabels
formatTCRDat <- function(dataset.id,sample.info.file)
{

    sample.info.tb <- fread(sample.info.file)
    sample.info.tb <- sample.info.tb[libType=="TCR" & cancerType==dataset.id,]
    sample.info.tb[,sample.uid:=as.integer(factor(sampleID))]

    tcr.tb <- as.data.table(ldply(seq_len(nrow(sample.info.tb)),function(i){
        sampleID <- sample.info.tb$sampleID[i]
        sample.uid <- sample.info.tb$sample.uid[i]
        in.file <- sprintf("OUT.TCR/filtered_contig_annotations.%s.csv.gz",sampleID)
        in.tb <- fread(in.file)
        a.name <- gsub("-\\d+","",in.tb$barcode)
        in.tb$barcode <- sprintf("%s-%d",a.name,sample.uid)
        b.name <- gsub(".+_contig","contig",in.tb$contig_id)
        in.tb$contig_id <- sprintf("%s_%s",in.tb$barcode,b.name)
        in.tb$sampleID <- sampleID
        return(in.tb)
    }))
    tcr.detail.tb <- merge(tcr.tb,sample.info.tb[,c("sampleID","patient"),with=F],by="sampleID")

}

assignCloneID <- function(tcr.detail.tb,pid="")
{
    tcr.slim.tb <- tcr.detail.tb[productive=="True" & high_confidence==TRUE & chain %in% c("TRA","TRB"),
                                 ][order(sampleID,barcode,chain),]

    ### filter out cell with more than 2 productive TRA or 2 productive TRB
    print(tcr.slim.tb[,.N,by=c("barcode","chain")][,table(chain,N)])
    stat.chain <- tcr.slim.tb[,.(.N),by=c("barcode","chain")]
    f.stat.chain <- stat.chain[N>2,]
    tcr.slim.tb <- tcr.slim.tb[!(barcode %in% f.stat.chain$barcode),]
    print(tcr.slim.tb[,.N,by=c("barcode","chain")][,table(chain,N)])

    ### assign clonetype id
    tcr.slim.tb[,clonetype.str:=sprintf("%s:%s:%s:%s:%s:%s",chain,v_gene,d_gene,j_gene,c_gene,cdr3_nt)]

    ####
    clonetypeInfo <- list()
    #contig2clonetype <- list()
    cell2clonetype <- list()

    tcr.slim.list <- split(tcr.slim.tb,list(tcr.slim.tb$barcode))
    ##tcr.slim.list <- split(tcr.slim.tb,list(tcr.slim.tb$patient,tcr.slim.tb$barcode))

    g.clone.uid <- 1

    id.test <- c()

    for(i in seq_along(tcr.slim.list)){
        x <- tcr.slim.list[[i]]
        cell.id <- names(tcr.slim.list)[i]
        if(length(unique(x$chain))==2){
            x.a <- x[chain=="TRA",]
            x.b <- x[chain=="TRB",]
            x.comb <- expand.grid(x.a$clonetype.str,x.b$clonetype.str,stringsAsFactors=F)
            x.comb$Var3 <- sprintf("%s.%s",x.comb$Var1,x.comb$Var2)
            for(j in seq_len(nrow(x.comb))){
                if(!(cell.id %in% names(cell2clonetype))){
                    if(x.comb$Var3[j] %in% names(clonetypeInfo)){
                        clonetypeInfo[[x.comb$Var3[j]]][["cells"]] <- c(clonetypeInfo[[x.comb$Var3[j]]][["cells"]],cell.id)
                    }else{
                        clonetypeInfo[[x.comb$Var3[j]]] <- list(cells=c(cell.id),uid=g.clone.uid)
                        g.clone.uid <- g.clone.uid+1
                        #cat(sprintf("******* %d\n",clone.uid))
                        id.test <- c(id.test,g.clone.uid)
                    }
                    cell2clonetype[[cell.id]] <- clonetypeInfo[[x.comb$Var3[j]]][["uid"]]
                }else{
                    if(x.comb$Var3[j] %in% names(clonetypeInfo)){
                        clonetypeInfo[[x.comb$Var3[j]]][["cells"]] <- c(clonetypeInfo[[x.comb$Var3[j]]][["cells"]],cell.id)
                    }else{
                        clonetypeInfo[[x.comb$Var3[j]]] <- list(cells=c(cell.id),
                                                                uid=cell2clonetype[[cell.id]])
                    }
                }
            }
        }
    }

    #tcr.slim.tb

    tcr.out.tb <- tcr.slim.tb[,cbind(.SD,data.table(chain.id=sprintf("%s%d",chain,seq_len(nrow(.SD))))),
                               by=c("sampleID","barcode","chain")]
    tcr.out.tb <- tcr.out.tb[,-c("chain","is_cell","high_confidence","productive","full_length",
                                 ##"clonetype.str", "patient","library.id")]
                                 "clonetype.str")]
    o.dat <- melt(tcr.out.tb,id.vars=c("patient","library.id","sampleID","barcode","chain.id"))
    vv <- expand.grid(levels(o.dat$variable),sort(unique(o.dat$chain.id)))
    o.dat[,variable:=factor(sprintf("%s.%s",chain.id,variable), levels=sprintf("%s.%s",vv$Var2,vv$Var1))]
    o.dat <- o.dat[,-c("chain.id")]
    o.dat <- dcast(o.dat,patient+library.id+sampleID+barcode~variable)
	tcr.out.tb <- o.dat
    #tcr.out.tb <- merge(o.dat,sample.info.tb[,-c("filepath")],by="sampleID")
    ####tcr.out.tb[,cloneID:=sprintf("C%04d",unlist(cell2clonetype)[sprintf("%s.%s",patient,barcode)])]
    tcr.out.tb[,cloneID:=sprintf("C%06d",unlist(cell2clonetype)[sprintf("%s",barcode)])]
    tcr.out.tb <- tcr.out.tb[cloneID!="C    NA",]
	tcr.out.tb[,cloneID:=sprintf("%s_%s",pid,cloneID)]
    tcr.out.tb <- tcr.out.tb[order(cloneID,sampleID,barcode),]
    stat.clone.freq <- tcr.out.tb[,.N,by="cloneID"]
    tcr.out.tb$clone.N <- stat.clone.freq[match(tcr.out.tb$cloneID,stat.clone.freq$cloneID),][["N"]]
	tcr.out.tb[,cloneID:=sprintf("%s:%d",cloneID,clone.N)]

	return(tcr.out.tb)
#    conn <- gzfile(sprintf("OUT.TCR/TCR.slim.%s.txt.gz",dataset.id),"w")
#    write.table(tcr.out.tb,file=conn,row.names=F,sep="\t",quote=F)
#    close(conn)
#
#    saveRDS(tcr.out.tb,file=sprintf("OUT.TCR/TCR.slim.%s.rds",dataset.id))

}


#######  dump info from json
cnsSeqFromJSON <- function(aid,in.json.file,out.file)
{
    ########
    in.dat <- fromJSON(file=in.file)

    ######### 
    out.df <- ldply(seq_along(in.dat),function(i){
      ext.cns <- get.ext.consensus(in.dat[[i]],species="human")
      data.frame(sampleID=aid,
                 raw_consensus_id=in.dat[[i]]$contig_name,
                 filtered=in.dat[[i]]$filtered,
                 high_confidence=in.dat[[i]]$high_confidence,
                 productive=in.dat[[i]]$productive,
                 VGene=get.annotation.feature(in.dat[[i]],"L-REGION+V-REGION"),
                 DGene=get.annotation.feature(in.dat[[i]],"D-REGION"),
                 JGene=get.annotation.feature(in.dat[[i]],"J-REGION"),
                 cdr3=in.dat[[i]]$cdr3,
                 cdr3_seq=in.dat[[i]]$cdr3_seq,
                 aa_squence=in.dat[[i]]$aa_sequence,
                 aa_ext_sequence=ext.cns[["AA"]],
                 neucletide_ext_sequence=ext.cns[["neucletide"]],
                 annotation_match_start_leftmost=get.annotation.match.start.leftmost(in.dat[[i]]),
                 cr_clonotype=in.dat[[i]]$clonotype,
                 cr_clonotype_freq=in.dat[[i]]$info$clonotype_freq,
                 cr_clonotype_prop=in.dat[[i]]$info$clonotype_prop
                )
    })

    write.table(out.df,file = out.file,quote = F,sep = "\t",row.names = F)

}











