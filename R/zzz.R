#' special gene list
#'
#' @docType data
#' @rdname g.geneOnUmap.list
#' @name g.geneOnUmap.list
#'
NULL

g.geneOnUmap.list <- list(
              "epithelial.01"=c("EPCAM", "SFN", "KRT8", "KRT10","KRT18",
                    "KRT19","KRTCAP3","KRT15",""),
              "endothelial"=c("PECAM1","VWF","CDH5","PLVAP","CLDN5",
                      "SLCO2A1","RAMP2","EGFL7","EMCN"),
              "CAF"=c("DCN", "COL1A2", "COL1A1", "COL3A1", "COL6A2", "COL6A1",
                  "PCOLCE", "C1S", "LUM"),
              "SMC"=c("MYH11", "ACTA2", "TAGLN", "CNN1", "MYL9",
                  "DES", "TPM2", "MYLK", "PLN"),
              "glia"=c("CRYAB", "CD9", "SPARC", "PMP22", "CDH19", "SEMA3B", "GPM6B", "NRXN1", "PLP1"),
              "B.01"=c("BLK", "CD19", "CD70", "CD79A","CD79B",
                   "MS4A1","FCER2","SDC1","JCHAIN"),
              "B.02"=c("TCL1A","RGS13","IL4R","IGKC","IGHG1","CD27","CD22","EBI3","STAG3"),
              "Mast"=c("ENPP3", "KIT", "PTPRC", "TPSB2","SLC18A2",
                   "FCER1A", "TPSAB1", "CPA3", "HPGDS"),
              "Neutrophil.01"=c("FUT4","CD3E","TLR2","CD14","CD33","FCGR1A","CSF3R","",""),
              "Neutrophil.02"=c("CD68","CSF3R","CXCR1","CXCR2","CXCR4","CXCL8","MMP9","VEGFA","ARG1"),
              "Monocyte"=c("CD14","VCAN", "S100A9","FCGR3A","FCEG1G",
                   "LST1","ITGAM", "HLA-DQA1","PTPRC"),
              "DC.01"=c("LILRA4","CD1C", "CLEC9A","CLEC10A","FCER1A",
                "RSG2","LAMP3", "THBD", "IDO1"),
              "DC.02"=c("ITGAX","IRF4", "IRF8", "CD274", "LTB",
              "CCL19","CCL21","IL3RA","NRP1","CLE4C"),
              "Macrophase.01"=c("CD68","THBS1","CD163", "C1QA","APOE",
                    "FCGR3A","MARCO", "S100A4","S100A9"),
              "Macrophase.02"=c("HLA-DQA1", "CST3", "LYZ", "ISGF6","FUT4",
                    "SEPP1","IFI30","FCGR1A","GPNMB"),
              "Macrophase.03"=c("CTSK", "MMP9", "SIGLEC15", "LYVE1","SPP1",
                    "FOLR2","C1QC","TREM2","LPL"),
              "Macrophase.04"=c("MT1G", "CCL18", "CD36", "CD52","FABP4",
                    "C3","CXCL8","PDK4","HK2"),
              "Macrophase.05"=c("FCN1", "IL1B", "NLRP3", "PLTP","CCR2",
                    "CCL2","CCL18","MMP9","SLC2A1"),
              "Macrophase.06"=c("GCHFR","CTSD","CTSL","CTSB", "APOE","APOC1",
                    "CCL2","CCL18","SLC2A1"),
              "Macrophase.07"=c("S100A8","S100A9", "CXCL10", "CXCL11", "OLR1",
                      "IDO1","CALHM6","FN1","SERPINA1"),
              "Macrophase.08"=c("CD14","FCGR3A", "CX3CR1", "MT1G", "CCL2",
                    "SLC2A1","FABP4","FOLR2","LYVE1"),
              "Macrophase.09"=c("FCN1","IL1B","GPNMB","C1QC","SPP1","TREM2","LYZ","CSF1R","CD163"),
              ##"Granulocyte"=c("CEACAM1","CEACAM6","CEACAM3","ENPP3","CD3E","","","",""),
              "NK"=c("FCGR3A", "NCR1","NCAM1","CD3E", "NKG7",
                 "KLRD1","KIR2DL3","KIR2DL4","KLRB1"),
              #"ILC"=c("IL4I1","KIT","KRT86","KRT81","ZFP36L1", "IL23R","CD83","NR4A1","EEF1A1"),
              "T.01"=c("CD4", "CD8A", "CD8B", "CD2", "CD3D",
                   "CD3G", "NKG7", "GNLY", "FGFBP2"),
              "T.02"=c("CX3CR1", "GZMK", "ZNF683","CD40LG","FOXP3",
                   "TNFRSF9","IL7R","IL2RA","TYROBP"),
              "T.03"=c("TRGC1","TRGC2","TRBC1","TRBC2","TRDC",
                   "TRAC","TRGV9","TRDV2","TRAV1-2"),
              "T.04"=c("CCR7","GZMK","LAYN","HAVCR2","CXCL13",
                   "CX3CR1","ZNF683","SLC4A10","CD160"),
              "T.05"=c("IL17A","RORC","IL23R","CCR6","BCL6",
                   "CXCR5", "CXCR6","XCL1","XCL2"),
              "T.06"=c("TCF7","MAL","CCL5","CCL4","CCL4L2",
                   "GZMA","GPR183","EOMES","KLRG1"),
              "T.07"=c("PDCD1","CTLA4","IFNG","TNF","GZMB",
                   "RGS1","S1PR1","ITGAE","MALAT1"),
              "cellCycle"=c("TYMSOS","PCLAF","TYMS","ZWINT","MCM2",
                    "MIR3917","UHRF1","TK1","MKI67"),
              "ISG"=c("ISG15","ISG20","IFIT1","IFIT2",
                  "IRF7", "OAS1","OAS2","STAT1"),
              "DIG"=c("DNAJB1","DUSP1","FOS","JUN","HSP90AA1",
                  "HSP90AB1","HSPA1A","HSPA6")
                        )


.onLoad <- function(libname , pkgname) {
    if(getRversion() >= "2.15.1"){
        utils::globalVariables(c(".",".N",".SD","variable","value","Group","x",
                                 "Dim1","Dim2",
                                 "vst.variance","vst.variance.standardized",
                                 "geneID","aid",
                                 "median.F.rank","hasSpeGene","fDataSets","medianRankSpeGene",
                                 "logFC","sig","comb.padj","comb.ES","sig.cate","meta.cluster",
                                 "freq.sig","geneSet.TF",
                                 "dataset",
                                 "dprime","vardprime","P.Value","adj.P.Val"))
    }
}




