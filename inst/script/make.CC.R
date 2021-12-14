

###### functions from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){

    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    humanx <- unique(genesV2[, 2])

    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
}


# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){

    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

    humanx <- unique(genesV2[, 2])

    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
}
#########

data("cc.genes",package="Seurat",envir=a.env)
a.env$cc.genes$s.genes
a.env$cc.genes$g2m.genes

cc.genes <- list("s.genes"=convertHumanGeneList(a.env$cc.genes$s.genes),
                       "g2m.genes"=convertHumanGeneList(a.env$cc.genes$g2m.genes))
save(cc.genes,file="cc.genes.mouse.rda")
