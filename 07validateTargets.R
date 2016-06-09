### libraries ###
library(biomaRt)
library(Vennerable)

### options ###
ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
chr <- c(1:22, "X", "Y", "MT")

### functions ###
calcEnrichmentHyper <- function(set1, set2, universe) {
    N <- length(universe)
    R <- length(set1)
    n <- length(set2)
    r <- sum(set1 %in% set2)
    res <- phyper(r-1, n, N-n, R, lower.tail=FALSE)
}

calcEnrichmentFisher <- function(set1, set2, universe) {
    a <- sum(set1 %in% set2)
    b <- length(set1) - a
    c <- length(set2) - a
    d <- length(universe) - a - b - c
    fisher <- fisher.test(matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE), alternative="two.sided")
}

### data ###
# protein coding genes
genes <- getBM(attributes="ensembl_gene_id", filters=c("chromosome_name", "biotype"), values=list(c(1:22, "X", "Y", "MT"), "protein_coding"), mart=ensembl)
genes <- unique(genes$ensembl_gene_id)

# prediction results
predres <- read.csv("../data/PredicitonResults.csv")
predres <- subset(predres, Prediction == 1)

# text mining results
gt <- read.delim("../data/docstore.gene_target.annotated.txt")
#gti <- read.delim("../data/docstore.gene_target_indication.annotated.txt")

# create sets
set1 <- genes[genes %in% unique(predres$Ensembl)]
set2 <- genes[genes %in% unique(gt$ensembl)]

# test for enrichment
res <- calcEnrichmentFisher(set1, set2, genes)
res$p.value
res$estimate

# Venn diagram
png(file.path("../data/Venn.png"), height=10*150, width=10*150, res=150)
plot(Venn(list("Predictions"=set1, "Text mining"=set2)))
dev.off()
