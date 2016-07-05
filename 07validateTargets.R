### libraries ###
library(biomaRt)
library(Vennerable)
library(ggplot2)

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

# training and test sets (to remove from text mining results)
dataset <- readRDS("../data/dataset.rds")
dataset <- rownames(dataset)

# text mining results
textmining <- read.delim("../data/docstore.gene_target.annotated.txt")
textmining <- textmining[!textmining$ensembl %in% dataset, ]

# create sets
set1 <- genes[genes %in% unique(predres$Ensembl)]
set2 <- genes[genes %in% unique(textmining$ensembl)]

# test for enrichment
res <- calcEnrichmentFisher(set1, set2, genes)
res$p.value
res$estimate

# Venn diagram
png(file.path("../data/Venn.png"), height=10*300, width=10*300, res=300)
plot(Venn(list("Predictions"=set1, "Text mining"=set2)))
dev.off()

# permutation test
n <- 10000
perm <- data.frame(pvalue=rep(1, n), oddsratio=rep(1, n))
for (i in 1:n) {
    set1 <- sample(genes, length(set1))
    res <- calcEnrichmentFisher(set1, set2, genes)
    perm$pvalue[i] <- res$p.value
    perm$oddsratio[i] <- res$estimate
}

# histograms
png("../data/HistogramPV.png", height=6*300, width=6*300, res=300)
print(
	  ggplot(perm, aes(pvalue)) +
		  geom_histogram(colour="black", fill="lightskyblue", bins=100) +
		  xlab("Log p-value") + 
		  ylab("Count") +
		  scale_x_log10() +
          theme_bw(14)
)
dev.off()
png("../data/HistogramOR.png", height=6*300, width=6*300, res=300)
print(
	  ggplot(perm, aes(oddsratio)) +
		  geom_histogram(colour="black", fill="lightseagreen", bins=100) +
		  xlab("Odds ratio") + 
		  ylab("Count") +
          theme_bw(14)
)
dev.off()
