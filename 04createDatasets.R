### libraries ###
library(biomaRt)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(gplots)
library(Rtsne)

### options ###
set.seed(16)
ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
chr <- c(1:22, "X", "Y", "MT")
type="protein_coding"

### data ###
# get all protein coding genes
genes <- getBM(
               attributes="ensembl_gene_id",
               filters=c("chromosome_name", "biotype"),
               values=list(chr, type),
               mart=ensembl)
genes <- genes[order(genes$ensembl_gene_id), , drop=FALSE]

## create numeric features from open targets
completeset <- read.csv("/GWD/bioinfo/projects/bix-analysis-stv/data/CTTV/v3.0/public.matrix.csv.gz")
# only use direct associations, remove known_drug and literature
completeset <- subset(completeset, Is.direct == TRUE, c(EnsemblId, OntologyId, affected_pathway, animal_model, genetic_association, rna_expression, somatic_mutation))
# remove lower confidence animal_model associations
completeset$animal_model[completeset$animal_model < 0.4] <- 0
# aggregate
completeset <- aggregate(completeset[3:ncol(completeset)], by=list(EnsemblId=completeset$EnsemblId), FUN=mean)
# merge
completeset <- merge(genes, completeset, by=1, all=FALSE)
saveRDS(completeset, file.path("../data/completeset.rds"))

# read pharmaprojects data for target information
pharmaprojects <- read.delim("/GWD/bioinfo/projects/bix-analysis-stv/data/pharmaceutical/pipeline/pipeline_triples.txt")
pharmaprojects <- subset(pharmaprojects,
                                GlobalStatus == "Clinical Trial" |
                                GlobalStatus == "Launched" |
                                GlobalStatus == "Phase I Clinical Trial" |
                                GlobalStatus == "Phase II Clinical Trial" |
                                GlobalStatus == "Phase III Clinical Trial" |
                                GlobalStatus == "Pre-registration" |
                                GlobalStatus == "Preclinical" |
                                GlobalStatus == "Registered")
pharmaprojects <- getBM(
                        attributes="ensembl_gene_id",
                        filters=c("entrezgene", "chromosome_name", "biotype"),
                        values=list(pharmaprojects$Target_EntrezGeneId, chr, type),
                        mart=ensembl)

# positive cases: these are targets according to targetpedia and/or pharmaprojects
positive <- unique(pharmaprojects)
positive <- completeset[completeset$ensembl_gene_id %in% positive$ensembl_gene_id, ]
positive$target <- 1

# unknown cases: it's not known whether these are targets or not
unknown <- completeset[!completeset$ensembl_gene_id %in% positive$ensembl_gene_id, ]
unknown <- unknown[sample(nrow(unknown), nrow(positive)), ]
unknown$target <- 0

# dataset is made of positive and unknown cases
# will be splitted into traing and test sets
dataset <- rbind(positive, unknown)
dataset$target <- as.factor(dataset$target)
rownames(dataset) <- dataset$ensembl_gene_id
dataset$ensembl_gene_id <- NULL
saveRDS(dataset, file.path("../data/dataset.rds"))

# prediction set will be kept for the actual prediction
predictionset <- completeset[!completeset$ensembl_gene_id %in% rownames(dataset), ]
rownames(predictionset) <- predictionset$ensembl_gene_id
predictionset$ensembl_gene_id <- NULL
saveRDS(predictionset, file.path("../data/predicitionset.rds"))

## data visualisation
# principal components analysis
pca <- prcomp(dataset[1:5], scale.=TRUE)
pve <- round(pca$sdev^2/sum(pca$sdev^2) * 100, 1)
pca <- cbind(pca$x, dataset["target"])
PC12 <- ggplot(pca, aes(PC1, PC2)) +
        geom_point(aes(fill=factor(target)), shape=21, size=3, alpha=0.4) +
        xlab(paste0("PC1 (", pve[1], "% variance)")) +
        ylab(paste0("PC2 (", pve[2], "% variance)")) +
        scale_fill_manual(values=c("forestgreen", "darkviolet"), name="Target", breaks=c(0, 1), labels=c("Yes", "Unknown")) +
        theme_bw(14)
PC13 <- ggplot(pca, aes(PC1, PC3)) +
        geom_point(aes(fill=factor(target)), shape=21, size=3, alpha=0.4) +
        xlab(paste0("PC1 (", pve[1], "% variance)")) +
        ylab(paste0("PC3 (", pve[3], "% variance)")) +
        scale_fill_manual(values=c("forestgreen", "darkviolet"), name="Target", breaks=c(0, 1), labels=c("Yes", "Unknown")) +
        theme_bw(14)
PC23 <- ggplot(pca, aes(PC2, PC3)) +
        geom_point(aes(fill=factor(target)), shape=21, size=3, alpha=0.4) +
        xlab(paste0("PC2 (", pve[2], "% variance)")) +
        ylab(paste0("PC3 (", pve[3], "% variance)")) +
        scale_fill_manual(values=c("forestgreen", "darkviolet"), name="Target", breaks=c(0, 1), labels=c("Yes", "Unknown")) +
        theme_bw(14)
png("../data/PCA.png", width=18*150, height=5*150, res=150)
grid.arrange(PC12, PC13, PC23, nrow=1)
dev.off()

# hierarchical  clustering
hcfun <- function(x) hclust(x, method="ward.D2")
hmcols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(255)
sidecols <- c("forestgreen", "darkviolet")[dataset$target]
png("../data/Heatmap.png", width=8*150, height=10*150, res=150)
heatmap.2(as.matrix(dataset[1:5]), hclustfun=hcfun, Rowv=TRUE, Colv=TRUE, dendrogram="both", scale="none", col=hmcols, RowSideColors=sidecols, density.info="none", trace="none", key=TRUE, srtCol=315, adjCol=c(0, 1), labRow="", lhei=c(2,8), margins=c(12,8))
dev.off()

# t-SNE
uniqueset <- dataset[!duplicated(dataset[1:5]),]
tsne <- Rtsne(as.matrix(uniqueset[1:5]))
tsne <- cbind(tsne$Y, uniqueset[6])
names(tsne) <- c("D1", "D2", "target")
png("../data/tSNE.png", width=10*150, height=10*150, res=150)
print(
      ggplot(tsne, aes(D1, D2)) +
          geom_point(aes(fill=factor(target)), shape=21, size=3, alpha=0.4) +
          scale_fill_manual(values=c("forestgreen", "darkviolet"), name="Target", breaks=c(0, 1), labels=c("Yes", "Unknown")) +
          theme_bw(14)
)
dev.off()

