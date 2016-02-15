### libraries ###
library(biomaRt)

### options ###
set.seed(16)
#ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host="www.ensembl.org")
chr <- c(1:22, "X", "Y", "MT")
type="protein_coding"

### functions ###
createNumericFeature <- function(genes, attribute, fun=median, mart=ensembl) {
    feat <- getBM(
                  attributes=c("ensembl_gene_id", attribute),
                  filters="ensembl_gene_id",
                  values=genes,
                  mart=mart)
    feat[feat==""] <- NA
    feat <- na.omit(feat)
    feat <- aggregate(feat[attribute], by=list(feat$ensembl_gene_id), FUN=fun)
    names(feat) <- c("ensembl_gene_id", attribute)
    feat <- merge(feat, genes, by="ensembl_gene_id", all=TRUE)
    feat <- feat[order(feat$ensembl_gene_id), , drop=FALSE]
    feat[is.na(feat)] <- 0
    feat[2:ncol(feat)] <- lapply(feat[2:ncol(feat)], as.numeric)
    return(feat)
}

createFactorFeature <- function(genes, attribute, mart=ensembl) {
    feat <- getBM(
                  attributes=c("ensembl_gene_id", attribute),
                  filters="ensembl_gene_id",
                  values=genes,
                  mart=mart)
    feat[feat==""] <- NA
    feat <- na.omit(feat)
    feat$feat <- 1
    feat <- reshape(feat, idvar="ensembl_gene_id", timevar=attribute, direction="wide")
    feat <- merge(feat, genes, by="ensembl_gene_id", all=TRUE)
    feat <- feat[order(feat$ensembl_gene_id), , drop=FALSE]
    feat[is.na(feat)] <- 0
    feat[2:ncol(feat)] <- lapply(feat[2:ncol(feat)], as.factor)
    return(feat)
}

### data ###
# get all protein coding genes
genes <- getBM(
               attributes="ensembl_gene_id",
               filters=c("chromosome_name", "biotype"),
               values=list(chr, type),
               mart=ensembl)
genes <- genes[order(genes$ensembl_gene_id), , drop=FALSE]

# create numeric and factor features
transcript_length <- createNumericFeature(genes, "transcript_length")
transcript_count <- createNumericFeature(genes, "transcript_count")
percentage_gc_content <- createNumericFeature(genes, "percentage_gc_content")
cfamiliaris_homolog_perc_id <- createNumericFeature(genes, "cfamiliaris_homolog_perc_id")
sscrofa_homolog_perc_id <- createNumericFeature(genes, "sscrofa_homolog_perc_id")
rnorvegicus_homolog_perc_id <- createNumericFeature(genes, "rnorvegicus_homolog_perc_id")
mmusculus_homolog_perc_id <- createNumericFeature(genes, "mmusculus_homolog_perc_id")
drerio_homolog_perc_id <- createNumericFeature(genes, "drerio_homolog_perc_id")
dmelanogaster_homolog_perc_id <- createNumericFeature(genes, "dmelanogaster_homolog_perc_id")
celegans_homolog_perc_id <- createNumericFeature(genes, "celegans_homolog_perc_id")
scerevisiae_homolog_perc_id <- createNumericFeature(genes, "scerevisiae_homolog_perc_id")
low_complexity <- createFactorFeature(genes, "low_complexity")
transmembrane_domain <- createFactorFeature(genes, "transmembrane_domain")
signal_domain <- createFactorFeature(genes, "signal_domain")
ncoils <- createFactorFeature(genes, "ncoils")
interpro <- createFactorFeature(genes, "interpro")
go_id <- createFactorFeature(genes, "go_id")
reactome <- createFactorFeature(genes, "reactome")

# generate complete set with all features
completeset <- cbind(genes,
                 transcript_length[2:ncol(transcript_length)],
                 transcript_count[2:ncol(transcript_count)],
                 percentage_gc_content[2:ncol(percentage_gc_content)],
                 mmusculus_homolog_perc_id[2:ncol(mmusculus_homolog_perc_id)],
                 drerio_homolog_perc_id[2:ncol(drerio_homolog_perc_id)],
                 dmelanogaster_homolog_perc_id[2:ncol(dmelanogaster_homolog_perc_id)],
                 celegans_homolog_perc_id[2:ncol(celegans_homolog_perc_id)],
                 scerevisiae_homolog_perc_id[2:ncol(scerevisiae_homolog_perc_id)],
                 low_complexity[2:ncol(low_complexity)],
                 transmembrane_domain[2:ncol(transmembrane_domain)],
                 signal_domain[2:ncol(signal_domain)],
                 ncoils[2:ncol(ncoils)],
                 interpro[2:ncol(interpro)],
                 go_id[2:ncol(go_id)],
                 reactome[2:ncol(reactome)]
                 )
names(completeset) <- gsub("[-:]", "_", names(completeset))
saveRDS(completeset, file.path("../data/completeset.rds"))

# set agenttype to small_molecule
agenttype="small_molecule"

# read tractability file
tractability <- read.delim("/GWD/bioinfo/projects/bix-analysis-stv/data/target_tractability/CombinedOutput/target_tractability_output_without_bucket5_incomplete_patents.txt", as.is=TRUE, na.strings=c("NA", ""))
tractability <- subset(tractability, Tractable.ranking.bucket <= 10, Ensembl.Gene.ID)
tractability <- unique(na.omit(tractability))

# positive cases: these are tractable targets
positive <- completeset[completeset$ensembl_gene_id %in% tractability$Ensembl.Gene.ID, ]
positive$target <- 1

# unknown cases: it's not known whether these are tractable targets or not
unknown <- completeset[!completeset$ensembl_gene_id %in% tractability$Ensembl.Gene.ID, ]
unknown <- unknown[sample(nrow(unknown), nrow(positive)), ]
unknown$target <- 0

# dataset is made of positive and unknown cases
# will be splitted into traing and test sets
dataset <- rbind(positive, unknown)
dataset$target <- as.factor(dataset$target)
rownames(dataset) <- dataset$ensembl_gene_id
dataset$ensembl_gene_id <- NULL
saveRDS(dataset, file.path(paste0("../data/dataset.", agenttype, ".rds")))

# prediction set will be kept for the actual prediction
predictionset <- completeset[!completeset$ensembl_gene_id %in% rownames(dataset), ]
rownames(predictionset) <- predictionset$ensembl_gene_id
predictionset$ensembl_gene_id <- NULL
saveRDS(predictionset, file.path(paste0("../data/predicitionset.", agenttype, ".rds")))
