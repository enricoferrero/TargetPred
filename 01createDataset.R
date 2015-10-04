### libraries ###
library(biomaRt)

### options ###
set.seed(16)
human <- useMart("ensembl", "hsapiens_gene_ensembl")
chr <- c(1:22, "X", "Y", "MT")
#chr <- 21
type="protein_coding"

### functions ###
createBooleanFeature <- function(genes, attribute, mart=human) {
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
}

createNumericFeature <- function(genes, attribute, fun=median, mart=human) {
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
}

### data ###

# get all protein coding genes
genes <- getBM(
               attributes="ensembl_gene_id",
               filters=c("chromosome_name", "biotype"),
               values=list(chr, type),
               mart=human)
genes <- genes[order(genes$ensembl_gene_id), , drop=FALSE]

# create numeric and boolean features
transcript_length <- createNumericFeature(genes, "transcript_length")
transcript_count <- createNumericFeature(genes, "transcript_count")
percentage_gc_content <- createNumericFeature(genes, "percentage_gc_content")
low_complexity <- createBooleanFeature(genes, "low_complexity")
transmembrane_domain <- createBooleanFeature(genes, "transmembrane_domain")
signal_domain <- createBooleanFeature(genes, "signal_domain")
ncoils <- createBooleanFeature(genes, "ncoils")
interpro <- createBooleanFeature(genes, "interpro")
go_id <- createBooleanFeature(genes, "go_id")
reactome <- createBooleanFeature(genes, "reactome")
mmusculus_homolog_perc_id <- createNumericFeature(genes, "mmusculus_homolog_perc_id")
drerio_homolog_perc_id <- createNumericFeature(genes, "drerio_homolog_perc_id")
dmelanogaster_homolog_perc_id <- createNumericFeature(genes, "dmelanogaster_homolog_perc_id")
celegans_homolog_perc_id <- createNumericFeature(genes, "celegans_homolog_perc_id")
scerevisiae_homolog_perc_id <- createNumericFeature(genes, "scerevisiae_homolog_perc_id")

# generate complete set with all features
completeset <- cbind(genes,
                 transcript_length[2:ncol(transcript_length)],
                 transcript_count[2:ncol(transcript_count)],
                 percentage_gc_content[2:ncol(percentage_gc_content)],
                 low_complexity[2:ncol(low_complexity)],
                 transmembrane_domain[2:ncol(transmembrane_domain)],
                 signal_domain[2:ncol(signal_domain)],
                 ncoils[2:ncol(ncoils)],
                 interpro[2:ncol(interpro)],
                 go_id[2:ncol(go_id)],
                 reactome[2:ncol(reactome)],
                 mmusculus_homolog_perc_id[2:ncol(mmusculus_homolog_perc_id)],
                 drerio_homolog_perc_id[2:ncol(drerio_homolog_perc_id)],
                 dmelanogaster_homolog_perc_id[2:ncol(dmelanogaster_homolog_perc_id)],
                 celegans_homolog_perc_id[2:ncol(celegans_homolog_perc_id)],
                 scerevisiae_homolog_perc_id[2:ncol(scerevisiae_homolog_perc_id)]
                 )
completeset[is.na(completeset)] <- 0

# read targetpedia data for target information
targetpedia <- read.delim("/GWD/bioinfo/projects/bix-analysis-stv/data/pharmaceutical/targetpedia/targetpedia_triples.txt")
targetpedia <- getBM(
                 attributes="ensembl_gene_id",
                 filters=c("entrezgene", "chromosome_name", "biotype"),
                 values=list(targetpedia$Target_EntrezGeneId, chr, type),
                 mart=mart)

# read pharmaprojects data for target information
pharmaprojects <- read.delim("/GWD/bioinfo/projects/bix-analysis-stv/data/pharmaceutical/pipeline/pipeline_triples.txt")
pharmaprojects <- getBM(
                        attributes="ensembl_gene_id",
                        filters=c("entrezgene", "chromosome_name", "biotype"),
                        values=list(pharmaprojects$Target_EntrezGeneId, chr, type),
                        mart=mart)

# positive cases: these are targets according to targetpedia and/or pharmaprojects
positive <- unique(rbind(targetpedia, pharmaprojects))
positive <- completeset[completeset$ensembl_gene_id %in% positive$ensembl_gene_id, ]
positive$target <- 1

# unknown cases: it's not known whether these are targets or not
unknown <- completeset[!completeset$ensembl_gene_id %in% positive$ensembl_gene_id, ]
unknown <- unknown[sample(nrow(unknown), 2*nrow(positive)), ]
unknown$target <- 0

dataset <- rbind(positive, unknown)
saveRDS(dataset, file.path("../data/dataset.rds"))

predictionset <- completeset[!completeset$ensembl_gene_id %in% datasetset$ensembl_gene_id, ]
saveRDS(predictionset, file.path("../data/predicitionset.rds")

