### libraries ###
library(biomaRt)
library(splitstackshape)

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

# read cttv master file
completeset <- read.csv("/GWD/bioinfo/projects/bix-analysis-stv/data/CTTV/v2.0/matrix.csv.gz")
# only use direct associations, remove known_drug and literature
completeset <- subset(completeset, Is.direct == "True", c(EnsemblId, OntologyId, genetic_association, somatic_mutation, rna_expression, affected_pathway, animal_model))
# remove lower confidence animal_model associations
completeset$animal_model[completeset$animal_model < 0.75] <- 0

# read therapeutic area mappings
tas <- read.csv("/GWD/bioinfo/projects/bix-analysis-stv/data/CTTV/v2.0/disease_tas.csv")
# split comma-separated therapeutic areas
tas <- as.data.frame(cSplit(tas, "Therapeutic.areas", sep=",", direction="long"))
# remove useless columns
tas <- tas[c("OntologyId", "Therapeutic.areas")]
# split into list
tas <- split(tas, tas$Therapeutic.areas)
tas <- lapply(tas, function(x) unique(x["OntologyId"]))
saveRDS(tas, file.path("../data/tas.rds"))

# read targetpedia data for target information
targetpedia <- read.delim("/GWD/bioinfo/projects/bix-analysis-stv/data/pharmaceutical/targetpedia/targetpedia_triples.txt")
targetpedia <- subset(targetpedia,
                            PROJECT_STATUS == "Active" |
                            PROJECT_STATUS == "Completed" |
                            PROJECT_STATUS == "Progressing" |
                            PROJECT_STATUS == "Proposed" |
                            PROJECT_STATUS == "Under Review-On Hold")
targetpedia <- getBM(
                attributes="ensembl_gene_id",
                filters=c("entrezgene", "chromosome_name", "biotype"),
                values=list(targetpedia$Target_EntrezGeneId, chr, type),
                mart=ensembl)

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


# divide workflow by therapeutic area
for (ta in names(tas)) {

    # create directory structure
    dir.create(file.path("../data", ta), showWarnings=FALSE)
    # merge with master file
    completeset.ta <- merge(completeset, tas[[ta]], by="OntologyId", all=FALSE)
    # aggregate
    completeset.ta <- aggregate(completeset.ta[3:ncol(completeset.ta)], by=list(EnsemblId=completeset.ta$EnsemblId), FUN=mean)
    # only keep protein coding genes
    completeset.ta <- merge(genes, completeset.ta, by=1, all=FALSE)
    saveRDS(completeset.ta, file.path("../data", ta, "completeset.rds"))

    # positive cases: these are targets according to targetpedia and/or pharmaprojects
    positive <- unique(rbind(targetpedia, pharmaprojects))
    positive <- completeset.ta[completeset.ta$ensembl_gene_id %in% positive$ensembl_gene_id, ]
    positive$target <- 1

    # unknown cases: it's not known whether these are targets or not
    unknown <- completeset.ta[!completeset.ta$ensembl_gene_id %in% positive$ensembl_gene_id, ]
    unknown <- unknown[sample(nrow(unknown), nrow(positive)), ]
    unknown$target <- 0

    # dataset is made of positive and unknown cases
    # will be splitted into traing and test sets
    dataset <- rbind(positive, unknown)
    dataset$target <- as.factor(dataset$target)
    rownames(dataset) <- dataset$ensembl_gene_id
    dataset$ensembl_gene_id <- NULL
    saveRDS(dataset, file.path("../data", ta, "dataset.rds"))

    # prediction set will be kept for the actual prediction
    predictionset <- completeset.ta[!completeset.ta$ensembl_gene_id %in% rownames(dataset), ]
    rownames(predictionset) <- predictionset$ensembl_gene_id
    predictionset$ensembl_gene_id <- NULL
    saveRDS(predictionset, file.path("../data", ta, "predictionset.rds"))

}
