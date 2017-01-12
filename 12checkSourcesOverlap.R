### libraries ###
library(biomaRt)
library(Vennerable)

### options ###
set.seed(986, kind="L'Ecuyer-CMRG")
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host="mar2016.archive.ensembl.org")
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

## Open Targets data
completeset <- read.csv("../data/open_targets_association_data.csv.gz")
opentargets <- unique(subset(completeset, Is.direct == TRUE & known_drug == 1, EnsemblId))
opentargets <- merge(genes, opentargets, by=1, all=FALSE)
opentargets <- opentargets$ensembl_gene_id

## Pharmaprojects data
pharmaprojects <- read.delim("../data/pipeline_triples.txt")
pharmaprojects <- getBM(
                        attributes="ensembl_gene_id",
                        filters=c("entrezgene", "chromosome_name", "biotype"),
                        values=list(pharmaprojects$Target_EntrezGeneId, chr, type),
                        mart=ensembl)
pharmaprojects <- pharmaprojects$ensembl_gene_id

# Venn diagram
png(file.path("../data/VennSources.png"), height=10*300, width=10*300, res=300)
plot(Venn(list("Open Targets"=opentargets, "Pharmaprojects"=pharmaprojects)))
dev.off()
