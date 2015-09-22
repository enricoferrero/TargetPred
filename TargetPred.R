library(biomaRt)
mart <- useMart("ensembl", "hsapiens_gene_ensembl")

chr <- c(1:22, "X", "Y", "MT")
type="protein_coding"

features <- getBM(
                 attributes=c("ensembl_gene_id",
                              "entrezgene",
                              "gene_biotype",
                              "transcript_length",
                              "transcript_count",
                              "percentage_gc_content"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)

domains <- getBM(
                 attributes=c("ensembl_gene_id",
                              "interpro"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)

structures <- getBM(
                 attributes=c("ensembl_gene_id",
                              "low_complexity",
                              "transmembrane_domain",
                              "signal_domain",
                              "ncoils"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)

functions <- getBM(
                 attributes=c("ensembl_gene_id",
                              "go_id"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)

pathways <- getBM(
                 attributes=c("ensembl_gene_id",
                              "reactome"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)

homologs <- getBM(
                 attributes=c("ensembl_gene_id",
                              "mmusculus_homolog_perc_id",
                              "drerio_homolog_perc_id",
                              "dmelanogaster_homolog_perc_id",
                              "celegans_homolog_perc_id",
                              "scerevisiae_homolog_perc_id"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)

ensembl <- Reduce(function(x, y) merge(x, y, by="ensembl_gene_id", all=TRUE), list(features, domains, structures, functions, pathways, homologs))
saveRDS(ensembl, "../data/ensembl.rds")
write.table(ensembl, "../data/ensembl.txt", sep="\t", quote=FALSE, rownames=FALSE, colnames=TRUE)

targetpedia <- read.delim("/GWD/bioinfo/projects/bix-analysis-stv/data/pharmaceutical/targetpedia/targetpedia_triples.txt")
pharmaprojects <- read.delim("/GWD/bioinfo/projects/bix-analysis-stv/data/pharmaceutical/pipeline/pipeline_triples.txt")
