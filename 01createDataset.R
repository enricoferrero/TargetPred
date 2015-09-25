library(RSQLite)
library(biomaRt)

db <- dbConnect(SQLite(), dbname="../data/targetpred.db")

mart <- useMart("ensembl", "hsapiens_gene_ensembl")
chr <- c(1:22, "X", "Y", "MT")
#chr <- 21
type="protein_coding"

# genes
genes <- getBM(
               attributes="ensembl_gene_id",
               filters=c("chromosome_name", "biotype"),
               values=list(chr, type),
               mart=mart)
genes[genes==""] <- NA
dbWriteTable(db, "genes", genes, overwrite=TRUE)

# features
features <- getBM(
                 attributes=c("ensembl_gene_id",
                              "gene_biotype",
                              "transcript_length",
                              "transcript_count",
                              "percentage_gc_content"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)
features[features==""] <- NA
dbWriteTable(db, "features", features, overwrite=TRUE)

# domains
domains <- getBM(
                 attributes=c("ensembl_gene_id",
                              "interpro"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)
domains[domains==""] <- NA
dbWriteTable(db, "domains", domains, overwrite=TRUE)

# structures
structures <- getBM(
                 attributes=c("ensembl_gene_id",
                              "low_complexity",
                              "transmembrane_domain",
                              "signal_domain",
                              "ncoils"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)
structures[structures==""] <- NA
dbWriteTable(db, "structures", structures, overwrite=TRUE)

# functions
functions <- getBM(
                 attributes=c("ensembl_gene_id",
                              "go_id"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)
functions[functions==""] <- NA
dbWriteTable(db, "functions", functions, overwrite=TRUE)

# pathways
pathways <- getBM(
                 attributes=c("ensembl_gene_id",
                              "reactome"),
                 filters=c("chromosome_name", "biotype"),
                 values=list(chr, type),
                 mart=mart)
pathways[pathways==""] <- NA
dbWriteTable(db, "pathways", pathways, overwrite=TRUE)

# homologs
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
homologs[homologs==""] <- NA
dbWriteTable(db, "homologs", homologs, overwrite=TRUE)

# targets
targetpedia <- read.delim("/GWD/bioinfo/projects/bix-analysis-stv/data/pharmaceutical/targetpedia/targetpedia_triples.txt")
targetpedia <- getBM(
                 attributes="ensembl_gene_id",
                 filters=c("entrezgene", "chromosome_name", "biotype"),
                 values=list(targetpedia$Target_EntrezGeneId, chr, type),
                 mart=mart)

pharmaprojects <- read.delim("/GWD/bioinfo/projects/bix-analysis-stv/data/pharmaceutical/pipeline/pipeline_triples.txt")
pharmaprojects <- getBM(
                 attributes="ensembl_gene_id",
                 filters=c("entrezgene", "chromosome_name", "biotype"),
                 values=list(pharmaprojects$Target_EntrezGeneId, chr, type),
                 mart=mart)

targets <- unique(rbind(targetpedia, pharmaprojects))
targets$target <- 1
targets <- merge(targets, genes, by="ensembl_gene_id", all=TRUE)
targets[is.na(targets)] <- 0
dbWriteTable(db, "targets", targets, overwrite=TRUE)

dataset <- dbGetQuery(db, "
                      select * from genes g
                      left join features f1 on g.ensembl_gene_id = f1.ensembl_gene_id
                      left join domains d on g.ensembl_gene_id = d.ensembl_gene_id
                      left join structures s on g.ensembl_gene_id = s.ensembl_gene_id
                      left join functions f2 on g.ensembl_gene_id = f2.ensembl_gene_id
                      left join pathways p on g.ensembl_gene_id = p.ensembl_gene_id
                      left join homologs h on g.ensembl_gene_id = h.ensembl_gene_id
                      left join targets t on g.ensembl_gene_id = t.ensembl_gene_id
                      ;")
# remove duplicated columns
dataset <- dataset[!duplicated(lapply(dataset, summary))]
dbWriteTable(db, "dataset", dataset, overwrite=TRUE)
