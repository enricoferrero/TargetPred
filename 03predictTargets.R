### libraries ###
library(mlr)
library(parallel)
library(parallelMap)
library(biomaRt)

### options ###
set.seed(16)
parallelStart("multicore", detectCores())
#ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host="www.ensembl.org")
chr <- c(1:22, "X", "Y", "MT")
type="protein_coding"

### data ###
# set agenttype to small_molecule
agenntype="small_molecule"

predictionset <- readRDS(file.path(paste0("../data/predicitionset.", agenttype, ".rds")))
mod <- readRDS(file.path(paste0("../data/mod.", agenttype, ".rds")))

### predict ###
pred <- predict(mod, newdata=predictionset)
saveRDS(pred, file.path(paste0("../data/pred.", agenttype, ".rds")))

### annotate ###
ann <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "external_gene_name"),
            filters=c("ensembl_gene_id", "chromosome_name", "biotype"),
            values=list(rownames(pred$data), chr, type),
            mart=ensembl)
predres <- merge(pred$data, ann, by.x="row.names", by.y="ensembl_gene_id", all=TRUE)
names(predres) <- c("Ensembl", "UnknownProb", "TargetProb", "Prediction", "Entrez", "Symbol")
predres <- predres[c("Ensembl", "Entrez", "Symbol", "Prediction", "TargetProb", "UnknownProb")]
write.csv(predres, file.path(paste0("../data/PredicitonResults.", agenttype, ".csv")), quote=FALSE, row.names=FALSE)

parallelStop()
