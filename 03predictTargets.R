### libraries ###
library(mlr)
library(parallel)
library(parallelMap)
library(biomaRt)

### options ###
set.seed(16)
parallelStartMulticore(detectCores())
ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
chr <- c(1:22, "X", "Y", "MT")
type="protein_coding"

### data ###
dataset <- readRDS(file.path("../data/dataset.rds"))
predictionset <- readRDS(file.path("../data/predicitionset.rds"))
mod <- readRDS(file.path("../data/mod.rds"))

### predict ###
pred <- predict(mod, newdata=predictionset)
saveRDS(pred, file.path("../data/pred.rds"))

### annotate ###
ann <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "external_gene_name"),
            filters=c("ensembl_gene_id", "chromosome_name", "biotype"),
            values=list(rownames(pred$data), chr, type),
            mart=ensembl)
predres <- merge(ann, pred$data, by.x="ensembl_gene_id", by.y="row.names", all=TRUE)
predres <- merge(predres, dataset, by.x="ensembl_gene_id", by.y="row.names", all=TRUE)
predres <- predres[c("ensembl_gene_id", "entrezgene", "external_gene_name", "target", "response", "prob.0", "prob.1")]
names(predres) <- c("Ensembl", "Entrez", "Symbol", "Value", "Prediction", "Prob0", "Prob1")
write.csv(predres, file.path("../data/PredicitonResults.csv"), quote=FALSE, row.names=FALSE)

parallelStop()
