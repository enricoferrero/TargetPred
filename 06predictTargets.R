### libraries ###
library(mlr)
library(parallel)
library(parallelMap)
library(biomaRt)

### options ###
set.seed(986, kind="L'Ecuyer-CMRG")
parallelStartMulticore(detectCores())
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host="mar2016.archive.ensembl.org")
chr <- c(1:22, "X", "Y", "MT")
type="protein_coding"

### data ###
predictionset <- readRDS(file.path("../data/predictionset.rds"))
gbm.mod <- readRDS(file.path("../data/gbm.mod.rds"))

### predict ###
pred <- predict(gbm.mod, newdata=predictionset)
pred <- setThreshold(pred, 0.9)
saveRDS(pred, file.path("../data/pred.rds"))

### annotate ###
ann <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "external_gene_name"),
            filters=c("ensembl_gene_id", "chromosome_name", "biotype"),
            values=list(rownames(pred$data), chr, type),
            mart=ensembl)
predres <- merge(pred$data, ann, by.x="row.names", by.y="ensembl_gene_id", all=TRUE)
names(predres) <- c("Ensembl", "UnknownProb", "TargetProb", "Prediction", "Entrez", "Symbol")
predres <- predres[c("Ensembl", "Entrez", "Symbol", "Prediction", "TargetProb", "UnknownProb")]
write.csv(predres, file.path("../data/PredictionResults.csv"), quote=FALSE, row.names=FALSE)

parallelStop()
