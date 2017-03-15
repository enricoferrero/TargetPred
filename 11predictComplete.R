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
completeset <- readRDS(file.path("../data/completeset.rds"))
nn.mod <- readRDS(file.path("../data/nn.mod.rds"))

### predict ###
complete.pred <- predict(nn.mod, newdata=completeset)
complete.pred <- setThreshold(complete.pred, 0.75)
rownames(complete.pred$data) <- completeset$ensembl_gene_id
saveRDS(complete.pred, file.path("../data/complete.pred.rds"))

### annotate ###
ann <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "external_gene_name"),
            filters=c("ensembl_gene_id", "chromosome_name", "biotype"),
            values=list(rownames(complete.pred$data), chr, type),
            mart=ensembl)
complete.predres <- merge(complete.pred$data, ann, by.x="row.names", by.y="ensembl_gene_id", all=TRUE)
names(complete.predres) <- c("Ensembl", "UnknownProb", "TargetProb", "Prediction", "Entrez", "Symbol")
complete.predres <- complete.predres[c("Ensembl", "Entrez", "Symbol", "Prediction", "TargetProb", "UnknownProb")]
write.csv(complete.predres, file.path("../data/CompleteResults.csv"), quote=FALSE, row.names=FALSE)

parallelStop()

