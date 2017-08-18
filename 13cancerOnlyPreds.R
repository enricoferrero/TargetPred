# libraries
library(biomaRt)
library(mlr)
library(parallel)
library(parallelMap)
# options
options(mc.preschedule = FALSE)
parallelStartMulticore(detectCores() - 1)
set.seed(986, kind="L'Ecuyer-CMRG")
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host="mar2016.archive.ensembl.org")
chr <- c(1:22, "X", "Y", "MT")
type="protein_coding"
cv.inner.n <- 4
cv.outer.n <- 4
bag.n <- 100

# get all protein coding genes
genes <- getBM(
               attributes="ensembl_gene_id",
               filters=c("chromosome_name", "biotype"),
               values=list(chr, type),
               mart=ensembl)
genes <- genes[order(genes$ensembl_gene_id), , drop=FALSE]


# read complete dataset
completeset <- read.csv("../data/open_targets_association_data.csv.gz")
# read therapeutic area mappings
tas <- read.delim("../data/opentargets.therapeuticareas.txt")
# merge
completeset <- merge(completeset, tas, by.x = "OntologyId", by.y = "ID")
# select neoplasm and somatic_mutation
completeset <- subset(completeset, TA_LABEL == "neoplasm" & somatic_mutation > 0)
# only use direct associations, remove known_drug and literature
completeset <- subset(completeset, Is.direct == TRUE, c(EnsemblId, OntologyId, affected_pathway, animal_model, genetic_association, rna_expression, somatic_mutation))
# remove lower confidence animal_model associations
completeset$animal_model[completeset$animal_model < 0.4] <- 0
# aggregate
completeset <- aggregate(completeset[3:ncol(completeset)], by=list(EnsemblId=completeset$EnsemblId), FUN=mean)
# merge
completeset <- merge(genes, completeset, by=1, all=FALSE)
rownames(completeset) <- completeset$ensembl_gene_id

# read pharmaprojects data for target information
pharmaprojects <- read.delim("../data/pipeline_triples.txt")
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

# positive cases: these are targets according to targetpedia and/or pharmaprojects
positive <- unique(pharmaprojects)
positive <- completeset[completeset$ensembl_gene_id %in% positive$ensembl_gene_id, ]
positive$target <- 1

# unknown cases: it's not known whether these are targets or not
unknown <- completeset[!completeset$ensembl_gene_id %in% positive$ensembl_gene_id, ]
unknown <- unknown[sample(nrow(unknown), nrow(positive)), ]
unknown$target <- 0

# dataset is made of positive and unknown cases
# will be splitted into traing and test sets
dataset <- rbind(positive, unknown)
dataset$target <- as.factor(dataset$target)
rownames(dataset) <- dataset$ensembl_gene_id
dataset$ensembl_gene_id <- NULL

# prediction set will be kept for the actual prediction
predictionset <- completeset[!completeset$ensembl_gene_id %in% rownames(dataset), ]
rownames(predictionset) <- predictionset$ensembl_gene_id
predictionset$ensembl_gene_id <- NULL

## task
classif.task <- makeClassifTask(id="TargetPred", data=dataset, target="target", positive="1")
## resampling strategy
rdesc.inner <- makeResampleDesc("CV", iters=cv.inner.n, stratify = TRUE)
rdesc.outer <- makeResampleDesc("CV", iters=cv.outer.n, stratify = TRUE)
## tuning control
ctrl <- makeTuneControlGrid()
## number of features and observations
nf <- getTaskNFeats(classif.task)
no <- getTaskSize(classif.task)
## training and test set
train.set <- sample(no, size = round(0.8*no))
test.set <- setdiff(seq(no), train.set)

## features
fv <- generateFilterValuesData(classif.task)
print(fv)

# neural network
nn.lrn <- makeLearner("classif.nnet", MaxNWts=5000, trace=FALSE)
nn.lrn <- makeBaggingWrapper(nn.lrn, bw.iters=bag.n)
nn.ps <- makeParamSet(
                      makeDiscreteParam("size", values = c(2, 3, 5, 7, 10)),
                      makeDiscreteParam("decay", values = c(0.5, 0.25, 0.1, 0))
)
nn.lrn <- makeTuneWrapper(nn.lrn, rdesc.inner, par.set=nn.ps, control=ctrl)
nn.lrn <- setPredictType(nn.lrn, predict.type="prob")
nn.lrn <- setId(nn.lrn, "Neural Network")

## benchmark/crossvalidation
bmrk <- benchmark(nn.lrn, subsetTask(classif.task, subset=train.set), rdesc.outer, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1))
print(bmrk)

## train model
nn.mod <- train(nn.lrn, classif.task, subset=train.set)

## evaluate performance on test set
nn.test.pred <- predict(nn.mod, task=classif.task, subset=test.set)
performance(nn.test.pred, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1))

