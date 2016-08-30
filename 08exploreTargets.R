library(mlr)
library(biomaRt)
library(ggplot2)
set.seed(986, kind="L'Ecuyer-CMRG")

# load data
classif.task <- readRDS("../data/classif.task.rds")
res <- readRDS("../data/res.rds")
test.pred <- readRDS("../data/test.pred.rds")

## number of features and observations
nf <- getTaskNFeats(classif.task)
no <- getTaskSize(classif.task)

# training and test set
train.set <- sample(no, size = round(0.8*no))
test.set <- setdiff(seq(no), train.set)

# annotate dataset
dataset <- getTaskData(classif.task)
dataset$id <- 1:nrow(dataset)
dataset$ensembl <- rownames(dataset)

# annotate resampling results
res <- res$pred$data
res <- res[order(res$id), ]
res$id <- train.set

# annotate test results
test.pred <- test.pred$data
test.pred <- test.pred[order(test.pred$id), ]
test.pred$id <- test.set

# merge and clean
dataset.train <- merge(dataset, res, all=FALSE)
dataset.test <- merge(dataset, test.pred, all=FALSE)
dataset <- merge(dataset.train, dataset.test, all=TRUE)
dataset <- subset(dataset, truth == 1, c(ensembl, response))

# get and process pharmaprojects data
ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
chr <- c(1:22, "X", "Y", "MT")
type="protein_coding"
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
pharmaprojects.id <- getBM(
                           attributes=c("ensembl_gene_id", "entrezgene"),
                           filters=c("entrezgene", "chromosome_name", "biotype"),
                           values=list(pharmaprojects$Target_EntrezGeneId, chr, type),
                           mart=ensembl)
pharmaprojects <- merge(pharmaprojects.id, pharmaprojects, by.x="entrezgene", by.y="Target_EntrezGeneId", all=FALSE)
pharmaprojects <- unique(pharmaprojects[c("ensembl_gene_id", "GlobalStatus")])

# annotate dataset
dataset <- merge(dataset, pharmaprojects, by.x="ensembl", by.y="ensembl_gene_id", all=FALSE)

# only consider latest stage
dataset$GlobalStatus <- factor(dataset$GlobalStatus, levels=c("Preclinical", "Clinical Trial", "Phase I Clinical Trial", "Phase II Clinical Trial", "Phase III Clinical Trial", "Pre-registration", "Registered", "Launched"), ordered=TRUE)
dataset <- split(dataset, dataset$ensembl)
dataset <- lapply(dataset, transform, Stage=max(GlobalStatus))
dataset <- do.call(rbind, dataset)
dataset <- unique(dataset[c("ensembl", "Stage", "response")])

# labels for plot
levels(dataset$response) <- c("Predicted non-target", "Predicted target")

# plot
png(file.path("../data/TargetStage.png"), height=6*300, width=8*300, res=300)
print(
      ggplot(dataset, aes(Stage)) +
          geom_bar(aes(fill=response), colour="black") +
          facet_wrap(~ response, ncol=2) +
          ylab("Number of targets") +
          theme_bw(base_size=14) +
          theme(axis.text.x = element_text(angle=45, hjust=1)) +
          theme(legend.position="none") +
          scale_fill_manual(values=c("darkviolet", "forestgreen"))
)
dev.off()

# Logistic regression: are differences significatives?
logit <- summary(glm(response ~ Stage - 1, data=dataset, family="binomial"))
print(logit)
