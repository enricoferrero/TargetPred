library(mlr)
library(biomaRt)
library(ggplot2)
set.seed(986, kind="L'Ecuyer-CMRG")

# load data
classif.task <- readRDS("../data/classif.task.rds")
nn.res <- readRDS("../data/nn.res.rds")
nn.test.pred <- readRDS("../data/nn.test.pred.rds")

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
dataset.train <- getTaskData(subsetTask(classif.task, subset = train.set))
dataset.train$id <- 1:nrow(dataset.train)
dataset.train$ensembl <- rownames(dataset.train)

# resampling results
nn.res <- nn.res$pred$data

# test results
nn.test.pred <- nn.test.pred$data

# merge and clean
dataset.train <- merge(dataset.train, nn.res, all=FALSE)
dataset.test <- merge(dataset, nn.test.pred, all=FALSE)
dataset <- merge(dataset.train, dataset.test, all=TRUE)
dataset <- subset(dataset, truth == 1, c(ensembl, response))

# get and process pharmaprojects data
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host="mar2016.archive.ensembl.org")
chr <- c(1:22, "X", "Y", "MT")
type="protein_coding"
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
png(file.path("../data/TargetStage.png"), height=8*300, width=10*300, res=300)
print(
      ggplot(dataset, aes(Stage)) +
          geom_bar(aes(fill=response), position=position_dodge(), colour="black") +
          #facet_wrap(~ response, ncol=2) +
          ylab("Number of targets") +
          theme_bw(base_size=24) +
          theme(axis.text.x = element_text(angle=45, hjust=1)) +
          #theme(legend.position="none") +
          theme(legend.title=element_blank()) +
          scale_fill_manual(values=c("darkviolet", "forestgreen"))
)
dev.off()

# Logistic regression: are differences significatives?
logit <- summary(glm(response ~ Stage - 1, data=dataset, family="binomial"))
print(logit)
