### libraries ###
library(mlr)
library(parallel)
library(parallelMap)

### options ###
set.seed(16)
cv.n <- 10
bag.n <- 100


### model selection
## data
dataset <- readRDS(file.path("../data/dataset.rds"))

## task
classif.task <- makeClassifTask(id="TargetPred", data=dataset, target="target", positive="1")
saveRDS(classif.task, file.path("../data/classif.task.rds"))

## resampling strategy
rdesc <- makeResampleDesc("CV", iters=cv.n)

## tuning control
ctrl <- makeTuneControlGrid()

## feature selection
# first, remove constant features and those that differ less than 1% from the mode (most frequent number) of the data
filtered.task <- removeConstantFeatures(classif.task, perc=0.01)
saveRDS(filtered.task, file.path("../data/filtered.task.rds"))
# then, perform feature selection using method of choice and keep top 250
filtered.task <- filterFeatures(filtered.task, method="mrmr", abs=250)
saveRDS(filtered.task, file.path("../data/filtered.task.rds"))
# filtered features
fv <- generateFilterValuesData(filtered.task, method="mrmr")
png(file.path("../data/FilteredFeatures.png"), height=10*150, width=10*150, res=150)
print(
    ggplot(data=fv$data, aes(x=reorder(name, -information.gain), y=information.gain)) +
        geom_bar(stat="identity", colour="black", fill="#FF6600") +
        xlab("") +
        ylab("Information gain") +
        theme_bw(base_size=16) +
        theme(axis.text.x=element_text(angle=315, hjust=0))
)
dev.off()
saveRDS(fv, file.path("../data/fv.rds"))
# number of observations
no <- getTaskSize(filtered.task)
# number of features
nf <- getTaskNFeats(filtered.task)

## training and test set
train.set <- sample(no, size = round(0.8*no))
test.set <- setdiff(seq(no), train.set)

## tuning
parallelStartMulticore(detectCores())

# decision tree
dt.lrn <- makeLearner("classif.rpart")
dt.ps <- makeParamSet(
                      makeDiscreteParam("cp", values = c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05))
                      )
dt.tun <- tuneParams(dt.lrn, filtered.task, rdesc, par.set=dt.ps, control=ctrl)
dt.lrn <- makeLearner("classif.rpart", par.vals = list(cp=dt.tun$x$cp))
dt.lrn <- makeBaggingWrapper(dt.lrn, bw.iters=bag.n)
dt.lrn <- setPredictType(dt.lrn, predict.type="prob")
dt.lrn$id <- "Bagged Decision Tree"

# random forest
rf.lrn <- makeLearner("classif.randomForest")
rf.ps <- makeParamSet(
                      makeDiscreteParam("ntree", values = c(250, 500, 1000, 2500, 5000)),
                      makeDiscreteParam("mtry", values = c(round(sqrt(nf)), round(sqrt(nf)*2), round(sqrt(nf)*5), round(sqrt(nf)*10)))
)
rf.tun <- tuneParams(rf.lrn, filtered.task, rdesc, par.set=rf.ps, control=ctrl)
saveRDS(rf.tun, file.path("../data/rf.tun.rds"))
rf.lrn <- makeLearner("classif.randomForest", par.vals = list(ntree=rf.tun$x$ntree, mtry=rf.tun$x$mtry))
rf.lrn <- setPredictType(rf.lrn, predict.type="prob")
rf.lrn$id <- "Random Forest"

# neural network
nn.lrn <- makeLearner("classif.nnet", par.vals = list(MaxNWts=5000, trace=FALSE))
nn.ps <- makeParamSet(
                      makeDiscreteParam("size", values = c(2, 3, 5, 7, 10)),
                      makeDiscreteParam("decay", values = c(0.5, 0.25, 0.1, 0))
)
nn.tun <- tuneParams(nn.lrn, filtered.task, rdesc, par.set=nn.ps, control=ctrl)
saveRDS(nn.tun, file.path("../data/nn.tun.rds"))
nn.lrn <- makeLearner("classif.nnet", par.vals = list(MaxNWts=5000, trace=FALSE, size=nn.tun$x$size, decay=nn.tun$x$decay))
nn.lrn <- makeBaggingWrapper(nn.lrn, bw.iters=bag.n)
nn.lrn <- setPredictType(nn.lrn, predict.type="prob")
nn.lrn$id <- "Bagged Neural Network"

# support vector machine
svm.lrn <- makeLearner("classif.svm")
svm.ps <- makeParamSet(
                        makeDiscreteParam("gamma", values = 2^(-2:2)),
                        makeDiscreteParam("cost", values = 2^(-2:2))
)
svm.tun <- tuneParams(svm.lrn, filtered.task, rdesc, par.set=svm.ps, control=ctrl)
saveRDS(svm.tun, file.path("../data/svm.tun.rds"))
svm.lrn <- makeLearner("classif.svm", par.vals = list(cost=svm.tun$x$cost, gamma=svm.tun$x$gamma))
svm.lrn <- makeBaggingWrapper(svm.lrn, bw.iters=bag.n)
svm.lrn <- setPredictType(svm.lrn, predict.type="prob")
svm.lrn$id <- "Bagged Support Vector Machine"

parallelStop()

## benchmark
parallelStartMulticore(detectCores(), level="mlr.resample")
lrns <- list(rf.lrn, nn.lrn, svm.lrn)
bmrk <- benchmark(lrns, filtered.task, rdesc, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1))
xlsx::write.xlsx(getBMRAggrPerformances(bmrk, as.df=TRUE), file.path("../data/Results.xlsx"), sheetName="Benchmark", row.names=FALSE, col.names=TRUE, append=FALSE)
parallelStop()

# boxplots of mean misclassification error
perf <- getBMRPerformances(bmrk, as.df=TRUE)
perf <- perf[c("learner.id", "mmce")]
png(file.path("../data/BenchmarkMmceBoxplots.png"), height=10*150, width=10*150, res=150)
print(
      ggplot(data=perf, aes(x=learner.id, y=mmce)) +
        geom_boxplot(aes(fill=learner.id)) +
        xlab("") +
        ylab("Misclassification error") +
        scale_fill_brewer(palette="Set1", name="Algorithm") +
        theme_bw(base_size=14) +
        theme(
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()
        )
)
dev.off()

# boxplots of other performance measures
perf <- getBMRPerformances(bmrk, as.df=TRUE)
names(perf)[5:ncol(perf)] <- c("Accuracy", "AUC", "Recall/Sensitivity", "Specificity", "Precision", "F1")
perf <- reshape(perf, varying=names(perf)[5:ncol(perf)], v.names="value", timevar="measure", times=names(perf)[5:ncol(perf)], direction="long")
png(file.path("../data/BenchmarkOtherBoxplots.png"), height=10*150, width=10*150, res=150)
print(
      ggplot(data=perf, aes(x=learner.id, y=value)) +
          geom_boxplot(aes(fill=learner.id)) +
          facet_wrap(~ measure, nrow=2) +
          xlab("Measure") +
          ylab("Performance") +
          scale_fill_brewer(palette="Set1", name="Algorithm") +
          theme_bw(base_size=14) +
          theme(
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()
        )
)
dev.off()

# ROC curves
roc <- generateThreshVsPerfData(bmrk, measures=list(fpr, tpr))
png(file.path("../data/BenchmarkROC.png"), height=10*150, width=10*150, res=150)
print(
      ggplot(data=roc$data, aes(x=fpr, y=tpr)) +
          geom_path(aes(color=learner), size=1.5) +
          xlab("False positive rate") +
          ylab("True positive rate") +
          theme_bw(base_size=14) +
          scale_colour_brewer(palette="Set1")
)
dev.off()

# PR curves
pr <- generateThreshVsPerfData(bmrk, measures=list(tpr, ppv))
png(file.path("../data/BenchmarkPR.png"), height=10*150, width=10*150, res=150)
print(
      ggplot(data=pr$data, aes(x=tpr, y=ppv)) +
          geom_path(aes(color=learner), size=1.5) +
          xlab("Recall") +
          ylab("Precision") +
          theme_bw(base_size=14) +
          scale_colour_brewer(palette="Set1")
)
dev.off()

parallelStop()

#### model testing
## cross-validation
parallelStartMulticore(detectCores())
res <- resample(rf.lrn, filtered.task, rdesc, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1))
saveRDS(res, file.path("../data/res.rds"))
parallelStop()

# export
xlsx::write.xlsx(res$aggr, file.path("../data/Results.xlsx"), sheetName="CV Performance Measures", row.names=TRUE, col.names=FALSE, append=TRUE)
xlsx::write.xlsx(as.data.frame(getConfMatrix(res$pred)), file.path("../data/Results.xlsx"), sheetName="CV Confusion Matrix", row.names=TRUE, col.names=TRUE, append=TRUE)

## train model
mod <- train(rf.lrn, filtered.task, subset=train.set)
saveRDS(mod, file.path("../data/mod.rds"))

## evaluate performance on test set
test.pred <- predict(mod, task=filtered.task, subset=test.set)
saveRDS(test.pred, file.path("../data/test.pred.rds"))

# export
xlsx::write.xlsx(performance(test.pred, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1)), file.path("../data/Results.xlsx"), sheetName="Test Performance Measures", row.names=TRUE, col.names=FALSE, append=TRUE)
xlsx::write.xlsx(as.data.frame(getConfMatrix(test.pred)), file.path("../data/Results.xlsx"), sheetName="Test Confusion Matrix", row.names=TRUE, col.names=TRUE, append=TRUE)

parallelStop()
