### libraries ###
library(mlr)
library(parallel)
library(parallelMap)
parallelStartMulticore(detectCores(), mc.preschdule=FALSE)

### options ###
set.seed(16)
cv.n <- 10
bag.n <- 50


### model selection
## data
dataset <- readRDS(file.path(paste0("../data/dataset.", ".rds")))

## task
classif.task <- makeClassifTask(id="TargetPred", data=dataset, target="target", positive="1")
saveRDS(classif.task, file.path(paste0("../data/classif.task.", ".rds")))

## resampling strategy
rdesc <- makeResampleDesc("CV", iters=cv.n)

## tuning control
ctrl <- makeTuneControlGrid()

## feature selection
# first, remove constant features and those that differ less than 5% from the mode (most frequent number) of the data
filtered.task <- removeConstantFeatures(classif.task, perc=0.05)
# then, perform feature selection using method of choice and keep top 250 
filtered.task <- filterFeatures(filtered.task, method="mrmr", abs=250)
saveRDS(filtered.task, file.path(paste0("../data/filtered.task.", ".rds")))
# filtered features
fv <- generateFilterValuesData(filtered.task, method="mrmr")
png(file.path(paste0("../data/FilteredFeatures.", ".png")), height=10*150, width=10*150, res=150)
print(
    plotFilterValues(fv)
)
dev.off()
saveRDS(fv, file.path(paste0("../data/fv.", ".rds")))
# number of observations
no <- getTaskSize(filtered.task)
# number of features
nf <- getTaskNFeats(filtered.task)

## training and test set
train.set <- sample(no, size = round(0.8*no))
test.set <- setdiff(seq(no), train.set)

# random forest
rf.lrn <- makeLearner("classif.randomForest")
rf.ps <- makeParamSet(
                        makeDiscreteParam("ntree", values = c(250, 500, 1000, 2500, 5000)),
                        makeDiscreteParam("mtry", values = c(round(sqrt(nf)), round(sqrt(nf)*2), round(sqrt(nf)*5), round(sqrt(nf)*10)))
)
rf.tun <- tuneParams(rf.lrn, filtered.task, rdesc, par.set=rf.ps, control=ctrl)
rf.lrn <- makeLearner("classif.randomForest", par.vals = list(ntree=rf.tun$x$ntree, mtry=rf.tun$x$mtry))
#rf.lrn <- makeBaggingWrapper(rf.lrn, bw.iters=bag.n)
rf.lrn <- setPredictType(rf.lrn, predict.type="prob")
rf.lrn$id <- "Random Forest"

# neural network
nn.lrn <- makeLearner("classif.nnet", par.vals = list(MaxNWts=5000))
nn.ps <- makeParamSet(
                        makeDiscreteParam("size", values = c(2, 3, 5, 7, 10)),
                        makeDiscreteParam("decay", values = c(0.5, 0.25, 0.1, 0))
)
nn.tun <- tuneParams(nn.lrn, filtered.task, rdesc, par.set=nn.ps, control=ctrl)
nn.lrn <- makeLearner("classif.nnet", par.vals = list(MaxNWts=5000, size=nn.tun$x$size, decay=nn.tun$x$decay))
nn.lrn <- makeBaggingWrapper(nn.lrn, bw.iters=bag.n)
nn.lrn <- setPredictType(nn.lrn, predict.type="prob")
nn.lrn$id <- "Neural Network"

# support vector machine
svm.lrn <- makeLearner("classif.svm")
svm.ps <- makeParamSet(
                        makeDiscreteParam("gamma", values = 2^(-2:2)),
                        makeDiscreteParam("cost", values = 2^(-2:2))
)
svm.tun <- tuneParams(svm.lrn, filtered.task, rdesc, par.set=svm.ps, control=ctrl)
svm.lrn <- makeLearner("classif.svm", par.vals = list(cost=svm.tun$x$cost, gamma=svm.tun$x$gamma))
svm.lrn <- makeBaggingWrapper(svm.lrn, bw.iters=bag.n)
svm.lrn <- setPredictType(svm.lrn, predict.type="prob")
svm.lrn$id <- "Support vector Machine"

## benchmark
lrns <- list(rf.lrn, nn.lrn, svm.lrn)
bmrk <- benchmark(lrns, filtered.task, rdesc)
xlsx::write.xlsx(print(bmrk), file.path(paste0("../data/Results.", ".xlsx")), sheetName="Benchmark", row.names=FALSE, col.names=TRUE, append=FALSE)

# boxplots
png(file.path(paste0("../data/BenchmarkBoxplots.", ".png")), height=10*150, width=10*150, res=150)
print(
    plotBMRBoxplots(bmrk, measure=mmce) +
        aes(colour=learner.id)
)
dev.off()

# ROC curves
roc <- generateThreshVsPerfData(bmrk, measures=list(fpr, tpr))
png(file.path(paste0("../data/BenchmarkROC.", ".png")), height=10*150, width=10*150, res=150)
print(
    #plotROCCurves(generateThreshVsPerfData(bmrk, measures=list(fpr, tpr)), diagonal=TRUE) # faceted plot
    qplot(x=fpr, y=tpr, color=learner, data=roc$data, geom="path", xlab="False positive rate", ylab="True positive rate")
)
dev.off()

# PR curves
pr <- generateThreshVsPerfData(bmrk, measures=list(tpr, ppv))
png(file.path(paste0("../data/BenchmarkPR.", ".png")), height=10*150, width=10*150, res=150)
print(
    #plotROCCurves(generateThreshVsPerfData(bmrk, measures=list(ppv, tpr)), diagonal=FALSE) # faceted plot
    qplot(x=ppv, y=tpr, color=learner, data=pr$data, geom="path", xlab="Precision", ylab="Recall")
)
dev.off()


#### model testing
## cross-validation
res <- resample(rf.lrn, filtered.task, rdesc)
saveRDS(res, file.path(paste0("../data/res.", ".rds")))

# export
xlsx::write.xlsx(res$aggr, file.path(paste0("../data/Results.", ".xlsx")), sheetName="CV Mean misclassification error", row.names=TRUE, col.names=FALSE, append=TRUE)
xlsx::write.xlsx(as.data.frame(getConfMatrix(res$pred)), file.path(paste0("../data/Results.", ".xlsx")), sheetName="CV Confusion matrix", row.names=TRUE, col.names=TRUE, append=TRUE)

## train model
mod <- train(rf.lrn, filtered.task, subset=train.set)
saveRDS(mod, file.path(paste0("../data/mod.", ".rds")))

## evaluate performance on test set
test.pred <- predict(mod, task=filtered.task, subset=test.set)
saveRDS(test.pred, file.path(paste0("../data/test.pred.", ".rds")))

# export
xlsx::write.xlsx(performance(test.pred), file.path(paste0("../data/Results.", ".xlsx")), sheetName="Test Mean misclassification error", row.names=TRUE, col.names=FALSE, append=TRUE)
xlsx::write.xlsx(as.data.frame(getConfMatrix(test.pred)), file.path(paste0("../data/Results.", ".xlsx")), sheetName="Test Confusion matrix", row.names=TRUE, col.names=TRUE, append=TRUE)

parallelStop()
