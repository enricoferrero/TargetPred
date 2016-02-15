### libraries ###
library(parallel)
library(mlr)
library(parallelMap)
library(xlsx)

### options ###
set.seed(16)
parallelStart("multicore", detectCores())
filter.method="kruskal.test"
filter.perc=0.01
cv.n <- 3
bag.n <- 20

# set agenttype to small_molecule
agenntype="small_molecule"

## data
dataset <- readRDS(file.path(paste0("../data/dataset.", agenttype, ".rds")))

## task
classif.task <- makeClassifTask(id="TargetPred", data=dataset, target="target", positive="1")

## resampling strategy
rdesc <- makeResampleDesc("CV", iters=cv.n)

## parameter set for tuning SVM
ps <- makeParamSet(
                makeDiscreteParam("cost", values = 2^(-2:2)),
                makeDiscreteParam("gamma", values = 2^(-2:2))
)

## feature selection
filtered.task <- filterFeatures(classif.task, method=filter.method, perc=filter.perc)
fv <- generateFilterValuesData(filtered.task, method=filter.method)
png(file.path(paste0("../data/FilteredFeatures.", agenttype, ".png")), height=10*150, width=10*150, res=150)
print(
    plotFilterValues(fv)
)
dev.off()
saveRDS(fv, file.path(paste0("../data/fv.", agenttype, ".rds")))

### test different algorithms ###
## learners
# k-nearest neighbour
knn.lrn <- makeLearner("classif.kknn", predict.type="prob")
knn.lrn$id <- "K-Nearest Neighbour"

# decision tree
dt.lrn <- makeLearner("classif.rpart", predict.type="prob")
dt.lrn$id <- "Decision Tree"

# random forest
rf.lrn <- makeLearner("classif.randomForest", predict.type="prob")
rf.lrn$id <- "Random Forest"

# neural network
nn.lrn <- makeLearner("classif.nnet", predict.type="prob")
nn.lrn$id <- "Neural Network"

# svm
svm.lrn <- makeLearner("classif.svm", predict.type="prob")
svm.lrn$id <- "SVM"

## tuned svm
#tun.svm.lrn <- makeLearner("classif.svm", predict.type="prob")
#tun.svm.lrn <- tuneParams(tun.svm.lrn, filtered.task, rdesc, par.set=ps, control=makeTuneControlGrid())
#tun.svm.lrn <- makeLearner("classif.svm", predict.type="prob", par.vals = list(cost=tun.svm.lrn$x$cost, gamma=tun.svm.lrn$x$gamma))
#tun.svm.lrn$id <- "Tuned SVM"

## bagged svm
#bag.svm.lrn <- makeLearner("classif.svm", predict.type="response")
#bag.svm.lrn <- makeBaggingWrapper(bag.svm.lrn, bw.iters=bag.n)
#bag.svm.lrn <- setPredictType(bag.svm.lrn, predict.type="prob")
#bag.svm.lrn$id <- "Bagged SVM"

## tuned bagged svm
#tun.bag.svm.lrn <- makeLearner("classif.svm", predict.type="response")
#tun.bag.svm.lrn <- makeBaggingWrapper(tun.bag.svm.lrn, bw.iters=bag.n)
#tun.bag.svm.lrn <- setPredictType(tun.bag.svm.lrn, predict.type="prob")
#tun.bag.svm.lrn <- tuneParams(tun.bag.svm.lrn, filtered.task, rdesc, par.set=ps, control=makeTuneControlGrid())
#tun.bag.svm.lrn <- makeLearner("classif.svm", predict.type="response", par.vals = list(cost=tun.bag.svm.lrn$x$cost, gamma=tun.bag.svm.lrn$x$gamma))
#tun.bag.svm.lrn <- makeBaggingWrapper(tun.bag.svm.lrn, bw.iters=bag.n)
#tun.bag.svm.lrn <- setPredictType(tun.bag.svm.lrn, predict.type="prob")
#tun.bag.svm.lrn$id <- "Tuned Bagged SVM"

## benchmark
lrns <- list(knn.lrn, dt.lrn, rf.lrn, nn.lrn, svm.lrn)
bmrk <- benchmark(lrns, filtered.task, rdesc)
write.xlsx(print(bmrk), file.path(paste0("../data/Results.", agenttype, ".xlsx")), sheetName="Benchmark", row.names=FALSE, col.names=TRUE, append=FALSE)
png(file.path(paste0("../data/BenchmarkBoxplots.", agenttype, ".png")), height=10*150, width=10*150, res=150)
print(
    plotBenchmarkResult(bmrk, measure=mmce) +
        aes(colour=learner.id)
)
dev.off()
png(file.path(paste0("../data/BenchmarkROC.", agenttype, ".png")), height=10*150, width=10*150, res=150)
print(
    plotROCRCurves(generateROCRCurvesData(bmrk), diagonal=TRUE)
)
dev.off()
png(file.path(paste0("../data/BenchmarkPR.", agenttype, ".png")), height=10*150, width=10*150, res=150)
print(
    plotROCRCurves(generateROCRCurvesData(bmrk, meas1="prec", meas2="rec"), diagonal=FALSE)
)
dev.off()

### use best algorithm ###
## resampling
res <- resample(rf.lrn, filtered.task, rdesc)
write.xlsx(res$aggr, file.path(paste0("../data/Results.", agenttype, ".xlsx")), sheetName="Mean misclassification error", row.names=TRUE, col.names=FALSE, append=TRUE)
getConfMatrix(res$pred)
write.xlsx(as.data.frame(getConfMatrix(res$pred)), file.path(paste0("../data/Results.", agenttype, ".xlsx")), sheetName="Confusion matrix", row.names=TRUE, col.names=TRUE, append=TRUE)
saveRDS(res, file.path(paste0("../data/res.", agenttype, ".rds")))

## train and test model
mod <- train(rf.lrn, filtered.task)
saveRDS(mod, file.path(paste0("../data/mod.", agenttype, ".rds")))

parallelStop()
