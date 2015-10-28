### libraries ###
library(mlr)
library(parallelMap)

### options ###
set.seed(16)
parallelStart("multicore", 16)
filter.method="kruskal.test"
cv.n <- 3
bag.n <- 20

### data ###
dataset <- readRDS(file.path("../data/dataset.rds"))

### setup ###
## task
classif.task <- makeClassifTask(id="TargetPred", data=dataset, target="target", positive="1")

## feature selection
filtered.task <- filterFeatures(classif.task, method=filter.method, threshold=1)

## resampling strategy
rdesc <- makeResampleDesc("CV", iters=cv.n)

## parameter set for tuning
ps <- makeParamSet(
                   makeDiscreteParam("cost", values = 2^(-2:2)),
                   makeDiscreteParam("gamma", values = 2^(-2:2)),
                   makeDiscreteParam("fw.perc", values = c(0.1, 0.05, 0.01))
)

## learners
# svm
svm.lrn <- makeLearner("classif.svm", predict.type="prob")
svm.lrn <- makeFilterWrapper(svm.lrn, fw.method=filter.method, fw.perc=0.1)
svm.lrn$id <- "SVM"

# tuned svm
lrn <- makeLearner("classif.svm", predict.type="prob")
lrn <- makeFilterWrapper(lrn, fw.method=filter.method)
tun <- tuneParams(lrn, filtered.task, rdesc, par.set=ps, control=makeTuneControlGrid())
tun.svm.lrn <- makeLearner("classif.svm", predict.type="prob", par.vals = list(cost=tun$x$cost, gamma=tun$x$gamma))
tun.svm.lrn <- makeFilterWrapper(tun.svm.lrn, fw.method=filter.method, fw.perc=tun$x$fw.perc)
tun.svm.lrn$id <- "Tuned SVM"

# bagged svm
lrn <- makeLearner("classif.svm", predict.type="response")
bag.svm.lrn <- makeBaggingWrapper(lrn, bw.iters=bag.n)
bag.svm.lrn <- setPredictType(bag.svm.lrn, predict.type="prob")
bag.svm.lrn <- makeFilterWrapper(bag.svm.lrn, fw.method=filter.method, fw.perc=0.1)
bag.svm.lrn$id <- "Bagged SVM"

# tuned bagged svm
lrn <- makeLearner("classif.svm", predict.type="response")
lrn <- makeBaggingWrapper(lrn, bw.iters=bag.n)
lrn <- setPredictType(lrn, predict.type="prob")
lrn <- makeFilterWrapper(lrn, fw.method=filter.method)
tun <- tuneParams(lrn, filtered.task, rdesc, par.set=ps, control=makeTuneControlGrid())
tun.bag.svm.lrn <- makeLearner("classif.svm", predict.type="response", par.vals = list(cost=tun$x$cost, gamma=tun$x$gamma))
tun.bag.svm.lrn <- makeBaggingWrapper(tun.bag.svm.lrn, bw.iters=bag.n)
tun.bag.svm.lrn <- setPredictType(tun.bag.svm.lrn, predict.type="prob")
tun.bag.svm.lrn <- makeFilterWrapper(tun.bag.svm.lrn, fw.method=filter.method, fw.perc=tun$x$fw.perc)
tun.bag.svm.lrn$id <- "Tuned Bagged SVM"

## benchmark
lrns <- list(svm.lrn, tun.svm.lrn, bag.svm.lrn, tun.bag.svm.lrn)
bmrk <- benchmark(lrns, filtered.task, rdesc)
print(bmrk)
png(file.path("../data/benchmarkBoxplots.png"), height=10*150, width=10*150, res=150)
plotBenchmarkResult(bmrk, measure=mmce) +
    aes(colour=learner.id)
dev.off()
png(file.path("../data/benchmarkROC.png"), height=10*150, width=10*150, res=150)
plotROCRCurves(generateROCRCurvesData(bmrk), diagonal=TRUE)
dev.off()

## resampling
res <- resample(tun.bag.svm.lrn, filtered.task, rdesc)
print(res$aggr)
getConfMatrix(res$pred)
saveRDS(res, file.path("../data/res.rds"))

## train and test model
mod <- train(tun.bag.svm.lrn, filtered.task)
saveRDS(mod, file.path("../data/mod.rds"))
