### libraries ###
library(mlr)
library(parallelMap)

### options ###
set.seed(16)
parallelStart("multicore", 16)


### data ###
dataset <- readRDS(file.path("../data/dataset.rds"))

### setup ###
## task
classif.task <- makeClassifTask(id="TargetPred", data=dataset, target="target")

## learner
svm.lrn <- makeLearner("classif.svm")
bag.svm.lrn <- makeBaggingWrapper(svm.lrn, bw.iters=10)
bag.svm.lrn <- setPredictType(bag.svm.lrn, predict.type="prob")

## resampling description
rdesc <- makeResampleDesc("CV", iters=3)

## feature selection tuning
lrn <- makeFilterWrapper(bagged.svm.learner, fw.method="information.gain")
ps <- makeParamSet(makeDiscreteParam("fw.perc", values = c(0.5, 0.25, 0.1, 0.05, 0.01, 0.001)))
res <- tuneParams(lrn, classif.task, rdesc, par.set=ps, control=makeTuneControlGrid())

## feature selection
bag.svm.lrn <- makeFilterWrapper(bag.svm.lrn, fw.method="information.gain", fw.perc=res$x$fw.perc)

## model with resampling
mod <- resample(bag.svm.lrn, classif.task, rdesc)
