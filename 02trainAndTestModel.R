### libraries ###
library(parallel)
library(mlr)
library(parallelMap)

### options ###
set.seed(16)
parallelStart("multicore", detectCores())
filter.method="kruskal.test"
filter.perc=0.05
cv.n <- 3
bag.n <- 20

# separate small molecules and antibodies
for (agenttype in c("small_molecule", "antibody")) {

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
    plotFilterValues(fv)
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

    # svm
    svm.lrn <- makeLearner("classif.svm", predict.type="prob")
    svm.lrn$id <- "SVM"

    ## benchmark
    all.lrns <- list(knn.lrn, dt.lrn, rf.lrn, svm.lrn)
    all.bmrk <- benchmark(all.lrns, filtered.task, rdesc)
    print(all.bmrk)
    png(file.path(paste0("../data/AllBenchmarkBoxplots.", agenttype, ".png")), height=10*150, width=10*150, res=150)
    plotBenchmarkResult(all.bmrk, measure=mmce) +
        aes(colour=learner.id)
    dev.off()
    png(file.path(paste0("../data/AllBenchmarkROC.", agenttype, ".png")), height=10*150, width=10*150, res=150)
    plotROCRCurves(generateROCRCurvesData(all.bmrk), diagonal=TRUE)
    dev.off()

    ### tune best algorithm ###
    ## learners
    # svm
    svm.lrn <- makeLearner("classif.svm", predict.type="prob")
    svm.lrn$id <- "SVM"

    # tuned svm
    tun.svm.lrn <- makeLearner("classif.svm", predict.type="prob")
    tun.svm.lrn <- tuneParams(tun.svm.lrn, filtered.task, rdesc, par.set=ps, control=makeTuneControlGrid())
    tun.svm.lrn <- makeLearner("classif.svm", predict.type="prob", par.vals = list(cost=tun.svm.lrn$x$cost, gamma=tun.svm.lrn$x$gamma))
    tun.svm.lrn$id <- "Tuned SVM"

    # bagged svm
    bag.svm.lrn <- makeLearner("classif.svm", predict.type="response")
    bag.svm.lrn <- makeBaggingWrapper(bag.svm.lrn, bw.iters=bag.n)
    bag.svm.lrn <- setPredictType(bag.svm.lrn, predict.type="prob")
    bag.svm.lrn$id <- "Bagged SVM"

    # tuned bagged svm
    tun.bag.svm.lrn <- makeLearner("classif.svm", predict.type="response")
    tun.bag.svm.lrn <- makeBaggingWrapper(tun.bag.svm.lrn, bw.iters=bag.n)
    tun.bag.svm.lrn <- setPredictType(tun.bag.svm.lrn, predict.type="prob")
    tun.bag.svm.lrn <- tuneParams(tun.bag.svm.lrn, filtered.task, rdesc, par.set=ps, control=makeTuneControlGrid())
    tun.bag.svm.lrn <- makeLearner("classif.svm", predict.type="response", par.vals = list(cost=tun.bag.svm.lrn$x$cost, gamma=tun.bag.svm.lrn$x$gamma))
    tun.bag.svm.lrn <- makeBaggingWrapper(tun.bag.svm.lrn, bw.iters=bag.n)
    tun.bag.svm.lrn <- setPredictType(tun.bag.svm.lrn, predict.type="prob")
    tun.bag.svm.lrn$id <- "Tuned Bagged SVM"

    ## benchmark
    best.lrns <- list(svm.lrn, tun.svm.lrn, bag.svm.lrn, tun.bag.svm.lrn)
    best.bmrk <- benchmark(best.lrns, filtered.task, rdesc)
    print(best.bmrk)
    png(file.path(paste0("../data/BestBenchmarkBoxplots.", agenttype, ".png")), height=10*150, width=10*150, res=150)
    plotBenchmarkResult(best.bmrk, measure=mmce) +
        aes(colour=learner.id)
    dev.off()
    png(file.path(paste0("../data/BestBenchmarkROC.", agenttype, ".png")), height=10*150, width=10*150, res=150)
    plotROCRCurves(generateROCRCurvesData(best.bmrk), diagonal=TRUE)
    dev.off()

    ### use best algorithm ###
    ## resampling
    res <- resample(tun.bag.svm.lrn, filtered.task, rdesc)
    print(res$aggr)
    getConfMatrix(res$pred)
    saveRDS(res, file.path(paste0("../data/res.", agenttype, ".rds")))

    ## train and test model
    mod <- train(tun.bag.svm.lrn, filtered.task)
    saveRDS(mod, file.path(paste0("../data/mod.", agenttype, ".rds")))

}

parallelStop()
