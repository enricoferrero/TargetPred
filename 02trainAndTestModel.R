### libraries ###
library(parallel)
library(mlr)
library(parallelMap)
library(xlsx)

### options ###
set.seed(16)
parallelStart("multicore", detectCores())
filter.method="mrmr"
filter.perc=0.01
cv.n <- 10
bag.n <- 50

# separate small molecules and antibodies
for (agenttype in c("small_molecule", "antibody")) {

    ## data
    dataset <- readRDS(file.path(paste0("../data/dataset.", agenttype, ".rds")))

    ## task
    classif.task <- makeClassifTask(id="TargetPred", data=dataset, target="target", positive="1")
    saveRDS(classif.task, file.path(paste0("../data/classif.task.", agenttype, ".rds")))

    ## resampling strategy
    rdesc <- makeResampleDesc("CV", iters=cv.n)

    ## parameter set for tuning SVM
    ps <- makeParamSet(
                    makeDiscreteParam("cost", values = 2^(-2:2)),
                    makeDiscreteParam("gamma", values = 2^(-2:2))
    )

    ## feature selection
    filtered.task <- filterFeatures(classif.task, method=filter.method, perc=filter.perc)
    saveRDS(filtered.task, file.path(paste0("../data/filtered.task.", agenttype, ".rds")))
    fv <- generateFilterValuesData(filtered.task, method=filter.method)
    png(file.path(paste0("../data/FilteredFeatures.", agenttype, ".png")), height=10*150, width=10*150, res=150)
    print(
        plotFilterValues(fv)
    )
    dev.off()
    saveRDS(fv, file.path(paste0("../data/fv.", agenttype, ".rds")))

    ## training and test set
    n <- getTaskSize(filtered.task)
    train.set <- sample(n, size = round(0.8*n))
    test.set <- setdiff(seq(n), train.set)

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

    # boxplots
    png(file.path(paste0("../data/BenchmarkBoxplots.", agenttype, ".png")), height=10*150, width=10*150, res=150)
    print(
        plotBMRBoxplots(bmrk, measure=mmce) +
            aes(colour=learner.id)
    )
    dev.off()

    # ROC curves
    roc <- generateThreshVsPerfData(bmrk, measures=list(fpr, tpr))
    png(file.path(paste0("../data/BenchmarkROC.", agenttype, ".png")), height=10*150, width=10*150, res=150)
    print(
        #plotROCCurves(generateThreshVsPerfData(bmrk, measures=list(fpr, tpr)), diagonal=TRUE) # faceted plot
        qplot(x=fpr, y=tpr, color=learner, data=roc$data, geom="path", xlab="False positive rate", ylab="True positive rate")
    )
    dev.off()

    # PR curves
    pr <- generateThreshVsPerfData(bmrk, measures=list(tpr, ppv))
    png(file.path(paste0("../data/BenchmarkPR.", agenttype, ".png")), height=10*150, width=10*150, res=150)
    print(
        #plotROCCurves(generateThreshVsPerfData(bmrk, measures=list(ppv, tpr)), diagonal=FALSE) # faceted plot
        qplot(x=ppv, y=tpr, color=learner, data=pr$data, geom="path", xlab="Precision", ylab="Recall")
    )
    dev.off()


    ### use best algorithm

    ## cross-validation
    res <- resample(rf.lrn, filtered.task, rdesc)
    saveRDS(res, file.path(paste0("../data/res.", agenttype, ".rds")))

    # export
    write.xlsx(res$aggr, file.path(paste0("../data/Results.", agenttype, ".xlsx")), sheetName="Mean misclassification error", row.names=TRUE, col.names=FALSE, append=TRUE)
    write.xlsx(as.data.frame(getConfMatrix(res$pred)), file.path(paste0("../data/Results.", agenttype, ".xlsx")), sheetName="Confusion matrix", row.names=TRUE, col.names=TRUE, append=TRUE)

    ## train model
    mod <- train(rf.lrn, filtered.task, subset=train.set)
    saveRDS(mod, file.path(paste0("../data/mod.", agenttype, ".rds")))

    ## evaluate performance on test set
    test.pred <- predict(mod, task=filtered.task, subset=test.set)
    saveRDS(test.pred, file.path(paste0("../data/test.pred.", agenttype, ".rds")))

    # export
    write.xlsx(performance(test.pred), file.path(paste0("../data/Results.", agenttype, ".xlsx")), sheetName="Mean misclassification error", row.names=TRUE, col.names=FALSE, append=TRUE)
    write.xlsx(as.data.frame(getConfMatrix(test.pred)), file.path(paste0("../data/Results.", agenttype, ".xlsx")), sheetName="Confusion matrix", row.names=TRUE, col.names=TRUE, append=TRUE)

}

parallelStop()
