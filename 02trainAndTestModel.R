### libraries ###
library(mlr)
library(parallel)
library(parallelMap)

### options ###
set.seed(16)
cv.n <- 10
bag.n <- 100

## resampling strategy
rdesc <- makeResampleDesc("CV", iters=cv.n)

## tuning control
ctrl <- makeTuneControlGrid()

## therapeutic areas
tas <- readRDS(file.path("../data/tas.rds"))

# divide workflow by therapeutic area
for (ta in names(tas)) {

    cat("Now working on", ta)

    ### model selection
    ## data
    dataset <- readRDS(file.path("../data", ta, "dataset.rds"))

    ## task
    classif.task <- makeClassifTask(id="TargetPred", data=dataset, target="target", positive="1")
    # simply remove constant features (if any)
    classif.task <- removeConstantFeatures(classif.task)
    saveRDS(classif.task, file.path("../data", ta, "classif.task.rds"))

    # features
    fv <- generateFilterValuesData(classif.task, method=c("variance", "kruskal.test", "chi.squared", "information.gain"))
    saveRDS(fv, file.path("../data/fv.rds"))
    png(file.path("../data/", ta, "Features.png"), height=10*150, width=10*150, res=150)
    print(
        plotFilterValues(fv) +
        geom_bar(aes(fill=method), stat="identity", colour="black") +
        ggtitle("") +
        theme_bw(base_size=14) +
        theme(legend.position="none") +
        theme(axis.text.x = element_text(angle=45, hjust=1))
    )
    dev.off()
    saveRDS(fv, file.path("../data", ta, "fv.rds"))
    # number of observations
    no <- getTaskSize(classif.task)
    # number of features
    nf <- getTaskNFeats(classif.task)

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
    dt.tun <- tuneParams(dt.lrn, classif.task, rdesc, par.set=dt.ps, control=ctrl)
    saveRDS(dt.tun, file.path("../data", ta, "dt.tun.rds"))
    dt.lrn <- makeLearner("classif.rpart", par.vals = list(cp=dt.tun$x$cp))
    dt.lrn <- setPredictType(dt.lrn, predict.type="prob")
    dt.lrn$id <- "Decision Tree"

    # random forest
    rf.lrn <- makeLearner("classif.randomForest")
    rf.ps <- makeParamSet(
                        makeDiscreteParam("ntree", values = c(250, 500, 1000, 2500, 5000)),
                        makeDiscreteParam("mtry", values = c(2, 3, 4, 5))
    )
    rf.tun <- tuneParams(rf.lrn, classif.task, rdesc, par.set=rf.ps, control=ctrl)
    saveRDS(rf.tun, file.path("../data", ta, "rf.tun.rds"))
    rf.lrn <- makeLearner("classif.randomForest", par.vals = list(ntree=rf.tun$x$ntree, mtry=rf.tun$x$mtry))
    rf.lrn <- setPredictType(rf.lrn, predict.type="prob")
    rf.lrn$id <- "Random Forest"

    # neural network
    nn.lrn <- makeLearner("classif.nnet", par.vals = list(MaxNWts=5000, trace=FALSE))
    nn.ps <- makeParamSet(
                        makeDiscreteParam("size", values = c(2, 3, 5, 7, 10)),
                        makeDiscreteParam("decay", values = c(0.5, 0.25, 0.1, 0))
    )
    nn.tun <- tuneParams(nn.lrn, classif.task, rdesc, par.set=nn.ps, control=ctrl)
    saveRDS(nn.tun, file.path("../data", ta, "nn.tun.rds"))
    nn.lrn <- makeLearner("classif.nnet", par.vals = list(MaxNWts=5000, trace=FALSE, size=nn.tun$x$size, decay=nn.tun$x$decay))
    nn.lrn <- makeBaggingWrapper(nn.lrn, bw.iters=bag.n)
    nn.lrn <- setPredictType(nn.lrn, predict.type="prob")
    nn.lrn$id <- "Neural Network"

    # support vector machine
    svm.lrn <- makeLearner("classif.svm")
    svm.ps <- makeParamSet(
                            makeDiscreteParam("gamma", values = 2^(-2:2)),
                            makeDiscreteParam("cost", values = 2^(-2:2))
    )
    svm.tun <- tuneParams(svm.lrn, classif.task, rdesc, par.set=svm.ps, control=ctrl)
    saveRDS(svm.tun, file.path("../data", ta, "svm.tun.rds"))
    svm.lrn <- makeLearner("classif.svm", par.vals = list(cost=svm.tun$x$cost, gamma=svm.tun$x$gamma))
    svm.lrn <- makeBaggingWrapper(svm.lrn, bw.iters=bag.n)
    svm.lrn <- setPredictType(svm.lrn, predict.type="prob")
    svm.lrn$id <- "Support Vector Machine"

    parallelStop()

    ## benchmark
    parallelStartMulticore(detectCores(), level="mlr.resample")
    lrns <- list(dt.lrn, rf.lrn, nn.lrn, svm.lrn)
    bmrk <- benchmark(lrns, classif.task, rdesc, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1))
    xlsx::write.xlsx(getBMRAggrPerformances(bmrk, as.df=TRUE), file.path("../data", ta, "Results.xlsx"), sheetName="Benchmark", row.names=FALSE, col.names=TRUE, append=FALSE)
    parallelStop()

    # boxplots of mean misclassification error
    perf <- getBMRPerformances(bmrk, as.df=TRUE)
    perf <- perf[c("learner.id", "mmce")]
    png(file.path("../data", ta, "BenchmarkMmceBoxplots.png"), height=10*150, width=10*150, res=150)
    print(
        ggplot(data=perf, aes(x=learner.id, y=mmce)) +
            geom_boxplot(aes(fill=learner.id)) +
            xlab("") +
            ylab("Misclassification error") +
            scale_fill_brewer(palette="Set1", name="Classifier") +
            theme_bw(base_size=14) +
            theme(axis.text.x=element_blank()) +
            theme(axis.ticks.x=element_blank())
    )
    dev.off()

    # boxplots of other performance measures
    perf <- getBMRPerformances(bmrk, as.df=TRUE)
    names(perf)[5:ncol(perf)] <- c("Accuracy", "AUC", "Recall/Sensitivity", "Specificity", "Precision", "F1")
    perf <- reshape(perf, varying=names(perf)[5:ncol(perf)], v.names="value", timevar="measure", times=names(perf)[5:ncol(perf)], direction="long")
    png(file.path("../data", ta, "BenchmarkOtherBoxplots.png"), height=10*150, width=10*150, res=150)
    print(
        ggplot(data=perf, aes(x=learner.id, y=value)) +
            geom_boxplot(aes(fill=learner.id)) +
            facet_wrap(~ measure, nrow=2) +
            xlab("Measure") +
            ylab("Performance") +
            scale_fill_brewer(palette="Set1", name="Classifier") +
            theme_bw(base_size=14) +
            theme(axis.text.x=element_blank()) +
            theme(axis.ticks.x=element_blank())
    )
    dev.off()

    # ROC curves
    roc <- generateThreshVsPerfData(bmrk, measures=list(fpr, tpr))
    png(file.path("../data", ta, "BenchmarkROC.png"), height=10*150, width=10*150, res=150)
    print(
        ggplot(data=roc$data, aes(x=fpr, y=tpr)) +
            geom_path(aes(colour=learner), size=1.5) +
            xlab("False positive rate") +
            ylab("True positive rate") +
            scale_colour_brewer(palette="Set1", name="Classifier") +
            theme_bw(base_size=14)
    )
    dev.off()

    # PR curves
    pr <- generateThreshVsPerfData(bmrk, measures=list(tpr, ppv))
    png(file.path("../data", ta, "BenchmarkPR.png"), height=10*150, width=10*150, res=150)
    print(
        ggplot(data=pr$data, aes(x=tpr, y=ppv)) +
            geom_path(aes(colour=learner), size=1.5) +
            xlab("Recall") +
            ylab("Precision") +
            scale_colour_brewer(palette="Set1", name="Classifier") +
            theme_bw(base_size=14)
    )
    dev.off()

    parallelStop()

    #### model testing
    bst.lrn <- rf.lrn

    ## cross-validation
    parallelStartMulticore(detectCores())
    res <- resample(bst.lrn, classif.task, rdesc, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1))
    saveRDS(res, file.path("../data", ta, "res.rds"))
    parallelStop()

    # export
    xlsx::write.xlsx(res$aggr, file.path("../data", ta, "Results.xlsx"), sheetName="CV Performance Measures", row.names=TRUE, col.names=FALSE, append=TRUE)
    xlsx::write.xlsx(as.data.frame(getConfMatrix(res$pred)), file.path("../data", ta, "Results.xlsx"), sheetName="CV Confusion Matrix", row.names=TRUE, col.names=TRUE, append=TRUE)

    ## train model
    mod <- train(bst.lrn, classif.task, subset=train.set)
    saveRDS(mod, file.path("../data", ta, "mod.rds"))

    ## evaluate performance on test set
    test.pred <- predict(mod, task=classif.task, subset=test.set)
    saveRDS(test.pred, file.path("../data", ta, "test.pred.rds"))

    # export
    xlsx::write.xlsx(performance(test.pred, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1)), file.path("../data", ta, "Results.xlsx"), sheetName="Test Performance Measures", row.names=TRUE, col.names=FALSE, append=TRUE)
    xlsx::write.xlsx(as.data.frame(getConfMatrix(test.pred)), file.path("../data", ta, "Results.xlsx"), sheetName="Test Confusion Matrix", row.names=TRUE, col.names=TRUE, append=TRUE)

    parallelStop()

    ### inference
    inf.lrn <- dt.lrn

    ## train model
    inf.mod <- train(inf.lrn, classif.task, subset=train.set)
    saveRDS(inf.mod, file.path("../data", ta, "inf.mod.rds"))
    inf.mod <- inf.mod$learner.model

    ## plot tree
    png(file.path("../data", ta, "DecisionTree.png"), height=10*300, width=10*300, res=300)
    rpart.plot::prp(inf.mod, type=2, extra=101, varlen=0, box.col=c("lightblue", "lightgreen")[inf.mod$frame$yval])
    dev.off()

}
