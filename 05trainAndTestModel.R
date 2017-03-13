options(mc.preschedule = FALSE)
### libraries ###
library(mlr)
library(parallel)
library(parallelMap)
parallelStartMulticore(detectCores() - 1)
library(Vennerable)

### options ###
set.seed(986, kind="L'Ecuyer-CMRG")
cv.inner.n <- 4
cv.outer.n <- 4
bag.n <- 100

## data
dataset <- readRDS(file.path("../data/dataset.rds"))

## task
classif.task <- makeClassifTask(id="TargetPred", data=dataset, target="target", positive="1")
# simply remove constatn features (if any)
classif.task <- removeConstantFeatures(classif.task)
saveRDS(classif.task, file.path("../data/classif.task.rds"))
## resampling strategy
rdesc.inner <- makeResampleDesc("CV", iters=cv.inner.n, stratify = TRUE)
rdesc.outer <- makeResampleDesc("CV", iters=cv.outer.n, stratify = TRUE)
## tuning control
ctrl <- makeTuneControlGrid()
## number of features and observations
nf <- getTaskNFeats(classif.task)
no <- getTaskSize(classif.task)
## training and test set
train.set <- sample(no, size = round(0.8*no))
test.set <- setdiff(seq(no), train.set)

## features
fv <- generateFilterValuesData(classif.task, method=c("chi.squared", "information.gain"))
saveRDS(fv, file.path("../data/fv.rds"))
png(file.path("../data/Features.png"), height=8*300, width=8*300, res=300)
print(
    plotFilterValues(fv, sort="none") +
    geom_bar(aes(fill=method), stat="identity", colour="black") +
    ggtitle("") +
    theme_bw(base_size=24) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    scale_fill_manual(values=c("tomato", "skyblue"))
)
dev.off()

## classifiers
## these are created with makeTuneWrapper() so that during the training a nested cross-validation procedure is used: parameters are tuned in an inner CV loop and the performance is estimated in the outer loop.

# decision tree
dt.lrn <- makeLearner("classif.rpart")
dt.ps <- makeParamSet(
                      makeDiscreteParam("cp", values = c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05))
                      )
dt.lrn <- makeTuneWrapper(dt.lrn, rdesc.inner, par.set=dt.ps, control=ctrl)
dt.lrn <- setPredictType(dt.lrn, predict.type="prob")
dt.lrn <- setId(dt.lrn, "Decision Tree")

# random forest
rf.lrn <- makeLearner("classif.randomForest")
rf.ps <- makeParamSet(
                      makeDiscreteParam("ntree", values = c(250, 500, 1000, 2500, 5000)),
                      makeDiscreteParam("mtry", values = c(2, 3, 4, 5))
)
rf.lrn <- makeTuneWrapper(rf.lrn, rdesc.inner, par.set=rf.ps, control=ctrl)
rf.lrn <- setPredictType(rf.lrn, predict.type="prob")
rf.lrn <- setId(rf.lrn, "Random Forest")

# neural network
nn.lrn <- makeLearner("classif.nnet", MaxNWts=5000, trace=FALSE)
nn.lrn <- makeBaggingWrapper(nn.lrn, bw.iters=bag.n)
nn.ps <- makeParamSet(
                      makeDiscreteParam("size", values = c(2, 3, 5, 7, 10)),
                      makeDiscreteParam("decay", values = c(0.5, 0.25, 0.1, 0))
)
nn.lrn <- makeTuneWrapper(nn.lrn, rdesc.inner, par.set=nn.ps, control=ctrl)
nn.lrn <- setPredictType(nn.lrn, predict.type="prob")
nn.lrn <- setId(nn.lrn, "Neural Network")

# support vector machine
svm.lrn <- makeLearner("classif.svm", kernel="radial")
svm.lrn <- makeBaggingWrapper(svm.lrn, bw.iters=bag.n)
svm.ps <- makeParamSet(
                        makeDiscreteParam("gamma", values = 2^(-2:2)),
                        makeDiscreteParam("cost", values = 2^(-2:2))
)
svm.lrn <- makeTuneWrapper(svm.lrn, rdesc.inner, par.set=svm.ps, control=ctrl)
svm.lrn <- setPredictType(svm.lrn, predict.type="prob")
svm.lrn <- setId(svm.lrn, "Support Vector Machine")

# gradient boosting machine
gbm.lrn <- makeLearner("classif.gbm", distribution="adaboost")
gbm.lrn <- makeBaggingWrapper(gbm.lrn, bw.iters=bag.n)
gbm.ps <- makeParamSet(
                      makeDiscreteParam("n.trees", values = c(250, 500, 1000, 2500, 5000)),
                      makeDiscreteParam("interaction.depth", values = c(2, 3, 4, 5))
)
gbm.lrn <- makeTuneWrapper(gbm.lrn, rdesc.inner, par.set=gbm.ps, control=ctrl)
gbm.lrn <- setPredictType(gbm.lrn, predict.type="prob")
gbm.lrn <- setId(gbm.lrn, "Gradient Boosting Machine")

## inference
dt.mod <- train(dt.lrn, classif.task, subset=train.set)
saveRDS(dt.mod, file.path("../data/dt.mod.rds"))
dt.mod <- dt.mod$learner.model$next.model$learner.model
# plot tree
png(file.path("../data/DecisionTree.png"), height=8*300, width=10*300, res=300)
rpart.plot::prp(dt.mod, type=2, extra=101, varlen=0, box.col=c(adjustcolor("darkviolet", alpha.f=0.4), adjustcolor("forestgreen", alpha.f=0.4))[dt.mod$frame$yval], cex = 1.2)
dev.off()

## benchmark
lrns <- list(rf.lrn, nn.lrn, svm.lrn, gbm.lrn)
bmrk <- benchmark(lrns, subsetTask(classif.task, subset=train.set), rdesc.outer, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1))

# boxplots of mean misclassification error
perf <- getBMRPerformances(bmrk, as.df=TRUE)
perf <- perf[c("learner.id", "mmce")]
png(file.path("../data/BenchmarkMmceBoxplots.png"), height=8*300, width=10*300, res=300)
print(
      ggplot(data=perf, aes(x=learner.id, y=mmce)) +
          geom_boxplot(aes(fill=learner.id)) +
          xlab("") +
          ylab("Misclassification error") +
          scale_fill_brewer(palette="Set1", name="Classifier") +
          theme_bw(base_size=24) +
          theme(axis.text.x=element_blank()) +
          theme(axis.ticks.x=element_blank())
)
dev.off()

# boxplots of other performance measures
perf <- getBMRPerformances(bmrk, as.df=TRUE)
names(perf)[5:ncol(perf)] <- c("Accuracy", "AUC", "Recall/Sensitivity", "Specificity", "Precision", "F1")
perf <- reshape(perf, varying=names(perf)[5:ncol(perf)], v.names="value", timevar="measure", times=names(perf)[5:ncol(perf)], direction="long")
png(file.path("../data/BenchmarkOtherBoxplots.png"), height=10*300, width=13*300, res=300)
print(
      ggplot(data=perf, aes(x=learner.id, y=value)) +
          geom_boxplot(aes(fill=learner.id)) +
          facet_wrap(~ measure, nrow=2) +
          xlab("Measure") +
          ylab("Performance") +
          scale_fill_brewer(palette="Set1", name="Classifier") +
          theme_bw(base_size=24) +
          theme(axis.text.x=element_blank()) +
          theme(axis.ticks.x=element_blank())
)
dev.off()

# ROC curves
roc <- generateThreshVsPerfData(bmrk, measures=list(fpr, tpr))
png(file.path("../data/BenchmarkROC.png"), height=8*300, width=12*300, res=300)
print(
      ggplot(data=roc$data, aes(x=fpr, y=tpr)) +
          geom_path(aes(colour=learner), size=1.5) +
          geom_abline(intercept=0, slope=1, linetype="dashed", colour="darkgrey") +
          xlab("False positive rate") +
          ylab("True positive rate") +
          scale_colour_brewer(palette="Set1", name="Classifier") +
          theme_bw(base_size=24)
)
dev.off()

# PR curves
pr <- generateThreshVsPerfData(bmrk, measures=list(tpr, ppv))
png(file.path("../data/BenchmarkPR.png"), height=8*300, width=12*300, res=300)
print(
      ggplot(data=pr$data, aes(x=tpr, y=ppv)) +
          geom_path(aes(colour=learner), size=1.5) +
          xlab("Recall") +
          ylab("Precision") +
          xlim(0, 1) + 
          ylim(0, 1) +
          scale_colour_brewer(palette="Set1", name="Classifier") +
          theme_bw(base_size=24)
)
dev.off()

## compare classifiers' predictions
rf.res <- bmrk$results$TargetPred$`Random Forest`
saveRDS(rf.res, file.path("../data/rf.res.rds"))
svm.res <- bmrk$results$TargetPred$`Support Vector Machine`
saveRDS(svm.res, file.path("../data/svm.res.rds"))
nn.res <- bmrk$results$TargetPred$`Neural Network`
saveRDS(nn.res, file.path("../data/nn.res.rds"))
gbm.res <- bmrk$results$TargetPred$`Gradient Boosting Machine`
saveRDS(gbm.res, file.path("../data/gbm.res.rds"))
# positive predictions
truth.1 <- subset(rf.res$pred$data, truth == 1, id, drop = TRUE)
rf.1 <- subset(rf.res$pred$data, response == 1, id, drop = TRUE)
svm.1 <- subset(svm.res$pred$data, response == 1, id, drop = TRUE)
nn.1 <- subset(nn.res$pred$data, response == 1, id, drop = TRUE)
gbm.1 <- subset(gbm.res$pred$data, response == 1, id, drop = TRUE)
# negative predictions drop = TRUE)
truth.0 <- subset(rf.res$pred$data, truth == 0, id, drop = TRUE)
rf.0 <- subset(rf.res$pred$data, response == 0, id, drop = TRUE)
svm.0 <- subset(svm.res$pred$data, response == 0, id, drop = TRUE)
nn.0 <- subset(nn.res$pred$data, response == 0, id, drop = TRUE)
gbm.0 <- subset(gbm.res$pred$data, response == 0, id, drop = TRUE)
# Venn diagrams
png(file.path("../data/VennTargets.png"), height=6*300, width=6*300, res=300)
plot(Venn(list("Random Forest" = rf.1, "Neural Network" = nn.1, "Support Vector Machine" = svm.1, "Gradient Boosting Machine" = gbm.1)), doWeights = FALSE)
dev.off()
png(file.path("../data/VennNontargets.png"), height=6*300, width=6*300, res=300)
plot(Venn(list("Random Forest" = rf.0, "Neural Network" = nn.0, "Support Vector Machine" = svm.0, "Gradient Boosting Machine" = gbm.0)), doWeights = FALSE)
dev.off()

## train model
rf.mod <- train(rf.lrn, classif.task, subset=train.set)
saveRDS(rf.mod, file.path("../data/rf.mod.rds"))
svm.mod <- train(svm.lrn, classif.task, subset=train.set)
saveRDS(svm.mod, file.path("../data/svm.mod.rds"))
nn.mod <- train(nn.lrn, classif.task, subset=train.set)
saveRDS(nn.mod, file.path("../data/nn.mod.rds"))
gbm.mod <- train(gbm.lrn, classif.task, subset=train.set)
saveRDS(gbm.mod, file.path("../data/gbm.mod.rds"))

## evaluate performance on test set
rf.test.pred <- predict(rf.mod, task=classif.task, subset=test.set)
saveRDS(rf.test.pred, file.path("../data/rf.test.pred.rds"))
svm.test.pred <- predict(svm.mod, task=classif.task, subset=test.set)
saveRDS(svm.test.pred, file.path("../data/svm.test.pred.rds"))
nn.test.pred <- predict(nn.mod, task=classif.task, subset=test.set)
saveRDS(nn.test.pred, file.path("../data/nn.test.pred.rds"))
gbm.test.pred <- predict(gbm.mod, task=classif.task, subset=test.set)
saveRDS(gbm.test.pred, file.path("../data/gbm.test.pred.rds"))

## export performance measures
# cross-validation
cv.perf <- getBMRAggrPerformances(bmrk, as.df=TRUE)[-1]
write.csv(cv.perf, "../data/crossvalidation.performance.csv", quote = FALSE, row.names=FALSE)
# test
test.perf <- cbind.data.frame(learner.id = c(rf.lrn$id, svm.lrn$id, nn.lrn$id, gbm.lrn$id), t(data.frame(rf = performance(rf.test.pred, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1)), svm = performance(svm.test.pred, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1)), nn = performance(nn.test.pred, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1)), gbm = performance(gbm.test.pred, measures=list(mmce, acc, auc, tpr, tnr, ppv, f1)))))
write.csv(test.perf, "../data/test.performance.csv", quote = FALSE, row.names=FALSE)
