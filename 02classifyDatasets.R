library(mlr)
set.seed(16)

# data
dataset <- readRDS(file.path("../data/dataset.rds"))

# mlr setup
classif.task <- makeClassifTask(id="TargetPred", data=dataset[-1], target="target")
svm.learner <- makeLearner("classif.svm", predict.type="response")
bagged.svm.learner <- makeBaggingWrapper(svm.learner, bw.iters=10)
resampl.desc <- makeResampleDesc("CV", iters=3)

# model
mod <- resample(bagged.svm.learner, classif.task, resampl.desc)
