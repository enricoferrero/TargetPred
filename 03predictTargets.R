### libraries ###
library(mlr)
library(parallelMap)

### options ###
set.seed(16)
parallelStart("multicore", 16)

### data ###
predicitonset <- readRDS(file.path("../data/predicitionset.rds"))
mod <- readRDS(file.path("../data/mod.rds"))

### predict ###
pred <- predict(mod, newdata=predicitonset)
