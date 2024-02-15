## ---- echo =  FALSE-----------------------------------------------------------
knitr::opts_chunk$set(eval = TRUE, 
                      echo = TRUE,
                      fig.align = "center",
                      warning = FALSE,
                      message = FALSE)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(timeOmics)

## ---- message=F---------------------------------------------------------------
library(tidyverse)

## -----------------------------------------------------------------------------
data("timeOmics.simdata")
sim.data <- timeOmics.simdata$sim

dim(sim.data) 
head(sim.data[,1:6])

## -----------------------------------------------------------------------------
remove.low.cv <- function(X, cutoff = 0.5){
  # var.coef
  cv <- unlist(lapply(as.data.frame(X), 
                      function(x) abs(sd(x)/mean(x))))
  return(X[,cv > cutoff])
}

data.filtered <- remove.low.cv(sim.data, 0.5)

## ---- message=FALSE-----------------------------------------------------------
# numeric vector containing the sample time point information
time <- timeOmics.simdata$time
head(time)

## ----eval=FALSE---------------------------------------------------------------
#  # example of lmms
#  lmms.output <- lmms::lmmSpline(data = data.filtered, time = time,
#                           sampleID = rownames(data.filtered), deri = FALSE,
#                           basis = "p-spline", numCores = 4, timePredict = 1:9,
#                           keepModels = TRUE)
#  modelled.data <- t(slot(lmms.output, 'predSpline'))

## ---- warning=FALSE, message=FALSE--------------------------------------------
lmms.output <- timeOmics.simdata$lmms.output
modelled.data <- timeOmics.simdata$modelled

## -----------------------------------------------------------------------------
# gather data
data.gathered <- modelled.data %>% as.data.frame() %>% 
  rownames_to_column("time") %>%
  mutate(time = as.numeric(time)) %>%
  pivot_longer(names_to="feature", values_to = 'value', -time)

# plot profiles
ggplot(data.gathered, aes(x = time, y = value, color = feature)) + geom_line() +
  theme_bw() + ggtitle("`lmms` profiles") + ylab("Feature expression") +
  xlab("Time")

## -----------------------------------------------------------------------------
filter.res <- lmms.filter.lines(data = data.filtered, 
                                lmms.obj = lmms.output, time = time)
profile.filtered <- filter.res$filtered

## -----------------------------------------------------------------------------
# run pca
pca.res <- pca(X = profile.filtered, ncomp = 5, scale=FALSE, center=FALSE)

# tuning ncomp
pca.ncomp <- getNcomp(pca.res, max.ncomp = 5, X = profile.filtered, 
                      scale = FALSE, center=FALSE)

pca.ncomp$choice.ncomp
plot(pca.ncomp)

## -----------------------------------------------------------------------------
# final model
pca.res <- pca(X = profile.filtered, ncomp = 2, scale = FALSE, center=FALSE)

## -----------------------------------------------------------------------------
# extract cluster
pca.cluster <- getCluster(pca.res)
head(pca.cluster)

## -----------------------------------------------------------------------------
plotIndiv(pca.res)

## -----------------------------------------------------------------------------
plotVar(pca.res)

## -----------------------------------------------------------------------------
plotLoadings(pca.res)

## -----------------------------------------------------------------------------
plotLong(pca.res, scale = FALSE, center = FALSE, 
         title = "PCA longitudinal clustering")

## -----------------------------------------------------------------------------
tune.spca.res <- tuneCluster.spca(X = profile.filtered, ncomp = 2, 
                                  test.keepX = c(2:10))
# selected features in each component
tune.spca.res$choice.keepX
plot(tune.spca.res)

## -----------------------------------------------------------------------------
# final model
spca.res <- spca(X = profile.filtered, ncomp = 2, 
                 keepX = tune.spca.res$choice.keepX, scale = FALSE)
plotLong(spca.res)

## -----------------------------------------------------------------------------
X <- profile.filtered
Y <- timeOmics.simdata$Y

pls.res <- pls(X,Y, ncomp = 5, scale = FALSE)
pls.ncomp <- getNcomp(pls.res, max.ncomp = 5, X=X, Y=Y, scale = FALSE)
pls.ncomp$choice.ncomp
plot(pls.ncomp)

## -----------------------------------------------------------------------------
# final model
pls.res <- pls(X,Y, ncomp = 2, scale = FALSE)

# info cluster
head(getCluster(pls.res))
# plot clusters
plotLong(pls.res, title = "PLS longitudinal clustering", legend = TRUE)

## -----------------------------------------------------------------------------
tune.spls <- tuneCluster.spls(X, Y, ncomp = 2, test.keepX = c(4:10), test.keepY <- c(2,4,6))

# selected features in each component on block X
tune.spls$choice.keepX
# selected features in each component on block Y
tune.spls$choice.keepY

# final model
spls.res <- spls(X,Y, ncomp = 2, scale = FALSE, 
                 keepX = tune.spls$choice.keepX, keepY = tune.spls$choice.keepY)

# spls cluster
spls.cluster <- getCluster(spls.res)

# longitudinal cluster plot
plotLong(spls.res, title = "sPLS clustering")

## -----------------------------------------------------------------------------
X <- list("X" = profile.filtered, "Z" = timeOmics.simdata$Z)
Y <- as.matrix(timeOmics.simdata$Y)

block.pls.res <- block.pls(X=X, Y=Y, ncomp = 5, 
                           scale = FALSE, mode = "canonical")
block.ncomp <- getNcomp(block.pls.res,X=X, Y=Y, 
                        scale = FALSE, mode = "canonical")
block.ncomp$choice.ncomp
plot(block.ncomp)

## -----------------------------------------------------------------------------
# final model
block.pls.res <- block.pls(X=X, Y=Y, ncomp = 1, scale = FALSE, mode = "canonical")
# block.pls cluster
block.pls.cluster <- getCluster(block.pls.res)

# longitudinal cluster plot
plotLong(block.pls.res)

## -----------------------------------------------------------------------------
test.list.keepX <- list("X" = 4:10, "Z" = c(2,4,6,8))
test.keepY <- c(2,4,6)

tune.block.res <- tuneCluster.block.spls(X= X, Y= Y, 
                                         test.list.keepX=test.list.keepX, 
                                         test.keepY= test.keepY, 
                                         scale=FALSE, 
                                         mode = "canonical", ncomp = 1)
# ncomp = 1 given by the getNcomp() function

# selected features in each component on block X
tune.block.res$choice.keepX
# selected features in each component on block Y
tune.block.res$choice.keepY

# final model
block.pls.res <- block.spls(X=X, Y=Y, 
                            ncomp = 1, 
                            scale = FALSE, 
                            mode = "canonical", 
                            keepX = tune.block.res$choice.keepX, 
                            keepY = tune.block.res$choice.keepY)

head(getCluster(block.pls.res))
plotLong(block.pls.res)

