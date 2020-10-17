library(lmms)

data(timeOmics.simdata)
data <- timeOmics.simdata$sim
lmms.output <- timeOmics.simdata$lmms.output
time <- timeOmics.simdata$time

tmp_ok <- lmms::lmmSpline(data = data[,c(2,3)], time = time, sampleID = rownames(data), keepModels = TRUE)
tmp_lmms.bad1 <- lmms::lmmSpline(data = data, time = time, sampleID = rownames(data))
tmp_lmms.bad2 <- lmms::lmmSpline(data = data, time = time, sampleID = rownames(data), timePredict = c(1:3), keepModels = TRUE)

setClass("lmmspline",slots= c(predSpline="data.frame", modelsUsed="numeric",models="list",derivative='logical', basis="character", knots="numeric",errorMolecules="character"))

ok <- new("lmmspline", 
          predSpline = tmp_ok@predSpline,
          modelsUsed = tmp_ok@modelsUsed,
          models = tmp_ok@models,
          derivative = tmp_ok@derivative, 
          basis = tmp_ok@basis, 
          knots = tmp_ok@knots, 
          errorMolecules = tmp_ok@errorMolecules)

lmms.bad1 <- new("lmmspline", 
          predSpline = tmp_lmms.bad1@predSpline,
          modelsUsed = tmp_lmms.bad1@modelsUsed,
          models = tmp_lmms.bad1@models,
          derivative = tmp_lmms.bad1@derivative, 
          basis = tmp_lmms.bad1@basis, 
          knots = tmp_lmms.bad1@knots, 
          errorMolecules = tmp_lmms.bad1@errorMolecules)

lmms.bad2 <- new("lmmspline", 
                     predSpline = tmp_lmms.bad2@predSpline,
                     modelsUsed = tmp_lmms.bad2@modelsUsed,
                     models = tmp_lmms.bad2@models,
                     derivative = tmp_lmms.bad2@derivative, 
                     basis = tmp_lmms.bad2@basis, 
                     knots = tmp_lmms.bad2@knots, 
                     errorMolecules = tmp_lmms.bad2@errorMolecules)

lmms.test.data <- list(ok = ok, lmms.bad1 = lmms.bad1, lmms.bad2 = lmms.bad2)
save(lmms.test.data, file = "~/Downloads/lmms.test.data.rda")
