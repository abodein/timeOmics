context('lmms.filter.lines')

 # data and lmms output
data(timeOmics.simdata)
data <- timeOmics.simdata$sim
lmms.output <- timeOmics.simdata$lmms.output
time <- timeOmics.simdata$time

ok <- lmms::lmmSpline(data = data[,c(2,3)], time = time, sampleID = rownames(data), keepModels = TRUE)

lmms.bad1 <- lmms::lmmSpline(data = data, time = time, sampleID = rownames(data))
lmms.bad2 <- lmms::lmmSpline(data = data, time = time, sampleID = rownames(data), timePredict = c(1:3), keepModels = TRUE)

test_that("lmms.filter.lines failed on invalid input - data",{
    expect_error(lmms.filter.lines(data = ""), "X must be a numeric matrix/data.frame")
    expect_error(lmms.filter.lines(data = 1), "X must be a numeric matrix/data.frame")
    expect_error(lmms.filter.lines(data = NULL), "X must be a numeric matrix/data.frame")
    expect_error(lmms.filter.lines(data = NA), "X must be a numeric matrix/data.frame")
})

test_that("lmms.filter.lines failed on invalid input - lmms.obj",{
    # invalide type
    expect_error(lmms.filter.lines(data = data, lmms.obj = "" ), "'lmms.obj' should be a 'lmms' object.")
    expect_error(lmms.filter.lines(data = data, lmms.obj = c()), "'lmms.obj' should be a 'lmms' object.")
    # no models
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.bad1), "No models found in 'lmms.obj', please use 'keepModels = TRUE' in 'lmmSpline()'", fixed=TRUE)
    # no correct time
    expect_error(lmms.filter.lines(data = data, time = time, lmms.obj = lmms.bad2), "wrong time between 'lmms.obj', and 'time'", fixed = TRUE)
})

test_that("lmms.filter.lines failed on invalid input - time",{
    # invalide type
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = ""),  "'time' should be a numeric vector with the same length as 'nrow(data)'", fixed = TRUE)
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = list()),  "'time' should be a numeric vector with the same length as 'nrow(data)'", fixed = TRUE)
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = NA),  "'time' should be a numeric vector with the same length as 'nrow(data)'", fixed = TRUE)
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = NULL),  "'time' should be a numeric vector with the same length as 'nrow(data)'", fixed = TRUE)
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = data.frame()),  "'time' should be a numeric vector with the same length as 'nrow(data)'", fixed = TRUE)
    
    # wrong length
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = c(1,2,3)),  "'time' should be a numeric vector with the same length as 'nrow(data)'", fixed = TRUE)
    
    # not identical
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = sample(x=9:20,size = 45,replace = TRUE)),  "wrong time between 'lmms.obj', and 'time'", fixed = TRUE)
})

test_that("lmms.filter.lines failed on invalid input - homoskedasticity.cutoff",{
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity.cutoff = ""), "'homoskedasticity.cutoff' should be a numeric between 0 and 1")
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity.cutoff = NULL), "'homoskedasticity.cutoff' should be a numeric between 0 and 1")
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity.cutoff = NA), "'homoskedasticity.cutoff' should be a numeric between 0 and 1")
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity.cutoff = list()), "'homoskedasticity.cutoff' should be a numeric between 0 and 1")
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity.cutoff = matrix()), "'homoskedasticity.cutoff' should be a numeric between 0 and 1")
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity.cutoff = -1), "'homoskedasticity.cutoff' should be a numeric between 0 and 1")
    expect_error(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity.cutoff = 1.4), "'homoskedasticity.cutoff' should be a numeric between 0 and 1")
})

test_that("lmms.filter.lines works",{
    expect_is(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time), "list")
    expect_is(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity = ""), "list")
    expect_is(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity = list()), "list")
    expect_is(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, MSE.filter = list()), "list")
    expect_is(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, MSE.filter = ""), "list")
    expect_is(lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time, homoskedasticity.cutoff = 0.01), "list")
    expect_is(lmms.filter.lines(data = data[,c(2,3)], lmms.obj = ok, time = time), "list")
})
    