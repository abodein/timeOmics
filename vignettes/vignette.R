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

lmms.output <- lmmSpline(data = data.filtered, time = time,
                         sampleID = rownames(data.filtered), deri = FALSE,
                         basis = "p-spline", numCores = 4, timePredict = 1:9,
                         keepModels = TRUE)
modelled.data <- t(slot(lmms.output, 'predSpline'))

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

