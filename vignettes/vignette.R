## ---- echo =  FALSE--------------------------------------------------------
knitr::opts_chunk$set(eval = TRUE, 
                      echo = TRUE,
                      fig.align = "center",
                      warning = FALSE,
                      message = FALSE)

## ---- message=FALSE, warning=FALSE-----------------------------------------
library(timeOmics)

## --------------------------------------------------------------------------
library(lmms)

data("timeOmics.simdata")
sim.data <- timeOmics.simdata$sim
sim.data[1:15,1:6]

