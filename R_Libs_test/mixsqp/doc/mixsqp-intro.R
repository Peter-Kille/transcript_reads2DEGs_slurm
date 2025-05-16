## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,results = "hold",comment = "#",
                      fig.align = "center",eval = FALSE)

## ----load-pkgs, eval=TRUE, message=FALSE--------------------------------------
library(mixsqp)

## ----set-seed, eval=TRUE------------------------------------------------------
set.seed(1)

## ----sim-data-small, eval=TRUE------------------------------------------------
L <- simulatemixdata(1000,20)$L
dim(L)

## ----fit-model-mixsqp-small, eval=TRUE----------------------------------------
fit.sqp <- mixsqp(L)

## ----plot-sqp-progress, eval=TRUE, fig.height=5, fig.width=7------------------
numiter <- nrow(fit.sqp$progress)
plot(1:numiter,fit.sqp$progress$objective,type = "b",
     pch = 20,lwd = 2,xlab = "SQP iteration",
     ylab = "objective",xaxp = c(1,numiter,numiter - 1))

## ----session-info, eval=TRUE--------------------------------------------------
sessionInfo()

