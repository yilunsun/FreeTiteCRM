PI <- c(0.10, 0.20, 0.40, 0.50, 0.60, 0.65)
prior <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70)
target <- 0.2
x0 <- 3
n <- 24
nsim <- 2
restrict <- TRUE
obswin <- 6
tgrp <- obswin
rate <- 4
accrual <- "poisson"
surv <- "uniform"
scheme <- "linear"
count <- TRUE
method <- "bayes"
model <- "empiric"
intcpt <- 3
scale <- sqrt(1.34)
seed <- 1009

source("R/util.R")
source("R/titecrm.R")
source("R/mtite.R")

