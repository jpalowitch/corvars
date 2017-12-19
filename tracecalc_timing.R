library(rbenchmark)
library(doParallel)
library(doRNG)
library(FNN)
library(lineprof)
source("mvrnormR.R")
source("tracecalcs.R")

m <- 1000
mY <- 1
rho <- 0.9
nsims <- 100
ndata <- 100
Beta <- 0
s2 <- 1

Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)
pvals1 <- pvals2 <- numeric(nsims)

for (i in 1:nsims) {
  
  cat("sim", i, "\n")
  
  # Data generation
  X <- mvrnormR(ndata, rep(0, m), Sig)
  Y <-
    Beta * rowSums(X) + matrix(rnorm(ndata * mY, sd = sqrt(s2)), ncol = mY)
  X <- as.matrix(scale(X))
  Y <- as.matrix(scale(Y))
  
  # General Calcs
  n <- nrow(X)
  
  X2 <- X ^ 2
  X3 <- X ^ 3
  X4 <- X ^ 4
  Y2 <- Y ^ 2
  Y3 <- Y ^ 3
  Y4 <- Y ^ 4
  
  
  tX <- t(X)
  tY <- t(Y)
  tX2 <- t(X2)
  tX3 <- t(X3)
  XXt <- tcrossprod(X)
  XXt2 <- XXt ^ 2
  
  Y4ColSums <- colSums(Y4)
  X4ColSums <- colSums(X4)
  X2RowSums <- rowSums(X2)
  
  allr <- crossprod(Y, X)
  allrSums <- rowSums(allr)
  allr22 <- crossprod(Y2, X2)
  allr31 <- crossprod(Y3, X)
  
  #trs.result <- trace_large_x_indx_timer(1)
  #trs <- trs.result$traces

  a1 <- trs1[2] / trs1[1]
  b1 <- trs1[1] ^ 2 / trs1[2]
  a2 <- trs2[2] / trs2[1]
  b2 <- trs2[1] ^ 2 / trs2[2]
  Tstat <- rowSums(allr ^ 2 / (ndata - 1) ^ 2)
  pvals1[i] <- pchisq(ndata * Tstat / a1, df = b1, lower.tail = FALSE)
  pvals2[i] <- pchisq(ndata * Tstat / a2, df = b2, lower.tail = FALSE)
  
}

library(ggplot2)
p <- qplot(pvals1, pvals2) + geom_abline(slope=1, intercept=0) + 
  xlab("slow version") + ylab("fast version") + 
  ggtitle("trace calc comparison")
ggsave("approx_compare.png", p)

# Save a single benchmark
time.df <- benchmark(trs1 <- trace_large_x_indx(1),
                     trs2 <- trace_large_x_indx_faster(1))
write.table(time.df, quote=FALSE, row.names=FALSE, file="approx_time.txt",
            sep="\t")

