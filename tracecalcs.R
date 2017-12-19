library(emulator)

trace_indx <- function(i) {
  if(ndata > 1.2*m){
    trace_uni_fast(i)
  } else {
    trace_large_x_indx(i)
  }
}

trace_uni <- function (i) {

  # General calcs
  xyCors <- as.vector(allr[i, ]) / (n - 1)
  r <- allrSums[i]
  rX2 <- crossprod(tX2, xyCors)

  # star 1
  S1 <- crossprod(as.vector(Y[ , i]) * X)

  # star 2
  S2 <- Y4ColSums[i] * tcrossprod(xyCors)

  # star 3
  S3 <- tcrossprod(xyCors, xyCors * colSums(X2 * as.vector(Y2[ , i])))

  # star 4
  S4 <- tcrossprod(xyCors * t(X2))

  # dagger 1
  D1 <- tcrossprod(xyCors, colSums(X * as.vector(Y3[ , i])))

  # dagger 2
  D2 <- t(crossprod(as.vector(Y[ , i]) * X, X2)) * xyCors

  bigmat <- (S1 + (S2 + S3 + t(S3) + S4) / 4 - (D1 + t(D1) + D2 + t(D2)) / 2) / (n - 1)

  trSig <- sum(diag(bigmat))
  trSig2 <- sum(diag(crossprod(bigmat)))

  return(c(trSig, trSig2))

}

trace_uni_fast <- function(i){
  
  r <- allr[i, ] / (n - 1)
  r2 <- allr22[i, ] / (n - 1)
  r3 <- allr31[i, ] / (n - 1)
  
  yX <- Y[ , i] * X
  tX2r <- tX2 * r
  # star 1
  S1 <- crossprod(yX)
  
  # star 2
  S2 <- Y4ColSums[i] * tcrossprod(r)
  
  # star 3
  S3 <- tcrossprod(r, r * r2) * (n-1)
  
  # star 4
  S4 <- tcrossprod(tX2r)
  
  # dagger 1
  D1 <- tcrossprod(r, r3) * (n-1)
  
  # dagger 2
  D2 <- tX2r %*% yX
  
  bigmat <- (S1 + (S2 + S3 + t(S3) + S4) / 4 - (D1 + t(D1) + D2 + t(D2)) / 2) / (n - 1)
  
  trSig <- sum(diag(bigmat))
  trSig2 <- sum(bigmat^2)
  
  return(c(trSig, trSig2))
}

trace_uni_mlist <- function (Yvec, Xmat) {

  Xmat <- as.matrix(scale(Xmat))
  Yvec <- as.vector(scale(Yvec))

  n <- nrow(Xmat)

  if (n != length(Yvec))
    stop('X and Y variables must be of same dimension\n')

  # General calcs
  xyCors <- as.vector(cor(Yvec, Xmat))
  y4 <- sum(Yvec^4)
  xRowSum <- rowSums(Xmat)
  xRowSum2 <- tcrossprod(xyCors, Xmat^2)

  # star 1
  S1 <- crossprod(Yvec * Xmat)

  # star 2
  S2 <- y4 * tcrossprod(xyCors)

  # star 3
  S3 <- tcrossprod(xyCors, xyCors * colSums(Xmat^2 * Yvec^2))

  # star 4
  S4 <- tcrossprod(xyCors * t(Xmat^2))

  # dagger 1
  D1 <- tcrossprod(xyCors, colSums(Xmat * Yvec^3))

  # dagger 2
  D2 <- t(crossprod(Yvec * Xmat, Xmat^2)) * xyCors

  # traces
  bigmat <- (S1 + (S2 + S3 + t(S3) + S4) / 4 - (D1 + t(D1) + D2 + t(D2)) / 2) / (n - 1)
  trSig <- sum(diag(bigmat))
  trSig2 <- sum(diag(crossprod(bigmat)))

  return(list(A = S1, B1 = S2, B2 = S3, B3 = t(S3), B4 = S4,
              C1 = t(D1), C2 = t(D2), C3 = D1, C4 = D2,
              tr1 = trSig, tr2 = trSig2))

}

# trace_uni_fast <- function () {
# 
#   # General calcs
#   xyCors <- crossprod(X, Y) / (n - 1)
#   xyCors2 <- xyCors^2
#   rX <- crossprod(tX, xyCors)
#   r2X <- crossprod(tX, xyCors2)
#   r2X2 <- crossprod(tX2, xyCors2)
#   rX3 <- crossprod(tX3, xyCors)
# 
# 
#   #AA
#   AA <- crossprod(Y2, crossprod(XXt2, Y2))
# 
# 
# 
#   #AB1
#   AB1 <- Y4ColSums * crossprod(Y2, rX^2) / 4
# 
# 
# 
#   #AB2
#   AB2array <- do.call(cbind, lapply(1:n, function (i) tX[ , i] * tX2))
#   AB2mat <- matrix(crossprod(AB2array, xyCors), n, byrow = TRUE)
#   AB2list <- lapply(1:ncol(AB2mat), function (i) matrix(AB2mat[ , i], nrow = n, byrow = TRUE))
#   AB2 <- crossprod(Y2, crossprod(AB2mat, Y2 * rX)) / 4
# 
#   #AB3
#   AB3 <- AB2
# 
# 
# 
#   #AB4
#   AB4 <- crossprod(rep(1, n), crossprod(AB2mat^2, Y2)) / 4
# 
# 
# 
#   #AC1
#   AC1vec1 <- crossprod(Y2 * rX, X)
#   AC1vec2 <- crossprod(X, Y3)
#   AC1 <- -AC1vec1 %*% AC1vec2 / 2
# 
#   #AC2
#   AC2 <- AC1
# 
# 
# 
#   #AC3
#   AC3 <- -crossprod(Y, crossprod(AB2mat * XXt, Y2)) / 2
# 
# 
# 
#   #AC4
#   AC4 <- AC3
# 
# 
# 
#   #B1^2
#   B1B1 <- Y4ColSums^2 * sum(xyCors2)^2 / 16
# 
# 
# 
#   #B1B2
#   B1B2 <- Y4ColSums * crossprod(Y2, r2X2) * sum(xyCors2) / 16
# 
# 
# 
#   #B1B3
#   B1B3 <- B1B2
# 
# 
# 
#   #B1B4
#   B1B4 <- Y4ColSums * crossprod(r2X2) / 16
# 
# 
# 
#   #B1C1
#   B1C1 <- -Y4ColSums * sum(xyCors2) * crossprod(Y3, rX) / 8
# 
# 
# 
#   #B1C2
#   B1C2 <- B1C1
# 
# 
# 
#   #B1C3
#   B1C3 <- -Y4ColSums * crossprod(Y, rX * r2X2) / 8
# 
# 
# 
#   #B1C4
#   B1C4 <- B1C3
# 
# 
# 
#   #B2^2
#   B2B2array <- do.call(cbind, lapply(1:n, function (i) tX2[ , i] * tX2))
#   B2B2mat <- matrix(crossprod(B2B2array, xyCors2), n, byrow = TRUE)
#   B2B2 <- sum(xyCors2) * crossprod(Y2, crossprod(B2B2mat, Y2)) / 16
# 
#   #B2B3
#   B2B3 <- crossprod(Y2, r2X2)^2 / 16
# 
# 
# 
#   #B2B4
#   B2B4 <- crossprod(Y2, crossprod(B2B2mat, r2X2)) / 16
# 
# 
# 
#   #B2C1
#   B2C1 <- -crossprod(Y3, rX) * crossprod(Y2, r2X2) / 8
# 
# 
#   #B2C2
#   B2C2mat <- matrix(crossprod(AB2array, xyCors), n, byrow = TRUE)
#   B2C2 <- -sum(xyCors2) * crossprod(Y2, crossprod(B2C2mat, Y3)) / 8
# 
# 
#   #B2C3
#   B2C3 <- -crossprod(Y2, crossprod(B2B2mat, Y * rX)) / 8
# 
# 
# 
#   #B2C4
#   B2C4 <- -crossprod(Y2, crossprod(AB2mat, Y * r2X2)) / 8
# 
# 
# 
#   #B3^2
#   B3B3 <- B2B2
# 
# 
# 
#   #B3B4
#   B3B4 <- B2B4
# 
# 
# 
#   #B3C1
#   B3C1 <- B2C2
# 
# 
# 
#   #B3C2
#   B3C2 <- B2C1
# 
# 
# 
#   #B3C3
#   B3C3 <- B2C4
# 
# 
# 
#   #B3C4
#   B3C4 <- B2C3
# 
# 
# 
#   #B4^2
#   B4B4 <- sum(B2B2mat^2) / 16
# 
# 
# 
#   #B4C1
#   B4C1 <- -crossprod(r2X2, crossprod(B2C2mat, Y3)) / 8
# 
# 
# 
#   #B4C2
#   B4C2 <- B4C1
# 
# 
# 
#   #B4C3
#   B4C3 <- -sum(crossprod(Y, B2B2mat * B2C2mat)) / 8
# 
# 
# 
#   #B4C4
#   B4C4 <- B4C3
# 
# 
# 
#   #C1^2
#   C1C1 <- sum(xyCors2) * crossprod(crossprod(X, Y3)) / 4
# 
# 
# 
#   #C1C2
#   C1C2 <- crossprod(Y3, rX)^2 / 4
# 
# 
# 
#   #C1C3
#   C1C3 <- crossprod(crossprod(X, Y3), crossprod(X, Y * r2X2)) / 4
# 
# 
# 
#   #C1C4
#   C1C4 <- crossprod(Y3, crossprod(t(AB2mat), Y * rX)) / 4
# 
# 
# 
#   #C2^2
#   C2C2 <- C1C1
# 
# 
# 
#   #C2C3
#   C2C3 <- C1C4
# 
# 
# 
#   #C2C4
#   C2C4 <- C1C3
# 
# 
# 
#   #C3^2
#   C3C3 <- crossprod(Y, crossprod(B2B2mat * XXt, Y)) / 4
# 
# 
# 
#   #C3C4
#   C3C4 <- crossprod(Y, crossprod(AB2mat * t(AB2mat), Y)) / 4
# 
# 
# 
#   #C4^2
#   C4C4 <- C3C3
# 
#   # Collect results
#   resvec <- c(AA,  AB1,  AB2,  AB3,  AB4,  AC1,  AC2,  AC3,  AC4,
#                   B1B1, B1B2, B1B3, B1B4, B1C1, B1C2, B1C3, B1C4,
#                         B2B2, B2B3, B2B4, B2C1, B2C2, B2C3, B2C4,
#                               B3B3, B3B4, B3C1, B3C2, B3C3, B3C4,
#                                     B4B4, B4C1, B4C2, B4C3, B4C4,
#                                           C1C1, C1C2, C1C3, C1C4,
#                                                 C2C2, C2C3, C2C4,
#                                                       C3C3, C3C4,
#                                                             C4C4)
# 
#   resmat <- matrix(0, 9, 9)
#   resmat[row(resmat) >= col(resmat)] <- resvec
#   trprod <- sum(diag(resmat)) + 2 * sum(resmat[row(resmat) > col(resmat)])
# 
#   ### Now calculating trace of the straight-up varmat
#   star1 <- crossprod(Y2, X2RowSums)
#   star2 <- Y4ColSums * sum(xyCors2)
#   star3 <- crossprod(Y2, r2X2)
#   star4 <- crossprod(xyCors2, X4ColSums)
#   dagger1 <- crossprod(Y3, rX)
#   dagger2 <- crossprod(Y, rX3)
#   trmat <- star1 + (star2 + 2 * star3 + star4) / 4 - dagger1 - dagger2
# 
#   return(c(trmat / (n - 1), trprod / (n - 1)^2))
# 
# }
# 
# trace_fast <- function () {
# 
#   # General calcs
#   xyCors <- as.vector(crossprod(Y, X)) / (n - 1)
#   xyCors2 <- xyCors^2
#   rX <- crossprod(tX, xyCors)
#   r2X <- crossprod(tX, xyCors2)
#   r2X2 <- crossprod(tX2, xyCors2)
#   rX3 <- crossprod(tX3, xyCors)
# 
# 
# 
#   #AA
#   AA <- crossprod(Y2, crossprod(XXt2, Y2))
# 
# 
# 
#   #AB1
#   AB1 <- Y4ColSums * crossprod(Y2, rX^2) / 4
# 
# 
# 
#   #AB2
#   AB2array <- do.call(cbind, lapply(1:n, function (i) tX[ , i] * tX2))
#   AB2mat <- matrix(crossprod(AB2array, xyCors), n, byrow = TRUE)
#   AB2 <- crossprod(Y2, crossprod(AB2mat, Y2 * rX)) / 4
# 
#   #AB3
#   AB3 <- AB2
# 
# 
# 
#   #AB4
#   AB4 <- sum(Y2 * AB2mat^2) / 4
# 
# 
# 
#   #AC1
#   AC1vec1 <- crossprod(Y2 * rX, X)
#   AC1vec2 <- crossprod(X, Y3)
#   AC1 <- -AC1vec1 %*% AC1vec2 / 2
# 
#   #AC2
#   AC2 <- AC1
# 
# 
# 
#   #AC3
#   AC3 <- -crossprod(Y, crossprod(AB2mat * XXt, Y2)) / 2
# 
# 
#   #AC4
#   AC4 <- AC3
# 
# 
# 
#   #B1^2
#   B1B1 <- Y4ColSums^2 * sum(xyCors2)^2 / 16
# 
# 
# 
#   #B1B2
#   B1B2 <- Y4ColSums * crossprod(Y2, r2X2) * sum(xyCors2) / 16
# 
# 
# 
#   #B1B3
#   B1B3 <- B1B2
# 
# 
# 
#   #B1B4
#   B1B4 <- Y4ColSums * crossprod(r2X2) / 16
# 
# 
# 
#   #B1C1
#   B1C1 <- -Y4ColSums * sum(xyCors2) * crossprod(Y3, rX) / 8
# 
# 
# 
#   #B1C2
#   B1C2 <- B1C1
# 
# 
# 
#   #B1C3
#   B1C3 <- -Y4ColSums * crossprod(Y, rX * r2X2) / 8
# 
# 
# 
#   #B1C4
#   B1C4 <- B1C3
# 
# 
# 
#   #B2^2
#   B2B2array <- do.call(cbind, lapply(1:n, function (i) tX2[ , i] * tX2))
#   B2B2mat <- matrix(crossprod(B2B2array, xyCors2), n, byrow = TRUE)
#   B2B2 <- sum(xyCors2) * crossprod(Y2, crossprod(B2B2mat, Y2)) / 16
# 
#   #B2B3
#   B2B3 <- crossprod(Y2, r2X2)^2 / 16
# 
# 
# 
#   #B2B4
#   B2B4 <- crossprod(Y2, crossprod(B2B2mat, r2X2)) / 16
# 
# 
# 
#   #B2C1
#   B2C1 <- -crossprod(Y3, rX) * crossprod(Y2, r2X2) / 8
# 
# 
#   #B2C2
#   B2C2mat <- matrix(crossprod(AB2array, xyCors), n, byrow = TRUE)
#   B2C2 <- -sum(xyCors2) * crossprod(Y2, crossprod(B2C2mat, Y3)) / 8
# 
# 
#   #B2C3
#   B2C3 <- -crossprod(Y2, crossprod(B2B2mat, Y * rX)) / 8
# 
# 
# 
#   #B2C4
#   B2C4 <- -crossprod(Y2, crossprod(AB2mat, Y * r2X2)) / 8
# 
# 
# 
#   #B3^2
#   B3B3 <- B2B2
# 
# 
# 
#   #B3B4
#   B3B4 <- B2B4
# 
# 
# 
#   #B3C1
#   B3C1 <- B2C2
# 
# 
# 
#   #B3C2
#   B3C2 <- B2C1
# 
# 
# 
#   #B3C3
#   B3C3 <- B2C4
# 
# 
# 
#   #B3C4
#   B3C4 <- B2C3
# 
# 
# 
#   #B4^2
#   B4B4 <- sum(B2B2mat^2) / 16
# 
# 
# 
#   #B4C1
#   B4C1 <- -crossprod(r2X2, crossprod(B2C2mat, Y3)) / 8
# 
# 
# 
#   #B4C2
#   B4C2 <- B4C1
# 
# 
# 
#   #B4C3
#   B4C3 <- -sum(crossprod(Y, B2B2mat * B2C2mat)) / 8
# 
# 
# 
#   #B4C4
#   B4C4 <- B4C3
# 
# 
# 
#   #C1^2
#   C1C1 <- sum(xyCors2) * crossprod(crossprod(X, Y3)) / 4
# 
# 
# 
#   #C1C2
#   C1C2 <- crossprod(Y3, rX)^2 / 4
# 
# 
# 
#   #C1C3
#   C1C3 <- crossprod(crossprod(X, Y3), crossprod(X, Y * r2X2)) / 4
# 
# 
# 
#   #C1C4
#   C1C4 <- crossprod(Y3, crossprod(t(AB2mat), Y * rX)) / 4
# 
# 
# 
#   #C2^2
#   C2C2 <- C1C1
# 
# 
# 
#   #C2C3
#   C2C3 <- C1C4
# 
# 
# 
#   #C2C4
#   C2C4 <- C1C3
# 
# 
# 
#   #C3^2
#   C3C3 <- crossprod(Y, crossprod(B2B2mat * XXt, Y)) / 4
# 
# 
# 
#   #C3C4
#   C3C4 <- crossprod(Y, crossprod(AB2mat * t(AB2mat), Y)) / 4
# 
# 
# 
#   #C4^2
#   C4C4 <- C3C3
# 
#   # Collect results
#   resvec <- c(AA,  AB1,  AB2,  AB3,  AB4,  AC1,  AC2,  AC3,  AC4,
#                   B1B1, B1B2, B1B3, B1B4, B1C1, B1C2, B1C3, B1C4,
#                         B2B2, B2B3, B2B4, B2C1, B2C2, B2C3, B2C4,
#                               B3B3, B3B4, B3C1, B3C2, B3C3, B3C4,
#                                     B4B4, B4C1, B4C2, B4C3, B4C4,
#                                           C1C1, C1C2, C1C3, C1C4,
#                                                 C2C2, C2C3, C2C4,
#                                                       C3C3, C3C4,
#                                                             C4C4)
# 
#   resmat <- matrix(0, 9, 9)
#   resmat[row(resmat) >= col(resmat)] <- resvec
#   trprod <- sum(diag(resmat)) + 2 * sum(resmat[row(resmat) > col(resmat)])
# 
#   ### Now calculating trace of the straight-up varmat
#   star1 <- crossprod(Y2, X2RowSums)
#   star2 <- Y4ColSums * sum(xyCors2)
#   star3 <- crossprod(Y2, r2X2)
#   star4 <- crossprod(xyCors2, X4ColSums)
#   dagger1 <- crossprod(Y3, rX)
#   dagger2 <- crossprod(Y, rX3)
#   trmat <- star1 + (star2 + 2 * star3 + star4) / 4 - dagger1 - dagger2
# 
#   return(c(trmat / (n - 1), trprod / (n - 1)^2))
# 
# }

trace_large_x <- function(
                       y, #m x 1 vector
                       X, #m x dx matrix

                       #Optional arguments that can be precomputed:
                       m = nrow(X),
                       XXt = tcrossprod(X),
                       XXt2 = XXt^2,
                       X2 = X^2,
                       X3 = X^3,
                       X4ColSums=colSums(X^4),
                       X2RowSums=rowSums(X2),
                       tX = t(X)) {

  e <- matrix(1, nrow=m, ncol=1)
  y2 <- y^2 # m x 1
  y3 <- y^3

  M <- m-1
  r <- crossprod(y, X) / M # 1 x dx
  r2 <- crossprod(y2, X2) / M # 1 x dx
  r3 <- crossprod(y3,X) / M # 1 x dx
  rX <- tcrossprod(X, r) # m x 1
  r3X <- tcrossprod(X, r3) # m x 1
  rX3 <- tcrossprod(X3, r) # m x 1


  rSq <- r^2
  rSqX2 <- tcrossprod(X2,rSq) # m x 1
  r2rSqX2 <- tcrossprod(X2, r2 * rSq) # m x 1
  r3rX2 <- tcrossprod(X2, r3*r) # m x 1

  tXr <- tX * as.vector(r) # dx x m

  y4Sum <- sum(y^4) # scalar
  rSqSum <- sum(rSq) #scalar

  Z1 <- X2 %*% tXr # m x m
  Z2 <- X2 %*% tXr^2 # m x m

  #All of these are scalars
  #TODO: Normalize them by (m-1)^2 after they are working correctly.

  AA <- quad.form(XXt2, y2)

  AB1 <- y4Sum * crossprod(y2, rX^2)

  AB2 <- quad.3form(XXt, y, y2 * tcrossprod(X, r2*r))

  AB3 <- AB2

  AB4 <- quad.3form(Z1^2, e, y2)

  AC1 <- crossprod(rX * r3X, y2) * M

  AC2 <- quad.3form(Z1 * XXt, y, y2)

  AC3 <- AC1

  AC4 <- AC2

  B1B1 <- (y4Sum*rSqSum)^2

  B1B2 <- rSqSum * y4Sum * tcrossprod(r2,rSq) * M

  B1B3 <- B1B2

  B1B4 <- y4Sum * sum(rSqX2^2)

  B1C1 <- y4Sum * rSqSum * tcrossprod(r,r3) * M

  B1C2 <- y4Sum * crossprod(y, rX * rSqX2)

  B1C3 <- B1C1

  B1C4 <- B1C2

  B2B2 <- tcrossprod(rSq, r2^2) * rSqSum * M^2

  B2B3 <- tcrossprod(rSq, r2)^2 * M^2

  B2B4 <- sum(rSqX2 * r2rSqX2) * M

  B2C1 <- tcrossprod(rSq, r2) * tcrossprod(r,r3) * M^2

  B2C2 <- crossprod(y, r2rSqX2 * rX) * M

  B2C3 <- sum(r * r2 * r3) * rSqSum * M^2

  B2C4 <- crossprod(y, tcrossprod(X, r * r2) * rSqX2) * M

  B3B3 <- B2B2

  B3B4 <- B2B4

  B3C1 <- B2C3

  B3C2 <- B2C4

  B3C3 <- B2C1

  B3C4 <- B2C2



  B4B4 <- sum(Z2^2)

  B4C1 <- sum(rSqX2 * r3rX2) * M

  B4C2 <- quad.3form(Z2 * Z1, e, y)

  B4C3 <- B4C1

  B4C4 <- B4C2

  C1C1 <- rSqSum * sum(r3^2) * M^2

  C1C2 <- crossprod(y, rSqX2 * r3X) * M

  C1C3 <- sum(r * r3)^2 * M^2

  C1C4 <- crossprod(y, rX * r3rX2) * M

  C2C2 <- quad.form(Z2 * XXt, y)

  C2C3 <- C1C4

  C2C4 <- quad.form(Z1*t(Z1), y)

  C3C3 <- C1C1

  C3C4 <- C1C2

  C4C4 <- C2C2

  resvec <- c(AA, c(AB1, AB2,  AB3,  AB4)/4, -1/2*c( AC1,  AC2,  AC3,  AC4),
               c(B1B1, B1B2, B1B3, B1B4)/16, -1/8*c(B1C1, B1C2, B1C3, B1C4),
                     c(B2B2, B2B3, B2B4)/16, -1/8*c(B2C1, B2C2, B2C3, B2C4),
                           c(B3B3, B3B4)/16, -1/8*c(B3C1, B3C2, B3C3, B3C4),
                                    B4B4/16, -1/8*c(B4C1, B4C2, B4C3, B4C4),
                                          c(C1C1, C1C2, C1C3, C1C4,
                                                C2C2, C2C3, C2C4,
                                                      C3C3, C3C4,
                                                            C4C4)/4)
  resmat <- matrix(0, 9, 9)
  resmat[row(resmat) >= col(resmat)] <- resvec
  trprod <- sum(diag(resmat)) + 2 * sum(resmat[row(resmat) > col(resmat)])

  #Trace of straight up varmat
  star1 <- crossprod(y2, X2RowSums) #O(m)
  star2 <- y4Sum * rSqSum
  star3 <- crossprod(y2, rSqX2) #O(m)
  star4 <- tcrossprod(rSq, X4ColSums) #O(n)
  dagger1 <- crossprod(y3, rX) #O(m)
  dagger2 <- crossprod(y, rX3) #O(m)
  trmat <- star1 + (star2 + 2 * star3 + star4) / 4 - dagger1 - dagger2

  return(c(trmat/M, trprod/M^2))
}

trace_kosher <- function(y, X){
  n <- ncol(X)
  y2 <- y^2
  y3 <- y^3
  M <- nrow(X) - 1

  u4 <- sum(y^4)/M

  s <- function(j,k){
    xj <- X[,j]
    xk <- X[,k]

    xj2 <- xj^2
    xk2 <- xk^2

    u2jk <- sum(y2*xk*xj)/M
    uj <- sum(y*xj)/M
    uk <- sum(y*xk)/M
    u2j2 <- sum(y2*xj2)/M
    u2k2 <- sum(y2*xk2)/M
    j2k2 <- sum(xj2*xk2)/M
    u3k <- sum(y3*xk)/M
    uj2k <- sum(y*xj2*xk)/M
    u3j <- sum(y3*xj)/M
    uk2j <- sum(y*xk2*xj)/M

    u2jk  + uj*uk*(u4 + u2j2 + u2k2 + j2k2)/4 - uj*(u3k + uj2k)/2 - uk*(u3j + uk2j)/2
  }

  S <- outer(seq(n), seq(n), Vectorize(s))
  c(sum(diag(S)), sum(S^2))
}


trace_large_x_indx <- function (i) {
  m <- ndata
  y <- Y[ , i]
  e <- matrix(1, nrow=m, ncol=1)
  y2 <- Y2[ , i] # m x 1  ## GET THIS FROM Y2
  y3 <- Y3[ , i] ## GET THIS FROM Y3

  M <- m - 1
  r <- allr[i, , drop = FALSE] / M # 1 x dx ## GET THIS FROM allr
  r2 <- allr22[i, , drop = FALSE] / M # 1 x dx
  r3 <- allr31[i, , drop = FALSE] / M # 1 x dx
  rX <- tcrossprod(X, r) # m x 1
  r3X <- tcrossprod(X, r3) # m x 1
  rX3 <- tcrossprod(X3, r) # m x 1


  rSq <- r^2
  rSqX2 <- tcrossprod(X2,rSq) # m x 1
  r2rSqX2 <- tcrossprod(X2, r2 * rSq) # m x 1
  r3rX2 <- tcrossprod(X2, r3*r) # m x 1

  tXr <- tX * as.vector(r) # dx x m

  y4Sum <- Y4ColSums[i] # scalar ## GET THIS FROM Y4ColSums
  rSqSum <- sum(rSq) #scalar

  Z1 <- X2 %*% tXr # m x m
  Z2 <- X2 %*% tXr^2 # m x m

  #All of these are scalars
  #TODO: Normalize them by (m-1)^2 after they are working correctly.

  AA <- quad.form(XXt2, y2)

  AB1 <- y4Sum * crossprod(y2, rX^2)

  AB2 <- quad.3form(XXt, y, y2 * tcrossprod(X, r2*r))

  AB3 <- AB2

  AB4 <- quad.3form(Z1^2, e, y2)

  AC1 <- crossprod(rX * r3X, y2) * M

  AC2 <- quad.3form(Z1 * XXt, y, y2)

  AC3 <- AC1

  AC4 <- AC2

  B1B1 <- (y4Sum*rSqSum)^2

  B1B2 <- rSqSum * y4Sum * tcrossprod(r2,rSq) * M

  B1B3 <- B1B2

  B1B4 <- y4Sum * sum(rSqX2^2)

  B1C1 <- y4Sum * rSqSum * tcrossprod(r,r3) * M

  B1C2 <- y4Sum * crossprod(y, rX * rSqX2)

  B1C3 <- B1C1

  B1C4 <- B1C2

  B2B2 <- tcrossprod(rSq, r2^2) * rSqSum * M^2

  B2B3 <- tcrossprod(rSq, r2)^2 * M^2

  B2B4 <- sum(rSqX2 * r2rSqX2) * M

  B2C1 <- tcrossprod(rSq, r2) * tcrossprod(r,r3) * M^2

  B2C2 <- crossprod(y, r2rSqX2 * rX) * M

  B2C3 <- sum(r * r2 * r3) * rSqSum * M^2

  B2C4 <- crossprod(y, tcrossprod(X, r * r2) * rSqX2) * M

  B3B3 <- B2B2

  B3B4 <- B2B4

  B3C1 <- B2C3

  B3C2 <- B2C4

  B3C3 <- B2C1

  B3C4 <- B2C2



  B4B4 <- sum(Z2^2)

  B4C1 <- sum(rSqX2 * r3rX2) * M

  B4C2 <- quad.3form(Z2 * Z1, e, y)

  B4C3 <- B4C1

  B4C4 <- B4C2

  C1C1 <- rSqSum * sum(r3^2) * M^2

  C1C2 <- crossprod(y, rSqX2 * r3X) * M

  C1C3 <- sum(r * r3)^2 * M^2

  C1C4 <- crossprod(y, rX * r3rX2) * M

  C2C2 <- quad.form(Z2 * XXt, y)

  C2C3 <- C1C4

  C2C4 <- quad.form(Z1*t(Z1), y)

  C3C3 <- C1C1

  C3C4 <- C1C2

  C4C4 <- C2C2

  resvec <- c(AA, c(AB1, AB2,  AB3,  AB4)/4, -1/2*c( AC1,  AC2,  AC3,  AC4),
              c(B1B1, B1B2, B1B3, B1B4)/16, -1/8*c(B1C1, B1C2, B1C3, B1C4),
              c(B2B2, B2B3, B2B4)/16, -1/8*c(B2C1, B2C2, B2C3, B2C4),
              c(B3B3, B3B4)/16, -1/8*c(B3C1, B3C2, B3C3, B3C4),
              B4B4/16, -1/8*c(B4C1, B4C2, B4C3, B4C4),
              c(C1C1, C1C2, C1C3, C1C4,
                C2C2, C2C3, C2C4,
                C3C3, C3C4,
                C4C4)/4)
  resmat <- matrix(0, 9, 9)
  resmat[row(resmat) >= col(resmat)] <- resvec
  trprod <- sum(diag(resmat)) + 2 * sum(resmat[row(resmat) > col(resmat)])

  #Trace of straight up varmat
  star1 <- crossprod(y2, X2RowSums) #O(m)
  star2 <- y4Sum * rSqSum
  star3 <- crossprod(y2, rSqX2) #O(m)
  star4 <- tcrossprod(rSq, X4ColSums) #O(n)
  dagger1 <- crossprod(y3, rX) #O(m)
  dagger2 <- crossprod(y, rX3) #O(m)
  trmat <- star1 + (star2 + 2 * star3 + star4) / 4 - dagger1 - dagger2

  return(c(trmat/M, trprod/M^2))
}

trace_large_x_indx_timer <- function (i) {
  
  timers <- list()
  
  m <- ndata
  y <- Y[ , i]
  e <- matrix(1, nrow=m, ncol=1)
  y2 <- Y2[ , i] # m x 1  ## GET THIS FROM Y2
  y3 <- Y3[ , i] ## GET THIS FROM Y3

  M <- m - 1
  r <- allr[i, , drop = FALSE] / M # 1 x dx ## GET THIS FROM allr
  r2 <- allr22[i, , drop = FALSE] / M # 1 x dx
  r3 <- allr31[i, , drop = FALSE] / M # 1 x dx
  rX <- tcrossprod(X, r) # m x 1
  r3X <- tcrossprod(X, r3) # m x 1
  rX3 <- tcrossprod(X3, r) # m x 1


  rSq <- r^2
  rSqX2 <- tcrossprod(X2,rSq) # m x 1
  r2rSqX2 <- tcrossprod(X2, r2 * rSq) # m x 1
  r3rX2 <- tcrossprod(X2, r3*r) # m x 1

  timers$tXr <- benchmark(
      tXr <- tX * as.vector(r), # dx x m
      replications=1)$elapsed

  y4Sum <- Y4ColSums[i] # scalar ## GET THIS FROM Y4ColSums
  rSqSum <- sum(rSq) #scalar

  timers$Z1 <- benchmark(
      Z1 <- X2 %*% tXr, # m x m
      replications=1)$elapsed
  timers$Z2 <- benchmark(
      Z2 <- X2 %*% tXr^2, # m x m
      replications=1)$elapsed

  #All of these are scalars
  #TODO: Normalize them by (m-1)^2 after they are working correctly.

  timers$AA <- benchmark(
      AA <- quad.form(XXt2, y2),
      replications=1)$elapsed

  timers$AB1 <- benchmark(
      AB1 <- y4Sum * crossprod(y2, rX^2),
      replications=1)$elapsed

  timers$AB2 <- benchmark(
      AB2 <- quad.3form(XXt, y, y2 * tcrossprod(X, r2*r)),
      replications=1)$elapsed

  timers$AB3 <- benchmark(
      AB3 <- AB2,
      replications=1)$elapsed

  timers$AB4 <- benchmark(
      AB4 <- quad.3form(Z1^2, e, y2),
      replications=1)$elapsed

  timers$AC1 <- benchmark(
      AC1 <- crossprod(rX * r3X, y2) * M,
      replications=1)$elapsed

  timers$AC2 <- benchmark(
      AC2 <- quad.3form(Z1 * XXt, y, y2),
      replications=1)$elapsed

  AC3 <- AC1

  AC4 <- AC2

  B1B1 <- (y4Sum*rSqSum)^2

  timers$B1B2 <- benchmark(
      B1B2 <- rSqSum * y4Sum * tcrossprod(r2,rSq) * M,
      replications=1)$elapsed

  B1B3 <- B1B2

  B1B4 <- y4Sum * sum(rSqX2^2)

  timers$B1C1 <- benchmark(
      B1C1 <- y4Sum * rSqSum * tcrossprod(r,r3) * M,
      replications=1)$elapsed

  timers$B1C2 <- benchmark(
      B1C2 <- y4Sum * crossprod(y, rX * rSqX2),
      replications=1)$elapsed

  B1C3 <- B1C1

  B1C4 <- B1C2

  timers$B2B2 <- benchmark(
      B2B2 <- tcrossprod(rSq, r2^2) * rSqSum * M^2,
      replications=1)$elapsed

  timers$B2B3 <- benchmark(
      B2B3 <- tcrossprod(rSq, r2)^2 * M^2,
      replications=1)$elapsed

  B2B4 <- sum(rSqX2 * r2rSqX2) * M

  timers$B2C1 <- benchmark(
      B2C1 <- tcrossprod(rSq, r2) * tcrossprod(r,r3) * M^2,
      replications=1)$elapsed

  timers$B2C2 <- benchmark(
      B2C2 <- crossprod(y, r2rSqX2 * rX) * M,
      replications=1)$elapsed

  B2C3 <- sum(r * r2 * r3) * rSqSum * M^2

  timers$B2C3 <- benchmark(
      B2C4 <- crossprod(y, tcrossprod(X, r * r2) * rSqX2) * M,
      replications=1)$elapsed

  B3B3 <- B2B2

  B3B4 <- B2B4

  B3C1 <- B2C3

  B3C2 <- B2C4

  B3C3 <- B2C1

  B3C4 <- B2C2



  B4B4 <- sum(Z2^2)

  B4C1 <- sum(rSqX2 * r3rX2) * M

  timers$B4C2 <- benchmark(
      B4C2 <- quad.3form(Z2 * Z1, e, y),
      replications=1)$elapsed

  B4C3 <- B4C1

  B4C4 <- B4C2

  C1C1 <- rSqSum * sum(r3^2) * M^2

  timers$C1C2 <- benchmark(
      C1C2 <- crossprod(y, rSqX2 * r3X) * M,
      replications=1)$elapsed

  C1C3 <- sum(r * r3)^2 * M^2

  timers$C1C4 <- benchmark(
      C1C4 <- crossprod(y, rX * r3rX2) * M,
      replications=1)$elapsed

  timers$C2C2 <- benchmark(
      C2C2 <- quad.form(Z2 * XXt, y),
      replications=1)$elapsed

  C2C3 <- C1C4

  timers$C2C4 <- benchmark(
      C2C4 <- quad.form(Z1*t(Z1), y),
      replications=1)$elapsed

  C3C3 <- C1C1

  C3C4 <- C1C2

  C4C4 <- C2C2

  resvec <- c(AA, c(AB1, AB2,  AB3,  AB4)/4, -1/2*c( AC1,  AC2,  AC3,  AC4),
              c(B1B1, B1B2, B1B3, B1B4)/16, -1/8*c(B1C1, B1C2, B1C3, B1C4),
              c(B2B2, B2B3, B2B4)/16, -1/8*c(B2C1, B2C2, B2C3, B2C4),
              c(B3B3, B3B4)/16, -1/8*c(B3C1, B3C2, B3C3, B3C4),
              B4B4/16, -1/8*c(B4C1, B4C2, B4C3, B4C4),
              c(C1C1, C1C2, C1C3, C1C4,
                C2C2, C2C3, C2C4,
                C3C3, C3C4,
                C4C4)/4)
  resmat <- matrix(0, 9, 9)
  resmat[row(resmat) >= col(resmat)] <- resvec
  trprod <- sum(diag(resmat)) + 2 * sum(resmat[row(resmat) > col(resmat)])

  #Trace of straight up varmat
  star1 <- crossprod(y2, X2RowSums) #O(m)
  star2 <- y4Sum * rSqSum
  star3 <- crossprod(y2, rSqX2) #O(m)
  star4 <- tcrossprod(rSq, X4ColSums) #O(n)
  dagger1 <- crossprod(y3, rX) #O(m)
  dagger2 <- crossprod(y, rX3) #O(m)
  trmat <- star1 + (star2 + 2 * star3 + star4) / 4 - dagger1 - dagger2

  return(list(traces=c(trmat/M, trprod/M^2),
              resmat=resmat,
              timers=timers))
}

trace_large_x_indx_faster <- function (i) {
  m <- ndata
  y <- Y[ , i]
  e <- matrix(1, nrow=m, ncol=1)
  y2 <- Y2[ , i] # m x 1  ## GET THIS FROM Y2
  y3 <- Y3[ , i] ## GET THIS FROM Y3

  M <- m - 1
  r <- allr[i, , drop = FALSE] / M # 1 x dx ## GET THIS FROM allr
  r2 <- allr22[i, , drop = FALSE] / M # 1 x dx
  r3 <- allr31[i, , drop = FALSE] / M # 1 x dx
  rX <- tcrossprod(X, r) # m x 1
  r3X <- tcrossprod(X, r3) # m x 1
  rX3 <- tcrossprod(X3, r) # m x 1


  rSq <- r^2
  rSqX2 <- tcrossprod(X2,rSq) # m x 1
  r2rSqX2 <- tcrossprod(X2, r2 * rSq) # m x 1
  r3rX2 <- tcrossprod(X2, r3*r) # m x 1

  tXr <- tX * as.vector(r) # dx x m

  y4Sum <- Y4ColSums[i] # scalar ## GET THIS FROM Y4ColSums
  rSqSum <- sum(rSq) #scalar

  #Z1 <- X2 %*% tXr # m x m
  #Z2 <- X2 %*% tXr^2 # m x m

  #All of these are scalars
  #TODO: Normalize them by (m-1)^2 after they are working correctly.

  AA <- quad.form(XXt2, y2)

  AB1 <- y4Sum * crossprod(y2, rX^2)

  AB2 <- quad.3form(XXt, y, y2 * tcrossprod(X, r2*r))

  AB3 <- AB2

  #AB4 <- quad.3form(Z1^2, e, y2)
  AB4 <- 0

  AC1 <- crossprod(rX * r3X, y2) * M

  #AC2 <- quad.3form(Z1 * XXt, y, y2)
  AC2 <- 0

  AC3 <- AC1

  AC4 <- AC2

  B1B1 <- (y4Sum*rSqSum)^2

  B1B2 <- rSqSum * y4Sum * tcrossprod(r2,rSq) * M

  B1B3 <- B1B2

  B1B4 <- y4Sum * sum(rSqX2^2)

  B1C1 <- y4Sum * rSqSum * tcrossprod(r,r3) * M

  B1C2 <- y4Sum * crossprod(y, rX * rSqX2)

  B1C3 <- B1C1

  B1C4 <- B1C2

  B2B2 <- tcrossprod(rSq, r2^2) * rSqSum * M^2

  B2B3 <- tcrossprod(rSq, r2)^2 * M^2

  B2B4 <- sum(rSqX2 * r2rSqX2) * M

  B2C1 <- tcrossprod(rSq, r2) * tcrossprod(r,r3) * M^2

  B2C2 <- crossprod(y, r2rSqX2 * rX) * M

  B2C3 <- sum(r * r2 * r3) * rSqSum * M^2

  B2C4 <- crossprod(y, tcrossprod(X, r * r2) * rSqX2) * M

  B3B3 <- B2B2

  B3B4 <- B2B4

  B3C1 <- B2C3

  B3C2 <- B2C4

  B3C3 <- B2C1

  B3C4 <- B2C2



  #B4B4 <- sum(Z2^2)
  B4B4 <- 0

  B4C1 <- sum(rSqX2 * r3rX2) * M

  #B4C2 <- quad.3form(Z2 * Z1, e, y)
  B4C2 <- 0

  B4C3 <- B4C1

  B4C4 <- B4C2

  C1C1 <- rSqSum * sum(r3^2) * M^2

  C1C2 <- crossprod(y, rSqX2 * r3X) * M

  C1C3 <- sum(r * r3)^2 * M^2

  C1C4 <- crossprod(y, rX * r3rX2) * M

  #C2C2 <- quad.form(Z2 * XXt, y)
  C2C2 <- 0

  C2C3 <- C1C4

  #C2C4 <- quad.form(Z1*t(Z1), y)
  C2C4 <- 0

  C3C3 <- C1C1

  C3C4 <- C1C2

  C4C4 <- C2C2

  resvec <- c(AA, c(AB1, AB2,  AB3,  AB4)/4, -1/2*c( AC1,  AC2,  AC3,  AC4),
              c(B1B1, B1B2, B1B3, B1B4)/16, -1/8*c(B1C1, B1C2, B1C3, B1C4),
              c(B2B2, B2B3, B2B4)/16, -1/8*c(B2C1, B2C2, B2C3, B2C4),
              c(B3B3, B3B4)/16, -1/8*c(B3C1, B3C2, B3C3, B3C4),
              B4B4/16, -1/8*c(B4C1, B4C2, B4C3, B4C4),
              c(C1C1, C1C2, C1C3, C1C4,
                C2C2, C2C3, C2C4,
                C3C3, C3C4,
                C4C4)/4)
  resmat <- matrix(0, 9, 9)
  resmat[row(resmat) >= col(resmat)] <- resvec
  trprod <- sum(diag(resmat)) + 2 * sum(resmat[row(resmat) > col(resmat)])

  #Trace of straight up varmat
  star1 <- crossprod(y2, X2RowSums) #O(m)
  star2 <- y4Sum * rSqSum
  star3 <- crossprod(y2, rSqX2) #O(m)
  star4 <- tcrossprod(rSq, X4ColSums) #O(n)
  dagger1 <- crossprod(y3, rX) #O(m)
  dagger2 <- crossprod(y, rX3) #O(m)
  trmat <- star1 + (star2 + 2 * star3 + star4) / 4 - dagger1 - dagger2

  return(c(trmat/M, trprod/M^2))
}
