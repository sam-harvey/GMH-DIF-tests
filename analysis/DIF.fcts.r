combined.test <- function(pvalues, L, CovL, alpha = 0.05, show = T) {
  m <- length(L)
  ind <- 1:m
  p <- 0
  sig.items <- NULL


  while (m > 0) {
    a <- as.matrix(L[ind])
    W <- t(a) %*% solve(CovL[ind, ind]) %*% a # df=c
    p <- 1 - pchisq(W, df = length(ind))

    if (p > alpha) {
      break
    } # if(p<=alpha){


    # determine smallest p-avlue
    j <- which.min(pvalues[ind])




    if (show) {
      cat("p-value", p, "chi2", W, "declared item", ind[j], "\n")
    }
    # cat(pvalues[ind])



    if (pvalues[ind][j] > alpha) {
      break
    }


    sig.items <- c(sig.items, ind[j])

    ind <- ind[-j]


    m <- length(L[ind])
  } # end while


  return(sig.items)
} # end



output1 <- function(mat, row.names = NULL, digits = 3, nsmall = 3) {
  hilf <- dim(mat)
  m <- hilf[1]
  n <- hilf[2]

  if (is.null(row.names)) {
    row.names <- rownames(mat)
  }
  output1 <- NULL
  for (i in 1:m) {
    output1 <- c(output1, row.names[i])
    for (j in 1:n) {
      if (!is.na(mat[i, j])) {
        if (j %% 2 == 1) {
          output1 <- c(output1, " & ($", format(mat[i, j], nsmall = nsmall, digits = digits))
        }
        if (j %% 2 == 0) {
          output1 <- c(output1, ",", format(mat[i, j], nsmall = nsmall, digits = digits), "$)")
        }
      } else {
        output1 <- c(output1, " & -- ")
      }
    }
    output1 <- c(output1, "\\\\\n")
  } # end
  cat(output1, "\n")
} # end


output <- function(mat, row.names = NULL, digits = 3, nsmall = 2) {
  hilf <- dim(mat)
  m <- hilf[1]
  n <- hilf[2]

  if (is.null(row.names)) {
    row.names <- rownames(mat)
  }
  output1 <- NULL
  for (i in 1:m) {
    output1 <- c(output1, row.names[i])
    for (j in 1:n) {
      if (!is.na(mat[i, j])) {
        output1 <- c(output1, " & ", format(mat[i, j], nsmall = nsmall, digits = digits))
      } else {
        output1 <- c(output1, " & -- ")
      }
    }
    output1 <- c(output1, "\\\\\n")
  } # end
  cat(output1, "\n")
} # end


obtain.counts <- function(y, c, pairwise = TRUE) {
  # c is number of items
  K <- dim(y)[1]

  # y is K by c matrix
  hilf <- apply(y, 2, sum, na.rm = TRUE)
  X <- hilf
  ns <- apply(!is.na(y), 2, sum)
  X1 <- ns - hilf # negative counts are simply n-X

  # print("hello")
  if (pairwise) {
    ind <- indexfct(c)
    L <- dim(ind)[2]
    X00 <- X10 <- X01 <- X11 <- rep(0, L)
    for (j in 1:L) {
      hilf <- matrix(y[, ind[, j]], K, 2)
      # print(hilf)
      X11[j] <- sum((hilf[, 1] == 1) * (hilf[, 2] == 1), na.rm = TRUE)
      X10[j] <- sum((hilf[, 1] == 1) * (hilf[, 2] == 0), na.rm = TRUE)
      X01[j] <- sum((hilf[, 1] == 0) * (hilf[, 2] == 1), na.rm = TRUE)
      X00[j] <- sum((hilf[, 1] == 0) * (hilf[, 2] == 0), na.rm = TRUE)
      # cat(X00[j],X10[j],X01[j],X11[j])
    } # end for j
    # print("hello1")
  } # end
  return(list(X = X, X1 = X1, X11 = X11, X10 = X10, X01 = X01, X00 = X00))
} #


############################################################################################################
# obtain marginal X,X1 and pairwise (X11,X10,X01,X00) counts directly and then submit to function MH.est()
# create K * 2 * 2^m data array
obtain.counts.y <- function(y, K, c) {
  X <- X1 <- matrix(NA, 2 * K, c)
  X00 <- X11 <- X01 <- X10 <- matrix(NA, 2 * K, c * (c - 1) / 2)
  nk <- rep(NA, 2 * K)

  for (k in 1:K) {
    for (i in 1:2) {
      # selects those rows with Total.score== i+4 and domestic or international
      ind <- y[, 1] == k & y[, 2] == i
      nk1 <- sum(ind)
      if (nk1 > 0) {

        # print(nk1)
        y.hilf <- matrix(y[ind, 3:(c + 2)], nk1, c)
        # print(y.hilf)
        hilf <- obtain.counts(y.hilf, c, pairwise = TRUE)
        X[2 * (k - 1) + i, ] <- hilf$X
        X1[2 * (k - 1) + i, ] <- hilf$X1

        X00[2 * (k - 1) + i, ] <- hilf$X00
        X01[2 * (k - 1) + i, ] <- hilf$X01
        X10[2 * (k - 1) + i, ] <- hilf$X10
        X11[2 * (k - 1) + i, ] <- hilf$X11
      } else {
        X[2 * (k - 1) + i, ] <- 0
        X1[2 * (k - 1) + i, ] <- 0

        X00[2 * (k - 1) + i, ] <- 0
        X01[2 * (k - 1) + i, ] <- 0
        X10[2 * (k - 1) + i, ] <- 0
        X11[2 * (k - 1) + i, ] <- 0
      }
      nk[2 * (k - 1) + i] <- nk1
    } # end i
  } # end k
  return(list(X = X, X1 = X1, X11 = X11, X10 = X10, X01 = X01, X00 = X00, nk = nk))
}
############################################################################################################
OR.pair <- function(X11, X10, X01, X00, K = dim(X11)[1] / 2, nk) {
  num <- den <- 0
  for (k in 1:K) {
    for (i in 1:2) {
      ind <- 2 * (k - 1) + i
      if (nk[ind] > 0) {
        num <- num + X11[, ind] * X00[, ind] / nk[ind]
        den <- den + X10[, ind] * X01[, ind] / nk[ind]
      }
    } #
  } #
  OR <- num / den
  OR[den == 0] <- NA
  OR[num == 0] <- NA
  return(OR)
} # end


obtain.bootstrap.samples.ind <- function(y, within.group = TRUE, K) {
  hilf <- dim(y)
  n <- hilf[1]
  m <- hilf[2] - 2

  y.new <- matrix(0, 0, m + 2)

  for (k in 1:K) {
    if (within.group) {
      for (i in 1:2) {


        # selects those rows with Total.score== i+4 and domestic or international
        ind <- y[, 1] == k & y[, 2] == i
        nk <- sum(ind)

        if (nk > 0) {
          hilf <- NULL
          for (j in 1:m) {
            s <- sample((1:n)[ind], nk, replace = TRUE)
            hilf.add <- y[s, 2 + j]



            hilf <- cbind(hilf, hilf.add)
          } # end

          # print(nk)
          # print(hilf)
          # print(dim(hilf))


          y.new <- rbind(y.new, cbind(k, i, hilf))
          # print(dim(y.new))
        } # end if nk>0
      } # end i
    } else {
      # whole strata

      ind <- y[, 1] == k
      nk <- sum(ind)

      ind1 <- y[, 1] == k & y[, 2] == 1
      nk1 <- sum(ind1)




      hilf <- NULL
      for (j in 1:m) {
        s <- sample((1:n)[ind], nk, replace = TRUE)
        hilf <- cbind(hilf, y[s, 2 + j])
      } # end


      # we sample from whole stratum, but we assign firt nk1 to group 1 and next nk2 to group2
      # nk2=nk- nk2
      # print(nk1)
      # print(nk-nk1)
      print(dim(hilf))

      hilf <- matrix(hilf, nk, m)

      hilf2 <- cbind(k, c(rep(1, nk1), rep(2, nk - nk1)), hilf)
      # print(dim(hilf2))

      y.new <- rbind(y.new, hilf2)
    } # end else
  } # end k
  return(y.new)
} # end function



obtain.bootstrap.samples <- function(y, within.group = TRUE, K = K) {
  hilf <- dim(y)
  n <- hilf[1]
  m <- hilf[2] - 2

  y.new <- matrix(0, 0, m + 2)

  for (k in 1:K) {
    if (within.group) {
      for (i in 1:2) {


        # selects those rows with Total.score== i+4 and domestic or international
        ind <- y[, 1] == k & y[, 2] == i
        nk <- sum(ind)

        s <- sample((1:n)[ind], nk, replace = TRUE)

        y.new <- rbind(y.new, cbind(k, i, y[s, 3:(m + 2)]))
      } # end i
    } else {
      # whole strata

      ind <- y[, 1] == k
      nk <- sum(ind)

      ind1 <- y[, 1] == k & y[, 2] == 1
      nk1 <- sum(ind1)


      # print(nk);print(nk1);print(k);print(dim(y.new))


      s <- NULL
      while (is.null(s)) {
        s <- sample((1:n)[ind], nk, replace = TRUE)
        # make sure in each group is at least 1 observation
        # if(sum(y[s,2]==1)>0 & sum(y[s,2]==2)>0){s<-s}else{s<-NULL}}

        # we sample from whole stratum, but we assign first nk1 to group 1 and next nk2 to group2
        # nk2=nk- nk2
        # print(s)
        # print(nk)
        y.add <- matrix(y[s, 3:(m + 2)], nk, m)

        h <- cbind(k, c(rep(1, nk1), rep(2, nk - nk1)), y.add)
        # print(dim(h))
        y.new <- rbind(y.new, h)
      }
    } # end else
  } # end k
  return(y.new)
}
