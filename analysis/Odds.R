# There are two sets of functions: the list functions and matrix functions
# The list functions are ones that were created for the language dataset
# The matrix functions are those for the contraception and child datasets.
#
#   NOTE For the .matrix functions
#
#   It is necessary that the data matrix be of the following form.
#   However it can be entered in the data-only form and the convert
#   function can then be used to add the item-rows to the matrix.
#
#   Examples of data (info) matrix
#   (1 <- Yes, 0 <- No)
#
#   Required form			   Form for convert function
#
#   Item 1   0 0 0 0 1 1 1 1
#   Item 2   0 0 1 1 0 0 1 1
#   Item 3   0 1 0 1 0 1 0 1
#   Sub 1 N  0 3 0 2 4 0 3 0               Sub 1 N  0 3 0 2 4 0 3 0
#   Sub 1 Y  1 2 0 0 0 0 3 0               Sub 1 Y  1 2 0 0 0 0 3 0
#   Sub 2 N  3 0 2 3 0 0 1 2               Sub 2 N  3 0 2 3 0 0 1 2
#   Sub 2 Y  1 1 1 0 0 1 1 1               Sub 2 Y  1 1 1 0 0 1 1 1
#
#
# FULL FUNCTION LIST:
#
# convert(info,n.items)
#
# b.s.matrix(info,n.items,n.ans,n.sub,n.sets,cov.mat<-NULL)
#
# b.s.list(info,n.items,n.ans,n.sub,n.sets,cov.mat<-NULL)
#
# m.h.matrix(info,n.items,n.ans,n.sub)
#
# m.h.list(info,n.items,n.ans,n.sub)
#
# subtable.matrix(info,item,row1,row2,sub,n.items,n.ans)
#
# subtable.list(info,item,row1,row2,sub,n.sub) #
##
## gen(n.items)
##
# bin.s(val,divs)
#
# cov.bs(bs.mat1,bs.mat2,n.sets,n.items)
#
#
#
# PRIMARY FUNCTIONS:
#
# convert(info,n.items)
# For use with data of the .matrix form. Converts a matrix in ordered, data-only form to the required
# form by adding a binary coding at the top. This should be used for .matrix functions.
#
# b.s.matrix(info,n.items,n.ans,n.sub,n.sets,cov.mat<-NULL)
# Based on the entered data, performs the bootstrap method by repeatedly
# resampling with replacement from the collected data. For datasets of
# the form of the contraception and child datasets
#
# b.s.list(info,n.items,n.ans,n.sub,n.sets,cov.mat<-NULL)
# Based on the entered data, performs the bootstrap method by repeatedly
# resampling with replacement from the collected data. For datasets of
# the form of the language dataset
#
# m.h.matrix(info,n.items,n.ans,n.sub)
# Calculates the generalised Mantel-Haenzsel statistic, its log and the variance of
# its log for datasets of the form of the contraception and child datasets.
#
# m.h.list(info,n.items,n.ans,n.sub)
# Calculates the generalised Mantel-Haenzsel statistic, its log and the variance of
# its log for datasets of the form of the language dataset.
#
#
# SECONDARY FUNCTIONS (should not be called):
#
# subtable.matrix(info,item,row1,row2,sub,n.items,n.ans)
# Calculates a subtable for the m.h.matrix function
#
# subtable.list(info,item,row1,row2,sub,n.sub)
# Calculates a subtable for the m.h.list function
#
# gen(n.items)
# Used in the convert function
#
# bin.s(val,divs)
# Used to generate a boostrap sample
#
# cov.bs(bs.mat1,bs.mat2,n.sets,n.items)
# Calculates the covariance in the bs functions
#


# PARAMETERS
# info	  <- matrix consisting of source data in form as above
# n.items  <- total number of binary-selectable items for data
# n.sub    <- number of subgroups required
# n.ans    <- number of possible responses that be given
# n.data   <- desired subgroup size for bootstrap method
# n.sets   <- required number of sets of data resampled for bootstrap method
# n.sims   <- required number of repetitions of bootstrap method
# theta    <- odds ratio for each item
# phi      <- probability ratio for 2 items
# previous <- list of filenames with previously simulated data to be loaded
# output   <- filename to save simulation data to
# sub      <- Numeric reference to a single subgroup, as needed for subtable
# val	  <- Value required to be sorted in bin.s
# divs	  <- Sorting boundarys, in ascending order, as required for bin.s
# row1	  <- row number for calculations
# row2     <- row number for calculations
# bs.mat1  <- Bootstrap data in matrix form, as required for cov.bs
# bs.mat2  <- Bootstrap data in matrix form, as required for cov.bs


#################################################

# FUNCTION CALLS

# For language dataset

# b.s.list(language,7,3,6,100)
## b.s.matrix(info,n.items,n.ans,n.sub,n.sets,cov.mat<-NULL)
# m.h.list(language,7,3,6)
## m.h.list(info,n.items,n.ans,n.sub)

# For child dataset

# b.s.matrix(child,5,3,10,100)
## b.s.list(info,n.items,n.ans,n.sub,n.sets,cov.mat<-NULL)
# m.h.matrix(child,5,3,10)
## m.h.matrix(info,n.items,n.ans,n.sub)

##################################################

## PRIMARY FUNCTIONS

##
convert <-
  function(info, n.items) {
    return(rbind(gen(n.items), info))
  }

##
b.s.matrix <-
  function(info, n.items, n.ans, n.sub, n.sets, cov.mat = NULL) {
    sub.samp <- NULL
    val <- NULL
    for (i in 1:n.sub) {
      val <- c(val, sum(info[n.items + n.ans * (i - 1) + (1:n.ans), ]))
    }
    sec <- 1 / val
    sub.samp <- NULL
    # Create probabilities for generating bootstrap samples
    for (k in 1:n.sub) {
      pop <- NULL
      p.val <- 0
      for (i in 1:(2^n.items)) {
        for (j in 1:n.ans) {
          if (info[n.items + (k - 1) * n.ans + j, i] != 0) {
            temp <- matrix(nrow <- 4, ncol <- 1, data <- c(p.val, i, j, 0))
            pop <- cbind(pop, temp)
            p.val <- p.val + info[n.items + (k - 1) * n.ans + j, i] * sec[k]
          }
        }
      }
      pop <- cbind(pop, matrix(nrow <- 4, ncol <- 1, data <- c(1, 0, 0, 0)))
      sub.samp <- c(sub.samp, list(pop))
    }
    ests <- list()
    len <- n.ans * (n.ans - 1) / 2
    length(ests) <- len
    l.ests <- list()
    length(l.ests) <- len
    ests.v <- list()
    length(ests.v) <- len
    rownam <- NULL
    # Generate bootstrap samples
    for (b in 1:n.sets) {
      res <- info[1:n.items, ]
      # Generate sample
      for (k in 1:n.sub) {
        curr <- sub.samp[[k]]
        for (h in 1:val[k]) {
          unit <- runif(1)
          pos <- bin.s(unit, curr[1, ])
          curr[4, pos] <- curr[4, pos] + 1
        }
        temp <- matrix(nrow <- n.ans, ncol <- 2^n.items, data <- 0)
        for (h in 1:ncol(curr)) {
          temp[curr[3, h], curr[2, h]] <- curr[4, h]
        }
        res <- rbind(res, temp)
      }
      # Calculate M-H estimates
      temp <- m.h.matrix(res, n.items, n.ans, n.sub)
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
          ests[[pos]] <-
            rbind(ests[[pos]], matrix(nrow <- 1, ncol <- n.items, data <- c(temp[[pos]][[2]])))
          l.ests[[pos]] <-
            rbind(l.ests[[pos]], matrix(nrow <- 1, ncol <- n.items, data <- c(temp[[pos]][[3]])))
          ests.v[[pos]] <-
            rbind(ests.v[[pos]], matrix(nrow <- 1, ncol <- n.items, data <- c(temp[[pos]][[4]])))
        }
      }
      rownam <- c(rownam, list(paste("Data Set ", b, sep <- "")))
    }
    colnam <- NULL
    l.colnam <- NULL
    cov.row <- NULL
    cov.col <- NULL
    se <- list()
    length(se) <- len
    cov <- list()
    out <- list()
    # Calculate covariance matricies
    if (is.null(cov.mat)) {
      length(cov) <- len
      for (i in 1:len) {
        cov[[i]] <- cov.bs(l.ests[[i]], l.ests[[i]], n.sets, n.items)
        se[[i]] <- matrix(nrow <- 1, ncol <- n.items)
        for (k in 1:n.items) {
          se[[i]][1, k] <- sqrt(cov[[i]][k, k])
        }
      }
      for (i in 1:n.items) {
        colnam <- c(colnam, list(paste("Theta-hat ", i, sep <- "")))
        l.colnam <- c(l.colnam, list(paste("Log-theta ", i, sep <- "")))
        cov.col <- c(cov.col, list(paste("Cov i, ", i, sep <- "")))
        cov.row <- c(cov.row, list(paste("Cov ", i, ", j", sep <- "")))
      }
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
          dimnames(se[[pos]]) <- list("SE Bootstrap", l.colnam)
          Parameters <- paste("Row ", i, " & Row ", j)
          out <- c(out, list(list(
            Parameters <- Parameters,
            "Standard errors" <- se[[pos]]
          )))
        }
      }
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * (n.ans) - i) + j - i
          Parameters <- paste("Row ", i, " & Row ", j)
          dimnames(cov[[pos]]) <- list(cov.row, cov.col)
          out <- c(out, list(list(Parameters <- Parameters, "Covariance Matrix" <- cov[[pos]])))
        }
      }
    } else {
      for (i in 1:len) {
        temp <- cov.bs(l.ests[[i]], l.ests[[i]], n.sets, n.items)
        se[[i]] <- matrix(nrow <- 1, ncol <- n.items)
        for (k in 1:n.items) {
          se[[i]][1, k] <- sqrt(temp[k, k])
        }
      }
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
          dimnames(se[[pos]]) <- list("SE Bootstrap", l.colnam)
          Parameters <- paste("Row ", i, " & Row ", j)
          out <- c(out, list(list(
            Parameters <- Parameters,
            "Standard errors" <- se[[pos]]
          )))
        }
      }
      covscol <- rep(0, nrow(cov.mat))
      for (z in 1:nrow(cov.mat)) {
        tot <- 0
        n <- cov.mat[z, 5]
        o <- cov.mat[z, 6]
        pos1 <- (cov.mat[z, 1] - 1) / 2 * (2 * (n.ans) - cov.mat[z, 1]) + cov.mat[z, 2] - cov.mat[z, 1]
        val.bar1 <- mean(l.ests[[pos1]][, n])
        pos2 <- (cov.mat[z, 3] - 1) / 2 * (2 * (n.ans) - cov.mat[z, 3]) + cov.mat[z, 4] - cov.mat[z, 3]
        val.bar2 <- mean(l.ests[[pos2]][, o])
        for (m in 1:n.sets) {
          curr <- (l.ests[[pos1]][m, n] - val.bar1) * (l.ests[[pos2]][m, o] - val.bar2)
          tot <- tot + curr
        }
        covscol[z] <- tot / (n.sets - 1)
      }
      cov.mat <- cbind(cov.mat, covscol)
      dimnames(cov.mat) <- list(NULL, list("i", "h", "i'", "h' ", "j", "j'", "Cov"))
      out <- c(out, list(list("Requested Covariances" <- cov.mat)))
    }
    return(out)
  }


b.s.list <-
  function(info, n.items, n.ans, n.sub, n.sets, cov.mat = NULL) {
    sub.samp <- NULL
    val <- nrow(info) / n.sub
    ests <- list()
    len <- n.ans * (n.ans - 1) / 2
    length(ests) <- len
    l.ests <- list()
    length(l.ests) <- len
    ests.v <- list()
    length(ests.v) <- len
    rownam <- NULL
    res <- info
    # print("hello")
    # Generate bootstrap samples
    for (b in 1:n.sets) {
      # Generate single bootstrap sample
      for (i in 0:(val - 1)) {
        unit <- floor(runif(1) * val)
        res[n.sub * i + 1:n.sub, ] <- info[unit * n.sub + 1:n.sub, ]
      }
      # Calculate m-h estimates
      temp <- m.h.list(res, n.items, n.ans, n.sub)
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
          ests[[pos]] <-
            rbind(ests[[pos]], matrix(nrow <- 1, ncol <- n.items, data <- c(temp[[pos]][[2]])))
          l.ests[[pos]] <-
            rbind(l.ests[[pos]], matrix(nrow <- 1, ncol <- n.items, data <- c(temp[[pos]][[3]])))
          ests.v[[pos]] <-
            rbind(ests.v[[pos]], matrix(nrow <- 1, ncol <- n.items, data <- c(temp[[pos]][[4]])))
        }
      }
      rownam <- c(rownam, list(paste("Data Set ", b, sep <- "")))
    }
    # print("hello")
    colnam <- NULL
    l.colnam <- NULL
    cov.row <- NULL
    cov.col <- NULL
    se <- list()
    length(se) <- len
    cov <- list()
    out <- list()
    # Calculate covariance matricies
    if (is.null(cov.mat)) {
      length(cov) <- len
      for (i in 1:len) {
        cov[[i]] <- cov.bs(l.ests[[i]], l.ests[[i]], n.sets, n.items)
        se[[i]] <- matrix(nrow <- 1, ncol <- n.items)
        for (k in 1:n.items) {
          se[[i]][1, k] <- sqrt(cov[[i]][k, k])
        }
      }
      for (i in 1:n.items) {
        colnam <- c(colnam, list(paste("Theta-hat ", i, sep <- "")))
        l.colnam <- c(l.colnam, list(paste("Log-theta ", i, sep <- "")))
        cov.col <- c(cov.col, list(paste("Cov i, ", i, sep <- "")))
        cov.row <- c(cov.row, list(paste("Cov ", i, ", j", sep <- "")))
      }
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
          dimnames(se[[pos]]) <- list("SE Bootstrap", l.colnam)
          Parameters <- paste("Row ", i, " & Row ", j)
          out <- c(out, list(list(
            Parameters <- Parameters,
            "Standard errors" <- se[[pos]]
          )))
        }
      }
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * (n.ans) - i) + j - i
          Parameters <- paste("Row ", i, " & Row ", j)
          dimnames(cov[[pos]]) <- list(cov.row, cov.col)
          out <- c(out, list(list(Parameters <- Parameters, "Covariance Matrix" <- cov[[pos]])))
        }
      }
    } else {
      for (i in 1:len) {
        temp <- cov.bs(l.ests[[i]], l.ests[[i]], n.sets, n.items)
        se[[i]] <- matrix(nrow <- 1, ncol <- n.items)
        for (k in 1:n.items) {
          se[[i]][1, k] <- sqrt(temp[k, k])
        }
      }
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
          dimnames(se[[pos]]) <- list("SE Bootstrap", l.colnam)
          Parameters <- paste("Row ", i, " & Row ", j)
          out <- c(out, list(list(
            Parameters <- Parameters,
            "Standard errors" <- se[[pos]]
          )))
        }
      }
      covscol <- rep(0, nrow(cov.mat))
      for (z in 1:nrow(cov.mat)) {
        tot <- 0
        n <- cov.mat[z, 5]
        o <- cov.mat[z, 6]
        pos1 <- (cov.mat[z, 1] - 1) / 2 * (2 * (n.ans) - cov.mat[z, 1]) + cov.mat[z, 2] - cov.mat[z, 1]
        val.bar1 <- mean(l.ests[[pos1]][, n])
        pos2 <- (cov.mat[z, 3] - 1) / 2 * (2 * (n.ans) - cov.mat[z, 3]) + cov.mat[z, 4] - cov.mat[z, 3]
        val.bar2 <- mean(l.ests[[pos2]][, o])
        for (m in 1:n.sets) {
          curr <- (l.ests[[pos1]][m, n] - val.bar1) * (l.ests[[pos2]][m, o] - val.bar2)
          tot <- tot + curr
        }
        covscol[z] <- tot / (n.sets - 1)
      }
      cov.mat <- cbind(cov.mat, covscol)
      dimnames(cov.mat) <- list(NULL, list("i", "h", "i'", "h'", "j", "j'", "Cov"))
      out <- c(out, list(list("Requested Covariances" <- cov.mat)))
    }
    return(out)
  }


##
m.h.matrix <-
  function(info, n.items, n.ans, n.sub) {
    out <- NULL
    colnam <- NULL
    for (i in 1:n.items) {
      colnam <- c(colnam, list(paste("Theta-hat ", i, sep = "")))
    }
    # Select 2 rows to calculate estimates for
    for (row1 in 1:(n.ans - 1)) {
      for (row2 in (row1 + 1):n.ans) {
        pos <- (row2 - 1) / 2 * (2 * n.ans - row1) + row2 - row1
        est <- matrix(nrow = 1, ncol = n.items)
        l.est <- matrix(nrow = 1, ncol = n.items)
        var <- matrix(nrow = 1, ncol = n.items)
        # Calculate estimate for each item
        for (i in 1:n.items) {
          top <- 0
          bot <- 0
          sum1T <- 0
          sum2T <- 0
          sum3T <- 0
          sumaB <- 0
          sumbB <- 0
          # Calculate contribution for each subgroup
          for (k in 1:n.sub) {
            # cat(i,row1,row2,k,n.items,n.ans,'\n')
            temp <- subtable.matrix(info, i, row1, row2, k, n.items, n.ans)
            total <- sum(info[n.items + n.ans * (k - 1) + (1:n.ans), ])
            # print(total)
            top <- top + (temp[1, 1] * temp[2, 2]) / total
            bot <- bot + (temp[1, 2] * temp[2, 1]) / total
            p1122 <- temp[1, 1] * temp[2, 2]
            p1221 <- temp[1, 2] * temp[2, 1]
            sum1T <- sum1T + (temp[1, 1] + temp[2, 2]) * (p1122) / (total^2)
            sumaB <- sumaB + (p1122 / total)
            int <- (temp[1, 1] + temp[2, 2]) * (p1221)
            int <- int + (temp[1, 2] + temp[2, 1]) * (p1122)
            sum2T <- sum2T + int / (total^2)
            sum3T <- sum3T + (temp[1, 2] + temp[2, 1]) * (p1221) / (total^2)
            sumbB <- sumbB + (p1221 / total)
          }
          if (top != 0 && bot != 0) {
            est[1, i] <- top / bot
            l.est[1, i] <- log(est[1, i])
            var[1, i] <- sum1T / 2 / (sumaB)^2 + sum2T / 2 / sumaB / sumbB + sum3T / 2 / (sumbB)^2
          }
          # Adjust largest cell if 0 on top or bottom of estimate
          else {
            largest <- 0
            for (k in 1:n.sub) {
              temp <- subtable.matrix(info, i, row1, row2, k, n.items, n.ans)
              val <- sum(temp)
              if (val > largest) {
                pos <- k
                largest <- val
              }
            }


            # temp <- subtable.matrix(info,i,row1,row2,k,n.items,n.ans)
            # total <- sum(info[n.items+n.ans*(i-1)+1:n.ans,]) #original

            temp <- subtable.matrix(info, i, row1, row2, pos, n.items, n.ans)
            print(temp)
            total <- sum(info[n.items + n.ans * (pos - 1) + (1:n.ans), ])


            # print(total)
            top <- top - (temp[1, 1] * temp[2, 2]) / total
            bot <- bot - (temp[1, 2] * temp[2, 1]) / total
            p1122 <- temp[1, 1] * temp[2, 2]
            p1221 <- temp[1, 2] * temp[2, 1]
            sum1T <- sum1T - (temp[1, 1] + temp[2, 2]) * (p1122) / (total^2)
            sumaB <- sumaB - (p1122 / total)
            int <- (temp[1, 1] + temp[2, 2]) * (p1221)
            int <- int - (temp[1, 2] + temp[2, 1]) * (p1122)
            sum2T <- sum2T - int / (total^2)
            sum3T <- sum3T - (temp[1, 2] + temp[2, 1]) * (p1221) / (total^2)
            sumbB <- sumbB - (p1221 / total)
            temp <- temp + 0.5
            total <- total + 2
            top <- top + (temp[1, 1] * temp[2, 2]) / total
            bot <- bot + (temp[1, 2] * temp[2, 1]) / total
            est[1, i] <- top / bot
            l.est[1, i] <- log(est[1, i])
            p1122 <- temp[1, 1] * temp[2, 2]
            p1221 <- temp[1, 2] * temp[2, 1]
            sum1T <- sum1T + (temp[1, 1] + temp[2, 2]) * (p1122) / (total^2)
            sumaB <- sumaB + (p1122 / total)
            int <- (temp[1, 1] + temp[2, 2]) * (p1221)
            int <- int + (temp[1, 2] + temp[2, 1]) * (p1122)
            sum2T <- sum2T + int / (total^2)
            sum3T <- sum3T + (temp[1, 2] + temp[2, 1]) * (p1221) / (total^2)
            sumbB <- sumbB + (p1221 / total)
            var[1, i] <- sum1T / 2 / (sumaB)^2 + sum2T / 2 / sumaB / sumbB + sum3T / 2 / (sumbB)^2
          }
        }
        dimnames(est) <- list("Theta", colnam)
        dimnames(l.est) <- list("Log-theta", colnam)
        dimnames(var) <- list("Variance", colnam)
        Parameters <- paste("Row ", row1, " & Row ", row2)
        out <- c(out, list(list(
          Parameters = Parameters, "Thetas" = est, "Log thetas" = l.est,
          "Log theta variances" = var
        )))
      }
    }
    inter <- list()
    length(inter) <- n.ans * (n.ans - 1) / 2
    for (i in 1:length(inter)) {
      inter[[i]] <- matrix(nrow = 1, ncol = n.items, data = 0)
    }
    # Calculate generalised estimates
    for (k in 1:n.items) {
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
          for (h in 1:n.ans) {
            if (i < h) {
              inter[[pos]][1, k] <- inter[[pos]][1, k] + out[[(i - 1) / 2 * (2 * n.ans - i) + h - i]][[3]][1, k]
            } else if (i > h) {
              inter[[pos]][1, k] <- inter[[pos]][1, k] - out[[(h - 1) / 2 * (2 * n.ans - h) + i - h]][[3]][1, k]
            }
            if (j < h) {
              inter[[pos]][1, k] <- inter[[pos]][1, k] - out[[(j - 1) / 2 * (2 * n.ans - j) + h - j]][[3]][1, k]
            } else if (j > h) {
              inter[[pos]][1, k] <- inter[[pos]][1, k] + out[[(h - 1) / 2 * (2 * n.ans - h) + j - h]][[3]][1, k]
            }
          }
        }
      }
    }
    for (i in 1:(n.ans - 1)) {
      for (j in (i + 1):n.ans) {
        pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
        out[[pos]][[3]][] <- inter[[pos]][] / n.ans
      }
    }
    return(out)
  }

##
m.h.list <-
  function(info, n.items, n.ans, n.sub) {
    # m.h.list(data,c,r,K)
    # n.ans   ... number of rows r
    # n.sub       K
    # n.c         c

    out <- NULL
    colnam <- NULL
    for (i in 1:n.items) {
      colnam <- c(colnam, list(paste("Theta-hat ", i, sep = "")))
    }
    # print("hello")
    # Calculate for each pair of rowa
    for (row1 in 1:(n.ans - 1)) {
      for (row2 in (row1 + 1):n.ans) {
        pos <- (row2 - 1) / 2 * (2 * n.ans - row1) + row2 - row1
        est <- matrix(nrow = 1, ncol = n.items)
        l.est <- matrix(nrow = 1, ncol = n.items)
        var <- matrix(nrow = 1, ncol = n.items)
        # Calculate for each item
        for (i in 1:n.items) {
          top <- 0
          bot <- 0
          sum1T <- 0
          sum2T <- 0
          sum3T <- 0
          sumaB <- 0
          sumbB <- 0
          # Calculate contribution for each subgroup
          for (k in 1:n.sub) {
            temp <- subtable.list(info, i, row1, row2, k, n.sub)
            total <- nrow(info) / n.sub
            top <- top + (temp[1, 1] * temp[2, 2]) / total
            bot <- bot + (temp[1, 2] * temp[2, 1]) / total
            p1122 <- temp[1, 1] * temp[2, 2]
            p1221 <- temp[1, 2] * temp[2, 1]
            sum1T <- sum1T + (temp[1, 1] + temp[2, 2]) * (p1122) / (total^2)
            sumaB <- sumaB + (p1122 / total)
            int <- (temp[1, 1] + temp[2, 2]) * (p1221)
            int <- int + (temp[1, 2] + temp[2, 1]) * (p1122)
            sum2T <- sum2T + int / (total^2)
            sum3T <- sum3T + (temp[1, 2] + temp[2, 1]) * (p1221) / (total^2)
            sumbB <- sumbB + (p1221 / total)
          }
          if (top != 0 && bot != 0) {
            est[1, i] <- top / bot
            l.est[1, i] <- log(est[1, i])
            var[1, i] <- sum1T / 2 / (sumaB)^2 + sum2T / 2 / sumaB / sumbB + sum3T / 2 / (sumbB)^2
          }
          # Adjust largest cell if 0 on top or bottom of estimate
          else {
            largest <- 0
            for (k in 1:n.sub) {
              temp <- subtable.list(info, i, row1, row2, k, n.sub)
              val <- sum(temp)
              if (val > largest) {
                pos <- k
                largest <- val
              }
            }
            temp <- subtable.list(info, i, row1, row2, pos, n.sub)
            total <- nrow(info) / n.sub
            top <- top - (temp[1, 1] * temp[2, 2]) / total
            bot <- bot - (temp[1, 2] * temp[2, 1]) / total
            p1122 <- temp[1, 1] * temp[2, 2]
            p1221 <- temp[1, 2] * temp[2, 1]
            sum1T <- sum1T - (temp[1, 1] + temp[2, 2]) * (p1122) / (total^2)
            sumaB <- sumaB - (p1122 / total)
            int <- (temp[1, 1] + temp[2, 2]) * (p1221)
            int <- int - (temp[1, 2] + temp[2, 1]) * (p1122)
            sum2T <- sum2T - int / (total^2)
            sum3T <- sum3T - (temp[1, 2] + temp[2, 1]) * (p1221) / (total^2)
            sumbB <- sumbB - (p1221 / total)
            temp <- temp + 0.5
            total <- total + 2
            top <- top + (temp[1, 1] * temp[2, 2]) / total
            bot <- bot + (temp[1, 2] * temp[2, 1]) / total
            est[1, i] <- top / bot
            l.est[1, i] <- log(est[1, i])
            p1122 <- temp[1, 1] * temp[2, 2]
            p1221 <- temp[1, 2] * temp[2, 1]
            sum1T <- sum1T + (temp[1, 1] + temp[2, 2]) * (p1122) / (total^2)
            sumaB <- sumaB + (p1122 / total)
            int <- (temp[1, 1] + temp[2, 2]) * (p1221)
            int <- int + (temp[1, 2] + temp[2, 1]) * (p1122)
            sum2T <- sum2T + int / (total^2)
            sum3T <- sum3T + (temp[1, 2] + temp[2, 1]) * (p1221) / (total^2)
            sumbB <- sumbB + (p1221 / total)
            var[1, i] <- sum1T / 2 / (sumaB)^2 + sum2T / 2 / sumaB / sumbB + sum3T / 2 / (sumbB)^2
          }
        }
        colnam <- NULL
        for (i in 1:n.items) {
          colnam <- c(colnam, list(paste("Theta-hat ", i, sep <- "")))
        }
        dimnames(est) <- list("Theta", colnam)
        dimnames(l.est) <- list("Log-theta", colnam)
        dimnames(var) <- list("Variance", colnam)
        Parameters <- paste("Row ", row1, " & Row ", row2)
        out <- c(out, list(list(
          Parameters <- Parameters, "Thetas" <- est, "Log thetas" <- l.est,
          "Log theta variances" <- var
        )))
      }
    }
    inter <- list()
    length(inter) <- n.ans * (n.ans - 1) / 2
    for (i in 1:length(inter)) {
      inter[[i]] <- matrix(nrow = 1, ncol = n.items, data = 0)
    }
    # Calculate generalised estimates
    for (k in 1:n.items) {
      for (i in 1:(n.ans - 1)) {
        for (j in (i + 1):n.ans) {
          pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
          for (h in 1:n.ans) {
            if (i < h) {
              inter[[pos]][1, k] <- inter[[pos]][1, k] + out[[(i - 1) / 2 * (2 * n.ans - i) + h - i]][[3]][1, k]
            } else if (i > h) {
              inter[[pos]][1, k] <- inter[[pos]][1, k] - out[[(h - 1) / 2 * (2 * n.ans - h) + i - h]][[3]][1, k]
            }
            if (j < h) {
              inter[[pos]][1, k] <- inter[[pos]][1, k] - out[[(j - 1) / 2 * (2 * n.ans - j) + h - j]][[3]][1, k]
            } else if (j > h) {
              inter[[pos]][1, k] <- inter[[pos]][1, k] + out[[(h - 1) / 2 * (2 * n.ans - h) + j - h]][[3]][1, k]
            }
          }
        }
      }
    }
    for (i in 1:(n.ans - 1)) {
      for (j in (i + 1):n.ans) {
        pos <- (i - 1) / 2 * (2 * n.ans - i) + j - i
        out[[pos]][[3]][] <- inter[[pos]][] / n.ans
      }
    }
    return(out)
  }

## SECONDARY FUNCTIONS

subtable.matrix <-
  function(info, item, row1, row2, sub, n.items, n.ans) {
    # cat('\n',item,row1,row2,sub,n.items,n.ans,'\n')
    tab <- matrix(nrow = 2, ncol = 2, data = 0)
    for (z in 1:ncol(info)) {
      if (info[item, z] == 1) {
        tab[1, 1] <- tab[1, 1] + info[n.items + (sub - 1) * n.ans + row1, z]
        tab[2, 1] <- tab[2, 1] + info[n.items + (sub - 1) * n.ans + row2, z]
      } else {
        tab[1, 2] <- tab[1, 2] + info[n.items + (sub - 1) * n.ans + row1, z]
        tab[2, 2] <- tab[2, 2] + info[n.items + (sub - 1) * n.ans + row2, z]
      }
    }
    return(tab)
  }

subtable.list <-
  function(info, item, row1, row2, sub, n.sub) {
    tab <- matrix(nrow = 2, ncol = 2, data = 0)
    for (z in 0:(nrow(info) / n.sub - 1)) {
      if (info[z * n.sub + sub, 3] == row1) {
        if (info[z * n.sub + sub, item + 3] == 1) {
          tab[1, 1] <- tab[1, 1] + 1
        } else {
          tab[1, 2] <- tab[1, 2] + 1
        }
      } else if (info[z * n.sub + sub, 3] == row2) {
        if (info[z * n.sub + sub, item + 3] == 1) {
          tab[2, 1] <- tab[2, 1] + 1
        } else {
          tab[2, 2] <- tab[2, 2] + 1
        }
      }
    }
    return(tab)
  }


## SECONDARY FUNCTIONS: (Should not need to be called)


##
gen <-
  function(n.items) {
    n <- 2^(n.items)
    ret <- matrix(nrow = n.items, ncol = n, data = 1)
    for (i in 1:n.items) {
      ret[i, ] <- c(rep(0, n / 2), rep(1, n / 2))
      n <- n / 2
    }
    return(ret)
  }

##
bin.s <-
  function(val, divs) {
    n <- floor(length(divs) / 2)
    l <- length(divs) - 1
    while (TRUE) {
      if (divs[[n + 1]] <<- val) {
        l <- ceiling((l - 1) / 2)
        n <- ceiling(n + l / 2)
      } else if (divs[[n]] > val) {
        l <- floor((l - 1) / 2)
        n <- ceiling(n - l / 2) - 1
      } else {
        return(n)
      }
    }
  }

##
cov.bs <-
  function(bs.mat1, bs.mat2, n.sets, n.items) {
    val.bar1 <- NULL
    val.bar2 <- NULL
    # Calculate means for each item
    for (i in 1:n.items) {
      val.bar1 <- c(val.bar1, mean(bs.mat1[, i]))
      val.bar2 <- c(val.bar2, mean(bs.mat2[, i]))
    }
    cov.mat <- matrix(nrow <- n.items, ncol <- n.items)
    # Calculate covariances for each pair of items
    for (i in 1:n.items) {
      for (j in 1:n.items) {
        tot <- 0
        for (k in 1:n.sets) {
          curr <- (bs.mat1[k, i] - val.bar1[i]) * (bs.mat2[k, j] - val.bar2[j])
          tot <- tot + curr
        }
        cov.mat[i, j] <- tot / (n.sets - 1)
      }
    }
    return(cov.mat)
  }

contraception <- matrix(nrow = 4, ncol = 32, data = c(
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 14, 0, 0, 0, 0, 0, 1, 1, 3, 0, 1, 0, 0, 2, 1,
  8, 10, 0, 0, 0, 0, 4, 0, 18, 12, 0, 1, 0, 1, 14, 5, 42, 44, 0, 0, 0, 3, 1, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 1, 15,
  0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 5, 7, 0, 0, 0, 0, 2, 0, 6, 3, 0,
  0, 0, 0
))


"child" <-
  structure(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 3, 5, 3, 3, 1, 1,
    0, 3, 5, 0, 2, 0, 0, 1, 0, 0, 0, 4, 0, 2, 0, 1, 2, 1, 0, 0, 0,
    0, 1, 0, 0, 1, 3, 5, 1, 2, 1, 1, 1, 0, 0, 1, 0, 4, 0, 2, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 1, 0, 1, 0,
    0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 1, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    2, 0, 2, 0, 0, 0, 2, 0, 0, 2, 0, 3, 0, 1, 4, 0, 0, 3, 4, 1, 0,
    4, 1, 2, 4, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 2, 1, 1, 0, 0, 1, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
    2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  ), .Dim = c(35, 32))

"language" <-
  structure(c(
    1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
    3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7,
    7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10,
    10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13,
    13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15,
    16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18,
    18, 18, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 21,
    21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23,
    24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26,
    26, 26, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 29, 29,
    29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31,
    32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 34, 34, 34, 34,
    34, 34, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36, 37, 37,
    37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39,
    40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42,
    42, 42, 43, 43, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 45, 45,
    45, 45, 45, 45, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47,
    48, 48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50,
    50, 50, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
    1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3,
    4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
    1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3,
    4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
    1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3,
    4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
    1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3,
    4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
    1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3,
    4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
    1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3,
    4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
    1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3,
    4, 5, 6, 1, 2, 3, 4, 5, 6, 2, 1, 2, 3, 3, 2, 1, 3, 3, 1, 3, 1,
    2, 1, 1, 2, 2, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 3, 2,
    3, 2, 1, 1, 2, 1, 3, 2, 1, 2, 3, 3, 2, 3, 2, 2, 2, 3, 3, 3, 2,
    3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 1, 2,
    2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 2, 2, 2, 2, 3, 2, 2,
    2, 1, 2, 2, 2, 1, 2, 3, 2, 3, 3, 3, 2, 2, 3, 3, 3, 3, 2, 2, 3,
    3, 3, 3, 2, 1, 2, 2, 2, 1, 2, 2, 3, 3, 2, 1, 3, 2, 3, 3, 3, 3,
    2, 2, 2, 3, 2, 1, 1, 1, 1, 2, 2, 1, 2, 3, 3, 3, 3, 3, 1, 2, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 2, 3, 2, 1, 1, 1, 1, 3, 2, 1,
    1, 2, 2, 2, 2, 1, 3, 3, 3, 3, 3, 3, 2, 3, 1, 2, 2, 1, 2, 3, 3,
    3, 3, 2, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3,
    2, 3, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 2, 1, 3, 2, 3, 1, 3, 3, 3,
    2, 3, 3, 2, 1, 1, 2, 2, 1, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 2,
    2, 3, 3, 3, 3, 3, 2, 2, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3,
    3, 3, 1, 1, 2, 3, 3, 2, 1, 3, 2, 3, 3, 3, 2, 1, 1, 1, 1, 1, 0,
    1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1,
    0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1,
    1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1,
    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
    1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1,
    0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
    1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
    1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0,
    1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1,
    1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
    1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1,
    0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1,
    0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0,
    0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1,
    1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1,
    1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1,
    1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1,
    1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1,
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1,
    1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1,
    0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1,
    1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1,
    0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
    0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1,
    1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1,
    0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1,
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1,
    0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0,
    0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0,
    0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0,
    1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0,
    0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1,
    0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,
    1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1,
    0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,
    0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
    0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
    0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  ), .Dim = c(300, 10), .Dimnames = list(NULL, c(
    "Utter", "Rater", "Rating", "A",
    "B", "C", "D", "E", "F", "G"
  )))

# 1 item
farmers00000 <- c(0, 0, 0, 0, 0)
farmers00001 <- c(8, 1, 6, 16, 1)
farmers00010 <- c(17, 3, 6, 30, 1)
farmers00100 <- c(5, 4, 2, 17, 2)
farmers01000 <- c(9, 3, 4, 11, 1)
farmers10000 <- c(9, 0, 0, 11, 1)
# 2 items
farmers11000 <- c(2, 0, 0, 0, 0)
farmers10100 <- c(0, 0, 1, 0, 0)
farmers10010 <- c(0, 0, 0, 0, 1)
farmers10001 <- c(0, 0, 0, 0, 0)
farmers01100 <- c(0, 0, 0, 1, 1)
farmers01010 <- c(1, 0, 0, 1, 0)
farmers01001 <- c(7, 0, 1, 2, 0)
farmers00110 <- c(3, 0, 1, 6, 1)
farmers00101 <- c(1, 0, 0, 1, 1)
farmers00011 <- c(4, 1, 2, 2, 1)
# 3 items
farmers11100 <- c(0, 0, 0, 1, 0)
farmers11010 <- c(0, 0, 0, 0, 0)
farmers11001 <- c(0, 0, 0, 0, 0)
farmers10110 <- c(0, 1, 0, 1, 0)
farmers10101 <- c(0, 0, 0, 0, 0)
farmers10011 <- c(0, 0, 0, 0, 0)
farmers01110 <- c(1, 0, 3, 3, 0)
farmers01101 <- c(0, 0, 0, 0, 1)
farmers01011 <- c(2, 0, 2, 0, 0)
farmers00111 <- c(3, 0, 0, 0, 1)
# 4 items
farmers11110 <- c(1, 1, 0, 2, 0)
farmers11101 <- c(0, 0, 0, 0, 0)
farmers11011 <- c(0, 0, 0, 0, 0)
farmers10111 <- c(0, 0, 0, 0, 0)
farmers01111 <- c(8, 2, 3, 4, 0)
# 5 items
farmers11111 <- c(7, 0, 0, 4, 1)
farmers <- NULL
farmers <- cbind(farmers, farmers00000, farmers00001, farmers00010, farmers00011)
farmers <- cbind(farmers, farmers00100, farmers00101, farmers00110, farmers00111)
farmers <- cbind(farmers, farmers01000, farmers01001, farmers01010, farmers01011)
farmers <- cbind(farmers, farmers01100, farmers01101, farmers01110, farmers01111)
farmers <- cbind(farmers, farmers10000, farmers10001, farmers10010, farmers10011)
farmers <- cbind(farmers, farmers10100, farmers10101, farmers10110, farmers10111)
farmers <- cbind(farmers, farmers11000, farmers11001, farmers11010, farmers11011)
farmers <- cbind(farmers, farmers11100, farmers11101, farmers11110, farmers11111)
colnames(farmers) <- NULL
