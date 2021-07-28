
# main function list:

# getjoint<-function(odds,margp,pair=1,method=1)
# dec_to_bin<-function(dec,c)
# binary_generate<-function(p_joint)
# pairwise<-function(odds,margp)
# perm_generate<-function(N)


# supplementary

# sum_rc<-function(X,method=1)
# prod_rc<-function(X,method=1)
# min_row_col<-function(x,y)
# max_row_col<-function(x,y)
# norm<-function(x)
# indexfct<-function(c)
# get_A
# get_X
# nchooseK<-function(n,K)
# count_joint<-function(X,c,k)
# biggerfct<-function(x,y){
# eqbigger<-function(x,y)
# count_joint_inverse<-function(y,c)
# obs_to_gee<-function(obs,c,k,X)
#############################################

#        sum, max, min

########################################

# sum function
sum_rc <- function(X, method = 1) {
  size <- dim(X)
  if (length(size) == 2) {
    if (method == 1) {
      matrix(1, 1, size[1]) %*% X
    } else {
      X %*% matrix(1, size[2], 1)
    }
  } else {
    tensor(rep(1, size[method]), X, alongA = 1, alongB = method)
  }
} # end function

# prod funtion
prod_rc <- function(X, method = 1) {
  size <- dim(X)
  # if(length(size)==2){
  if (method == 1) {
    h <- matrix(1, 1, size[2])
    for (i in 1:size[1]) {
      h <- h * X[i, ]
    }
  } else {
    h <- matrix(1, size[1], 1)
    for (i in 1:size[2]) {
      h <- h * X[, i]
    }
  }
  # }else{
  # h<-matrix
  # for(i in 1:size[method]){
  # }}
  h
}

# max-function of vectors
max_rc_vec <- function(x, y) {
  if (is.vector(x)) {
    lx <- length(x)
    x <- matrix(x, lx, 1)
  }
  if (is.vector(y)) {
    ly <- length(y)
    y <- matrix(y, ly, 1)
  }
  diff <- y - x
  ones <- matrix(1, lx, 1)
  ones[diff < 0] <- 0
  hilf <- x + (y - x) * ones
  hilf
}

# min-function of vectors
min_rc_vec <- function(x, y) {
  if (is.vector(x)) {
    lx <- length(x)
    x <- matrix(x, lx, 1)
  }
  if (is.vector(y)) {
    ly <- length(y)
    y <- matrix(y, ly, 1)
  }
  diff <- y - x
  ones <- matrix(1, lx, 1)
  ones[diff > 0] <- 0
  hilf <- x + (y - x) * ones
  hilf
}

# max-function Matrix
max_rc_mat <- function(X, method = 1) {
  size <- dim(X)
  hilf <- NULL
  if (method == 1) {
    for (i in 1:size[2]) {
      hilf <- cbind(hilf, max(X[, i]))
    }
  } else {
    for (i in 1:size[1]) {
      hilf <- rbind(hilf, max(X[i, ]))
    }
  }
  hilf
}

# min-function Matrix
min_rc_mat <- function(X, method = 1) {
  size <- dim(X)
  hilf <- NULL
  if (method == 1) {
    for (i in 1:size[2]) {
      hilf <- cbind(hilf, min(X[, i]))
    }
  } else {
    for (i in 1:size[1]) {
      hilf <- rbind(hilf, min(X[i, ]))
    }
  }
  hilf
}

#############################################

#        sum, max, min

########################################






# l2-norm function
norm <- function(x) {
  sqrt(sum(x * x))
}

# function for indexing odds
indexfct <- function(c) {
  LODDS <- c * (c - 1) / 2
  hilf <- NULL
  for (i in 1:(c - 1)) {
    for (j in (i + 1):c) {
      hilf <- cbind(hilf, matrix(c(i, j), 2, 1))
    }
  }
  hilf
} # end indexfct


# function for decimal to binary numbers
# return bin: 0>000,7>111
dec_to_bin <- function(dec, c) {
  c <- max(c, ceiling(log2(dec + 1)))
  bin <- matrix(0, c, 1)

  for (i in 1:c) {
    hilf <- dec - 2^(c - i)
    f <- (hilf >= 0)
    if (f) {
      dec <- hilf
      bin[i] <- 1
    } else {
      dec <- dec
    }
  } # end for
  bin
} # end function

# function for binary to decimal numbers
# return dec: 000>0,111>7
bin_to_dec <- function(bin) {
  c <- length(bin)
  dec <- 0
  for (i in 1:c) {
    dec <- dec + bin[i] * 2^(c - i)
  } # end for
  dec
} # end function


# generating random binary vector from joint distribution
binary_generate <- function(p_joint) {
  lc <- length(p_joint)
  # c<-log2(lc);
  ps <- cumsum(p_joint)
  ps <- c(0, ps)
  r <- runif(1)
  ind <- min(max((1:(lc + 1))[ps <= r]), lc)
  # b<-dec_to_bin(ind-1,c);
  # b
  ind
} # end  function


# counting joint counts
counts <- function(X) {
  size <- dim(X)
  c2 <- 2^size[2]
  count <- matrix(0, 1, c2)
  twos <- 2^(seq(c - 1, 0, by = -1))
  for (i in 1:size[1]) {
    number <- sum(X[i, ] * twos) + 1
    count[number] <- count[number] + 1
  } # end for
  count
} # end  function


###################################
## begin function pairwise
###################################
pairwise <- function(odds, margp) {
  pairpp <- NULL
  c <- length(margp)
  c2 <- 2^c
  LODDS <- c * (c - 1) / 2
  hilfodds <- seq(c - 1, 1, by = -1)
  hilfodds <- c(0, cumsum(hilfodds))
  for (i in 1:LODDS) {
    hilf <- (1:c)[i > hilfodds]
    row <- max(hilf)
    col <- i - hilfodds[row] + row
    # i,row,col
    coeff <- c(odds[i] - 1, -(1 + (odds[i] - 1) * (margp[row] + margp[col])), odds[i] * margp[row] * margp[col])
    coeff <- coeff[3:1]

    hilf <- Re(polyroot(coeff))

    s1 <- min(margp[row], margp[col])
    s2 <- max(0, margp[row] + margp[col] - 1)

    cond <- (hilf <= s1 & hilf >= s2)
    # print(cond);print(hilf);print(s1);print(s2)
    if (cond[1]) {
      pairpp <- rbind(pairpp, hilf[1])
    } else {
      if (cond[2]) {
        pairpp <- rbind(pairpp, hilf[2])
      }
    }
    if (!max(cond)) {
      print("ERROR computation of roots")
    }
  } # for i=1:LODDS

  # print(pairpp)

  pairpp

  ###############################
} # end function pairwise
###############################


###############################
# begin function getjoint
###############################

# computation of joint probabilities from marginal and odds ratios

getjoint <- function(odds, margp, pair = 1, method = 1) {
  # method=1: iteration algorithm default
  # method=2: linear programming

  c <- length(margp)
  c2 <- 2^c
  LODDS <- c * (c - 1) / 2

  margp <- matrix(margp, c, 1)
  odds <- matrix(odds, LODDS, 1)

  indexhilf <- indexfct(c)
  # print(indexhilf)

  # either odds or pairwise probs given
  if (pair) {
    pairpp <- pairwise(odds, margp)
  } else {
    pairpp <- odds
  }
  # print(pairpp)
  B <- matrix(0, c2, c)
  for (i in 1:c2) {
    B[i, ] <- t(dec_to_bin(i - 1, c))
  } # for(i in 1:c2)
  # print(B)
  print(t(pairpp))

  # methods: either iteration method or linear programming
  if (method == 1) {
    p0 <- matrix(0, c2, 1) # 2^c vector
    p <- matrix(1, c2, 1) # 2^c vector
    iterations <- 0
    diff <- 1
    # l.choice<-(2*c+4*LODDS+1);
    # print(LODDS);print(c);print(l.choice)
    # print(norm(p-p0))
    all <- 1
    while (diff > 1e-12) {
      diff <- 0
      # if(iterations%%1e3==0){cat(" ",norm(p-p0),iterations,"\n")}
      p0 <- p
      iterations <- iterations + 1
      for (order in c(1, 2, 3)) {

        # marginal

        if (order == 1) {
          # pairwise
          for (i in 1:LODDS) {
            pairpp11 <- pairpp[i]
            row <- indexhilf[1, i]
            col <- indexhilf[2, i]

            index11 <- (1:c2)[t(B[, row]) == 1 & t(B[, col]) == 1]
            hilf11 <- sum(p[index11])
            p[index11] <- p[index11] * pairpp11 / hilf11
            diff <- max(diff, norm(p - p0))

            if (all) {
              index10 <- (1:c2)[t(B[, row]) == 1 & t(B[, col]) == 0]
              pairpp10 <- margp[row] - pairpp11
              hilf10 <- sum(p[index10])
              p[index10] <- p[index10] * pairpp10 / hilf10
              diff <- max(diff, norm(p - p0))

              index01 <- (1:c2)[t(B[, row]) == 0 & t(B[, col]) == 1]
              pairpp01 <- margp[col] - pairpp11
              hilf01 <- sum(p[index01])
              p[index01] <- p[index01] * pairpp01 / hilf01
              diff <- max(diff, norm(p - p0))

              index00 <- (1:c2)[t(B[, row]) == 0 & t(B[, col]) == 0]
              pairpp00 <- 1 - margp[row] - margp[col] + pairpp11
              hilf00 <- sum(p[index00])
              p[index00] <- p[index00] * pairpp00 / hilf00
              diff <- max(diff, norm(p - p0))
            } # if(all){
          } # for(i in 1:LODDS)
        } # if(order==1)

        if (order == 2) {
          # p
          for (i in 1:c) {
            # if(choice<=c){
            # i<-choice;
            index <- (1:c2)[t(B[, i]) == 1]
            hilf <- sum(p[index])
            # cat("1:marg-i:",i,"index:",index,"\n");
            # print(index);
            # print(hilf)
            p[index] <- p[index] * margp[i] / hilf
            diff <- max(diff, norm(p - p0))
            # }#if(choice<=c){
          } # for(i in 1:c)
          # 1-p
          if (all) {
            for (i in 1:c) {
              index <- (1:c2)[t(B[, i]) == 0]
              hilf <- sum(p[index])
              # cat("0:marg-i:",i,"index:",index,"\n");
              # print(index);
              # print(hilf)
              p[index] <- p[index] * (1 - margp[i]) / hilf
              diff <- max(diff, norm(p - p0))
            } # for(i in 1:c)
          } # if(all){
        } # order==2

        # all joint must sum to 1
        if (order == 3) {
          p <- p * 1 / sum(p)
          diff <- max(diff, norm(p - p0))
        }
      } # for(order in 1:3){
      if (iterations %% 1e2 == 0) {
        cat(" ", diff, iterations, "\n")
      }
    } # while(norm(p-p0))>1e-6)

    print("Number iterations")
    print(iterations)
    p_joint <- p
    # print(p)
  } else {

    # A is created: Ax=b


    A <- rbind(matrix(1, 1, c2), t(B))
    # print(A);print(B)
    for (i in 1:LODDS) {
      # print(indexhilf[,i])
      A <- rbind(A, t(min_rc_vec(B[, indexhilf[1, i]], B[, indexhilf[2, i]])))
    }
    # print("A,odds")
    # print(A);print(LODDS);
    # print("A,odds ende")

    # b is created: Ax=b
    b <- c(1, margp, pairpp)

    # print(matrix(0,c2,1));print(diag(c2));print(matrix(0,c2,1));print(A);print(b)
    p_joint <- simplex(a = matrix(0, c2, 1), A1 = -diag(c2), b1 = matrix(0, c2, 1), A3 = A, b3 = b)
    p_joint <- p_joint$soln
    # print(p_joint)
  }


  # Control

  # marginal probs
  cmarg <- NULL
  for (i in 1:c) {
    cmarg <- rbind(cmarg, sum(p_joint[((1:c2) * t(B[, i])) > 0]))
  }

  cmarg <- cbind(cmarg, margp)
  print(cmarg)

  # joint probs
  cpair <- NULL
  for (i in 1:LODDS) {
    index <- t(min_rc_vec(B[, indexhilf[1, i]], B[, indexhilf[2, i]]))
    cpair <- rbind(cpair, sum(p_joint[((1:c2) * index) > 0]))
  }

  cpair <- cbind(cpair, pairpp)
  print(cpair)

  L <- list(p_joint = matrix(p_joint, length(p_joint), 1), cmarg = cmarg, cpair = cpair)
  return(L)

  ############################################
} # end  function getjoint
############################################



###############################################
# function for generating A
# c blocks of length k, where each block represents 1 item
get_A <- function(c, k) {
  A <- NULL
  for (i in 1:c) {
    h <- rep(0:1, each = 2^(c - i), times = 2^(i - 1))
    # print(h)
    A <- rbind(A, matrix(h, 1, 2^c, byrow = TRUE))
    # print(A)
  }
  A <- kronecker(diag(k), A)
  A
}
##################################

###################################
# get_X Independence model
get_X <- function(c, r, q, model = 1) {
  switch(model,
    # Independence for ML: model=1
    {
      hilf <- kronecker(diag(q), matrix(1, r, 1))
    },
    # Independence for GEE: model=2
    {
      hilf <- diag(r)
      hilf[, 1] <- 1
      hilf <- kronecker(diag(q), hilf)
    }
    # parameters for Wald statistic=2:r,r+2:2*r,...,q*(c-1)*r+2:q*c*r
  ) # end switch

  X <- kronecker(diag(c), t(hilf[1, ]))
  for (i in 2:k) {
    X <- rbind(X, kronecker(diag(c), t(hilf[i, ])))
  } # end for

  X
}
#################################




################################
# generate random permutation ...without replacement
perm_generate <- function(N) {
  # N is vector with r entries - the group sizes
  # X is data matrix: n*c Matrix (observations are rows)
  # X;N
  r <- length(N)
  csN <- c(0, cumsum(N))
  n <- sum(N)

  # perm ist Permutation
  zahl <- 2^n - 1
  twos <- 2^seq(n - 1, 0, by = -1)

  while (zahl == 2^n - 1) {
    perm <- sample(n, n, replace = FALSE, prob = NULL)
    # perm<-c(3,2,1,4,6,5,10,9,7,8)
    v <- NULL
    for (i in 1:r) {
      hilf <- is.element(seq(csN[i] + 1, csN[i + 1]), perm[seq(csN[i] + 1, csN[i + 1])])
      v <- c(v, hilf)
    } # end for

    zahl <- sum(v * twos)
  } # end while

  L <- list(zahl, perm)
  L
} # end perm_generate<-function(N){
##################################


##### function nchoosek
nchooseK <- function(n, K) {
  NP <- choose(n, K[1])
  r <- length(K)
  if (r > 2) {
    sK <- cumsum(K)
    for (i in 2:(r - 1)) {
      NP <- NP * choose(n - sK[i - 1], K[i])
    } # end for;
  } # end if;
  NP
} # end function nchooseK
#################################################

#################################################
# count joint table - function
count_joint <- function(X, c, k) {
  y <- matrix(0, k * 2^c, 1)
  for (i in 1:length(X)) {
    y[X[i]] <- y[X[i]] + 1
  }
  y
}
##################################################

##############################################
# function compare if bigger
bigger <- function(x, y) {
  hilf <- matrix(0, dim(x))
  hilf[x > y] <- 1
  hilf
} # end comparefct<-function(x,y){
####################################################

##############################################
# function compare if bigger
eqbigger <- function(x, y) {
  hilf <- matrix(0, dim(x))
  hilf[x >= y] <- 1
  hilf
} # end comparefct<-function(x,y){
####################################################

#######################
L.fct <- function(m) {
  p <- diag(c(1 / (Z %*% t(Z) %*% m))) %*% m
  hilf <- (A %*% p)
  log(hilf) - log(1 - hilf)
}
##########################

####################################################
# Derivative function of L.fct
derLt.fct <- function(m) {
  Ninv <- diag(c(1 / (Z %*% t(Z) %*% m)))
  hilf <- diag(1 / as.vector(A %*% Ninv %*% m))
  Ninv %*% t((hilf - (1 - hilf)) %*% A)
}
#############################################


######################################################
count_joint_inverse <- function(y) {
  hilf <- NULL
  l <- length(y)
  for (i in 1:l) {
    if (y[i] > 0) {
      # hilf<-rbind(hilf,as.matrix(rep(i,each=y[i])));
      hilf <- c(hilf, rep(i, each = y[i]))
    } # end if
  } # end for(i in 1:l){
  hilf
} # end function count_joint_inverse
#################################

obs_to_gee <- function(obs, c, k, X, N) {
  c2 <- 2^c
  l.obs <- length(obs)
  # hilfvec<-rep(c2*0:(k-1),N);
  Y.gee <- NULL
  X.gee <- NULL
  for (i in 1:l.obs) {
    ind.k <- (obs[i] - 1) %/% c2 + 1 # print(ind.k);
    bin <- (obs[i] - 1) %% c2 + 1 # print(bin);
    Y.gee <- rbind(Y.gee, dec_to_bin(bin - 1, c))
    X.gee <- rbind(X.gee, X[seq((ind.k - 1) * c + 1, ind.k * c), ])
  } # end for(i in 1:l.obs){
  L <- list(X = X.gee, Y = Y.gee)
  L
}

check.pair_marg <- function(pair, marg) {
  c <- length(marg)

  LODDS <- c * (c - 1) / 2
  index <- indexfct(c)
  errors <- 0
  marg.ind <- NULL
  pair.ind <- NULL
  if (c > 2) {
    for (i in 1:max((LODDS - 1), 1)) {
      for (j in (i + 1):max(2, LODDS)) {
        # cat(i,j,"\n")
        ind <- NULL
        if (index[1, i] == index[1, j]) {
          ind <- c(index[1, i], index[2, i], index[2, j])
        } # if
        # not necessary because always covered by the case before
        # if(index[2,i]==index[2,j]){
        # ind<-c(index[1,i],index[1,j],index[2,j]);
        # };#if
        if (!is.null(ind)) {
          p.ind <- NULL
          # l is third index beside i and j
          hilf <- c(0, cumsum((c - 1):1))
          l <- hilf[ind[2]] + (ind[3] - ind[2])
          hilf <- sum(marg[ind]) - sum(pair[c(i, j, l)])
          if (hilf > 1 | hilf < 0) {
            cat("ERROR: Pair+Marg do not aggree:", hilf, " \n")
            errors <- errors + 1
            marg.ind <- rbind(marg.ind, ind)
            pair.ind <- rbind(pair.ind, c(i, j, l))
          }
        } # if
      } # end for
    } # end for
  } # if(c>2){
  cat("There were ", errors, " errors\n")
  list(errors, marg.ind, pair.ind)
} # end check.pair_marg


# get_Z_first
get_Z_first <- function(c) {
  Z <- matrix(0, 2 * c, 2^c)
  for (i in 1:c) {
    hilf01 <- rep(0:1, each = 2^(c - i), times = 2^(i - 1))
    hilf10 <- rep(1:0, each = 2^(c - i), times = 2^(i - 1))
    # print(hilf01);print(hilf10)
    Z[2 * (i - 1) + 1, ] <- hilf01
    Z[2 * (i - 1) + 2, ] <- hilf10
  } #
  Z
} # end

# get_Z_second
get_Z_second <- function(c) {
  h <- matrix(c(0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0), 4, 4, byrow = TRUE)
  if (c > 2) {
    Z1 <- get_Z_second(c - 1)
  } else {
    Z1 <- h
    return(Z1)
  }
  Z1 <- cbind(Z1, Z1)

  c1 <- c - 1
  Z <- matrix(0, c1 * 4, 2^c)
  for (i in 1:c1) {
    for (j in 1:4) {
      hilf <- c(rep(h[j, 1:2], each = 2^(c1 - i), times = 2^(i - 1)), rep(h[j, 3:4], each = 2^(c1 - i), times = 2^(i - 1)))
      # print(hilf);print(2^(c1-i));print(2^(i-1))
      Z[(i - 1) * 4 + j, ] <- hilf
    } # for(j in 1:4){
  } # end for(i in 1:(c-1)){

  # print(Z);print(Z1)
  Z <- rbind(Z, Z1)
  Z
} # end get_Z_second

################## MH function ##################
M.H.Est <- function(data2 = NULL, c, r, K, Z1 = NULL, Z2 = NULL, zeros = FALSE, i.Cab = 0, i.Cba = 0, X = NULL, X1 = NULL, h11 = NULL, h10 = NULL, h01 = NULL, h00 = NULL, n = NULL) {
  r <- 2
  if (!is.null(data2)) {
    X <- Z1[[1]] %*% t(data2) # p
    X1 <- Z1[[2]] %*% t(data2) # 1-p
    h11 <- Z2[[1]] %*% t(data2) # p11 pairwise prob
    h10 <- Z2[[2]] %*% t(data2) # p10 pairwise prob
    h01 <- Z2[[3]] %*% t(data2) # p01 pairwise prob
    h00 <- Z2[[4]] %*% t(data2) # p00 pairwise prob
    n <- matrix(1, 1, 2^c) %*% t(data2)
  }

  # cc2<-c*(c-1)/2;
  # transform Pij to arrays pij
  UT0 <- upper.tri(diag(c))
  # 2 is number of rows per strata
  UT <- array(UT0, c(c, c, K * r)) # make index for array from index from matrix
  X11 <- X10 <- X01 <- X00 <- array(0, c(c, c, K * r)) # create NULL arrays
  # assign values from matrices to uppertri arrays
  X11[UT] <- h11
  X10[UT] <- h10
  X01[UT] <- h01
  X00[UT] <- h00

  # if zeros==TRUE, then add 0.5, because previous computation was not succesful
  if (zeros) {
    # find stratum with lowest observation
    N <- kronecker(diag(K), matrix(1, 1, r)) %*% t(n)
    i1 <- which.max(N)
    i1 <- seq((i1 - 1) * r + 1, i1 * r)
    # find observations with zero
    i2 <- union(i.Cab, i.Cba)
    X[i2, i1] <- X[i2, i1] + 0.5
    X1[i2, i1] <- X1[i2, i1] + 0.5
    # n[i1]<-n[i1]+1;
    # print(i1);print(i2);print(X);print(X1)
  } # if(zeros){
  # X<-X+0.5;X1<-X1+0.5;X11<-X11+0.5;X10<-X10+0.5;X01<-X01+0.5;X00<-X00+0.5;}

  Cov <- tc1 <- tc2 <- tc3 <- tc4 <- matrix(0, c, c) # Cov has 4 numerators
  tv1 <- tv2 <- tv3 <- matrix(0, 1, c) # 4 numerators for var
  Cab <- Cba <- matrix(0, 1, c)

  for (i in 1:K) {
    # actually index for a=1
    a <- (i - 1) * r + 1
    # actually index for b=2
    b <- a + 1
    nk <- n[a] + n[b]
    # cab, cba both are c-dimensional
    # a=1,b=2
    cab <- X[, a] * X1[, b] / nk
    cba <- X[, b] * X1[, a] / nk
    dab <- (X[, a] + X1[, b]) / nk
    dba <- (X[, b] + X1[, a]) / nk
    Cab <- Cab + cab
    Cba <- Cba + cba

    tv1 <- tv1 + cab * dab
    tv2 <- tv2 + cab * dba + cba * dab
    tv3 <- tv3 + cba * dba



    # now covariance terms
    for (x in 1:(c - 1)) {
      for (y in (x + 1):c) {
        tc1[x, y] <- tc1[x, y] + (X[x, a] * X[y, a] * X00[x, y, b] + X11[x, y, a] * X1[x, b] * X1[y, b] - X11[x, y, a] * X00[x, y, b]) / (nk^2)

        tc2[x, y] <- tc2[x, y] + (X1[x, a] * X[y, a] * X10[x, y, b] + X01[x, y, a] * X[x, b] * X1[y, b] - X01[x, y, a] * X10[x, y, b]) / (nk^2)

        tc3[x, y] <- tc3[x, y] + (X[x, a] * X1[y, a] * X01[x, y, b] + X10[x, y, a] * X1[x, b] * X[y, b] - X10[x, y, a] * X01[x, y, b]) / (nk^2)

        tc4[x, y] <- tc4[x, y] + (X1[x, a] * X1[y, a] * X11[x, y, b] + X00[x, y, a] * X[x, b] * X[y, b] - X00[x, y, a] * X11[x, y, b]) / (nk^2)
      } # for(j in (c+1):c){
    } # for(i in 1:c){
  } # end for(i in 1:K){M.H.Est<-function(data2,c,r,K,Nk,Z1,Z2){

  # print(Cab)

  # check if denumerators are zero
  if (min(Cab) > 0 && min(Cba) > 0) {
    Psi <- Cab / Cba
    L <- log(Psi)
    # finally compute variance estmators
    Var <- tv1 / (2 * Cab^2) + tv2 / (2 * Cab * Cba) + tv3 / (2 * Cba^2)
    # and coariance estmators
    for (x in 1:(c - 1)) {
      for (y in (x + 1):c) {
        Cov[x, y] <- tc1[x, y] / (Cab[x] * Cab[y]) - tc2[x, y] / (Cba[x] * Cab[y]) - tc3[x, y] / (Cab[x] * Cba[y]) + tc4[x, y] / (Cba[x] * Cba[y])
      }
    }
    #            1 2  3  4  5 6  7  8  9  10
    return(list(Psi, L, Var, Cov, X, X1, X11, X10, X01, X00))
  } else {
    return(list(NULL, ind.min.Cab = which.min(Cab), ind.min.Cba = which.min(Cba), Cab = Cab, Cba = Cba))
  } # end if
} # end function M.H.Est

# MH estimator only, no covariance function
M.H.Est.L <- function(data2 = NULL, c, r, K, Z1 = NULL, zeros = FALSE, X = NULL, X1 = NULL, n = NULL) {
  r <- 2
  if (!is.null(data2)) {
    X <- Z1[[1]] %*% t(data2) # p
    X1 <- Z1[[2]] %*% t(data2) # 1-p
    n <- matrix(1, 1, 2^c) %*% t(data2)
  }

  if (zeros) {
    X <- X + 0.5
    X1 <- X1 - 0.5
  }
  Cab <- Cba <- matrix(0, 1, c)

  for (i in 1:K) {
    # actually index for a=1
    a <- (i - 1) * r + 1
    # actually index for b=2
    b <- a + 1
    nk <- n[a] + n[b]
    # cab, cba both are c-dimensional
    # a=1,b=2
    cab <- X[, a] * X1[, b] / nk
    cba <- X[, b] * X1[, a] / nk
    Cab <- Cab + cab
    Cba <- Cba + cba
  } # end for(i in 1:K){M.H.Est<-function(data2,c,r,K,Nk,Z1,Z2){

  # check if denumerators are zero
  if (min(Cab) > 0 && min(Cba) > 0) {
    Psi <- Cab / Cba
    L <- log(Psi)

    #            1 2  3  4  5 6  7  8  9  10
    return(list(Psi, L))
  } else {
    return(NULL)
  } # end if
} # end function M.H.Est

##########################################################################

M.H.Est.alt <- function(data2, c, r, K, Z1, Z2, VAR = TRUE, zeros = FALSE, i.Cab = 0, i.Cba = 0, COV = FALSE) {
  c2 <- 2^c
  X <- Z1[[1]] %*% t(data2) # p
  X1 <- Z1[[2]] %*% t(data2) # 1-p
  h11 <- Z2[[1]] %*% t(data2) # p11 pairwise prob
  h10 <- Z2[[2]] %*% t(data2) # p10 pairwise prob
  h01 <- Z2[[3]] %*% t(data2) # p01 pairwise prob
  h00 <- Z2[[4]] %*% t(data2) # p00 pairwise prob
  n <- matrix(1, 1, c2) %*% t(data2)
  # cc2<-c*(c-1)/2;
  # transform Pij to arrays pij
  UT0 <- lower.tri(diag(c))
  # 2 is number of rows per strata
  UT <- array(UT0, c(c, c, K * r)) # make index for array from index from matrix
  X11 <- X10 <- X01 <- X00 <- array(0, c(c, c, K * r)) # create NULL arrays
  # assign values from matrices to uppertri arrays
  X11[UT] <- h11
  X10[UT] <- h10
  X01[UT] <- h01
  X00[UT] <- h00
  X11 <- aperm(X11, c(2, 1, 3))
  X10 <- aperm(X10, c(2, 1, 3))
  X01 <- aperm(X01, c(2, 1, 3))
  X00 <- aperm(X00, c(2, 1, 3))

  # print(X)
  # print(X1)
  # print(X11)

  if (zeros) {
    # find stratum with lowest observation
    N <- kronecker(diag(K), matrix(1, 1, r)) %*% t(n)
    i1 <- which.max(N)
    i1 <- seq((i1 - 1) * r + 1, i1 * r)
    # find observations with zero
    i2 <- union(i.Cab, i.Cba)
    X[i2, i1] <- X[i2, i1] + 0.5
    X1[i2, i1] <- X1[i2, i1] + 0.5
    # n[i1]<-n[i1]+1;
    # print(i1);print(i2);print(X);print(X1)
  } # if(zeros){



  # h11 first (1,2), then (1,3), finally (2,3)
  ###########################################################################
  # r==1
  #############################################################################
  if (r == 1) {
    C <- U <- O <- matrix(0, c, c)
    Psi <- Psit <- matrix(1, c, c)

    for (x in 1:(c - 1)) {
      for (y in (x + 1):c) {
        Ctxy <- Ctyx <- Cxy <- Cyx <- tv1 <- tv2 <- tv3 <- tov1 <- tov2 <- tov3 <- 0

        for (i in 1:K) {
          n0 <- n[i]
          n1 <- n0 - 1
          n2 <- n0 - 2
          n02 <- n0^2
          n12 <- n1^2

          if (n0 > 0) {
            # print(n0)
            # computation of Cab and Cba
            cxy <- X[x, i] * X1[y, i] # not divided by n0
            cyx <- X1[x, i] * X[y, i] # not divided by n0
            Cxy <- Cxy + cxy / n0
            Cyx <- Cyx + cyx / n0 # but here
          }
          if (n1 > 1) {
            ctxy <- (cxy - h10[i]) / n1
            ctyx <- (cyx - h01[i]) / n1
            Ctxy <- Ctxy + ctxy
            Ctyx <- Ctyx + ctyx
          }

          if (VAR) {
            if (n0 > 0) {
              dxy <- (X[x, i] + X1[y, i]) / n02 # here divided by n02 because cxy was not divided by n0
              dyx <- (X1[x, i] + X[y, i]) / n02
              # old variance
              tov1 <- tov1 + cxy * dxy
              tov3 <- tov3 + cyx * dyx
              tov2 <- tov2 + cxy * dyx + cyx * dxy
            }
            if (n2 > 0) {
              tv1 <- tv1 + (h10[i]^2 - h10[i]) / n12
              tv3 <- tv3 + (h01[i]^2 - h01[i]) / n12
              tv2 <- tv2 + ((n2 * n1 + 2 * n1 - 1) * (h10[i] + h01[i]) + 2 * h10[i] * h01[i] - n2 * (h10[i] - h01[i])^2) / n12
            }
          } # if
        } # end for(i in 1:K)

        C[y, x] <- C[x, y] <- Cxy
        # print(C)


        if (min(Cxy) > 0 && min(Cyx) > 0 && min(Ctxy) > 0 && min(Ctyx) > 0) {
          Psi[x, y] <- as.matrix(Cxy / Cyx) # Psi[y,x]<- -Psi[x,y];
          Psit[x, y] <- as.matrix(Ctxy / Ctyx) # Psit[y,x]<- -Psit[x,y];
        } else { # if min
          # cat("not defined\n")
          return(list(NULL, which.min(Cxy), which.min(Cyx)))
        } # end if


        if (VAR) {
          z1 <- Cxy^2
          z2 <- Cxy * Cyx
          z3 <- Cyx^2
          v1 <- Ctxy^2
          v2 <- Ctxy * Ctyx
          v3 <- Ctyx^2


          U[x, y] <- tv1 / v1 + tv2 / v2 + tv3 / v3
          O[x, y] <- tov1 / (2 * z1) + tov2 / (2 * z2) + tov3 / (2 * z3)
        }
      } # for y
    } # for x

    # computation of
    if (COV) {
      O3 <- array(0, dim = c(c, c, c))
      for (x in 1:c) {
        for (y in 1:(c - 1)) {
          for (z in (y + 1):c) {
            if (x != y & x != z & y != z) {
              tc1 <- tc2 <- tc3 <- tc4 <- 0
              for (i in 1:K) {
                n0 <- n[i]
                n02 <- n0^2

                if (n0 > 0) {
                  tc1 <- tc1 + X[x, i] * X1[y, i] * X1[z, i] / n02
                  tc2 <- tc2 + (X[x, i] + X1[x, i]) * X1[y, i] * X[z, i] / n02
                  tc3 <- tc3 + (X[x, i] + X1[x, i]) * X[y, i] * X1[z, i] / n02
                  tc4 <- tc4 + X1[x, i] * X[y, i] * X[z, i] / n02
                }
              }
              z1 <- C[x, y] * C[x, z]
              z2 <- C[y, x] * C[x, z]
              z3 <- C[x, y] * C[z, x]
              z4 <- C[y, x] * C[z, x]

              # cat("x,y,z :",x,y,z,"\n")
              # cat("z1,z2,z3,z4 : ",z1,z2,z3,z4,"\n")
              O3[x, y, z] <- O3[x, z, y] <- tc1 / (3 * z1) + tc2 / (3 * z3) + tc3 / (3 * z2) + tc4 / (3 * z4)
            }
          }
        }
      }
    }



    L <- log(Psi)
    Lt <- log(Psit)
    if (VAR) {
      if (COV) {
        return(list(L = L, Lt = Lt, Psi = Psi, Psit = Psit, X = X, X1 = X1, X11 = X11, X10 = X10, X01 = X01, X00 = X00, U = U, O = O, O3 = O3))
      }

      return(list(L = L, Lt = Lt, Psi = Psi, Psit = Psit, X = X, X1 = X1, X11 = X11, X10 = X10, X01 = X01, X00 = X00, U = U, O = O))
    } else { # if VAR
      return(list(L = L, Lt = Lt, Psi = Psi, Psit = Psit, X = X, X1 = X1, X11 = X11, X10 = X10, X01 = X01, X00 = X00))
    }
  } else {


    ########################################################
    ## r==2
    #########################################################
    # item 1 is in row 1 of X and item 2 is in row 2 of X
    # c will be eitehr 2 or 3
    L <- matrix(0, c, c)

    comb2 <- combn(c, 2)
    n.comb2 <- dim(comb2)[2]

    if (c > 2) {
      comb3 <- combn(c, 3)
      n.comb3 <- dim(comb3)[2]
    } else {
      comb3 <- NULL
      n.comb3 <- 0
    }
    if (c > 3) {
      comb4 <- combn(c, 4)
      n.comb4 <- dim(comb4)[2]
    } else {
      comb4 <- NULL
      n.comb4 <- 0
    }

    C <- Psi <- L <- matrix(0, c, c)
    O <- U <- VA12 <- VA21 <- VB12 <- VB21 <- V <- array(0, dim = c(c, c, c))
    U4bar <- V4 <- U4 <- S <- array(0, dim = c(c, c, c, c))

    for (l in 1:n.comb2) {
      x <- comb2[1, l]
      y <- comb2[2, l]

      ## STEP 1: COMPUTATION of L, U and O

      ## Variance
      Cxy <- Cyx <- tv1 <- tv2 <- tv3 <- tv4 <- 0
      # cat("x,y: ",x,y,"\n")
      for (i in 1:K) {

        # computation of Cab and Cba
        # actually index for a=1
        a <- (i - 1) * r + 1
        # actually index for b=2
        b <- a + 1
        nk <- n[a] + n[b]
        nk <- nk
        nk2 <- nk^2
        # cab, cba both are c-dimensional
        # a=1,b=2
        if (nk > 0) {
          cxy <- X[x, a] * X[y, b] / nk
          cyx <- X[y, a] * X[x, b] / nk
          dxy <- (X[x, a] + X[y, b]) / nk
          dyx <- (X[x, b] + X[y, a]) / nk
          Cxy <- Cxy + cxy
          Cyx <- Cyx + cyx
          if (VAR) {
            tv1 <- tv1 + cxy * dxy
            tv3 <- tv3 + cyx * dyx
            tv2 <- tv2 + cxy * dyx + cyx * dxy
            tv4 <- tv4 - 4 * X[x, a] * X[y, a] * X11[x, y, b] / nk2 - 4 * X[x, b] * X[y, b] * X11[x, y, a] / nk2
            tv4 <- tv4 + 4 * X11[x, y, a] * X11[x, y, b] / nk2
            tv4 <- tv4 - X11[x, y, a] * (X[x, b] + X[y, b]) / nk2 - X11[x, y, b] * (X[x, a] + X[y, a]) / nk2
          } # end if nk>0
        } # end if
      } # end for(i in 1:K)

      # print(Cxy)
      # print(Cyx)

      if (Cxy > 0 && Cyx > 0) {
        Psi[x, y] <- Cxy / Cyx
        C[x, y] <- Cxy
        C[y, x] <- Cyx
        L[x, y] <- log(Psi[x, y])
        L[y, x] <- -L[x, y]
        # finally compute variance estmators

        if (VAR) {
          z1 <- Cxy^2
          z2 <- Cxy * Cyx
          z3 <- Cyx^2
          est <- tv1 / (2 * z1) + tv2 / (2 * z2) + tv3 / (2 * z3)
          O[y, x, x] <- O[x, y, y] <- est
          U[y, x, x] <- U[x, y, y] <- est + tv4 / (2 * z2)

          # print(U)
          # print(O)
        }
      } else {
        O[y, x, x] <- O[x, y, y] <- U[y, x, x] <- U[x, y, y] <- NA
        return(list(NULL, which.min(Cxy), which.min(Cyx)))
      } # if(min(Cxy)>0 && min(Cyx)>0)
    } # for(l in 1:n.comb){

    ## computation covariances
    if (VAR) {
      ## STEP2:  COMPUTATION of VA12, VA21, VB12, VB21, V_{xy,xz} and U[x,y,z] and O[x,y,z]

      Xh <- X11 + aperm(X11, c(2, 1, 3))

      if (c > 2) {
        set3 <- 1:c

        for (x in set3) {
          combxyz <- combn(set3[-x], 2)

          n.combxyz <- dim(combxyz)[2]

          for (l in 1:n.combxyz) {
            y <- combxyz[1, l]
            z <- combxyz[2, l]
            # cat("x,y,z: ",x,y,z,"\n")

            tc1 <- tc2 <- tc3 <- tc4 <- 0

            for (i in 1:K) {
              a <- (i - 1) * r + 1
              b <- a + 1

              nk <- n[a] + n[b]
              if (nk > 0) {
                nk2 <- nk^2
                # cat("a,b : ",a,b,"\n")

                # V_xyz,12^A, V_xyz,21^A
                VA12[x, y, z] <- VA12[x, y, z] + ((X[x, a]^2) * Xh[y, z, b]) / nk2 # OK
                VA21[x, y, z] <- VA21[x, y, z] + ((X[x, b]^2) * Xh[y, z, a]) / nk2 # OK
                # V_xyz,12^B, V_xyz,21^B
                VB12[x, y, z] <- VB12[x, y, z] - X[x, a] * Xh[y, z, b] / nk2 # OK
                VB21[x, y, z] <- VB21[x, y, z] - X[x, b] * Xh[y, z, a] / nk2 # OK
                # V_xy,xz and V_xz,xy
                V[x, y, z] <- V[x, y, z] + (X[x, a] * X[y, a] * Xh[x, z, b] + Xh[x, y, a] * X[x, b] * X[z, b] - Xh[x, y, a] * Xh[x, z, b]) / nk2 # OK
                V[x, z, y] <- V[x, z, y] + (X[x, a] * X[z, a] * Xh[x, y, b] + Xh[x, z, a] * X[x, b] * X[y, b] - Xh[x, z, a] * Xh[x, y, b]) / nk2 # OK

                tc1 <- tc1 + X[x, a] * X[y, b] * X[z, b] / nk2 # OK
                tc2 <- tc2 + (X[x, a] + X[x, b]) * X[y, b] * X[z, a] / nk2 # OK
                tc3 <- tc3 + (X[x, a] + X[x, b]) * X[y, a] * X[z, b] / nk2 # OK
                tc4 <- tc4 + X[x, b] * X[y, a] * X[z, a] / nk2 # OK
              } # end if nk>0
            } # end for i in 1:K

            # VA and VB are symmetric in y and z, and we have y != z
            VA12[x, z, y] <- VA12[x, y, z]
            VA21[x, z, y] <- VA21[x, y, z] # OK
            VB12[x, z, y] <- VB12[x, y, z]
            VB21[x, z, y] <- VB21[x, y, z] # OK

            z1 <- C[x, y] * C[x, z] # OK
            z2 <- C[y, x] * C[x, z] # OK
            z3 <- C[x, y] * C[z, x] # OK
            z4 <- C[y, x] * C[z, x] # OK

            O[x, y, z] <- O[x, z, y] <- tc1 / (3 * z1) + tc2 / (3 * z3) + tc3 / (3 * z2) + tc4 / (3 * z4) # OK
          } # for(l in 1:n.combxyz){
        } # end for x in 1:c
      } # if (c>2)

      if (c > 3) {
        ## STEP3:  COMPUTATION of V4 = v_{xy,wz}
        for (l in 1:n.comb4) {
          # ind4 consists of 4 indices
          ind4 <- comb4[, l]

          for (k in 1:n.comb2) {
            # now we take all combinations over these 4 indices to compute V4, hat{v}_{xw,yz}
            # here were go over all combinations of 2 over the fist 2 indices and the
            # last 2 are the remaining of the 4 indices
            x <- ind4[comb2[1, k]]
            y <- ind4[comb2[2, k]]
            hilf <- setdiff(ind4, comb2[, k])
            w <- hilf[1]
            z <- hilf[2]

            # cat("x,y,w,z:",x,y,w,z, " \n ")
            for (i in 1:K) {
              a <- (i - 1) * r + 1
              b <- a + 1
              nk <- n[a] + n[b]
              if (nk > 0) {
                nk2 <- nk^2
                # cat(nk2,"\n")

                # term V4 is v_{xw,yz}
                V4[x, y, w, z] <- V4[x, y, w, z] + (X[x, a] * X[y, a] * Xh[w, z, b] + X[w, b] * X[z, b] * Xh[x, y, a] - Xh[w, z, b] * Xh[x, y, a]) / nk2 # OK
              } # end if nk>0
            } # end i in 1:K)
            # symetric: (x and w)  and   (y and z)
            V4[y, x, w, z] <- V4[x, y, z, w] <- V4[y, x, z, w] <- V4[x, y, w, z] # OK
          } # comb2
        } # for(l in 1:n.comb){
      } # if(c>3){

      # Covariance U_xyz symmetric in y and z

      if (c > 2) {

        ## STEP 4: COMPUTATION of  U[x,y,z]

        for (x in set3) {
          combxyz <- combn(set3[-x], 2)
          n.combxyz <- dim(combxyz)[2]

          for (l in 1:n.combxyz) {
            y <- combxyz[1, l]
            z <- combxyz[2, l]
            # cat("x,y,z:",x,y,z, " \n ")
            # cat("x,y,z :",x,y,z,"\n")
            z1 <- C[x, y] * C[x, z] # OK
            z2 <- C[y, x] * C[x, z] # OK
            z3 <- C[x, y] * C[z, x] # OK
            z4 <- C[y, x] * C[z, x] # OK

            # 2nd line of formulae in paper
            d2 <- VB12[x, y, z] / (3 * z1) + VB21[x, y, z] / (3 * z4) # OK
            d2 <- d2 + (VB12[y, x, z] + VB21[z, x, y]) / (3 * z2) + (VB21[y, x, z] + VB12[z, x, y]) / (3 * z3) # OK

            # 1st line of formulae in paper
            # middle terms in first line
            d3 <- V[x, y, z] / z2 + V[x, z, y] / z3 # OK
            # 1st and 4th term in first line
            d4 <- VA12[x, y, z] / z1 + VA21[x, y, z] / z4 # OK

            # h1<-d4-d3+d1;
            h2 <- d4 - d3 + d2 # OK

            U[x, z, y] <- U[x, y, z] <- h2 + O[x, y, z] # OK
          } # for(x in set3){
        } # for(l in 1:n.combxyz){
      } # if(c>2){



      if (c > 3) {
        ## STEP 5 : COMPUTATION of U4
        for (l in 1:n.comb4) {
          x <- comb4[1, l]
          y <- comb4[2, l]
          w <- comb4[3, l]
          z <- comb4[4, l]
          # cat("x,y,w,z:",x,y,w,z, " \n ")
          z1 <- C[x, y] * C[w, z]
          z2 <- C[y, x] * C[w, z]
          z3 <- C[x, y] * C[z, w]
          z4 <- C[y, x] * C[z, w]

          # changing pairwise indices; changing just once changinges sign
          # changing twice doesn't change sign

          U4[x, y, w, z] <- U4[y, x, z, w] <- V4[x, w, y, z] / z1 - V4[y, w, x, z] / z2 - V4[x, z, y, w] / z3 +
            V4[y, z, x, w] / z4 # OK

          U4[x, y, z, w] <- U4[y, x, w, z] <- -U4[x, y, w, z] # OK

          # exchaning last 2 with first 2 indices
          U4[w, z, x, y] <- U4[z, w, y, x] <- U4[x, y, w, z] # OK

          U4[z, w, x, y] <- U4[w, z, y, x] <- -U4[x, y, w, z] # OK
        } # for(l in 1:n.comb4){
      } # if(c>3){
    } # end if VAR


    if (c > 2) { # otherwise no generalised estimator

      # now compute generalized estimators
      Lplus <- L %*% matrix(1, c, 1)

      ## STEP 6 : COMPUTATION of U+ and O+
      # compute first Uplus
      Lbar <- Oplus <- Uplus <- matrix(0, c, c)
      for (i in 1:c) {
        for (j in i:c) {
          if (i == j) {
            Uplus[i, i] <- sum(U[i, , ])
            Oplus[i, i] <- sum(O[i, , ]) # OK
          } else {
            Lbar[i, j] <- (Lplus[i] - Lplus[j]) / c
            Lbar[j, i] <- -Lbar[i, j] # OK
            if (VAR) {
              Uplus[j, i] <- Uplus[i, j] <- sum(U[, i, j]) - sum(U[i, j, ]) - sum(U[j, i, ]) + U[i, j, j] + sum(U4[j, , i, ]) # OK
              Oplus[j, i] <- Oplus[i, j] <- sum(O[, i, j]) - sum(O[i, j, ]) - sum(O[j, i, ]) + O[i, j, j] # OK
              # print(sum(U4[j,,i,]))
            } # if VAR
          } # if
        } # end j
      } # end i

      # x=1 is fixed

      Obar <- Ubar <- array(0, dim = c(c, c, c))
      if (VAR) {
        ## STEP 7: COMPUTATION of Ubar and Obar
        for (x in set3) {
          for (y in set3[-x]) {
            for (z in set3[-x]) {
              # cat("x,y,z",x,y,z,"\n")
              Ubar[x, z, y] <- Ubar[x, y, z] <- (Uplus[x, x] - Uplus[x, z] - Uplus[y, x] + Uplus[y, z]) / c^2
              Obar[x, z, y] <- Obar[x, y, z] <- (Oplus[x, x] - Oplus[x, z] - Oplus[y, x] + Oplus[y, z]) / c^2
            }
          }
        } # end for
      } else {
        Ubar <- U
        Obar <- O
        Lbar <- L
      }

      if (c > 3) {
        for (l in 1:n.comb4) {
          # ind4 consists of 4 indices
          ind4 <- comb4[, l]

          for (k in 1:n.comb2) {
            # now we take all combinations over these 4 indices to compute U4bar
            # here were go over all combinations of 2 over the fist 2 indices and the
            # last 2 are the remaining of the 4 indices
            x <- ind4[comb2[1, k]]
            y <- ind4[comb2[2, k]]
            hilf <- setdiff(ind4, comb2[, k])
            w <- hilf[1]
            z <- hilf[2]
            # cat("x,y,w,z",x,y,w,z,"\n")
            # changing pairwise indices; changing just once changinges sign
            # changing twice doesn't change sign
            U4bar[x, y, w, z] <- U4bar[y, x, z, w] <- (Uplus[x, w] - Uplus[x, z] - Uplus[y, w] + Uplus[y, z]) / c^2
            U4bar[x, y, z, w] <- U4bar[y, x, w, z] <- -U4bar[x, y, w, z]

            # exchaning last 2 with first 2 indices
            # U4bar[w,z,x,y] <- U4bar[z,w,y,x] <-   U4bar[x,y,w,z]#OK
            # U4bar[z,w,x,y] <- U4bar[w,z,y,x] <- - U4bar[x,y,w,z];#OK
          }
        }
      } # if(c>3){

      return(list(L = L, Lbar = Lbar, Psi = Psi, O = O, U = U, Ubar = Ubar, Obar = Obar, U4 = U4, U4bar = U4bar))
    } # end if(c>2)


    return(list(L = L, Psi = Psi, U = U, O = O))
  } # end if(r==1)

  # finally compute variance estmators
} # end function M.H.Est.alt

#########################################################################




U2Ubar <- function(L, U, U4, c) {
  # now compute generalized estimators
  Lplus <- L %*% matrix(1, c, 1)

  U4bar <- array(0, dim = c(c, c, c, c))
  # compute first Uplus
  Lbar <- Oplus <- Uplus <- matrix(0, c, c)
  for (i in 1:c) {
    for (j in i:c) {
      if (i == j) {
        Uplus[i, i] <- sum(U[i, , ]) # Oplus[i,i]<-sum(O[i,,]);
      } else {
        Lbar[i, j] <- (Lplus[i] - Lplus[j]) / c
        Lbar[j, i] <- Lbar[i, j]
        # if(VAR){
        Uplus[j, i] <- Uplus[i, j] <- sum(U[, i, j]) - sum(U[i, j, ]) - sum(U[j, i, ]) + U[i, j, j] + sum(U4[j, , i, ])
        # Oplus[j,i]<-Oplus[i,j]<-sum(O[,i,j])-sum(O[i,j,])-sum(O[j,i,])+O[i,j,j];
        # }#if VAR
      } # if
    } # end j
  } # end i

  # x=1 is fixed
  Ubar <- array(0, dim = c(c, c, c))

  set3 <- 1:c
  for (x in set3) {
    for (y in set3[-x]) {
      for (z in set3[-x]) {
        # cat("x,y,z",x,y,z,"\n")
        Ubar[x, z, y] <- Ubar[x, y, z] <- (Uplus[x, x] - Uplus[x, z] - Uplus[y, x] + Uplus[y, z]) / c^2
        # Obar[x,z,y]<-Obar[x,y,z]<-(Oplus[x,x]-Oplus[x,z]-Oplus[y,x]+Oplus[y,z])/c^2
      }
    }
  }

  if (c > 3) {
    comb2 <- combn(c, 2)
    n.comb2 <- dim(comb2)[2]
    comb4 <- combn(c, 4)
    n.comb4 <- dim(comb4)[2]

    for (l in 1:n.comb4) {
      # ind4 consists of 4 indices
      ind4 <- comb4[, l]

      for (k in 1:n.comb2) {
        # now we take all combinations over these 4 indices to compute U4bar
        # here were go over all combinations of 2 over the fist 2 indices and the
        # last 2 are the remaining of the 4 indices
        x <- ind4[comb2[1, k]]
        y <- ind4[comb2[2, k]]
        hilf <- setdiff(ind4, comb2[, k])
        w <- hilf[1]
        z <- hilf[2]
        # cat("x,y,w,z",x,y,w,z,"\n")
        # changing pairwise indices; changing just once changinges sign
        # changing twice doesn't change sign
        U4bar[x, y, w, z] <- U4bar[y, x, z, w] <- (Uplus[x, w] - Uplus[x, z] - Uplus[y, w] + Uplus[y, z]) / c^2
        U4bar[x, y, z, w] <- U4bar[y, x, w, z] <- -U4bar[x, y, w, z]

        # exchaning last 2 with first 2 indices
        # U4bar[w,z,x,y] <- U4bar[z,w,y,x] <-   U4bar[x,y,w,z]#OK
        # U4bar[z,w,x,y] <- U4bar[w,z,y,x] <- - U4bar[x,y,w,z];#OK
      }
    }
  } # if(c>3){

  return(list(Lbar = Lbar, Ubar = Ubar, U4bar = U4bar))
}



L.vector <- function(L) {
  c <- dim(L)[1]
  if (is.null(c)) {
    return(L)
  }
  if (dim(L)[2] == 1) {
    return(c(t(L)))
  }
  ind <- indexfct(c)
  l.ind <- dim(ind)[2]
  L.vec <- matrix(0, 1, l.ind)
  L.names <- matrix(0, 1, l.ind)
  for (i in 1:l.ind) {
    L.vec[i] <- L[ind[1, i], ind[2, i]]
    L.names[i] <- paste(ind[1, i], ind[2, i])
  }
  colnames(L.vec) <- L.names
  L.vec
}




transf.data <- function(data1, c, r, K, Nk) {
  # data1 is (K*Nk)*2 matrix
  data2 <- matrix(0, K * r, 2^c)
  # data is (Nk*K)*(c+3) matrix, 1st column observation, 2nd strata, 3rd rating (group)
  for (i in 1:Nk) {
    for (k in 1:K) {
      # compute index in data1 from current obs in data
      hilf <- data1[(i - 1) * K + k, ]
      rating <- hilf[1]
      dec <- hilf[2]
      index <- r * (k - 1) + rating
      data2[index, dec] <- data2[index, dec] + 1
    } # end K
  } # end for i
  data2
} # end transform

#################################

transformZ <- function(Z, k) {
  l.Z <- dim(Z)[1]
  newZ <- NULL
  d <- l.Z / k
  for (i in 1:k) {
    index <- seq(i, l.Z, by = k)
    newZ[[i]] <- Z[index, ]
  } # for(i in 1:k){
  newZ
} # end function transform


create.Z.Y <- function(data0, c, r, K, Nk) {
  Y <- NULL
  if (r > 1) {
    NkK <- Nk * K
    cNkK <- NkK * c
    Z <- cbind(matrix(0, cNkK, c), kronecker(matrix(1, Nk, 1), diag(K * c)))

    for (i in 1:NkK) {
      Y <- c(Y, data0[i, (3 + 1):(3 + c)])
      index <- seq((i - 1) * c + 1, i * c)
      if (data0[i, 3] == 1) {
        Z[index, 1:c] <- diag(c)
      }
    } # end for i
  }
  if (r == 1) {
    if (K > 1) {
      Z <- kronecker(matrix(1, Nk, 1), cbind(kronecker(matrix(1, K, 1), diag(c)), kronecker(diag(K), matrix(1, c, 1))[, 2:K]))
    } else { # if(K>1){

      Z <- kronecker(matrix(1, Nk, 1), cbind(kronecker(matrix(1, K, 1), diag(c))))
    }

    for (i in 1:NkK) {
      Y <- c(Y, data0[i, (3 + 1):(3 + c)])
    } # end for i

    return(list(as.matrix(Y), Z))
  } # if(r==1){
} # end create.Z.Y


# h<-create.Z.Y_from_data2(ff,c,5)
create.Y_from_data2 <- function(data2, c, K) {
  c2 <- 2^c
  # data2 has dimension (K*r)*2^c;
  # where r rows stand for one stratum

  # Y will have dimension sum(data2)*c
  Y <- NULL
  for (k in 1:K) {
    for (j in 1:c2) {
      no.obs <- data2[k, j]
      if (no.obs > 0) {
        bin <- t(dec_to_bin(j - 1, c))
        # cat("k,j,no.obs : ", k,j,no.obs, " \n");print(hilf);print(bin);
        for (i in 1:no.obs) {
          Y <- c(Y, bin)
        } # end for i
      } # end if(no.obs]>1)
    } # end j
  } # end k
  return(list(Y = as.matrix(Y)))
} # end function create.Z.Y_from_data2


show.elapsed.time <- function(start.time, string = "") {
  end.time <- proc.time()[3]
  time <- end.time - start.time
  hours <- floor(time / 3600)
  minutes <- floor(time / 60) - hours * 60
  seconds <- time - hours * 3600 - minutes * 60
  cat(string, "Time: ", hours, "h", minutes, "min", seconds, "sec", "\n")
  return(time)
} # End FUNCTION: show.elapsed.time


########################################
get.Cov <- function(V, V4, c) {
  # V is c*c*c array
  # V4 is c*c*c*c array
  ind <- indexfct(c)
  n <- c * (c - 1) / 2
  Cov <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      i1 <- ind[, i][1]
      j1 <- ind[, j][1]
      i2 <- ind[, i][2]
      j2 <- ind[, j][2]

      # V2_{i1,i2} or equivalently i==j
      if (i1 == j1 & i2 == j2) {
        Cov[i, j] <- V[i1, i2, i2]
      }
      # V3_{i1,i2,j2}
      if (i1 == j1 & i2 != j2) {
        Cov[i, j] <- V[i1, i2, j2]
      }
      # U_{i1,i2,j2}
      if (i1 != j1 & i2 == j2) {
        Cov[i, j] <- V[i2, i1, j1]
      }
      # if i1!=j1 and i2!=j2
      if (i1 != j1 & i2 != j2) {
        # i1==j2 which means i2!=j1
        if (i1 == j2) {
          Cov[i, j] <- -V[i1, j1, i2]
        }
        # i2==j1 which means i1!=j2
        if (i2 == j1) {
          Cov[i, j] <- -V[j1, j2, i1]
        }
        # now i1\neq j2 and i2 \neq j1
        if (i2 != j1 & i1 != j2) {
          Cov[i, j] <- V4[i1, i2, j1, j2]
        }
      }
    } # end i
  } # end j
  return(Cov = Cov)
} # end function
###############################################


BT2F <- function(Cov.BT, c, Ut = NULL, BT = TRUE) {
  ind <- indexfct(c)
  l.ind <- dim(ind)[2]

  if (BT) {
    # Cov.BT is n*n matrix
    U <- array(0, dim = c(c, c, c))
    U4 <- array(0, dim = c(c, c, c, c))
    for (i in 1:l.ind) {
      for (j in i:l.ind) {
        i1 <- ind[1, i]
        j1 <- ind[1, j]
        i2 <- ind[2, i]
        j2 <- ind[2, j]

        # note: i1<i2 and j1<j2

        # U_{i1,i2,i2}
        if (i1 == j1 & i2 == j2) {
          if (is.null(Ut)) {
            # use Bootstrap variance
            U[i2, i1, i1] <- U[i1, i2, i2] <- Cov.BT[i, j]
          } else {
            # use formula variance
            U[i2, i1, i1] <- U[i1, i2, i2] <- Ut[i1, i2]
          }
        }
        # U_{i1,i2,j2}
        if (i1 == j1 & i2 != j2) {
          U[i1, j2, i2] <- U[i1, i2, j2] <- Cov.BT[i, j]
        }
        # U_{i1,i2,j2}
        if (i1 != j1 & i2 == j2) {
          U[i2, i1, j1] <- U[i2, j1, i1] <- Cov.BT[i, j]
        }
        # if i1!=j1 and i2!=j2
        if (i1 != j1 & i2 != j2) {
          # i1==j2 which means i2!=j1
          if (i1 == j2) {
            U[i1, j1, i2] <- U[i1, i2, j1] <- -Cov.BT[i, j]
          }
          # i2==j1 which means i1!=j2
          if (i2 == j1) {
            U[i2, j2, i1] <- U[i2, i1, j2] <- -Cov.BT[i, j]
          }
          # now i1\neq j2 and i2 \neq j1
          if (i2 != j1 & i1 != j2) {
            U4[j2, j1, i2, i1] <- U4[i2, i1, j2, j1] <- U4[i1, i2, j1, j2] <- U4[j1, j2, i1, i2] <- Cov.BT[i, j]
            U4[j1, j2, i2, i1] <- U4[i1, i2, j2, j1] <- U4[i2, i1, j1, j2] <- U4[j2, j1, i1, i2] <- -Cov.BT[i, j]
          }
        }
      }
    }



    return(list(U = U, U4 = U4))
  } else {
    U <- Cov.BT
    for (x in 1:(c - 1)) {
      for (y in (x + 1):c) {
        cat("x,y: ", x, y, "\n")
        U[x, y, y] <- U[y, x, x] <- Ut[x, y]
      }
    }
    return(U)
  }
} # end function


get.gen <- function(Cov, c, L = NULL, LT = TRUE) {
  ind <- indexfct(c)
  l.ind <- dim(ind)[2]
  # print(LT)

  if (is.null(L)) {
    L <- rep(0, 1, l.ind)
  }

  if (!LT) {
    hilf1 <- BT2F(Cov, c)
    U <- hilf1$U
    U4 <- hilf1$U4

    check.U(U)
    check.U4(U4)

    Lnew <- matrix(0, c, c)
    Lnew[lower.tri(Lnew)] <- L
    Lnew <- t(Lnew)
    L1 <- Lnew + t(-Lnew)
    hilf2 <- U2Ubar(L1, U, U4, c)

    Ubar <- hilf2$Ubar
    U4bar <- hilf2$U4bar
    Lbar <- hilf2$Lbar

    check.U(Ubar)
    check.U4(U4bar)

    CovL <- get.Cov(U, U4, c)
    CovLbar <- get.Cov(Ubar, U4bar, c)
    Lbar <- L.vector(Lbar)$L.vec
    L <- L.vector(Lnew)$L.vec

    list(L = L, Lbar = Lbar, CovL = CovL, CovLbar = CovLbar)
  } else {
    # compute generalized estimator as a linear transformation


    # transform L to a vector
    L <- L.vector(L)
    # print(L)
    # make sure Cov is not a triangular matrix
    C <- t(Cov)
    P <- hilf <- matrix(0, l.ind, l.ind)
    hilf[lower.tri(hilf)] <- C[lower.tri(C)]
    CovL <- hilf + t(hilf) + diag(diag(Cov))
    # now Cov is symmetric matrix

    for (i in 1:l.ind) {
      i1 <- ind[1, i]
      i2 <- ind[2, i]
      plus1 <- (1:l.ind)[i1 == ind[1, ]]
      plus2 <- (1:l.ind)[i2 == ind[2, ]]

      minus1 <- (1:l.ind)[i1 == ind[2, ]]
      minus2 <- (1:l.ind)[i2 == ind[1, ]]

      P[i, plus1] <- P[i, plus1] + 1
      P[i, plus2] <- P[i, plus2] + 1
      P[i, minus1] <- P[i, minus1] - 1
      P[i, minus2] <- P[i, minus2] - 1
    } # end for(i in 1:l.ind) {
    P <- P / c

    Lbar <- c(P %*% as.matrix(c(L)))
    CovLbar <- P %*% CovL %*% t(P)
    list(L = L, Lbar = Lbar, CovL = CovL, CovLbar = CovLbar)
  }
}


check.U <- function(U) {
  c <- dim(U)[1] # is of dim c*c*c
  set3 <- 1:c
  for (x in set3) {
    set2 <- combn(set3[-x], 2)
    y <- set2[1]
    z <- set2[2]
    cat("x,y,z: ", x, y, z, "U[x,y,z]: ", U[x, y, z], "U[x,z,y]: ", U[x, z, y], "\n")
    cat("abs(U[x,y,z]-U[x,z,y]):", abs(U[x, y, z] - U[x, z, y]), "\n\n")
  }
}


check.U4 <- function(U4) {
  c <- dim(U4)[1] # is of dim c*c*c*c

  if (c > 3) {
    comb4 <- combn(c, 4)
    n.comb4 <- dim(comb4)[2]



    for (l in 1:n.comb4) {
      # ind4 consists of 4 indices
      ind4 <- comb4[, l]


      x <- ind4[1] # take first index
      for (y in ind4[2:4]) { # take one of the 3 other indices
        hilf <- setdiff(ind4, c(x, y))
        w <- hilf[1] # y and z are the remaining indices
        z <- hilf[2]
        cat("x,y,w,z", x, y, w, z, "\n")
        set.of.equal.U <- NULL
        # changing pairwise indices; changing just once changinges sign
        # changing twice doesn't change sign

        # exchanging first two and last indices
        set.of.equal.U <- c(set.of.equal.U, U4[x, y, w, z], U4[y, x, z, w], -U4[x, y, z, w], -U4[y, x, w, z])

        # exchaning last 2 with first 2 indices
        set.of.equal.U <- c(set.of.equal.U, U4[w, z, x, y], U4[z, w, y, x], -U4[z, w, x, y], -U4[w, z, y, x])
      } # endfor(y in ind4[2:4]){
      cat("x,y,w,z: ", x, y, w, z, "U[x,y,w,z]", set.of.equal.U, "\n")
      cat("abs(max-min):", abs(max(set.of.equal.U) - min(set.of.equal.U)), "\n\n")
    } # for(l in 1:n.comb4){
  } # if(c>3){
}


remove.item <- function(data, i) {
  d <- dim(data)
  c2 <- d[2]
  K <- d[1]
  c <- log2(c2)
  A <- kronecker(diag(2^(i - 1)), kronecker(matrix(1, 1, 2), diag(2^(c - i))))
  data %*% t(A)
}




output.variances <- function(l, ind, var = F) {
  l.list <- length(l)
  d <- dim(l[[1]])[2] # is d times d matrix
  i <- ind[1]
  j <- ind[2]
  a <- NULL
  cat("i,j: ", i, j, "\n")
  for (k in 1:l.list) {
    if (!var) {
      a <- c(a, paste("$", format(l[[k]][i, j], digits = 3, nsmall = 3), "$ & \n", sep = ""))
    } else {
      a <- c(a, paste("($", format(l[[k]][i, j], digits = 3, nsmall = 3), "$)  \n", sep = ""))
    }
  }
  cat(a)
}

fct.table <- function(a) {
  output <- NULL
  for (i in 1:length(a)) {
    if (i == 4) {
      output <- paste(output, a[i], " & & ", sep = "")
    } else {
      output <- paste(output, a[i], " & ", sep = "")
    }
  }
  cat(output)
}

# transf.data.alt<-function(data1,c,K,Nk){
# c2<-2^c;
# data2<-matrix(0,K,c2)
# Nk<-rep(Nk,K-length(Nk)+1)
# for(k in 1:K){
# data2[k,]<-tabulate(data1[seq(k,NkK,by=K),2],nbins=c2)
# }
# data2
# }#end function  trans.data.alt



# Function for bootstrapping, that gives L, Lt, U and Ut
stat.both <- function(datanew, index, c, r, K, Nknew, Z1, Z2, initial = NULL) {
  # index is from 1 to sum Nk=K*Nk=500 in most cases
  # data1 is vector of dimension Nk*K=500
  # data1.bt<-data1[index,];
  # data2.bt<-transf.data.alt(data1.bt,c,K,Nk);
  # print(index)

  c2 <- 2^c
  data.bt <- datanew[index]

  data2.bt <- matrix(0, r * K, c2)

  indk <- c(0, cumsum(Nknew))
  for (k in 1:(r * K)) {
    if (Nknew[k] > 0) {
      chosenind <- data.bt[(indk[k] + 1):indk[k + 1]]
      data2.bt[k, ] <- tabulate(chosenind, nbins = c2)
    } # end
  } # end for k

  # print(sum(data2.bt))
  # print(data2.bt)

  # usually zeros=FALSE
  # M.H.Est.alt(data2,c,r,K,Z1,Z2);
  hilf <- M.H.Est.alt(data2.bt, c, r, K, Z1, Z2) # M.H.Est.alt(data2,c,r,K,Z1,Z2);
  # print(hilf$L)
  # print(index)
  # print(hilf)
  if (is.null(hilf$L)) {
    # cat("NULL \n")
    initial
  } else {
    if (r == 1) {
      # cat("not NULL \n")
      L0 <- hilf$L[1, 2]
      Lt0 <- hilf$Lt[1, 2]
      Ut0 <- hilf$U[1, 2]
      U0 <- hilf$O[1, 2]
      c(L0, U0, Lt0, Ut0)
    } else {
      # cat("hello1\n")
      L0 <- c(hilf$L[1, 2])
      # Lbar0<-c(hilf$Lbar[1,2]);
      # cat("hello2\n")
      U0 <- c(hilf$U[1, 2, 2])
      # Ubar0<-c(hilf$Ubar[1,2,2],hilf$Ubar[1,3,3]);
      O0 <- c(hilf$O[1, 2, 2])
      # cat("hello3\n")
      # Obar0<-c(hilf$Obar[1,2,2],hilf$Obar[1,3,3]);
      # hilf<-c(L0,Lbar0,U0,Ubar0,O0,Obar0)
      hilf <- c(L0, U0, O0)
      names(hilf) <- c("L12", "U12", "O12")
      hilf
    }
  }
} # end function stat.both


# Function for bootstrapping, that gives L, Lt, U and Ut
stat.both.rr <- function(datanew, index, K, Nknew, Z1, Z2, initial = NULL) {
  # index is from 1 to sum Nk=K*Nk=500 in most cases
  # data1 is vector of dimension Nk*K=500
  # data1.bt<-data1[index,];
  # data2.bt<-transf.data.alt(data1.bt,c,K,Nk);
  # print(index)
  r <- 2
  c <- 2
  c2 <- 2^c
  data.bt <- datanew[index]

  data2.bt <- matrix(0, r * K, c2)

  indk <- c(0, cumsum(Nknew))
  for (k in 1:(r * K)) {
    if (Nknew[k] > 0) {
      chosenind <- data.bt[(indk[k] + 1):indk[k + 1]]
      data2.bt[k, ] <- tabulate(chosenind, nbins = c2)
    }
  } # end for k

  # usually zeros=FALSE
  # M.H.Est.alt(data2,c,r,K,Z1,Z2);
  hilf <- M.H.Est.rr(data2.bt, K, Z1, Z2)
  # print(hilf$L)
  # print(index)
  # print(hilf)
  if (is.null(hilf)) {
    # cat("NULL \n")
    initial
  } else {
    L0 <- hilf$logratio
    U0 <- hilf$Var.logratio
    O0 <- hilf$VarL1 + hilf$VarL2
    hilf <- c(L0, U0, O0)
    names(hilf) <- c("Lalt", "Ualt", "Oalt")
    hilf
  } # end if else
} # end function stat.both


# Begin FUNCTION: show.elapsed.time
show.elapsed.time <- function(start.time, string = "") {
  end.time <- proc.time()[3]
  time <- end.time - start.time
  hours <- floor(time / 3600)
  minutes <- floor(time / 60) - hours * 60
  seconds <- time - hours * 3600 - minutes * 60
  cat(string, "Time: ", hours, "h", minutes, "min", seconds, "sec", "\n")
  return(time)
} # End FUNCTION: show.elapsed.time

# logit function
logit <- function(pi) {
  log(pi / (1 - pi))
}
# expit function
expit <- function(eta) {
  exp(eta) / (1 + exp(eta))
}


# compute mse for lPsi-methods
compute.mse.odds <- function(est.stat, index, true.stat, names, VAR = T) {
  l <- dim(index)[2]
  if (length(true.stat) == 1) {
    true.stat <- matrix(true.stat, l, 1)
  }
  var0 <- m0 <- mse0 <- matrix(0, 1, l)
  for (i in 1:l) {
    mse0[i] <- mean((est.stat[, i][index[, i]] - true.stat[i])^2)
    m0[i] <- mean(est.stat[, i][index[, i]])
    if (VAR) {
      var0[i] <- var(est.stat[, i][index[, i]])
    }
  }
  colnames(mse0) <- names
  colnames(m0) <- names
  if (VAR) {
    colnames(var0) <- names
    return(list(mse = mse0, m = m0, var = var0))
  } else {
    return(list(mse = mse0, m = m0))
  }
} # end compute.mse.lPsi


##########################################################################################
##########################################################################################
M.H.Est.rr <- function(data2, K, Z1, Z2, VAR = TRUE) {
  # c<-2 and r=2, for simplicity
  r <- 2
  c <- 2
  c2 <- 2^c
  X <- Z1[[1]] %*% t(data2) # p
  X1 <- Z1[[2]] %*% t(data2) # 1-p
  h11 <- Z2[[1]] %*% t(data2) # p11 pairwise prob
  h10 <- Z2[[2]] %*% t(data2) # p10 pairwise prob
  h01 <- Z2[[3]] %*% t(data2) # p01 pairwise prob
  h00 <- Z2[[4]] %*% t(data2) # p00 pairwise prob
  n <- matrix(1, 1, c2) %*% t(data2)
  # cc2<-c*(c-1)/2;
  # transform Pij to arrays pij
  UT0 <- lower.tri(diag(c))
  # 2 is number of rows per strata
  UT <- array(UT0, c(c, c, K * r)) # make index for array from index from matrix
  X11 <- X10 <- X01 <- X00 <- array(0, c(c, c, K * r)) # create NULL arrays
  # assign values from matrices to uppertri arrays
  X11[UT] <- h11
  X10[UT] <- h10
  X01[UT] <- h01
  X00[UT] <- h00
  X11 <- aperm(X11, c(2, 1, 3))
  X10 <- aperm(X10, c(2, 1, 3))
  X01 <- aperm(X01, c(2, 1, 3))
  X00 <- aperm(X00, c(2, 1, 3))

  # c=2 items, therefor 2 rr
  U <- matrix(0, c, c) # variance and covariance estimator
  Psi <- matrix(1, c, c)

  C1ba <- C1ab <- 0 # item 1
  C2ba <- C2ab <- 0 # item 2

  V1t1 <- V1t2 <- V1t3 <- 0 # 3 variance terms  for item 1
  V2t1 <- V2t2 <- V2t3 <- 0 # 3 variance terms  for item 2

  c1 <- c2 <- 0 # 2 covariance terms  for cov between items 1 and 2

  for (k in 1:K) {
    a <- (k - 1) * r + 1 # first row
    # actually index for b=2
    b <- a + 1 # second row
    nk <- n[a] + n[b]
    nk <- nk
    nk2 <- nk^2



    if (nk > 0) {
      # print(n0)
      # computation of Cab and Cba
      c1ab <- X[1, a] * n[b] / nk # not divided by n0
      c2ab <- X[2, a] * n[b] / nk
      c1ba <- X[1, b] * n[a] / nk # not divided by n0
      c2ba <- X[2, b] * n[a] / nk



      # var item 1
      h1 <- c1ab * n[b] / nk
      h2 <- c1ba * n[b] / nk
      V1t1 <- V1t1 + h1
      V1t3 <- V1t3 + h2
      V1t2 <- V1t2 + h1 + h2

      # var item 2
      h1 <- c2ab * n[b] / nk
      h2 <- c2ba * n[b] / nk
      V2t1 <- V2t1 + h1
      V2t3 <- V2t3 + h2
      V2t2 <- V2t2 + h1 + h2

      # cov item 1 and 2
      if (n[a] > 1 & n[b] > 1) {
        c1 <- c1 + (n[b]^2 / nk^2) * (h11[a] - (X[1, a] * X[2, a] - h11[a]) / (n[a] - 1))
        c2 <- c2 + (n[a]^2 / nk^2) * (h11[b] - (X[1, b] * X[2, b] - h11[b]) / (n[b] - 1))
      } # end


      C1ab <- C1ab + c1ab
      C1ba <- C1ba + c1ba
      C2ab <- C2ab + c2ab
      C2ba <- C2ba + c2ba
    }
  } # end for K


  # print(C1ab);print(C2ab);print(C1ba);print(C2ba)


  #


  if (C1ba == 0 | C2ba == 0 | C1ab == 0 | C2ab == 0) {
    return(NULL)
  } else {
    logratio <- log(C1ab) + log(C2ba) - log(C1ba) - log(C2ab)
    theta1 <- C1ab / C1ba
    theta2 <- C2ab / C2ba
    L1 <- log(theta1)
    L2 <- log(theta2)
    Var1 <- V1t1 / (2 * C1ab^2) + V1t2 / (2 * C1ab * C1ba) + V1t3 / (2 * C1ba^2)
    Var2 <- V2t1 / (2 * C2ab^2) + V2t2 / (2 * C2ab * C2ba) + V2t3 / (2 * C2ba^2)
    Cov12 <- (c1 + c2 * theta1 * theta2) / (C1ab * C2ab)



    Var.logratio <- Var1 + Var2 - 2 * Cov12

    return(list(L1 = L1, L2 = L2, logratio = logratio, Var.logratio = Var.logratio, VarL1 = Var1, VarL2 = Var2, CovL1L2 = Cov12))
  } # end if
} # end function


# outout functions

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


output1 <- function(mat, row.names = NULL, digits = 3, nsmall = 2) {
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
          output1 <- c(output1, " & $", format(mat[i, j], nsmall = nsmall, digits = digits), "$")
        }
        if (j %% 2 == 0) {
          output1 <- c(output1, "  ($", format(mat[i, j], nsmall = nsmall, digits = digits), "$)")
        }
      } else {
        output1 <- c(output1, " & -- ")
      }
    }
    output1 <- c(output1, "\\\\\n")
  } # end
  cat(output1, "\n")
} # end
