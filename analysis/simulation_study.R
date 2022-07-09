
########################################################################
# SIMULATION
########################################################################

# source("sim.DIF.new.r")
# Rscript sim.DIF.new.r
# source("Lfct.r")

setup_sim_cluster = function(clusters = detectCores() - 1,
                             sim_environment = .GlobalEnv){
  sim_cluster <- makeCluster(clusters)
  env_objects = ls(envir = sim_environment)
  env_objects = env_objects[!(env_objects %in% c('clusters', "sim_cluster"))]
  
  clusterExport(sim_cluster, env_objects, envir = sim_environment)
  clusterEvalQ(sim_cluster, {
    source("libraries.R")
    source("analysis/Odds.R") # that's where the UTI data comes from
    source("analysis/Odds_suppl.R") # functions, e.g. to get MH estimators
    source("analysis/Lfct.r")
    source("analysis/DIF.fcts.r")
  })
  
  return(sim_cluster)
}

generate_joint_distribution = function(m = 2, 
                                       K = 20, 
                                       Gamma = 1,
                                       OR = 1.5,
                                       seed_val = 1
                                       ){
  set.seed(seed_val) # in order to replicate results
  
  #Parameters for 20 items and 2 rows
  #reference group Raju et al. (1995)
  not_conv_list<-vector("list",m+1)
  
  # beta_ij   m by 2 matrix
  # first refe
  beta_ij<- matrix(rnorm(m,-1,1),m,2,byrow=FALSE)
  
  # group 1 is reference group, group 2 is focal group (groups i= columns, items j =items)
  beta_ij[,1] <- beta_ij[,2] + log(OR)
  
  # tau_jk is m by K mnatrix
  X<-rnorm(K, sd=1)
  X<- X - mean(X)
  tau_jk<- tcrossprod(rep(0.5,m),X)
  
  #now obtain marginal probabilities
  # p is
  or <- matrix(Gamma, m, m)
  rownames(or) <- colnames(or) <- 1:m
  
  p.joint <- list()  # length will be 2*K
  # obtain joint distribution for each group i and stratum k
  for(FH in 0:m){
    
    p.joint[[FH+1]]<-list()
    for(k in 1:K){
      
      # marginal probabilities given by logistic regression model
      if(FH==0){
        p <- expit(tau_jk[,k]+beta_ij[,2]) # take all from H0
      }else{
        
        if(FH==m){
          hilf<-c(beta_ij[1:m,1])
        }else{
          hilf<-c(beta_ij[1:FH,1],beta_ij[(FH+1):m,2]) # first FH items are false (H1), last m-FH are true (H0)
        }
        
        p <- expit(tau_jk[,k]+hilf)
      }
      
      cat("FH=",FH,"  k=",k," p=",round(p,3),"\n")
      
      # estimating the joint-distribution
      hilf  <- ObtainMultBinaryDist(odds = or, 
                                    marg.probs = p,
                                    iter=5e4,
                                    print=F,
                                    tol=1e-10)
      p.joint[[FH+1]][[k]]<-hilf
      # simulating 100,000 draws from the obtained joint-distribution
      #only consider if p.joint[[FH+1]][[k]]$conv == TRUE
    }
    
    not.conv<-0
    not.conv.set<-NULL
    for(k in 1:K){
      if(!p.joint[[FH+1]][[k]]$conv){
        not.conv<-not.conv+1
        not.conv.set<-c(not.conv.set,k)
      }
    }
    
    not_conv_list[[FH+1]]<-not.conv.set
    
    # replace the non-converged strata with other converged strata
    conv.set<- setdiff(1:K, not.conv.set)
    ind.replace<-sample(conv.set, length(not.conv.set))
    
    p.joint[[FH+1]][not.conv.set]<- p.joint[[FH+1]][ind.replace]
    
  }
  
  file.name.gen <- paste("Sim_GenDif_m",m,"_K",K,"_Gamma",Gamma,"_OR",OR,".RData",sep="")
  file.name.gen = glue("data/joint-distributions/{file.name.gen}")
  
  save(p.joint,beta_ij,tau_jk,X,K,Gamma,m,X,OR,not_conv_list,
       file = file.name.gen)
}


f<-function(nk,N,K,X,fixed,nk.min){
  nk<-ceiling(abs(nk))
  #
  val.min<-sum((nk/N-1/K)^2)  # equivalent to (nk-N/K)^2
  cond<-  sum(nk/N*X)-fixed
  res<-val.min+1e6*cond^2
  attr(res,"cond")<-cond
  return(res)
}

gen.nk<-function(nk,N,K,X,fixed,nk.min){
  K<-length(nk)
  ind<-c(1,1)
  cond<-TRUE
  while(ind[1]==ind[2] | nk[ind[1]]>=(N-nk.min) | nk[ind[2]]<=nk.min){
    ind<-ceiling(K*runif(2))
  }
  
  nk[ind[1]]<-nk[ind[1]]+1 #  and add 1
  nk[ind[2]]<-nk[ind[2]]-1 # take 1 away
  
  return(nk)
}

adjust_reference_sample_size = function(X, K, N_ref, mu.delta){
  meanX<-mean(X[1:K])
  nk.min<- ceiling(N_ref/K*0.3)
  
  res <- optim(par=nk, 
               fn=f,
               gr=gen.nk, 
               method = "SANN",
               control = list(
                 maxit = 1e6, 
                 temp = 5000, 
                 trace = TRUE,
                 REPORT = 500),
               N=N_ref, 
               K=K, 
               X=X, 
               fixed=mu.delta+meanX, 
               nk.min=nk.min)
  
  cat(mu.delta)
  print(f(res$par,N=N_ref,K=K,X=X,fixed=mu.delta))
  nk1<-res$par
  
  return(nk1)
}


dif_simulate = function(x){
  # generate data with
  y.sim<-matrix(0,0,m+2)
  
  # nk1 already defined
  # our aim is to have a mean of +1, i.e. if mean is 0.05, then the mean should be 1.05
  # determine now the nks for the
  st.BT<-system.time({
    for(k in 1:K){
      # ROW 1 (reference)
      # generate nk binary vectors of length m for each strata by row
      
      # For large values of m we take independent samples from the joint dist used
      if(m>12){
        draws = ceiling(m/m0)
        
        y1 = map(1:draws,
            function(x){
              RMultBinary(n = nk1[k], mult.bin.dist = p.joint[[1]][[k]])$binary.sequence
            }) %>% 
          reduce(cbind)
        
        y1=y1[,1:m]
        
      } else{
        y1<-RMultBinary(n = nk1[k], mult.bin.dist = p.joint[[1]][[k]])$binary.sequence
      }
      hilf2<-cbind(k,1,y1[,1:m])
      rm(y1)
      y.sim <- rbind(y.sim,hilf2)
      
      # ROW 2 (focal with FH false hypothese)
      # reference group: i=2
      y2<-RMultBinary(n = nk2[k], mult.bin.dist = p.joint[[FH+1]][[k]])$binary.sequence
      
      hilf1<-cbind(k,2,y2[,1:m])
      rm(y2)
      y.sim <- rbind(y.sim,hilf1)
      
    }
  })
  
  # obtain data in particular format
  # items<-1:m
  
  hilf<- obtain.counts.y(y.sim,K,m)
  X<-t(hilf$X)
  X1<-t(hilf$X1)
  
  X00<-t(hilf$X00)
  X01<-t(hilf$X01)
  X10<-t(hilf$X10)
  X11<-t(hilf$X11)
  nk<-t(as.matrix(hilf$nk))
  
  # we don't need more than that
  # Now apply MH estimator
  
  hilf<-  M.H.Est(data2=NULL,c=m,r=2,K,X=X,X1=X1,h00=X00,h10=X10,h01=X01,h11=X11,zeros=zeros,n=nk)
  L<- -hilf[[2]]
  VarL<-hilf[[3]]  # variance of Lab for particular item
  CovL<-hilf[[4]]  # covariance between items, c=4, e.g. Cov(L(item1),L(item2)) which is in entru CovL[1,2]
  
  # Now obtain bootstrap samples
  cat("create Bootstrap samples and calculate Bootstrap Covariance matrix\n")
  
  st<-system.time({
    L.BT <- matrix(NA,BT,m)
    W.BT<-rep(NA,BT)
    Wind.BT<-rep(NA,BT)
    
    for(b in 1:BT){
      y.boot<- obtain.bootstrap.samples(y.sim,within.group=within.group)
      hilf<- obtain.counts.y(y.boot,K,m)
      
      X.BT<-t(hilf$X)
      X1.BT<-t(hilf$X1)
      nk.BT<-t(as.matrix(hilf$nk))
      X00.BT<-t(hilf$X00)
      X10.BT<-t(hilf$X10)
      X01.BT<-t(hilf$X01)
      X11.BT<-t(hilf$X11)
      
      # use formulae for (co)covariance
      hilfBT1<-  M.H.Est(data2=NULL,c=m,r=2,K,X=X.BT,X1=X1.BT,h00=X00.BT,h10=X10.BT,h01=X01.BT,h11=X11.BT,zeros=zeros,n=nk)
      L.BT[b,]<- -hilfBT1[[2]]
      
      # use estimated (co)variance to get test statistic
      hilfBT2<-try(construct.tests(L=hilfBT1[[2]],VarL=hilfBT1[[3]],CovL=hilfBT1[[4]],show=F))
      if(!inherits(hilfBT2,'try-error')){
        W.BT[b]<-c(hilfBT2$W)
        Wind.BT[b]<- c(hilfBT2$Wind)
      }
      
    }#end for
  })
  
  sim_results = list(L = L, VarL = VarL, CovL = CovL)
  
  return(sim_results)
}

#' Title
#'
#' @param Gamma # OR dependence between items
#' @param m0 #number items (max should be 12 as 2^12=4096)
#' @param m # same ability
#' @param mu.delta # 0 or 1
#' @param FH #Number of questions for which the hypothesis H_0 is false, i.e. number of questions with hypothesis of no DIF rejected
#' @param N_ref #Samples from reference group
#' @param N_foc #Samples from focal group
#' @param BT #Number of bootstrap simulations
#' @param sim #Number of simulations
#' @param zeros 
#' @param within.group 
#' @param K #Number of strata in the ability distribution, p.joint has 100 strata, but can be less  here
#' @param OR #Odds ratio for ref vs focal
#'
#' @return
#' @export
#'
#' @examples
dif_simulation = function(Gamma=1,
                          m0=10,
                          m=2,
                          mu.delta=0,
                          FH=1,
                          K=20,
                          N_ref=1000,
                          N_foc=200,
                          OR=1.5,
                          BT=5e2,
                          sim= 1e3,
                          zeros=FALSE,
                          within.group=F,
                          seed_val=1){
  set.seed(seed_val) # in order to replicate results
  upp.tri.ind<-upper.tri(diag(m))
  
  # K, m, Gamma and OR is given
  cat("K, m0, m, Gamma , OR, mufoc, i, FH ",K, m0, m, Gamma , OR, mu.delta,i,FH,"\n")
  
  # Load joint prod dist to use in simulation results from generate_joint_distribution()
  joint_dist_name<-paste("Sim_GenDif_m",m0,"_K",100,"_Gamma",Gamma,"_OR",OR,".RData",sep="")
  load(glue('data/joint-distributions/{joint_dist_name}'))
  
  # simulate sim data sets
  method.names<-c("combined","Bonf","BH","Holm","Hochberg","Hommel","BY")
  n.methods<-7
  
  X<-rnorm(K,sd=1)
  X<- X - mean(X)
  
  K1<-K
  K<-K1
  nk<- rep(N_ref/K,K)
  m<-m1
  
  # adjust nk's to get different abilities in ref and focal group
  if(mu.delta>0){
    nk1 <- adjust_reference_sample_size(X, K, N_ref, mu.delta)
  }else{
    nk1<- rep(N_ref/K,K)
  }
  
  nk2<- rep(N_ref/K,K)
  
  
  log.OR <- c(rep(log(OR),FH),rep(0,m-FH))
  log.OR.diff <- matrix(log.OR,m,m,byrow=FALSE) -matrix(log.OR,m,m,byrow=TRUE)
  log.OR.diff <- log.OR.diff[upp.tri.ind]
  
  m2<-m*(m-1)/2
  
  log.OR.mat       <- matrix(log.OR,sim,m,byrow=TRUE)
  log.OR.diff.mat  <- matrix(log.OR.diff,sim,m2,byrow=TRUE)
  
  K<-K1
  
  CI.L<-CI.L.BT<- array(NA,dim=c(sim,m,2))
  CIdiff.L<-CIdiff.L.BT<-array(NA,dim=c(sim,m2,2))
  Wind.pvalue.BT<- Wind.pvalue<- W.pvalue.BT<-W.pvalue.BT1<-W.pvalue.BT2<-W.pvalue.BT3<- W.pvalue <- Wold.pvalue <- rep(NA,sim)
  L.pvalue<-array(NA,dim=c(sim,m))
  L.BT.pvalue<- array(NA,dim=c(sim,m))
  
  W0 <- W0ind <- W0old <- rep(NA,sim)
  
  sim_cluster = setup_sim_cluster(sim_environment = environment())
  
  cat("FH=",FH,"  K=",K,"\n")
  cat("START SIMULATIONS\n")
  
  simulation_results = parLapply(
    sim_cluster,
    1:sim,
    dif_simulate)
  
  sim_name<-paste("Sim_GenDif_m",m,"_K",100,"_Gamma",Gamma,"_OR",OR,"_sim_results.RData",sep="")
  
  save(simulation_results,
       file = sim_name)
  
  return(simulation_results)
  
}
