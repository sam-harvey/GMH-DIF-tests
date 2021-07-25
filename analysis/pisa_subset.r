load("~/R/DIF/pisa_AUS_bookID_1.RData")


#my working directory
# setwd( "D:/Thomas/PhD-WLG/r/Odds")
# setwd("C:/Users/tsuesse/Dropbox/2row/DIF")
# setwd("~/R/DIF")

# source("sample.code.L.r")

#rm(list=ls());  # delete everything

library(MASS);
source("Odds.R")   # that's where the UTI data comes from
source("Odds_suppl.R")  # functions, e.g. to get MH estimators
#source("gen_mult.R");   # needed for generating artificial data
#source("mph.Rcode.R")   # fitting Lang's code
library(mltools)

# mean(apply(dat[,3:37]==8,1,mean,na.rm=TRUE)>0)


##DATA SETS
alpha<-0.05
within.group<-F
zeros<-F
m<-35

dat<-as.data.frame(pisa.AUS.book1)

dat[,3:(m+2)]<-as.matrix(dat[,3:(m+2)])

hilf<-dat[,3:(m+2)]
hilf[hilf==2 | hilf==1]<-1
hilf[hilf==8]<-0
ind<- apply(hilf==7,1,mean,na.rm=TRUE)>0

dat<-dat[!ind,]
dat[,3:(m+2)]<-hilf[!ind,]

dat[,"Total"]<- apply(dat[,3:(m+2)],1,sum)
hist(dat[,"Total"])
dat[, "group"] <- bin_data(dat[,"Total"], bins=c(-Inf,seq(5,30,by=5), Inf), binType = "explicit")

dat1<-dat[,3:(m+2)]
dat1<-cbind(as.integer(dat[, "group"]),as.integer(dat[,"Gender"]),dat1)

colnames(dat1)<- c("Stratum","Gender",colnames(dat1)[-(1:2)])

h<-dim(dat1)

n<-h[1]
K<-7
y.sim<-as.matrix(dat1)


source("Odds_suppl.R")
source("Lfct.r")
library("mipfp")
source("DIF.fcts.r")
zeros<-FALSE


#(1:m)[abs(Zs)<2.5]
#(1:m)[abs(Zs)>=2.5]
# sub.items<- c(3,5,10,18,23,32,34)
sub.items<- c(10,23,34)


n<-h[1]
m<-m-length(sub.items)
K<-7
y.sim<-y.sim[,-(sub.items+2)]

source("Odds_suppl.R")
source("Lfct.r")
library("mipfp")
source("DIF.fcts.r")
zeros<-FALSE

# y.sim : 1st column is k and 2nd column is row
hilf<- obtain.counts.y(y.sim,K,m)
X<-t(hilf$X)
X1<-t(hilf$X1)

X00<-t(hilf$X00)
X01<-t(hilf$X01)
X10<-t(hilf$X10)
X11<-t(hilf$X11)
nk<-hilf$nk


# we don't need more than that
# Now apply MH estimator

hilf<-  M.H.Est(data2=NULL,c=m,r=2,K,X=X,X1=X1,h00=X00,h10=X10,h01=X01,h11=X11,zeros=zeros,n=nk)
L<- hilf[[2]]
VarL<-hilf[[3]]  # variance of Lab for particular item
CovL<-hilf[[4]]  # covariance between items, c=4, e.g. Cov(L(item1),L(item2)) which is in entru CovL[1,2]

Zs<- c(L/sqrt(VarL))
# nc.1 <- nearPD(pr, conv.tol = 1e-7, corr = TRUE)

# Now obtain bootstrap samples
# dim(y.sim)
ORpair<-OR.pair(X11,X10,X01,X00,K=K,nk=nk)

summary(Zs)
summary(ORpair)
pval.Zs<-2*(1-pnorm(abs(Zs)))

m.pairs<-m*(m-1)/2





cat("create Bootstrap samples and calculate Bootstrap Covariance matrix\n")
BT<-10e3
st<-system.time({
  L.BT <- matrix(NA,0,m)
  W.BT<-NULL
  Wind.BT<-NULL
  ORpair.BT <- matrix(NA,0,m.pairs)


  while(dim(L.BT)[1]<BT){

    y.boot<- obtain.bootstrap.samples(y.sim,within.group=within.group,K=K)
    hilf<- obtain.counts.y(y.boot,K,m)

    X.BT<-t(hilf$X)
    X1.BT<-t(hilf$X1)
    nk.BT<-t(as.matrix(hilf$nk))
    X00.BT<-t(hilf$X00)
    X10.BT<-t(hilf$X10)
    X01.BT<-t(hilf$X01)
    X11.BT<-t(hilf$X11)
    hilfBT1<-  M.H.Est(data2=NULL,c=m,r=2,K,X=X.BT,X1=X1.BT,h00=X00.BT,h10=X10.BT,h01=X01.BT,h11=X11.BT,zeros=zeros,n=nk)

    y.boot.ind<- obtain.bootstrap.samples.ind(y.sim,within.group=TRUE,K=K)
    hilf<- obtain.counts.y(y.boot.ind,K,m)

    X.BT.ind<-t(hilf$X)
    X1.BT.ind<-t(hilf$X1)
    nk.BT.ind<-t(as.matrix(hilf$nk))
    X00.BT.ind<-t(hilf$X00)
    X10.BT.ind<-t(hilf$X10)
    X01.BT.ind<-t(hilf$X01)
    X11.BT.ind<-t(hilf$X11)
    nk<-hilf$nk

    ORpair.BT <- rbind(ORpair.BT,OR.pair(X11.BT.ind,X10.BT.ind,X01.BT.ind,X00.BT.ind,K=K,nk=nk))

    if(!is.null(hilfBT1[[1]])){

      L.BT<- rbind(L.BT,hilfBT1[[2]])

      #VarLBT<-hilfBT1[[3]]  # variance of Lab for particular item
      #CovLBT<-hilfBT1[[4]]  # covariance between items, c=4, e.g. Cov(L(item1),L(item2)) which is in entru CovL[1,2]

      if(dim(L.BT)[1]%%100==0){cat("Iteration",dim(L.BT)[1],"\n")}

      hilfBT2<-try(construct.tests(L=hilfBT1[[2]],VarL=hilfBT1[[3]],CovL=hilfBT1[[4]],show=F))




      if(!inherits(hilfBT2,'try-error')){
        W.BT<-c(W.BT,c(hilfBT2$W));
        Wind.BT<- c(Wind.BT,hilfBT2$Wind)
        }

    }#end

  }#end for
})
print(st)

#X.BT[1:5,1:4]


# hilf<-abs(log(ORpair.BT))>=abs(matrix(log(ORpair),dim(ORpair.BT)[1],741,byrow=TRUE))
hilf0<- ORpair.BT >= matrix(ORpair,dim(ORpair.BT)[1],m*(m-1)/2,byrow=TRUE)
pval.ORpair<-apply(hilf0,2,mean,na.rm=TRUE)
summary(pval.ORpair)




# estimate COV from BT sample
L.Cov.BT<-cov(L.BT)
# nc.1 <- nearPD(pr, conv.tol = 1e-7, corr = TRUE)

# compare Variances of both
# cbind(diag(L.Cov.BT),c(VarL))

#eigen(L.Cov.BT)$values
#eigen(diag(VarL)+CovL+t(CovL))$values

# source("Lfct.r")
upp.tri.ind<-upper.tri(diag(m))
#
cat("using MH Cov estimator\n")
tests<-try(construct.tests(L,VarL,CovL))
# CovL[1:4,1:4]                  # re-constructed estimates
# use corrected CovL
if(!inherits(tests,'try-error')){
  L.Cov <- tests$CovL  # use the corrected CovL which is now labelled L.Cov

  W0<- tests$W
  W0ind<-tests$Wind
  # CI's diffL
  CIdiff.L<-cbind(tests$lowerdiffL[upp.tri.ind],tests$upperdiffL[upp.tri.ind])
  # W
  W.pvalue<-tests$pvalueW  #using standard test
  Wind.pvalue<-tests$pvalueWind
  L.pvalue<- tests$pvalues
}#end


L.Cov.BT[1:4,1:4]  # BT covariance matrix
cat("using BT Cov estimator\n")
tests1<-try(construct.tests(L,VarL=NULL,CovL=L.Cov.BT))
if(!inherits(tests1,'try-error')){
  # hilf1$pvalues
  # hilf$pvalues

  # CI's diffL
  CIdiff.L.BT<-cbind(tests1$lowerdiffL[upp.tri.ind],tests1$upperdiffL[upp.tri.ind])
  W0.BT<- tests1$W
  # W
  W.pvalue.BT3<- tests1$pvalueW  # using bootstrap
  L.BT.pvalue<- tests1$pvalues
}#end


# doesn"t require tests and tests1

#CI's of L's
CI.L<-cbind(c(L-1.96*sqrt(VarL)),c(L+1.96*sqrt(VarL)))
CI.L.BT<-cbind(c(L-1.96*sqrt(diag(L.Cov.BT))),c(L+1.96*sqrt(diag(L.Cov.BT))))







# Ws are W's from BT
# W is just single test statistic from original data set


# hilf1$W
# hilf$W
W.BT1 <- diag(L.BT%*%solve(L.Cov)%*%t(L.BT)) # use same L.Cov for all BT samples
# plot(density(W.BT1));abline(v=hilf$W)
W.BT2 <- diag(L.BT%*%solve(L.Cov.BT)%*%t(L.BT))  # use L.Cov.BT for all BT samples
# plot(density(W.BT2));abline(v=hilf1$W)

Wind.pvalue.BT<- mean(abs(Wind.BT)>=abs(c(W0ind)) ,na.rm=TRUE)
Wind.pvalue<- tests$pvalueWind

W.pvalue.BT1 <- mean(abs(W.BT1)>=abs(c(W0)) ,na.rm=TRUE)
W.pvalue.BT2 <- mean(W.BT2>=c(W0.BT),na.rm=TRUE)
W.pvalue.BT <- mean(W.BT>=c(W0),na.rm=TRUE) # standard method calculate BT for each sample again

W.pvalue
W.pvalue.BT
W.pvalue.BT1
W.pvalue.BT2
W.pvalue.BT3

Wind.pvalue.BT
Wind.pvalue

Wpvals<-c(Wind.pvalue,Wind.pvalue.BT,W.pvalue,W.pvalue.BT,W.pvalue.BT1,W.pvalue.BT2,W.pvalue.BT3)
names(Wpvals)<-c("ind-asymp","ind-BT","asymp","BT","BT-CovL","BT- CovL.BT","asymp - CovL.BT")

pairs<-indexfct(m)

ind<-lower.tri(diag(m))
ORpair.matrix <- matrix(0,m,m)
pval.ORpair.matrix<-matrix(0,m,m)
ORpair.matrix[ind]<-ORpair
ORpair.matrix<-t(ORpair.matrix)
pval.ORpair.matrix[ind]<-pval.ORpair
pval.ORpair.matrix<-t(pval.ORpair.matrix)

save(L,CovL,VarL,L.Cov.BT,Wpvals,ORpair,ORpair.BT,pval.ORpair,pval.ORpair.matrix,ORpair.matrix,pairs,tests,
     file="pisa_subset_results.RData")
