
library(MASS);
source("Odds.R")   # that's where the UTI data comes from
source("Odds_suppl.R")  # functions, e.g. to get MH estimators
source("DIF.fcts.r")


obtain.counts<-function(y,c,pairwise=TRUE){
  # c is number of items
  #print(y)

  K<-dim(y)[1]

  #print(K)

  #y is K by c matrix
  hilf<- matrix(1,1,K)%*%y
  X<- hilf
  X1<- K-hilf # negative counts are simply 1-X

  if(pairwise){
    ind<-indexfct(c);
    L<-dim(ind)[2]
    X00<-X10<-X01<-X11<-rep(0,L)
    for(j in 1:L){
      #print(y[,ind[,j]])
      hilf<- matrix(y[,ind[,j]],K,2)
      X11[j]<- sum((hilf[,1]==1) * (hilf[,2]==1))
      X10[j]<- sum((hilf[,1]==1) * (hilf[,2]==0))
      X01[j]<- sum((hilf[,1]==0) * (hilf[,2]==1))
      X00[j]<- sum((hilf[,1]==0) * (hilf[,2]==0))
      # cat(X00[j],X10[j],X01[j],X11[j])
    }#end for j
    nk<-K
  }#end
  return(list(X=X,X1=X1,X11=X11,X10=X10,X01=X01,X00=X00,nk=nk))
}


obtain.counts.y<-function(y,K,c){

  nk<-rep(NA,2*K)
  X<-X1<-matrix(NA,2*K,c)
  X00<-X11<-X01<-X10<-matrix(NA,2*K,c*(c-1)/2)


  for(k in 1:K){
    for(i in 1:2){
      # selects those rows with Total.score== i+4 and domestic or international
      ind<- y[,1]==k & y[,2]==i

      l.ind<-sum(ind)



      hilf<- obtain.counts(matrix(y[ind,3:(c+2)],l.ind,c),c,pairwise=TRUE)

      X[2*(k-1)+i,]<-hilf$X
      X1[2*(k-1)+i,]<-hilf$X1

      X00[2*(k-1)+i,]<-hilf$X00
      X01[2*(k-1)+i,]<-hilf$X01
      X10[2*(k-1)+i,]<-hilf$X10
      X11[2*(k-1)+i,]<-hilf$X11
      nk[2*(k-1)+i]<-hilf$nk



    }#end i
  }#end k
  return(list(X=X,X1=X1,X11=X11,X10=X10,X01=X01,X00=X00,nk=nk))

}

obtain.bootstrap.samples<-function(y,within.group=TRUE){

hilf<-dim(y)
n<-hilf[1]
m<-hilf[2]-2

y.new<-matrix(0,0,m+2)

  for(k in 1:K){

    if(within.group){
    for(i in 1:2){


      # selects those rows with Total.score== i+4 and domestic or international
      ind<- y[,1]==k & y[,2]==i
      nk<-sum(ind)

      s<-sample((1:n)[ind],nk,replace=TRUE)

      y.new<-rbind(y.new,cbind(k,i,y[s,3:(m+2)]))

    }#end i

    }else{
      # whole strata

      ind<- y[,1]==k
      nk<-sum(ind)

      ind1<- y[,1]==k & y[,2]==1
      nk1<-sum(ind1)

      s<-NULL
      while(is.null(s)){
      s<-sample((1:n)[ind],nk,replace=TRUE)
      #make sure in each group is at least 1 observation
      #if(sum(y[s,2]==1)>0 & sum(y[s,2]==2)>0){s<-s}else{s<-NULL}}

      # we sample from whole stratum, but we assign first nk1 to group 1 and next nk2 to group2
      # nk2=nk- nk2
      y.new<-rbind(y.new,cbind(k,c(rep(1,nk1),rep(2,nk-nk1)),y[s,3:(m+2)]))

      }
    }#end else


  }#end k
return(y.new)
}# end function



construct.tests<-function(L,VarL,CovL,alpha=0.05,show=TRUE){
# CI
c<-length(L)
z<-qnorm(1-alpha/2)
#z<-qnorm(0.95)






# for Test 1: at 90%,  Q2 + Q6  are significant positive and Q12 is negative


# optional
# CovL<-CovL+t(CovL)
# diag(CovL)<-VarL

# generalised MH estimator only defined when r>2, i.e. more than 2 rows only, here we have 2 rows
# so generalsied and ordinary coincide

if(!is.null(VarL)){
  # when VarL is given then construct CovL
CovL<- CovL+t(CovL)
diag(CovL)<-VarL
}else{
  #otherwise CovL is already CovL
CovL<-CovL
VarL<-diag(CovL) #and extract diagonal as Var
}
CovLold<-CovL

CI<-rbind(L-z*sqrt(VarL),L+z*sqrt(VarL))
rownames(CI)<-c("lower","upper")
#colnames(CI)<- items


library(corpcor)
library(Matrix)

# now compare
# first adjustment




if(!is.positive.definite(CovL)){
  if(show){cat("Covariance matrix is not pos definite")}


is.positive.definite(CovL)

if(show){cat("eigenvalues of original matrix CovL\n")
print(eigen(CovL)$values)
}

CovLnew<- nearPD(CovL, corr = FALSE,doSym=TRUE,eig.tol=0.01,posd.tol=0.1,maxit = 1000,keepDiag=TRUE)$mat

#nc.1 <- nearPD(pr, conv.tol = 1e-7, corr = TRUE)

if(show){cat("eigenvalues of new matrix Corr\n")
print(eigen(CovLnew)$values)
}



if(show){round(CovL-CovLnew,3)}
CovL<-as.matrix(CovLnew)

# is.positive.definite(CovL)

}else{

  if(show){cat("Provided Covariance matrix is pos definite")}
}#end if(!is.positive.definite(CorrL)){

  hilf<-matrix(NA,c,c)
  #rownames(hilf)<-colnames(hilf)<-items
  lowerdiffL<-upperdiffL<-diffL<-hilf



  for(i in  1:(c-1)){
    for(j in (i+1):c){

  diffL[i,j] <- L[i]-L[j]
  lowerdiffL[i,j] <- L[i]-L[j] - z*sqrt(VarL[i]+VarL[j]-2*CovL[i,j])
  upperdiffL[i,j] <- L[i]-L[j] + z*sqrt(VarL[i]+VarL[j]-2*CovL[i,j])


    }
  }

#cat("estimate of difference\n")
#print(diffL)
#cat("CI lower endpoint\n")
#print(lower)
#cat("CI upper endpoint\n")
#print(upper)

# Which ones both signs of CI the same?
#cat("Plus 1 says that both CI endpoints have same sign, i.e. difference in L is significant\n")
#print(sign(upper) * sign(lower))


# multivariate Wald-statistic

W <- L%*%solve(CovL)%*%t(L)  #df=c
#cat("\n\nWald test W with all ",c," items/Questions  with alpha=",alpha," and W=",W,"\n")
#cat("p-value:",1-pchisq(W,df=c),"\n")
pvalueW<- 1-pchisq(W,df=c)

# print(c)

# apply now multiple testing
Z<- L/sqrt(VarL)
pvalues<-  2*(1 - pnorm(abs(Z)))

#library(FSA)

Wind<- sum(Z^2)
pvalueWind<- 1-pchisq(Wind,df=c)


pad.Bonf<-p.adjust(pvalues,
         method = "bonferroni")

pad.BH<-  p.adjust(pvalues,
           method = "BH")


pad.Holm<-  p.adjust(pvalues,
           method = "holm")


pad.Hochberg<-  p.adjust(pvalues,
           method = "hochberg")


pad.Hommel<-  p.adjust(pvalues,
           method = "hommel")


pad.BY<-  p.adjust(pvalues,
           method = "BY")


Bonf<-(1:c)[pad.Bonf<=alpha]
BH<-(1:c)[pad.BH<=alpha]
Holm<-(1:c)[pad.Holm<=alpha]
Hochberg<-(1:c)[pad.Hochberg<=alpha]
Hommel<-(1:c)[pad.Hommel<=alpha]
BY<-(1:c)[pad.BY<=alpha]

if(length(Bonf)==0){Bonf<-NULL}
if(length(BH)==0){BH<-NULL}
if(length(Holm)==0){Holm<-NULL}
if(length(Hochberg)==0){Hochberg<-NULL}
if(length(Hommel)==0){Hommel<-NULL}
if(length(BY)==0){BY<-NULL}

#print(pvalues)
#print(Z)
#print(L)
#print(VarL)

comb<-  combined.test(pvalues,L,CovL,alpha=alpha,show=show)



return(list(pvalues=pvalues,CovL=CovL,CovLold=CovLold,W=W,Wind=Wind,
            diffL=diffL,lowerdiffL=lowerdiffL,upperdiffL=upperdiffL,
            pvalueW=pvalueW,pvalueWind=pvalueWind,
            Bonf=Bonf,BH=BH,Holm=Holm,Hochberg=Hochberg,Hommel=Hommel,BY=BY,comb=comb,

            pad.Bonf=pad.Bonf,pad.BH=pad.BH,pad.Holm=pad.Holm,pad.Hochberg=pad.Hochberg,pad.Hommel=pad.Hommel,
            pad.BY=pad.BY

            ))




}#end function apply.L


Odds<-function(p){
  p/(1-p)
}#end

Oddsinv<-function(Odds){
  Odds/(1+Odds)
}
