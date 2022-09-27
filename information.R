### this is the basic file defining the s.Cij function.
### The inputs should be:
### (1) number of treatment "t". some check should be made here.
### (2) number of period "p". some check should be made here.
### (3) corvariance matrix "sigma".
### (4) model index, denoted by "model" with model=0 as NULL:
# "0" for crossover design      # "1" for interference model
### (5) the effect to estimate, denoted by "effect" with effect=0 as NULL:
# "0" for direct effect         # "1" for total effect
### (6) the flag of circular/non-circular, denoted by "flag_circular" with effect=0 as NULL:
# "0" for non-circular          # "1" for circular
library(nloptr)
library(partitions)
library(MASS)
library(Matrix)
library(psych)
library(plyr)
library(permute)
library(combinat)
library(lpSolve)
s.Cij1<-function(s,t,p,model=1,effect=1,flag_circular=0,sigma=diag(p)){
  if(model==1 || model==0)
  {
    ident=diag(t) # this is only used for constructional purpose
    sigma=sigma # covariance matrix
    #for (i in 1:p){
     # for(j in 1:p){
      #    sigma[i,j]=0.2^abs(i-j)
      #}
    #}
    sigma_inverse=solve(sigma)
    v=sigma_inverse
    tb=v-v%*%matrix(1,p,p)%*%v/sum(v)
    bt=diag(1,t)-matrix(1/t,t,t)
    B.tilde=sigma_inverse-sigma_inverse%*%((rep(1,p))%*%t(rep(1,p)))%*%
      sigma_inverse/as.numeric(t(rep(1,p))%*%sigma_inverse%*%rep(1,p))
    tu=ident[s,]
    if(model==0){
      if(model==0&&effect==0&&flag_circular==0){ld=rbind(rep(0,t),tu[1:(p-1),]);fu=ld}
      if(model==0&&effect==0&&flag_circular==1){ld=rbind(tu[p,],tu[1:(p-1),]);fu=ld}
      if(model==0&&effect==1&&flag_circular==0){ld=rbind(rep(0,t),tu[1:(p-1),])-tu;fu=ld}
      if(model==0&&effect==1&&flag_circular==1){ld=rbind(tu[p,],tu[1:(p-1),])-tu;fu=ld}
      coefb1=tr(bt%*%t(tu)%*%tb%*%tu%*%bt);coefb2=tr(bt%*%t(tu)%*%tb%*%fu%*%bt);coefb3=tr(bt%*%t(fu)%*%tb%*%fu%*%bt)
      C00=t(tu)%*%B.tilde%*%tu;C01=t(tu)%*%B.tilde%*%ld;C11=t(ld)%*%B.tilde%*%ld
      return(list(C00=C00,C01=C01,C11=C11,coefb1=coefb1,coefb2=coefb2,coefb3=coefb3))
    }
    if(model==1){
      if(model==1&&effect==0&&flag_circular==0){ld=rbind(rep(0,t),tu[1:(p-1),]);rd=rbind(tu[2:p,],rep(0,t));fu=ld+rd}
      if(model==1&&effect==0&&flag_circular==1){ld=rbind(tu[p,],tu[1:(p-1),]);rd=rbind(tu[2:p,],tu[1,]);fu=ld+rd}
      if(model==1&&effect==1&&flag_circular==0){ld=rbind(rep(0,t),tu[1:(p-1),])-tu;rd=rbind(tu[2:p,],rep(0,t))-tu;fu=ld+rd}
      if(model==1&&effect==1&&flag_circular==1){ld=rbind(tu[p,],tu[1:(p-1),])-tu;rd=rbind(tu[2:p,],tu[1,])-tu;fu=ld+rd}
      coefb1=tr(bt%*%t(tu)%*%tb%*%tu%*%bt);coefb2=tr(bt%*%t(tu)%*%tb%*%fu%*%bt);coefb3=tr(bt%*%t(fu)%*%tb%*%fu%*%bt)
      C00=t(tu)%*%B.tilde%*%tu;C01=t(tu)%*%B.tilde%*%ld;C02=t(tu)%*%B.tilde%*%rd
      C11=t(ld)%*%B.tilde%*%ld;C12=t(ld)%*%B.tilde%*%rd;C22=t(rd)%*%B.tilde%*%rd
      return(list(C00=C00,C01=C01,C02=C02,C11=C11,C12=C12,C22=C22,coefb1=coefb1,coefb2=coefb2,coefb3=coefb3))
    } 
  }
  else{print("A self-defined model is used. Make sure you define the information matrix of each given sequence correctly.")}
}
