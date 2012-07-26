GUM.validate <- function(var.name, x.i, u.i, nu.i,
		                 type,distribution,measurement.fnc,
						 correlation=diag(length(var.name)),
						 shared.u.i=var.name, cl=0.95, 
						 cov.factor="Student's t",
						 sig.digits.U=2)
{

# generate muitivariate correlated unif(-1,1)
mvrunif <-
function(n,Corr){ 
p <- ncol(as.matrix(Corr))
2*(pnorm(matrix(rnorm(n*p),ncol=p)%*%chol(Corr))-.5)
}

# generate multivariate chi-squares
r.u.obs <-
function(df, p = nrow(SqrtSigma), SqrtSigma = diag(p)) {
Z <- matrix(0, p, p)
diag(Z) <- sqrt(rchisq(p, df:(df-p+1)))
if(p > 1) {
pseq <- 1:(p-1)
Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
sqrt(diag(crossprod(Z %*% SqrtSigma)/(df)))
}
else Z * SqrtSigma/sqrt(df)
}


### Check for Numeric Input Values for xmat, umat, dfmat, correlation ###
if(!is.numeric(x.i) | !is.numeric(u.i) | !is.numeric(nu.i) | !is.numeric(correlation)) return(NA)
### Check for empty or invalid measurement function ### 
if(measurement.fnc=="" | is.na(measurement.fnc) | is.numeric(measurement.fnc)) return(NA)
### Check if variable names in meqs match namevec ###
measurement.fnc.names <- all.vars(measurement.fnc)
if(any(!(measurement.fnc.names %in% var.name))) return(NA)
### Check for Positive Definite Correlation Matrix ###
if(min(eigen(correlation)$values)<=.Machine$double.eps) return(NA)
#if(length(measurement.fnc.names)==length(var.name)) if(sort(meq.names)!=sort(namevec)) return(NA)
###measurement.fnc,x.i,u.i,nu.i,var.name,gp.ans,shared.u.i,corr.ans,correlation,cov.factor,type,distribution){
       nreps <- 1000
       #set.seed(411)
       n <- length(x.i)
       ### Generate True y and Simulate Input Data ###
       tm.y <- GUM(var.name, x.i, u.i, nu.i, measurement.fnc, correlation, shared.u.i, cl, cov.factor, sig.digits.U)$y
       if(n==1) correlation <- matrix(1)
       obs.xi <- matrix(NA,nrow=nreps,ncol=n)
       u.obs.xi <- matrix(NA,nrow=nreps,ncol=n)
       df.obs.xi <- matrix(rep(nu.i,rep(nreps,n)),nrow=nreps)
       if(is.na(sum(as.numeric(correlation)))) correlation <- diag(n)
       subset <- type=="A"
       if(sum(subset)!=0){
corrltn <- correlation[subset,subset]
tsd <- matrix(u.i[subset],nrow=1)
library(MASS)
obs.xi[,subset] <- mvrnorm(nreps,x.i[subset],t(tsd)%*%tsd*corrltn)
subsetG <- subset
while(sum(subset)!=0){ # Simulate a Block (corrsub  components) of u.obs's  At a Time
 corrsub.old <- rep(F,n)
 corrsub <- correlation[(1:n)[subset][1],]!=0 & subset
 while(sum(corrsub) != sum(corrsub.old)){
  corrsub.old <- corrsub
  for(i in 1:sum(corrsub.old)) corrsub<-corrsub | (correlation[(1:n)[corrsub.old][i],]!=0 & subset)
  }
 C<-chol(correlation[corrsub,corrsub])
 for(i in 1:nreps) u.obs.xi[i,corrsub]<-u.i[corrsub]*r.u.obs(floor(mean(nu.i[corrsub])),sum(corrsub),C)
 subset[corrsub] <- F
 }
if(sum(shared.u.i == var.name)!=n){ ### Reassign Grouped Type A u.obs.xi ###
 shared.u.i <- rank(shared.u.i)
 gname <- unique(shared.u.i)
 for(grp in gname){
  subg <- grp==shared.u.i & subsetG
  if(sum(subg)>=2) u.obs.xi[,subg]<-rep(u.obs.xi[,subg][,1],sum(subg))*rep(u.i[subg]/(u.i[subg][1]),rep(nreps,sum(subg)))
  }
 }
}
       subset <- type=="B"
       if(sum(subset)!=0){
subsubset <- subset&distribution=="Rectangular"
p <- sum(subsubset)
if(p!=0){
 corrltn <- correlation[subsubset,subsubset]
 obs.xi[,subsubset] <- rep(x.i[subsubset],rep(nreps,p))+sqrt(3)*rep(u.i[subsubset],rep(nreps,p))*mvrunif(nreps,corrltn)
 }
subsubset <- subset&distribution=="Triangular"
p <- sum(subsubset)
if(p!=0){
 corrltn <- correlation[subsubset,subsubset]
 obs.xi[,subsubset] <- rep(x.i[subsubset],rep(nreps,p))+sqrt(3/2)*rep(u.i[subsubset],rep(nreps,p))*(mvrunif(nreps,corrltn)+
           mvrunif(nreps,corrltn))
 }
subsubset <- subset&distribution=="Normal"
if(sum(subsubset)!=0){
 corrltn <- correlation[subsubset,subsubset]
 tsd <- matrix(u.i[subsubset],nrow=1)
 obs.xi[,subsubset] <- mvrnorm(nreps,x.i[subsubset],t(tsd)%*%tsd*corrltn)
 }
u.obs.xi[,subset] <- rep(u.i[subset],rep(nreps,length(u.i[subset])))
}
### Apply GUM to Simulated Data and Compute Coverage Prob ###
DP <- GUM(var.name,obs.xi,u.obs.xi,df.obs.xi,measurement.fnc,correlation,shared.u.i,cl,cov.factor,sig.digits.U)
        sum(DP$y-DP$U<tm.y & tm.y<DP$y+DP$U)/nreps
}

