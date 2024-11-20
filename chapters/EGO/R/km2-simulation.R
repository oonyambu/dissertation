## km2-simulation.R :  modified from km-simulation.R  7/10/24
## call test2-functions.R
## add uniform projection designs in sim1.mse()

##
## R code for JASA paper "A Projection Space-Filling Criterion and Related Optimality Results"
## simulations using kriging models; see Section 3.1 and supplementary material
## Note: It takes a long time to generate 64x15 and 128x31 maximin, maxpro and uniform designs. So these designs are read from saved files (*.txt).
## main function:  Run() to run simulations.
## (authors) 3/11/2023

library(DiceKriging)
require(parallel) # mclapply
require(SLHD)
require(MaxPro)
#require(UniDOE) # not available for R 4.2.1
require(LHD)

# source("km2-simulation.R") # source this file

source("test2-functions.R") # test functions
source("generate-soa.R") # get.soa.new() 

scale.ab <- function(x, a=-1, b=1) 
{ # scale x to [a,b]
	xmin= min(x); xmax=max(x)
	x=(x-xmin)/(xmax-xmin); 
	a+ (b-a)*x 
}

scale01=function(X, mid=F, rand=F)
{ # scale x to [0,1], modified on 7/15/22
	X=X-min(X)
	q = max(X)+1 # number of levels
	if(mid==T){
		if(rand==T){
			u=runif(nrow(X)*ncol(X),0, 1)
			(X+u)/q  # scale to a random point within each interval
		}
		else (X+0.5)/q  # scale to the middle point within each interval
	} 
	else X/(q-1) # scale to [0,1]
}
randomLHD=function(n,k)
{ # generate a random LHD
  x=matrix(0,n,k)
  for(j in 1:k) x[,j]= sample(1:n, n)
  x
}

gen.response = function(D, fun=sim.fun1)
{
   apply(D ,1, fun)  # apply fun to each row
}

randomPermuteSign=function(x)
{ # randomly permute levels
	if(max(x)>1) cat("error: max(x)>1")
  	n=ncol(x)
	sgn=sample(c(0,1), n, rep=T)
  	for(i in 1:n) if(sgn[i]==1) x[,i]=1-x[,i]  # assume x is within [0,1]
  	x
}
randomPermuteCol=function(x)
{ # randomly permute columns
  n=ncol(x)
  a=sample(1:n, n)
  xa=x[,a]
  dimnames(xa)[[2]]=dimnames(x)[[2]]  # change variable names
  xa
}

collapse= function(x0, snew) 
{ # collapse design x0 into snew levels
	x0=x0-min(x0)
  s0=max(x0)-min(x0)+1
  floor(x0/(s0/snew))
}


fit.km=function(D0, y0, nugget=1e-8,  covtype="matern5_2",upper=NULL)
{  	
  fit0 = km(~1, design=data.frame(D0), response=data.frame(y0), covtype=covtype,
            nugget=nugget, upper=upper, control = list(pop.size = 100, trace = FALSE))

  if( sum(coef(fit0)$range < 1e-5) ){ # some theta==0, fit again
    print(fit0)
    fit0 = km(~1, design=data.frame(D0), response=data.frame(y0), covtype=covtype,
              nugget=nugget, upper=upper, control = list(pop.size = 100, trace = FALSE))
    print(fit0)
  }
  
  fit0
}

test.mse=function(fit0, X, Y)
{ # X is the test data
  pre0=predict(fit0,data.frame(X),"UK")$mean;
  mean((pre0-Y)^2)
}

fit.mse=function(x0, nugget=1e-8, xTest, yTest, sim.fun1=sim.fun1, covtype="matern5_2", permuteCol=TRUE)
{
	 # return normalized mse as in Chen et al. (2016)
	if(min(x0)<0 || max(x0)>1) D0=scale01(x0) # scale to 0-1
	else D0=x0 # x0 is within [0,1] already, do not rescale
  	if(permuteCol){
	  	D0 = randomPermuteCol(D0)
  	  	D0 = randomPermuteSign(D0)
  	  	}
	dimnames(D0)[[2]]=paste("V", 1:ncol(x0), sep="")
	y0 = gen.response(D0, fun=sim.fun1)	
	mu0=mean(y0)

  	fit0=fit.km(D0, y0, nugget=nugget, covtype=covtype)  # fit with D0 and y0
   	mse= test.mse(fit0, xTest, yTest)/mean((yTest-mu0)^2)
  if(mse >20) {  # estimate of theta=0, so fit again
    print(fit0)  # 
    fit0=fit.km(x0, nugget=nugget, covtype=covtype)
    mse= test.mse(fit0, xTest, yTest)/mean((yTest-mu0)^2)
  }
  mse
}  

fit.mse1=function(dummy=1, x0, nugget=1e-8, xTest, yTest, sim.fun1=sim.fun1, covtype="matern5_2", permuteCol=TRUE)
{ # used by mclapply
	fit.mse(x0, nugget=nugget, xTest=xTest, yTest=yTest, sim.fun1=sim.fun1, covtype= covtype, permuteCol= permuteCol)
}


read.design=function(N=64, k=15, type)
{ # read design if the file exists
	filename=paste(type,N,"x",k,".txt", sep="")
	if(!file.exists(filename)) return(FALSE)
	x=read.table(filename, h=T)
	as.matrix(x)[,1:k]
}

get.design=function(N=64, k=15, type, readfile=T)
{
	if(readfile){
		a=read.design(N,k,type)
		if(is.matrix(a)) return(a) # successful in reading the design
	} 
		
	D15=get.soa.new(N=64,k=15)+1 # soa64x15
	D20=get.soa.new(N=64,k=20)+1 # soa64x20
	
	D31=get.soa.new(N=128,k=31)+1 # soa128x31
	D40=get.soa.new(N=128,k=40)+1 # soa128x40

	switch(type,
		 lhd=randomLHD(N,k), # randomLHD 
		 maximin=maximinSLHD(1,N,k)$Design, # maximin
		 maxpro=MaxProLHD(n=N, p=k)$Design*N+0.5, # maxpro
		 ud=GenUD(N, k, q=N, crit="CD2")$final, # uniform, extremely slow
		 
		 soa15=D15[,1:k], # soa64x15
		 soa15lhd=OA2LHD(D15)[, 1:k], # first k columns 	 
		 soa20=D20[,1:k], #  soa64x20
		 soa20lhd=OA2LHD(D20[,1:k]), # first k columns
		 	 
		 soa31=D31[,1:k], # soa128x31 
		 soa31lhd=OA2LHD(D31)[, 1:k], # first k columns 
		 soa40=D40[,1:k], # soa128x40 
		 soa40lhd=OA2LHD(D40)[, 1:k], # first k columns 
	)
}


gen.TestX=function(k = pFac, nTest=10000, method=c("LHD", "ghalton", "sobol")[1])
{ 
	if(method=="ghalton") X = ghalton(nTest, d=k)
	else if (method=="sobol") X = sobol(nTest, d=k, rand="Owen.Faure.Tezuka", seed=sample(1:100000,1))
	else X=randomLHD(nTest, k)
	
	X=as.data.frame(X)
 	dimnames(X)[[2]]=paste("V", 1:k, sep="")
	X
}


sim1.mse=function(dummy=1, N, k=pFac, nugget=NULL, X.test=X.test, Y.test=Y.test, sim.fun1=sim.fun1, covtype="matern5_2")
{  # dummy is used for mclapply
# N is number of runs 	
## setting for Shi and Xu (2023)
#	if(N==64) types=c("lhd","maximin", "maxpro",  "ud", "soa15lhd",  "soa20lhd") # 
#	else if(N==128) types=c("lhd","maximin", "maxpro",  "ud", "soa31lhd",  "soa40lhd") # 

	# add upd, upd16q and upd8q 7/10/24
	if(N==64 && k==15) types=c("lhd","maximin", "maxpro",  "ud", "upd", "upd16q","upd8q")  
	else if(N==128 && k==31) types=c("lhd","maximin", "maxpro",  "ud",  "upd", "upd16q") 
		
	q=length(types)
	mse=rep(0,q)
	 for(i in 1:q){ 
	 	D=get.design(N, k, types[i])
	  	mse[i] = fit.mse(D, nugget, xTest=X.test, yTest=Y.test, sim.fun1=sim.fun1, covtype=covtype, permuteCol=TRUE)
		}
	names(mse)=types
 	mse 	
}

Run.sim=function(id=1,N=32, k0=pFac, nRep=100, covtype="matern5_2", sim.mse=sim1.mse)
{
	nugget= 1e-8 # NULL # 1e-8; use a small nugget to avoid singular issue
	mc.cores = if(.Platform$OS.type == 'windows') 1 else detectCores() # number of cores. 
	pFac=fun.lib[[id]]$dim;  
	fun.name=fun.lib[[id]]$name; sim.f1 = fun.lib[[id]]$fun;
	sim.fun1=sim.f1; 
	cat("===> Simulating", pFac, "-dim", fun.name, "function with", k0, "dimensions\n")
	nTest=10000 # test sample size
	X=gen.TestX(k=k0, nTest) # use randomLHD
	X.test=scale01(X) # scale to [0,1]
	Y.test=gen.response(X.test, fun=sim.fun1) # X is required in [0,1]

	print(date()) # the first argument of mcapply is X
	res <- mclapply(X=1:nRep, sim.mse, N=N, k=k0, nugget=nugget, X.test=X.test, Y.test=Y.test,  sim.fun1=sim.fun1, covtype=covtype, mc.cores = mc.cores)
	print(date())
	
	out = matrix(0, length(res), length(res[[1]])) 
	for(i in 1:length(res)) out[i,]=res[[i]]
	colnames(out)=names(res[[1]]); 
	rmse = sqrt(out) # normalized rmse
	m3=apply(rmse,2,max, na.rm=T); round(m3,3)
	title = paste(fun.name, ":", covtype, "N=", N, "k=", k0, "Rep=", nRep)
	par(mfrow=c(1,1)); boxplot(rmse, main=title, ylab="Normalized RMSE")
	if(sum(rmse>0.7)>0 && sum(rmse>0.7)<10){ ## exclude a few cases with rmse > 0.7
		print(apply(rmse>0.7, 2, sum));  rmse[rmse>0.7]=NA
		par(mfrow=c(1,1)); boxplot(rmse, main=title, ylab="Normalized RMSE")
		}
	rmse
}

Run=function()
{
	source("km2-simulation.R") # source this file
	## show available test functions
	show.functions()
	
	## try Borehole function with 10 repetitions
	a=Run.sim(id=5, N=64, k=15, nRep=10, covtype="matern5_2") ## 
	a=Run.sim(id=5, N=128, k=31, nRep=10, covtype="matern5_2") ## 
	
	## try first 8 test functions with 100 repetitions
	pdf("Rplots64x15.pdf", width=8, height=6)
	a <- list()
	for(i in 1:8) 
	  a[[i]]=Run.sim(id=i, N=64, k=15, nRep=100, covtype="matern5_2") 
	dev.off()

	## try all 9 test functions with 100 repetitions
	pdf("Rplots128x31.pdf", width=8, height=6)
	b<- list()
	for(i in 1:9) b[[i]]=Run.sim(id=i, N=128, k=31, nRep=100, covtype="matern5_2")	
	dev.off()
	
}

