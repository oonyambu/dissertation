# ego2x.R: evaluate initial designs for fastEGO.nsteps(), 7/10/24
# modified from ego.R from Onyambu on 7/9/24; see EGO2.pdf
# update: replace UniPro with UniPro2, maximinLHS with maximinSLHD, which is better
# increase initial design size n=10*k, and nreps=100 or 1000
# 7/13/24: add examine_designs() to compare different types of designs in terms of different criteria

library(parallel)
library(Meta4Design)
library(egoOptim)
library(UniDOE)
library(LHD)	# AvgAbsCor

# source("ego2x.R")		# this file
  
UniPro2X=function(N, k, q=NULL, maxTime=10)
{	# run DE twice with 2 different settings and return the better one,  up to maxTime=10s
# try 1
	a=Meta4Design::myMetaX(N, k, q=q, NP=100, itermax=1e6, method="DE", 
		pMut = 0.1, pCR=0.9, pGBest=1, pSelf=0, beta=1.0,
		type=-10, maxTime= maxTime, maxNoImprove=1000)		# up to maxTime=10s
# try 2
#	b=Meta4Design::UniPro(N, k, q=q, NP = 100, itermax = 1500, pMut = 0.2, pCR = 0.5, pGBest = 0.95)
	b=Meta4Design::myMetaX(N, k, q=q, NP=100, itermax=1500, method="DE", 
		pMut = 0.2, pCR=0.5, pGBest=0.95, pSelf=(1-0.95)/2, beta=1.0,
		type=-10, maxTime= maxTime, maxNoImprove= 1500)		# 
	
	if(a$opt < b$opt)	a
	else b
}


get_design <- function(name, n, d)
{	# levels: 0, ..., n-1
	switch(name,
                MaxPro = MaxPro::MaxProLHD(n, d)$Design * n  - 0.5,
                randomLHD = replicate(d, sample(n)-1),
                randomLHS = floor(lhs::randomLHS(n, d) * n),
                UniPro = UniPro2X(n, d, q=n)$xbest - 1,		#      
                	ud=UniDOE::GenUD(n, d, q=n, crit="CD2")$final -1, # uniform, extremely slow        
                maximin = SLHD::maximinSLHD(1, n, d)$Design -1,		# better than maximinLHS
                maximinLHS = floor(lhs::maximinLHS(n, d) * n))
}

library(DiceOptim)
EGO <- function(i, name, fun_name, domain, n = NULL, nsteps = 50, result_dir = "data/compareEGO"){
  i <<- i
  fun <- match.fun(fun_name)
  d <- length(domain$lower)
  if(missing(n)||is.null(n)) n <- 10*d  # d*4
  X <- get_design(name, n, d)
   x <- t(t(X/(n - 1))*(domain$upper - domain$lower) + domain$lower)
   y <- apply(x, 1, fun)
   model <-km(design = x, response = y, covtype = "matern5_2", control = list(trace = 0))
   ys_file <- file.path(result_dir,sprintf('%s_%dx%d_ymin.txt', fun_name, n, d))
   res <- fastEGO.nsteps(model, fun, nsteps, domain$lower, domain$upper)
   sn <- seq_len(n)
   vals <- res$value
   vals[1] <- min(y, vals[1])
   vals <- cummin(vals)
   write.table(data.frame(iteration = i, design = name, nrow = n, ncol = d,nstep = seq_along(vals),
                          y = vals)%>%unname(),
                file = ys_file, append = TRUE, row.names = FALSE, quote = FALSE)
   vals
}

design_compare <- function(fun, domain, n = NULL, nsteps = 20, result_dir = "data/compareEGO",
                           overwrite = TRUE, nreps=100){

  d <- length(domain$lower)
  if(is.null(n)) n <- 4*d
  if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
  ys_file <- file.path(result_dir,sprintf('%s_%dx%d_ymin.txt', fun, n, d))
  if(overwrite) if(file.exists(ys_file)) file.remove(ys_file)
  cat('iteration design nrow ncol nstep y \n',  file = ys_file, append = TRUE)
  ncores <- if(.Platform$OS.type == 'windows') 1 else parallel::detectCores()
  if(d>8) ncores = max(ncores/2, 1) # use half cores to avoid memory issue when k=10
  designs <- c('UniPro', 'MaxPro', 'randomLHD', 'maximin')
#  reps <- function(name, nreps = 20){	
  reps <- function(name){	
    dx <- function(x)   EGO(x, name, fun, domain, n, nsteps = nsteps)
    unlist(parallel::mclapply(1:nreps, dx, mc.cores = ncores))
  }
#  mclapply(designs, reps, mc.cores = ncores)
  lapply(designs, reps)		# avoid nesting mclapply
}


library(tidyverse)
box_plot <- function(path, nstep_max = NULL){
  read.table(path, header =TRUE)%>%
    filter(nstep <= if(is.null(nstep_max)) max(nstep) else nstep_max) %>%
    ggplot(aes(design, y)) +
    labs(title = path) +
    geom_boxplot()
}

lineplot <- function(path, transformer = identity, nstep_max = NULL){
  suppressMessages(suppressWarnings(library(tidyverse)))
  read.table(path, header = TRUE) %>%
    filter(nstep <= if(is.null(nstep_max)) max(nstep) else nstep_max) %>%
    ggplot(aes(nstep, y, color = design)) +
    stat_summary() + 
    labs(title = path) +
    stat_summary(geom = 'line', fun = mean, linewidth=1) 
}

eval_design=function(x)
{
	require(Meta4Design)
	x=(x-min(x)+1)  # make sure the minimum level is 1, not 0, required for some criterion
	cd2 = 100*UniDOE::DesignEval(x,crit="CD2") # CD2 value, same as cd2(x)
	phip= myphip(x) # phip for L2 distance with pow=15
	maxpro=dmaxpro2(x, scale=T)  # maxpro criterion, same as bid(x, beta=0)
	upd=100*pcd(x) # upd criterion
	cor1 = AvgAbsCor(x)
	c(phip=phip, maxpro=maxpro, cd2=cd2, upd=upd, cor1=cor1)
}


cmp_design1=function(n ,d, plot.it=F)
{	# scatter plot of the first 2 dimensions
	designs <- c('maximin', 'MaxPro', "ud", 'UniPro', 'randomLHD'  )
	a = lapply(designs, function(name) get_design(name, n, d))
	if(plot.it){
		op=par(mfrow=c(2,3))
		for(i in 1:length(designs))		plot(a[[i]], main=designs[i])
		par(op)		
	}
	res=t(sapply(a, eval_design))
	res = as.data.frame(res)
	dimnames(res)[[1]]=designs
#	list(design=designs, res=res)
	res$design=designs
	res
}

cmp_design=function(n ,d, nreps=10)
{
	ncores <- if(.Platform$OS.type == 'windows') 1 else parallel::detectCores()
	cmp1 = function(i){
		a = cmp_design1(n, d)
		a$iteration = i
		a
	}
	  
	res = parallel::mclapply(1:nreps, cmp1, mc.cores = ncores)
#   unlist(res)
	df = NULL
	for(i in 1:nreps) df = rbind(df, res[[i]])
	df
}

examine_designs=function()
{
# make sure the designs are constructed properly,
# maximinSLHD is better than maximinLHS
	k=4; n=10*k; a=cmp_design1(n, k, plot.it=F); a 

	# generate data
	k=8; n=10*k; 
	df=df0=cmp_design(n, k, nreps=100); cor(df[,1:5])
	file=paste0("Res_", n, "x", k, ".txt")
	write.table(df0, file)

	# make plots
	k=8; n=10*k; 
	file=paste0("Res_", n, "x", k, ".txt")
	df=read.table(file, h=T)
	
	cols = as.numeric(factor(df$design))
	des = sort(unique(df$design))
	vars = dimnames(df)[[2]]
#	subs=combn(5,2); 
	pdffile=paste0("Rplot_", n, "x", k, ".pdf")	
	pdf(pdffile, w=12, h=5)
	op=par(mfrow=c(2,5))
	par(mar=c(5,4, 1,0.5)+0.1)		# default c(5,4,4,2)+0.1
	for(i in 1:4) for(j in (i+1):5){
		plot(df[,i], df[,j], col=cols, pch=cols, xlab=vars[i], ylab=vars[j])
	}	
	legend("topleft", legend=des, col=1:5, pch=1:5)
	par(op)
	dev.off()

	# boxplots	
	op=par(mfrow=c(2,5))
	for(i in 1:5)		boxplot(df[,i]~df$design, main=vars[i])

	df1=df
	df=df1[df1$design != "randomLHD",]
	for(i in 1:5)		boxplot(df[,i]~df$design, main=vars[i])
	par(op)

	for(i in 1:5)	{ cat("criterion:", vars[i], "\n"); print( xtabs(df[,i]~design, df)/xtabs(~design, df) ) }

	by(df[,1:5], ~df$design, summary)
#	by(df[,1:5], ~df$design, mean)

}

#	remotes::install_github('oonyambu/egoOptim')

good.examples = function()
{	
	pdf("ego2x_plots.pdf")

# branin: Maxpro and UniPro are better than Maximin and randomLHS with n=10
	k=2; n=5*k		# with n and nreps, 
	ptm = proc.time()
	branin <- egoOptim::branin # unnormalized branin
	a=design_compare('branin', egoOptim::domain('branin'), n=n,  nsteps=10, nreps=200)
	file=paste0('data/compareEGO/branin_', n, "x", k, '_ymin.txt')
	print( lineplot(file) )		
	ptm2 = proc.time() - ptm;  print(ptm2)


# camel3: use n=10*k, UniPro is good 
	k=2; n=10*k		# with n and nreps, 
	ptm = proc.time()
	camel3 <- egoOptim::camel3 # 
	a=design_compare('camel3', egoOptim::domain('camel3'), n=n,  nsteps=10, nreps=200)
	file=paste0('data/compareEGO/camel3_', n, "x", k, '_ymin.txt')
	print( lineplot(file) )		
	ptm2 = proc.time() - ptm;  print(ptm2)

# camel6: use n=10*k, UniPro is good
	k=2; n=10*k		# with n and nreps, 
	ptm = proc.time()
	camel6 <- egoOptim::camel6 # 
	a=design_compare('camel6', egoOptim::domain('camel6'), n=n,  nsteps=10, nreps=200)
	file=paste0('data/compareEGO/camel6_', n, "x", k, '_ymin.txt')
	print( lineplot(file) )		
	ptm2 = proc.time() - ptm;  print(ptm2)


## levy: try k=2, 4, 6, 8, UniPro is good
for(k in c(4,6,8)){
	n=10*k		# 
	ptm = proc.time()
	levy <- function(x) egoOptim:: levy(x) 
	a=design_compare('levy', egoOptim::domain('levy')(k), n=n, nsteps=10, nreps=200)
	file=paste0('data/compareEGO/levy_', n, "x", k, '_ymin.txt')
	print( lineplot(file) )
	ptm2 = proc.time() - ptm;  print(ptm2)
	}


## michal: try k=2, 4, 6, 8, UniPro is good 
for(k in c(4,6,8)){
	n=10*k		
	ptm = proc.time()
	michal <- function(x) egoOptim::michal(x, m=10) # default m=10
	a=design_compare('michal', egoOptim::domain('michal')(k), n=n, nsteps=10, nreps=200)
	file=paste0('data/compareEGO/michal_', n, "x", k, '_ymin.txt')
	print( lineplot(file) )
	ptm2 = proc.time() - ptm;  print(ptm2)
	}

	dev.off()
}

ackley.example=function()
{ # nsteps=30, nreps=200; n= 10*k
  pdf("ackley-r200.pdf")
# ackley: for k=6, 8 and n=10*k; UniPro is better than others;  maximin is worst; randomLHD is better than maxpro and maximin
for(k in c(2,4,6,8)){ # small difference for k=2, 4
  n = 10*k   # use k=6, n=5*k or 10*k,
  ptm = proc.time()
  ackley <- egoOptim::ackley
  a=design_compare('ackley', egoOptim::domain('ackley')(k), n=n, nsteps=30, nreps=200)
  file=paste0('data/compareEGO/ackley_', n, "x", k, '_ymin.txt')
  print( lineplot(file) )		
  ptm2 = proc.time() - ptm;  print(ptm2)
}
  dev.off()
}

tmp=function()
{
  pdf("ackley-r200.pdf")
  
  for(k in c(2, 4, 6, 8)){ # small difference for k=2, 4
    for (n in c(10*k)){   # use k=6, n=5*k or 10*k,
  file=paste0('data/compareEGO/ackley_', n, "x", k, '_ymin.txt')
#      file=paste0('data/ackley-r200/ackley_', n, "x", k, '_ymin.txt')
      print( lineplot(file) )	
    }
  }
  dev.off()
}

other.examples = function()
{

# hartman4: difference is small
	k=4; n=10*k		# with n and nreps, 
	ptm = proc.time()
	hartman4 <- egoOptim::hart4 # 
	a=design_compare('hartman4', egoOptim::domain('hart4'), n=n,  nsteps=10, nreps=100)
	file=paste0('data/compareEGO/hartman4_', n, "x", k, '_ymin.txt')
	print( lineplot(file) )		
	ptm2 = proc.time() - ptm;  print(ptm2)

# hartman6: difference is small
	k=6; n=10*k		# with n and nreps, 
	ptm = proc.time()
	hartman6 <- egoOptim::hart6 # 
	a=design_compare('hartman6', egoOptim::domain('hart6'), n=n,  nsteps=10, nreps=100)
	file=paste0('data/compareEGO/hartman6_', n, "x", k, '_ymin.txt')
	print( lineplot(file) )		
	ptm2 = proc.time() - ptm;  print(ptm2)

#	easom <- egoOptim::easom # error
#	a=design_compare('easom', egoOptim::domain('easom'), nsteps=20)
#	lineplot('data/compareEGO/easom_8x2_ymin.txt')



 for(k in c(2, 4,6)){ # small difference for k=2, 4. for k=6, maximin is worst
	n= 10*k
	ptm = proc.time()
	rosenbrock <- egoOptim::rosen
	a=design_compare('rosenbrock', egoOptim::domain('rosen')(k), n=n, nsteps=20, nreps=100)	# [-2.04, 2.048]
#	a=design_compare('rosenbrock', list(lower = rep(-5, k), upper = rep(10,k)) ) # [-5, 10]
	file=paste0('data/compareEGO/rosenbrock_', n, "x", k, '_ymin.txt')
	print( lineplot(file) )		
	ptm2 = proc.time() - ptm;  print(ptm2)
	}
	
	k=2; n=10*k
	ptm = proc.time()
	spheref <- egoOptim::spheref
	a=design_compare('spheref', egoOptim::domain('spheref')(k), n=n, nsteps=10, nreps=100)	
	file=paste0('data/compareEGO/spheref_', n, "x", k, '_ymin.txt')
	print( lineplot(file) )		
	ptm2 = proc.time() - ptm;  print(ptm2)


	## others:
	## factor-sparse
	ptm = proc.time()
	spheref2k <- function(x) egoOptim::spheref(x[1:2]) # 2 active factors
	k=3; n=10*k
	a=design_compare('spheref2k', egoOptim::domain('spheref')(2*k), nsteps=10)
	file=paste0('data/compareEGO/spheref2k_', n, "x", k*2, '_ymin.txt')
	print( lineplot(file) )		
	ptm2 = proc.time() - ptm;  print(ptm2)

	sumpow <- egoOptim::sumpow
	a=design_compare('sumpow', egoOptim::domain('sumpow')(10))
	
#	branin1 <- DiceKriging::branin #normalized branin
#	a=design_compare('branin1', list(lower = c(0,0), upper = c(1,1)))
#	lineplot('data/compareEGO/branin_8x2_ymin.txt')
	
	## no practical difference among designs
	branin_un <- egoOptim::branin # unnormalized branin
	a=design_compare('branin_un', egoOptim::domain('branin'), nsteps=10)
	lineplot('data/compareEGO/branin_un_8x2_ymin.txt')
	
	## factor-sparse
	ptm = proc.time()
	branin2k <- function(x) egoOptim::branin(x[1:2]) # unnormalized branin
	k=3; lower = rep(c(-5,0), k);  upper = rep(c(10,15), k)
	a=design_compare('branin2k', list(lower=lower, upper=upper), nsteps=10)
	file=paste0('data/compareEGO/branin2k_', 4*k*2, "x", k*2, '_ymin.txt')
	print( lineplot(file) )		
	ptm2 = proc.time() - ptm;  print(ptm2)
	
}
