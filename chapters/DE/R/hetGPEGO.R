
library(tidyr)
library(tidyverse)
library(parallel)
library(hetGP)

# change dir
sr=function()
{
	source("~/2021/onyambu/hetGPEGO.R")
}

setwd("~/2021/onyambu")
require(Meta4Design)	# version 0.1.5 has UniProX

version=0 #	0 or 1

dir.hetGP = paste0("hetGPv", version)
dir.data=paste0("dataV", version)

# UniProX matches myDEx1e6 in HPO-DE.R, 6/10/24
ftest.UniPro = Meta4Design::UniPro ## original version 0
if( version == 1 ) ftest.UniPro = Meta4Design::UniProX # version 1: maxNoImprove=itermax and maxTime=60


hetEGO <- function(path, nsteps = 5, nreps=10, save = TRUE){
  print(path)
  dat <- read.csv(path, row.names = 'X') %>%
    rename(m=k)%>%
    pivot_longer(starts_with('X'), names_to = NULL, values_to = 'y')

  bounds <- dat %>%  reframe(across(-y, range))
  size <- select(bounds, n, m)%>%slice(1)

  dat %>%
    mutate(across(-y, ~scales::rescale(.x, 0:1, unlist(bounds[cur_column()]))),
           .keep = 'used') %>%
    select(-n, -m) %>%
    as.matrix() -> X

  ftest <- function(X, point = FALSE){
    if (is.null(dim(X)))
      X <- matrix(X, 1)
    if(is.null(names(X)))
      colnames(X) <- c("NP", "itermax", "pMut", "pCR", "pGBest")
    if(is.array(X))
      X <- as.data.frame(X)
    X <- imap(X, ~scales::rescale(.x,bounds[[.y]], 0:1)) %>%
      merge(size)
    if(point) return(X)
#    res <- do.call(Meta4Design::UniPro, X)
    res <- do.call(ftest.UniPro, X)
    structure(res$opt.value, class="unipro", skeleton = res)
  }

	ftest.reps = function(x, nreps=10)
	{
		ftest1 = function(dummy, x)	ftest(x)
		res = simplify2array( mclapply(1:nreps, ftest1, x, mc.cores = min(10, nreps)) )
		res 
	}

  model <- hetGP::mleHetGP(X, dat$y, covtype = 'Matern5_2')
  x_y <- vector('list', nsteps)
  y_ord <- numeric(nsteps)

  parallel <-  TRUE # FALSE
   if(parallel) ncores <- detectCores() else ncores <- 1

  for(i in 1:nsteps){
    res <- crit_optim(model, crit = 'crit_EI',ncores = ncores,
                      control = list(multi.start = 50, maxit = 50))
    newX <- res$par
    newZ <- ftest.reps(newX, nreps=nreps)
    if(i < 4){
	   	cat("iteration", i , ": ")
	 	cat("x=(", round(newX,2), "); z=", round(newZ,4) ) 
		# If a replicate is selected
		if(!res$path[[1]]$new) cat(" (replicate)\n") else cat("\n")
	}
	else cat("iteration", i , "\r")
	
    x_y[[i]] <- list(par = newX, val = newZ)
    y_ord[i] <- mean( newZ )	# mean
    
 	if(nreps > 1) 	newX=matrix(newX, nrow=nreps, ncol=length(newX), byrow=T)	# 
   model <- update(object = model, Xnew = newX, Znew = newZ)
  }
  id <- which.min(y_ord)

  # get the results
  result <- list(par = ftest(x_y[[id]][[1]], TRUE),
                 value = mean( x_y[[id]][[2]] ),		# mean Y
                 id = id,
                 model = model, x_y = x_y)

  # save the results to a Rda file
  if(save){
    filename <- str_replace_all(path,
                              setNames(c(dir.hetGP, "Rda"),c(dir.data, "csv")))
    result_name <- chartr("/", "_", gsub("^.*?/|.Rda", "", filename))
    assign(result_name, result)
    if(!dir.exists(dirname(filename)))
      dir.create(dirname(filename), recursive = TRUE)
    save(list = result_name, file = filename)
  }

  # return the parameters and the values:
  return(result[1:2])
}


print.unipro <- function(x){
  x <- c(x)
  NextMethod()
}


test=function()
{
	result <- hetEGO(paste0(dir.data, "/up30x3/ccd3_43.csv"), 5)

# file stored in ~hetGP/file_dirname/file_basename.Rda

	load(paste0(dir.hetGP,"/up30x3/ccd3_43.Rda"))
	up30x3_ccd3_43 # all the data from the results:
	b=up30x3_ccd3_43
	res=rep(0,K)
	for(i in 1:K) res[i] = mean( b$x_y[[i]]$val )
	summary(res)
	plot(res)
}

run.hetGP=function(K=5, dir.wk="up50x5")
{
#	K = 200
	designs=c("ccd3_43", "oacd3_50", "lhd_50", "maximin_50", "maxpro_50")
	for(file in designs){
		print(date())
#		(result <- hetEGO("data/up30x3/ccd3_43.csv", 5))
		filename = paste0(dir.data, "/", dir.wk, "/", file, ".csv")
		a <- hetEGO(filename, K)
	}
	print(date())
}


## all results
plot.res=function(K=200, dir.wk="up30x3", quant=NULL)
{ # oldi=1,2,3,4,5,6 correspond to the summary of oldY
# K=200; n=50; m=5
	load.Rda(dir.wk)		# run this outside the function

	if(!is.na(pmatch("up30x3", dir.wk))) res.all=list(up30x3_ccd3_43, up30x3_oacd3_50, up30x3_lhd_50, up30x3_maximin_50, up30x3_maxpro_50)
	else if(!is.na(pmatch("up50x5", dir.wk))) res.all=list(up50x5_ccd3_43, up50x5_oacd3_50, up50x5_lhd_50, up50x5_maximin_50, up50x5_maxpro_50)
	else if(!is.na(pmatch("up70x7", dir.wk))) res.all=list(up70x7_ccd3_43, up70x7_oacd3_50, up70x7_lhd_50, up70x7_maximin_50, up70x7_maxpro_50)
	else if(!is.na(pmatch("up50x10", dir.wk))) res.all=list(up50x10_ccd3_43, up50x10_oacd3_50, up50x10_lhd_50, up50x10_maximin_50, up50x10_maxpro_50)
	else if(!is.na(pmatch("up60x20", dir.wk))) res.all=list(up60x20_ccd3_43, up60x20_oacd3_50, up60x20_lhd_50, up60x20_maximin_50, up60x20_maxpro_50)
	
	designs=c("ccd3_43", "oacd3_50", "lhd_50", "maximin_50", "maxpro_50")
	# exclude lhd_50
	designs=designs[-3]; res.all = res.all[-3]

	newY=matrix(0, K, length(res.all))
	oldY=matrix(0, 6, length(res.all))	# summary of old Y
	Y0 = rep(0, length(res.all))	# initial plot
	for(j in 1:length(res.all)){
		res = res.all[[j]]
		model = res$model	# final hetGP model
#		model$mult		# Z is re-arranged, need to use the intitial data
#		K0 = nrow(model$X) - length(res$x_y)	# total - additional points, 
#		oldY[ ,j] = summary(model$Z[1:K0], quant)		# initial mean value
#		if(!is.null(quant)) 
#				Y0[j] = quantile(model$Z[1:K0], quant)		# initial mean value
		for(i in 1:K) newY[i, j] = mean( res$x_y[[i]]$val )	 # mean
	}
	
	dimnames(newY)[[2]]=designs
	summary(newY)
	boxplot(newY, main=dir.wk)
#	matplot(t(apply(newY, 2, quantile)[c(1:4),]), type="l", main="Quantiles")	# w/o max
	
	matplot(newY, type="p", main=dir.wk); 
	legend("topright", legend=designs,  col=1:5, lty=1:5)
	
	minY=apply(newY, 2, cummin)	
	if(!is.null(quant)){
		minY0 = rbind(Y0, minY)		# 
		matplot(minY0, type="l", main=paste0(dir.wk, "-Q-", quant)); 
	}
	else 	matplot(minY, type="l", main=paste(dir.hetGP, "/", dir.wk)); 
	legend("topright", legend=designs,  col=1:5, lty=1:5)
	invisible(list(oldY=oldY, newY=newY)) 
}

load.Rda=function(dir.wk="up30x3")
{
	# dir.wk="up30x3";  # dir.wk="up50x5";  # dir.wk="up70x7"; #  dir.wk="up60x20"
	designs=c("ccd3_43", "oacd3_50", "lhd_50", "maximin_50", "maxpro_50")
	for(file in designs){
#		load("hetGP/up30x3/ccd3_43.Rda")
		filename = paste0(dir.hetGP, "/", dir.wk, "/", file, ".Rda")
		a <- load(filename)
	}
	
}

run.order=function()
{
	# difference is not significnt among designs after calling hetGPEGO
	# possible strategies: (1) try second-order model with EGO, (2) apply adaptive region with less nsteps
	
	sr()
	dir.wk="up30x3"; # dir.wk="up50x5"; # dir.wk="up70x7";
	run.hetGP(K=50, dir.wk=dir.wk)
	
	run.hetGP(K=50, dir.wk="up50x10")
	run.hetGP(K=50, dir.wk="up60x20")

	load.Rda(dir.wk)	# need to run manually
	
	# reps=10, K=100
	b=plot.res(K=50, dir.wk=dir.wk); 	summary(b$newY)
#	b$oldY # summary of old Y
	
	# run1: K=200
	b=plot.res(K=200, dir.wk=dir.wk, quant=0.25)
	
}