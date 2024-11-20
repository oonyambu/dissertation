UniPro2X <- function(N, k, q=NULL, maxTime=10)
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


get_design <- function(name, n, d){	# levels: 0, ..., n-1
  switch(name,
         UniPro = UniPro2X(n, d, q=n)$xbest - 1,		#
         MaxPro = MaxPro::MaxProLHD(n, d)$Design * n  - 0.5,
         #randomLHD = replicate(d, sample(n)-1),
         randomLHS = floor(lhs::randomLHS(n, d) * n),
         maximin = SLHD::maximinSLHD(1, n, d)$Design -1		# better than maximinLHS
         #maximinLHS = floor(lhs::maximinLHS(n, d) * n)
  )
}



get_model <- function(name, fun, n, d){
  count <- 1
  repeat{
    X <- get_design(name, n, d)/(n - 1)
    x <- rbind(X, X)
    y <- apply(x, 1, fun)
    if(!is.null(dim(y))) y <- colMeans(y)
    prdata <- hetGP::find_reps(x, y, rescale = FALSE, normalize = FALSE)
    model <- hetGP::mleHetGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult),
                             Z = prdata$Z, covtype = "Matern5_2",
                             settings = list(trace = -1))
    if(inherits(model, 'hetGP')) break
    if(count > 5) {warning('HomGP model fitted');break}
    count <- count + 1
  }
  model
}

library(hetGP)
hetEGO <- function(i, name, fun, n, d, nsteps, ys_file, ncores){
  model <- get_model(name, fun, n %/% 2, d)
  ymins <- numeric(nsteps + 1)
  ymin <- min(model$Z0)
  ymins[1] <- ymin
  for(j in seq_len(nsteps)){
    cat("\t", i, name, ' iteration: ',j,'\r')
    res <- crit_optim(model, 'crit_EI', h = 0, ncores = ncores)
    newX <- matrix(res$par, 1)
    newZ <- fun(newX)
    z <- mean(newZ)
    if(ymin > z){
      ymin <- z
      xmin <- newX
    }
    ymins[j + 1] <- ymin
    model <- hetGP::rebuild(update(object = model, Xnew = newX[rep(1,length(newZ)),, drop = FALSE], Znew = newZ), TRUE)
  }
  cat("\t\t\t\t\t\r\r")
  write.table(data.frame(i, name, n,d, seq(0,nsteps), ymins), file = ys_file,
              append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

design_compare <- function(fun_name, domain, n = NULL, nsteps = 50, result_dir = "data/compareHetEGO",
                           overwrite = TRUE, nreps = 100){

  d <- length(domain$lower)
  if(is.null(n)) n <- 10*d
  if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
  ys_file <- file.path(result_dir,sprintf('%s_%dx%d_ymin.txt', fun_name, n, d))
  if(!file.exists(ys_file) | overwrite)
    cat('iteration design nrow ncol nstep y \n',  file = ys_file)
  ncores <- if(.Platform$OS.type == 'windows') 1 else parallel::detectCores()%/% 2
  designs <- c('UniPro', 'MaxPro', 'randomLHS', 'maximin')
  fn <- function(x)match.fun(fun_name)(x * (domain$upper - domain$lower) + domain$lower)
  reps <- function(name)
    lapply(seq_len(nreps), hetEGO, name, fn, n, d, nsteps, ys_file, ceiling(ncores/2))
  invisible(parallel::mclapply(designs, reps, mc.cores = ncores))
}

# branin <- function(x)  rnorm(10, egoOptim::branin(x), runif(10))
# branin_mean10 <- function(x) rowMeans(matrix(rnorm(1000, egoOptim::branin(x), runif(1000)), 10))
# design_compare('branin', egoOptim::domain('branin'))
# design_compare('branin_mean10', egoOptim::domain('branin'))


camel3 <- function(x)  rnorm(10, egoOptim::camel3(x), runif(10))
camel3_mean10 <- function(x)  replicate(10,mean(rnorm(100, egoOptim::camel3(x), runif(100))))
design_compare('camel3', egoOptim::domain('camel3'))
design_compare('camel3_mean10', egoOptim::domain('camel3'))

# camel6 <- function(x)rnorm(10, egoOptim::camel6(x), runif(10))
# camel6_mean <- function(x) mean(camel6(x))
# camel6_mean10 <- function(x) replicate(10, camel6_mean(x))
# invisible(design_compare('camel6', egoOptim::domain('camel6')))
# invisible(design_compare('camel6_mean', egoOptim::domain('camel6')))
# invisible(design_compare('camel6_mean10', egoOptim::domain('camel6')))
#
# ackley_mean10 <- function(x)matrix(rnorm(1000, egoOptim::ackley(x), runif(1000)),1000)%>%colMeans()
#
# invisible(design_compare('ackley_mean10', egoOptim::domain('ackley')(4)))
# invisible(design_compare('ackley_mean10', egoOptim::domain('ackley')(6)))
# invisible(design_compare('ackley_mean10', egoOptim::domain('ackley')(8)))

lineplot <- function(path, transformer = identity, nstep_max = NULL){
  suppressMessages(suppressWarnings(library(ggplot2)))
  read.table(path, header = TRUE) %>%
    dplyr::filter(nstep <= if(is.null(nstep_max)) max(nstep) else nstep_max) %>%
    ggplot(aes(nstep, y, color = design)) +
    stat_summary() +
    labs(title = path) +
    stat_summary(geom = 'line', fun = mean, linewidth=1)
}










