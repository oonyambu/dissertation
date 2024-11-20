library(hetGP)
hetEGO <- function(i, name, fun_name, domain, n, nsteps = 50, ys_file = NULL){
  i <<- i
  fun <- match.fun(fun_name)
  ftest <- function(x, point = FALSE){
    X <- x * (domain$upper - domain$lower) + domain$lower
    if(point) X
    else fun(X)
  }
  d <- length(domain$lower)
  if(missing(n)||is.null(n)) n <- d*10
  n <-  n%/%2
  count <- 0
  repeat{
    count <- count + 1
    if(count > 10) {
      warning('HomGP model fitted')
      break
    }
    X <- switch(name,
                MaxPro = MaxPro::MaxProLHD(n, d)$Design * n  - 0.5,
                random = replicate(d, sample(n)-1),
                randomLHS = floor(lhs::randomLHS(n, d) * n),
                UniPro= Meta4Design::UniPro(n, d, NP = 100, itermax = 1500, pMut = 0.2,
                                            pCR = 0.5, pGBest = 0.95)$xbest - 1,
                maximin = floor(lhs::maximinLHS(n, d) * n))
    x <- rbind(X, X)
    x <- x/(n - 1)
    y <- apply(x, 1, ftest)
    if(!is.null(dim(y))) y <- colMeans(y)
    ymin <- min(y)
    xmin <- x[which.min(y),]

    prdata <- hetGP::find_reps(x, y, rescale = FALSE, normalize = FALSE)
    (model <- hetGP::mleHetGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult),
                      Z = prdata$Z, covtype = "Matern5_2", settings = list(trace = -1)))
    if(inherits(model, 'hetGP')) break
  }
  n <- nrow(x)
  ncores <- if(.Platform$OS.type == 'windows') 1 else min(5,parallel::detectCores())
  for(j in 1:nsteps){
    cat(i, name, fun_name, ' iteration:', sprintf("%02d", j), "\b\b")
    res <- crit_optim(model, 'crit_EI', h = 0, ncores = ncores)
    newX <- matrix(res$par, 1)
    newZ <- ftest(newX)
    z <- mean(newZ)
    if(ymin > z){
      ymin <- z
      xmin <- newX
    }
    cat(i, name, n, d, j, ymin,"\n", file = ys_file, append = TRUE)
    model <- hetGP::rebuild(update(object = model, Xnew = newX[rep(1,length(newZ)),, drop = FALSE], Znew = newZ), TRUE)
  }
  par <- ftest(xmin, TRUE)
  list(par = par, value = ymin, model = model)
}



design_compare <- function(fun, domain, n = NULL, nsteps = 50, nreps = 200,
                           overwrite = TRUE, dir = "data/compareHetEGO"){
  library(parallel)
  d <- length(domain$lower)
  if(is.null(n)) n <- 10*d
  ncores <- if(.Platform$OS.type == 'windows') 1 else min(8, parallel::detectCores())
  ys_file <- sprintf('%s/%s_%dx%d_ymin.txt', dir, fun, n, d)
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  if(overwrite && file.exists(ys_file)) file.remove(ys_file)
  cat('iteration design nrow ncol nstep y \n',  file = ys_file, append = TRUE)
  reps <- function(name){
    sapply(1:nreps, \(x)hetEGO(x, name, fun, domain, n,nsteps, ys_file)$value)
  }
  designs <- c('UniPro', 'MaxPro','maximin', 'randomLHS')
  mclapply(designs, reps, mc.cores = ncores)
}


branin <- function(x)  replicate(10, mean(rnorm(10, egoOptim::branin(x), runif(10))))
design_compare('branin', egoOptim::domain('branin'))

ackley <- function(x)  rnorm(10, egoOptim::ackley(x), runif(10))
design_compare('ackley', egoOptim::domain('ackley')(8))


sumpow <- function(x){
  set.seed(i)
  rnorm(10, egoOptim::sumpow(x), runif(10))
}
dom_sum5 <- egoOptim::domain('sumpow')(5)
dom_sum10 <- egoOptim::domain('sumpow')(10)

levy <- function(x){
  rnorm(10, egoOptim::levy(x), runif(10))
}
dom_lev <- egoOptim::domain('levy')(6)


hart6 <- function(x){
  set.seed(i)
  rnorm(10, egoOptim::hart6(x), runif(10))
}
dom_hart <- egoOptim::domain('hart6')

branin <- function(x){
  rnorm(10, egoOptim::branin(x), runif(10))
}

branin_unr <- function(x){
  egoOptim::branin(x)
}

dom_branin <- egoOptim::domain('branin')[1:2]






branin6 <- function(x){
  y <- egoOptim::branin(x[1:2])*egoOptim::branin(x[3:4])*egoOptim::branin(x[5:6])
  rnorm(10,y, runif(10))
}

a <- rep(egoOptim::domain('branin')[1:2], 3)
dom_branin6 <- lapply(split(a,names(a)), unlist)


# mclapply(designs, reps, 'levy', dom_lev, mc.cores = ncores)
# mclapply(designs, reps, 'sumpow', dom_sum5, mc.cores = ncores)
# mclapply(designs, reps, 'sumpow', dom_sum10, mc.cores = ncores)
# mclapply(designs, reps, 'hart6', dom_hart, mc.cores = ncores)
# mclapply(designs, reps, 'spheref', dom_spheref, mc.cores = ncores)

library(tidyverse)
box_plot <- function(path){
  read.table(path, header = grepl('ymin.txt', path))%>%
    filter(nstep == max(nstep)) %>%
    ggplot(aes(design, y)) +
    geom_boxplot()
}



# Plot for branin 8*2. Not random.

lineplot <- function(path, transformer = identity, nstep_max = NULL){
  read.table(path, header = TRUE) %>%
    dplyr::filter(nstep <= if(is.null(nstep_max)) max(nstep) else nstep_max) %>%
    ggplot(aes(nstep, y, color = design)) +
    stat_summary(geom = 'pointrange',fun.data = ~mean_se(.x, 2)) +
    stat_summary(geom = 'line', fun = mean, linewidth=1)
}



branin_r <- function(x){
  rnorm(10, egoOptim::branin(x), runif(10))
}
# mclapply(designs, reps, 'branin_r', dom_branin, mc.cores = ncores)
#
# lineplot("compare_random/branin_unr_8x2_ymin.txt")
# histplot("compare_random/branin_unr_8x2_ymin.txt")
# mclapply(designs, reps, 'branin_unr', dom_branin, mc.cores = ncores)

# mclapply(designs, reps, 'ackley', dom_ack, mc.cores = ncores)

# lineplot('compare_random/ackley_32x8_ymin.txt')


dom_lev <- egoOptim::domain('levy')(8)

spheref <- function(x){
  #set.seed(i)
  rnorm(10, egoOptim::spheref(x), runif(10))
}

dom_spheref <- egoOptim::domain('spheref')(6)


sumpow <- function(x){
  set.seed(i)
  rnorm(10, egoOptim::sumpow(x), runif(10))
}
dom_sum8 <- egoOptim::domain('sumpow')(8)
dom_sum10 <- egoOptim::domain('sumpow')(10)


michal <- function(x){
  rnorm(10, egoOptim::michal(x), runif(10))
}
dom_michal <- egoOptim::domain('michal')(10)


dejong <- function(x){
 rnorm(10, egoOptim::dejong5(x), runif(10))
}



