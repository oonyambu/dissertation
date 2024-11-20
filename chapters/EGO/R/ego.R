
EGO <- function(i, name, fun_name, domain, n = NULL, nsteps = 50, result_dir = "data/compareEGO"){
  i <<- i
  fun <- match.fun(fun_name)
  d <- length(domain$lower)
  if(missing(n)||is.null(n)) n <- d*4
  X <- switch(name,
                MaxPro = MaxPro::MaxProLHD(n, d)$Design * n  - 0.5,
                random = replicate(d, sample(n)-1),
                randomLHS = floor(lhs::randomLHS(n, d) * n),
                UniPro= Meta4Design::UniPro(n, d, NP = 100, itermax = 1500, pMut = 0.2,
                                            pCR = 0.5, pGBest = 0.95)$xbest - 1,
                maximin = floor(lhs::maximinLHS(n, d) * n))
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

design_compare <- function(fun, domain, n = NULL, nsteps = 30, result_dir = "data/compareEGO",
                           overwrite = TRUE, nreps = 20){
  library(parallel)
  library(DiceOptim)
  d <- length(domain$lower)
  if(is.null(n)) n <- 4*d
  if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
  ys_file <- file.path(result_dir,sprintf('%s_%dx%d_ymin.txt', fun, n, d))
  if(overwrite && file.exists(ys_file)) file.remove(ys_file)
  cat('iteration design nrow ncol nstep y \n',  file = ys_file, append = TRUE)
  ncores <- if(.Platform$OS.type == 'windows') 1 else parallel::detectCores()
  designs <- c('UniPro', 'MaxPro', 'randomLHS', 'maximin')
  reps <- function(name){
    dx <- function(x) EGO(x, name, fun, domain, n, nsteps = nsteps)
    unlist(lapply(1:nreps, dx))
  }
  mclapply(designs, reps, mc.cores = ncores)
}



# rosenbrock4d <- DiceOptim::rosenbrock4
# sphere6 <- DiceOptim::sphere6
# design_compare('rosenbrock4d', list(lower = numeric(4), upper = rep(1,4)))
# design_compare('sphere6', list(lower = numeric(6), upper = rep(1,6)))

# remotes::install_github('oonyambu/egoOptim')
# rosenbrock <- egoOptim::rosen
# design_compare('rosenbrock', egoOptim::domain('rosen')(4))
# spheref <- egoOptim::spheref
# design_compare('spheref', egoOptim::domain('spheref')(6))
#
# design_compare('rosenbrock', egoOptim::domain('rosen')(8))
#
#
# easom <- egoOptim::easom
# design_compare('easom', list(lower = c(-10,-10), upper = c(10,10)),nreps = 200)

# branin <- DiceKriging::branin #normalized branin
# design_compare('branin1', list(lower = c(0,0), upper = c(1,1)))
# lineplot('data/compareEGO/branin_8x2_ymin.txt')

# branin_un <- egoOptim::branin # unnormalized branin
# design_compare('branin_un', egoOptim::domain('branin'), 20,nsteps = 20)
# lineplot('data/compareEGO/branin_un_20x2_ymin.txt')

# 
# krige <- function(fun_name, domain, n = NULL){
#   fun <- match.fun(fun_name)
#   d <- length(domain$lower)
#   if(missing(n)||is.null(n)) n <- d*4
#   X <- list(maximin = floor(lhs::maximinLHS(n, d) * n),
#             MaxPro = MaxPro::MaxProLHD(n, d)$Design * n  - 0.5,
#             randomLHS = floor(lhs::randomLHS(n, d) * n),
#             UniPro= Meta4Design::UniPro(n, d, NP = 100, itermax = 1500, pMut = 0.2,
#                                           pCR = 0.5, pGBest = 0.95)$xbest - 1
#             )
#   X$UniPro_half <- X$UniPro %% (n%/% 2)
#   test <- data.frame(t(t(lhs::randomLHS(10000, d))*(domain$upper - domain$lower) + domain$lower))
#   
#   fn <- function(X){
#     x <- t(t(X/(n - 1))*(domain$upper - domain$lower) + domain$lower)
#     y <- apply(x, 1, fun)
#     model <- km(design = x, response = y, covtype = "matern5_2", control = list(trace = 0))
#     Metrics::mse(apply(test, 1, fun), predict(model, test, type='UK')$mean)
#   }
#   sapply(X, fn)
# }
# 
