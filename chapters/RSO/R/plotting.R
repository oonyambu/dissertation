ftest <- function(X, size){

  replicate(10,
            do.call(Meta4Design::UniPro, as.list(setNames(c(size, X), names(bounds))))$opt.value)%>%
    mean()
}



library(tidyverse)

path <- "../DE/data/up30x3/ccd3_43.csv"
dat <- read.csv(path, row.names = 'X') %>%
  rename(m=k)%>%
  pivot_longer(starts_with('X'), names_to = NULL, values_to = 'y')
bounds <- dat %>%  reframe(across(-y, range))
size <- select(bounds, n, m)%>%slice(1)

lower <- unlist(bounds[1,-(1:2)])



comp <- function (fun, lower, upper, ..., p = NULL, rho = 0.3, maximize = FALSE,
                  reps = 20L, file = NULL, control = list(), overwrite = FALSE)
{
  fun_name <- gsub("\\W", "", deparse1(substitute(fun)))
  time_file <- sub("\\.([^.]+$)", "_time.\\1", file)
  nsteps <- if (is.null(nsteps <- control$nsteps))
    5
  else nsteps
  if (is.null(control$budget))
    control$budget <- 50
  if (!is.null(file)) {
    if (file.exists(file) && overwrite) {
      file.remove(c(file, time_file))
    }
    cat("\n", strrep("\n#", 80), "\n", sep = "", file = time_file,
        append = TRUE)
    cat(strrep("#", 80), "\n# nsteps: ", nsteps, "\n# rho: ",
        rho, "\n# Lower Bound: (", toString(lower), ")",
        "\n# Upper Bound: (", toString(upper), ")", "\n# Budget: ",
        control$budget, "\n", file = file, append = TRUE,
        sep = "")
  }
  optimal <- control$trueglobal
  if (missing(lower)) {
    dom <- domain(fun)
    fun <- getFromNamespace(fun, "egoOptim")
    if (is.function(dom)) {
      if (is.null(p))
        stop("the dimension `p` must be provided for ",
             fun_name)
      else dom <- dom(p)
    }
    if (is.null(optimal))
      optimal <- dom$opt$f[1]
    lower <- if (!is.null(dom$lower))
      dom$lower
    else rep(0, p)
    upper <- if (!is.null(dom$upper))
      dom$upper
    else rep(1, p)
  }
  if (maximize & is.null(optimal))
    optimal <- -1
  control$trueglobal <- optimal
  res <- setNames(vector("list", 3), c("RSO", "EGO", "TREGO"))
  errors_list <- res
  RScontrol <- modifyList(control, list(basicEGO = FALSE))
  EGcontrol <- modifyList(control, list(basicEGO = TRUE))
  TRcontrol <- modifyList(control, list(basicEGO = TRUE, method = "TREGO"))
  time <- matrix(NA, nrow = reps, ncol = 3 + length(control$expansion_rate),
                 dimnames = list(NULL, unique(c(names(res), paste0("RSO",
                                                                   control$expansion_rate)))))
  parallel::mclapply(seq_len(reps), function(i) {
    X <- lhs::maximinLHS(5 * length(lower), length(lower))
    X1 <- mapply(rescale, data.frame(X), data.frame(rbind(lower,
                                                          upper)))
    cat("\n\nComputing f(x)...")
    y1 <- apply(X1, 1, function(x) (-1)^(maximize) * fun(x,
                                                         ...))
    cat("Done\nRSO ITERATION:", i, "\n")
    t1 <- proc.time()[["elapsed"]]
    res[["RSO"]][[i]] <<- optimize_fun(fun, lower, upper,
                                       ..., X = X1, y = y1, maximize = maximize, rho = rho,
                                       control = modifyList(RScontrol, list(expansion_rate = 0)))
    time[i, "RSO"] <<- proc.time()[["elapsed"]] - t1
    if (!is.null(file)) {
      cat("RSO", time[i, "RSO"], "\n", file = time_file,
          append = TRUE)
      cat("RSO", i, res[["RSO"]][[i]]$errors, "\n", file = file,
          append = TRUE)
    }
    cat("EGO ITERATION:", i, "\n")
    t2 <- proc.time()[["elapsed"]]
    res[["EGO"]][[i]] <<- optimize_fun(fun, lower, upper,
                                       ..., X = X1, y = y1, maximize = maximize, rho = rho,
                                       control = EGcontrol)
    time[i, "EGO"] <<- proc.time()[["elapsed"]] - t2
    if (!is.null(file)) {
      cat("EGO", time[i, "EGO"], "\n", file = time_file,
          append = TRUE)
      cat("EGO", i, res[["EGO"]][[i]]$errors, "\n", file = file,
          append = TRUE)
    }
    cat("TREGO ITERATION:", i, "\n")
    t3 <- proc.time()[["elapsed"]]
    res[["TREGO"]][[i]] <<- optimize_fun(fun, lower, upper,
                                         ..., X = X1, y = y1, maximize = maximize, rho = rho,
                                         control = TRcontrol)
    time[i, "TREGO"] <<- proc.time()[["elapsed"]] - t3
    if (!is.null(file)) {
      cat("TREGO", time[i, "TREGO"], "\n", file = time_file,
          append = TRUE)
      cat("TREGO", i, res[["TREGO"]][[i]]$errors, "\n",
          file = file, append = TRUE)
    }
    if (!is.null(control$expansion_rate) && any(control$expansion_rate >
                                                0)) {
      for (j in control$expansion_rate) {
        nms <- paste0("RSO", j)
        cat(nms, "ITERATION: ", i, "\n", sep = "")
        t4 <- proc.time()[["elapsed"]]
        res[[nms]][[i]] <<- optimize_fun(fun, lower,
                                         upper, ..., X = X1, y = y1, maximize = maximize,
                                         rho = rho, control = modifyList(RScontrol,
                                                                         list(expansion_rate = j)))
        time[i, nms] <<- proc.time()[["elapsed"]] -
          t4
        if (!is.null(file)) {
          cat(nms, time[i, nms], "\n", file = time_file,
              append = TRUE)
          cat(nms, i, res[[nms]][[i]]$errors, "\n",
              file = file, append = TRUE)
        }
      }
    }
  }, mc.cores = if (.Platform$OS.type == "windows")
    1
  else parallel::detectCores())
  len <- (control$budget - if (is.null(n <- control$n))
    10
    else n)/nsteps
  lapply(res, function(x) sapply(x, getElement, "errors"))
  res

}
environment(comp) <- asNamespace("egoOptim")



lineplot <- function(dat, transformer = identity, nstep_max = NULL, dir='..'){
  dat %>%
    mutate(y = match.fun(transformer)(y))%>%
    dplyr::filter(nstep <= if(is.null(nstep_max)) max(nstep) else nstep_max) %>%
    ggplot(aes(nstep, y, color = method)) +
    stat_summary(geom = 'pointrange',fun.data = ~mean_se(.x, 2)) +
    stat_summary(geom = 'line', fun = mean, linewidth=1)
}


boxplots <- function(data){
  data%>% slice_max(nstep)%>%select(y,method)%>%unstack()%>%boxplot()
}

get_data <- function(n, m, K=100){
  cbind(read.table(sprintf("data/up%dx%d_%d.txt", n, m, K),
                   header = TRUE),
        read.table(sprintf("data/up%dx%d_%dEGO.txt",  n, m, K),
                   col.names = "RSO"),
        read.table(sprintf("data/grid/1024_%dx%d.txt", n, m),
                   col.names = "Grid"),
        read.table(sprintf("data/random/1024_%dx%d.txt", n, m),
                   col.names = "Random"))
}
comparison <- function(n, m, K=100, force = FALSE){
  # a<- try(read.table(sprintf("data/up%dx%d_%d.txt",n,m, K), header = TRUE), silent = TRUE)%>% suppressWarnings()
  #
  # if(force || inherits(a, "try-error")) {
  #   s <- read.table(sprintf("data/results%dx%d_par.txt", n, m), header = TRUE)
  #   v <- read.table(sprintf("data/results%dx%d.txt", n, m), h=TRUE)
  #   mn <- slice_min(subset(v, method == 'RSO'), values,with_ties = FALSE)
  #   par <- s[parse_number(mn$ind),]
  #   opt <- replicate(K,
  #                    do.call(Meta4Design::UniPro, as.list(c(n=n, m=m, par)))[[1]])
  #   cat(opt, file=sprintf("data/up%dx%d_%dEGO.txt",n,m, K), sep="\n")
  # }
  boxplot(get_data(n, m, K))
}

corplot <- function(n, m, point=FALSE, return_s=FALSE){
  s <- read.table(sprintf("data/results%dx%d_par.txt", n, m), header = TRUE)
  s <- mapply(scales::rescale, s, from=bounds[-(1:2)])
  if(return_s) return(s)
  if(point) boxplot(s)
  else corrplot::corrplot(cor(s), type = 'lower', diag = FALSE)
}
#
#
# a<-comp(ftest, unlist(bounds[1,-(1:2)]),
#         unlist(bounds[2,-(1:2)]),
#         size = c(n=70, m=7))
#
#
# b <- lapply(a, \(x)stack(data.frame(sapply(x, \(x)fn(x$model@y, 25)))))%>%
#   bind_rows(.id='method')%>%
#   transform(nstep= ave(values, ind,method, FUN=seq),y=values)

#
# lineplot(b,transformer = 'log10')
#
# results70x7 = list(a=a, b=b, fn=fn)
# save(results70x7, file = "results70x7_1.Rda")
# write.table(b, "results70x7_1.txt", row.names = FALSE, quote = FALSE)

