## change some up50x5 to up30x3 for illustrations on 6/3/24
## modify optimal pMut and pCR used in comparison() on 8/24/24

lab <- 0
suppressWarnings(suppressMessages(library(tidyverse,quietly = TRUE,
                                          warn.conflicts = FALSE, verbose = FALSE)))
# READ RESULTS DATA
fn1 <- function(dir = "data/up30x3"){
  path <- paste0(dir, "/hetGP_", basename(dir), ".csv")
  fn <- function(dat){
    dat %>%
      select(contains("test")) %>%
      cbind(name = dat[,1], .) %>%
      na.omit() %>%
      set_names(c('Design', paste0(rep(c("correlation_", "RMSE_"), each = 3), c("lm", "km", "hetGP")))) %>%
#      filter(str_detect(Design, 'ccd|oacd|lhd_(50|243)|full|max')) %>%
      filter(str_detect(Design, 'ccd|oacd|lhd_(50)|max')) %>%
      type.convert(as.is = TRUE)
  }
  a <- if(grepl("csv", path)) read.csv(path, na.strings = "") else readxl::read_xlsx(path)
  lapply(split(a,cumsum(is.na(a[,1]))), fn)%>%
    setNames(c("(a) Testing on the $3^5$ FFD",
               "(b) Testing on the 243 LHD",
               "(c) Testing on the $3^5$ FFD+243 LHD",
               "(d) Testing on the $4^5$ FFD"))
}


tables <- function(dir, digits = 3){

  cap <- gsub(".*?(\\d+)x(\\d+)",
               "Comparison of designs and model evaluations with target size $\\1 \\\\times \\2$",
               basename(dir))
  lab <<- lab + 1
  tbs <- function(lst){
    nms <- names(lst)
    exec(left_join, !!!c(unname(lst), by = "Design")) %>%
      kableExtra::kbl(format = 'latex',caption = cap, digits = digits,
                      label = paste0('tab',lab)) %>%
      gsub("\\w+\\W_(\\w+).(x|y)", "\\1", .) %>%
      kableExtra::add_header_above(c(" " = 1,
                              rep(c(correlation = 3, RMSE = 3),2))) %>%
      kableExtra::add_header_above(c(" " = 1,
                              setNames(c(6,6), nms)),TRUE,escape = FALSE) %>%
      gsub("([|{][rlc])([}])", "\\1|\\2", .) %>%
      gsub("[{]([lrc][|}])", "{|\\1", .)

  }
  v <- lapply(split(fn1(dir), gl(2,2)), tbs)
  v[[1]] <- sub("(?s)\\\\((?!hline).)+$", "", v[[1]], perl=TRUE)
  v[[2]] <- sub("(?s).*?hline", "\n", v[[2]], perl=TRUE)
  structure(paste(v[[1]], v[[2]], "\n"), class = class(v[[1]]), format='latex')
}

# generate boxplots:
generate_boxplots <- function(dir = "data/up30x3"){
  tabs <- fn1(dir)
  plts <- list()
  #for (i in seq(1,4,1)){
  dat <- map_df(tabs, ~.x%>%
    pivot_longer(-Design, names_to = c('grp1', 'model type'), names_sep = '_') %>%
    filter(grp1=='RMSE', # `model type` != 'km',
           str_detect(Design, 'ccd3|oacd|lhd_(50)|max')), #%>%
    #sort_by(with(., ave(value, Design))),
      .id = 'test_on') %>%
    mutate(Design = str_remove(Design, "[0-9_].*"),
           Design = factor(Design, unique(Design)),
           test_on = factor(as.numeric(factor(test_on, unique(test_on)))))


  plts <- dat %>%
    ggplot(aes(Design, value, fill = `model type`)) +
    geom_col(position = position_dodge(0.5), width=0.4) +
    facet_wrap(~test_on,
               labeller = \(x)map(x,
                                ~map(names(tabs), ~latex2exp::TeX(.x)))
               ,shrink = TRUE) +
    theme_bw() +
    theme(#panel.border = element_blank(),
          panel.grid.major = element_blank(),
          #strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'top',
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     vjust = 0.5, hjust=1,
                                     #margin = margin(t=-8)
                                     ),
          axis.ticks.x = element_blank(),
          legend.justification = "center",
          legend.key = element_blank(),
          legend.key.size =  unit(10, "pt"),
          #legend.box.spacing = unit(-4,"pt"),
          #legend.margin=margin(0,0)
          #legend.box.margin=margin(0,0,0,0)
    ) +
    ylab('RMSE') +
    guides(fill=guide_legend(title=""))
  #}
  plts
}
#generate_boxplots()

generate_boxplots_base <- function(dir = "data/up30x3"){
  tabs <- fn1(dir)
  for (x in tabs){
   x %>%
    pivot_longer(-Design,
                  names_to = c('grp1', 'model type'), names_sep = '_') %>%
    filter(grp1=='RMSE',`model type` != 'km',
            str_detect(Design, 'ccd3|oacd|max')) %>%
    #sort_by(with(., ave(value, Design))) %>%
      mutate(Design = str_remove(Design, "[0-9_].*"),
             Design = factor(Design, unique(Design))) %>%
      with(barplot(value~`model type`+Design, beside = TRUE))
  }
}

fn_mse<- function(dir, return_models = FALSE){
  f_models <- function(path){
    read.csv(path) %>%
      {if ("NP" %in% names(.)){
        mutate(., y = rowMeans(pick(X1:last_col())), .keep = 'unused') %>%
          select(!X:k) %>%
   #       rename(N = NP)%>%
          #rename_with(~str_replace_all(.x, nms)) %>%
          select(order(names(.))) %>%
          mutate(across(-y, ~scales::rescale(.x, to=c(-1,1)))) %>%
          mutate(across(-y, ~.^2, .names = "{col}_q")) %>%
          lm(y~(NP + pMut + pGBest + pCR + itermax)^2+., .)
      }}
  }
  nms <- c(NP = "A", pMut = "B", pGBest = "C", pCR ="D", itermax = "E")
  nms <- setNames(names(nms), names(nms)) # Maintain the same names
  b <- list.files(dir, full.names = T, pattern = ".csv$")
  d <- lapply(b, f_models) %>%
    setNames(tools::file_path_sans_ext(basename(b))) %>%
    discard(is.null)
  if(return_models) return(d)
  pred <- list.files(dir, full.names = T, pattern = "(full|lhd_243).*.csv$")

  e <- map(pred, ~read.csv(.x) %>%
             mutate(., y = rowMeans(pick(X1:last_col())), .keep = 'unused')%>%
             select(all_of(names(nms)), y) %>%
        rename_with(~str_replace_all(.x, nms))%>%
          mutate(across(-y, ~scales::rescale(.x, to=c(-1,1)))) %>%
          mutate(across(-y, ~.^2, .names = "q{col}"))) %>%
    setNames(tools::file_path_sans_ext(basename(pred)))

  map_df(d, \(model)map_dbl(e, \(tdata)sqrt(Metrics::mse(tdata$y, predict(model, tdata)))*100), .id = 'design')%>%
    arrange(full_1024)
}

# fn_mse("data/up30x3/")%>%data.frame()
# fn_mse("data/up50x5/")%>%data.frame()
# fn_mse("data/up70x7/")%>%data.frame()


auto.lm <- function(dir, remove_sd = TRUE) {
  pick <- c("ccd3_43", "oacd3_50", "maximin_50", "maxpro_50")  # "full_243", "full_1024")
  models <- fn_mse(dir, TRUE)[pick]
  names(models) <- gsub("_", "\\\\_", names(models))
  models <- texreg::texreg(models, digits = 4)
  if(remove_sd) {
    models <- models %>%
      str_remove_all("(?<!hline)\n\\s+&.*(?=\n)")%>%
      structure(class = oldClass(models))
  }
  models
}





# DISTANCE AND DENSITY:

density.plot.d=function(x, main="distance density")
{	# x is a design
  dat= apply(x, 2, scale.ab, -1, 1) # scale to [-1,1]
  #		if(i==1) dat=round(dat,1) # FFD243, 3-levels: -1,0,1
  #		print(summary(dat))
  dd = sqrt(apply(dat^2,1,sum))
  print(summary(dd))
  hist(dd, freq=F, xlim=c(0, 2.5),  breaks=seq(0,2.5, .1), main=main)
  lines(density(dd, from=0, to=2.5), col="red")
}




read.dat=function(filename, dir="up30x3", path="data/", all=T)
{
  file=paste0(path, dir,"/", filename, ".csv")
  dat=read.csv(file)[,-1];  dim(dat)
  cols=ncol(dat)-c(9:0);
  yy = dat[,cols]
  dat$ybar=apply(yy, 1, mean)
  dat$sigma=apply(yy, 1, sd)
  #	plot(dat$ybar, dat$sigma);
  #	cor(dat$ybar, dat$sigma)
  dat$y2=apply(yy^2, 1, mean)		# use mean 11/27/23
  if(all) dat
  else 	dat[,-cols]
}


density.plot=function(cols=c(3:7),  reps=c(8:17), dir.wk="up50x5", hist = TRUE)
{
  training =c( "ccd3_43", "oacd3_50", "full_243", "lhd_50",  "maximin_50", "maxpro_50","lhd_243") #
  testing= c("lhd_243", "full_243+lhd_243", "full_1024") # "lhd_1024"

  datasets = unique(c(training, testing))

  if(hist){
    par(mfrow=c(3,3))
    par(mar=c(2,4,4,1)+0.1) #
    for(i in 1:length(datasets)){
      dat1 = read.dat(datasets[i], dir=dir.wk)
      density.plot.d(dat1[,cols], main=datasets[i])
    }

  }
  else{
    par(mfrow=c(3,3))
    par(mar=c(2,4,4,1)+0.1) #
    dat1=read.dat(testing[3], dir=dir.wk); ylim=range(dat1$ybar)

    for(i in 1:length(datasets)){
      dat1=read.dat(datasets[i], dir=dir.wk)
      plot(density(dat1$ybar, from=ylim[1], to=ylim[2]), main=datasets[i])
      abline(h=0)
    }
  }

}
scale.ab <- function(x, a=-1, b=1)
{ # scale x to [a,b]
  xmin= min(x); xmax=max(x)
  x=(x-xmin)/(xmax-xmin);
  a+ (b-a)*x
}
#density.plot()



# interaction plots

# # all interaction plots
# for(i in asplit(combn(names(nms), 2), 2))
#   do.call(interaction.plot, c(unname(d[c(i, 'value')]), xlab = i[1],
#                               trace.label=i[2], fixed = TRUE))


my_interaction_plot_ggplot <- function(dir = "data/up30x3", x = "pMut", y = "pCR"){
  read.csv(file.path(dir, "full_1024.csv"), row.names = 'X') %>%
    pivot_longer(starts_with('X')) %>%
    rename(x = all_of(x), y = all_of(y)) %>%
    select(x, y, value) %>%
    mutate(across(c(x,y), factor)) %>%
    ggplot(aes(x, value, color = y, group = y))+
    stat_summary(geom="line", fun = mean) +
    xlab(x) +
    ylab("") +
    labs(color = y)
    #guides(color=guide_legend(title=""))
}


interaction.plot <- function (x.factor, trace.factor, response, fun = mean,
                              type = c("l", "p", "b", "o", "c"), legend = TRUE,
                              trace.label = deparse(substitute(trace.factor)),
                              fixed = FALSE, xlab = deparse(substitute(x.factor)),
                              ylab = ylabel, ylim = range(cells, na.rm = TRUE),
                              lty = nc:1, col = 1, pch = c(1L:9, 0, letters),
                              xpd = NULL, leg.bg = par("bg"), leg.bty = "n",
                              xtick = FALSE, xaxt = par("xaxt"), axes = TRUE,
                              xleg=NULL, yleg=NULL, ...) {
  ylabel <- paste(deparse(substitute(fun)), "of ", deparse(substitute(response)))
  type <- match.arg(type)
  cells <- tapply(response, list(x.factor, trace.factor), fun)
  nr <- nrow(cells)
  nc <- ncol(cells)
  xvals <- 1L:nr
  if (is.ordered(x.factor)) {
    wn <- getOption("warn")
    options(warn = -1)
    xnm <- as.numeric(levels(x.factor))
    options(warn = wn)
    if (!anyNA(xnm))
      xvals <- xnm
  }
  xlabs <- rownames(cells)
  ylabs <- colnames(cells)
  nch <- max(sapply(ylabs, nchar, type = "width"))
  if (is.null(xlabs))
    xlabs <- as.character(xvals)
  if (is.null(ylabs))
    ylabs <- as.character(1L:nc)
  xlim <- range(xvals)
  if (is.null(xleg)) {
    xleg <- xlim[2L] + 0.05 * diff(xlim)
    xlim <- xlim + c(-0.2/nr, if (legend) 0.2 + 0.02 * nch else 0.2/nr) *
      diff(xlim)
  }
  dev.hold()
  on.exit(dev.flush())
  matplot(xvals, cells, ..., type = type, xlim = xlim, ylim = ylim,
          xlab = xlab, ylab = ylab, axes = axes, xaxt = "n",
          col = col, lty = lty, pch = pch)
  if (axes && xaxt != "n") {
    axisInt <- function(x, main, sub, lwd, bg, log, asp,
                        ...) axis(1, x, ...)
    mgp. <- par("mgp")
    if (!xtick)
      mgp.[2L] <- 0
    axisInt(1, at = xvals, labels = xlabs, tick = xtick,
            mgp = mgp., xaxt = xaxt, ...)
  }
  if (legend) {
    yrng <- diff(ylim)
    if (is.null(yleg))
      yleg <- ylim[2L] - 0.1 * yrng
    if (!is.null(xpd) || {
      xpd. <- par("xpd")
      !is.na(xpd.) && !xpd. && (xpd <- TRUE)
    }) {
      op <- par(xpd = xpd)
      on.exit(par(op), add = TRUE)
    }
    # text(xleg, ylim[2L] - 0.05 * yrng, paste("  ",
    #                                          trace.label), adj = 0)
    if (!fixed) {
      ord <- sort.list(cells[nr, ], decreasing = TRUE)
      ylabs <- ylabs[ord]
      lty <- lty[1 + (ord - 1)%%length(lty)]
      col <- col[1 + (ord - 1)%%length(col)]
      pch <- pch[ord]
    }
    legend(xleg, yleg, legend = ylabs, col = col,
           title = if (trace.label == "") NULL else trace.label,
           pch = if (type %in% c("p", "b"))
             pch, lty = if (type %in% c("l", "b"))
               lty, bty = leg.bty, bg = leg.bg)
  }
  invisible()
}

my_interaction_plot_base <- function(x, y, dir = "data/up30x3", height = 4, width = 8){
  v <- call('interaction.plot', as.name(substitute(x)),
       as.name(substitute(y)), substitute(value))
  png("interaction.png")
  read.csv(file.path(dir, "full_1024.csv"), row.names = 'X') %>%
    pivot_longer(starts_with('X'))%>%
    with(eval(v))
  dev.off()

  sprintf("\\begin{figure}
  \\centering
  \\includegraphics[height=%dcm, width=%dcm]{interaction.png}
    \\caption{Interaction plot between pMut and pCR with target size $30\\times 3$}
  \\label{fig:enter-label}
  \\end{figure}", height, width)

}

nms<- c("NP", "itermax", "pMut", "pCR", "pGBest")

my_interaction_plot_base2 <- function(x = 'pMut', dir = "data/up30x3", xleg = 1.5,
                                      yleg = .46){
  dat <- read.csv(file.path(dir, "full_1024.csv"), row.names = 'X') %>%
    pivot_longer(starts_with('X'))
  x <- as.name(substitute(x))
  nms <- setdiff(nms, as.character(x))
  vs <- list(c(1.8, .435), c(1.8, .435), c(1.8,.435), c(1.5, .46))
  par(mfrow = c(2,2))
  for (i in seq_along(nms)){
    v <- call('interaction.plot', x, as.name(substitute(y, list(y = nms[i]))),
              substitute(value), xleg = vs[[i]][1], yleg = vs[[i]][2])
    with(dat, eval(v))
  }
}

#my_interaction_plot_base2()



###############################
#contours
###############################
library(rsm)
code_fun <- function(from, var){
  halfwd <- diff(from)/2
  mid <- sum(from)/ 2
  substitute(x~(var - center) / halfwd,
             list(x = as.name(paste0(var, '_c')),
                  var = as.name(var), center = mid, halfwd = halfwd))%>%
    as.formula()
}




contour_path <- function(dir, plot = TRUE){

  dat <- read.csv(dir, row.names = 'X') %>%
    pivot_longer(starts_with('X')) %>%
    summarise(value = mean(value), .by = c(NP, itermax,pMut,pCR,pGBest, name))

  ranges <- dat %>%
    reframe(across(all_of(nms), range))

  codes <- ranges %>%
    imap(code_fun) %>%
    unname()

  d <- coded.data(dat, formulas = codes)
  mod <- rsm(value ~ SO(NP_c, itermax_c, pMut_c, pGBest_c, pCR_c), d)

  fn <- function(x, get_p = FALSE){
    p <- list(NP_c=1, itermax_c = 1, pGBest_c=1, pMut_c = x[1],pCR_c = x[2])
    if(get_p) return(p)
    predict(mod, p)
  }
  if(plot) contour(mod, ~ pMut_c + pCR_c, at = c(NP_c=1, itermax_c = 1, pGBest_c = 0.95), nlevels = 50)

  vals <- list(par = rsm::canonical(mod,threshold = 0)$xs[c('pMut_c', 'pCR_c')])

  for(i in 1:5) vals <- optim(vals$par, fn, lower = c(-1,-1), upper = 1, method = 'L-BFGS-B')
  c(list(code2val(fn(vals$par, TRUE), codes)), vals)
}

f_contours <- function(file){
  v <- c(3,5,7)
  par(mfrow = c(length(file),3))
  dr <- sprintf("data/up%d0x%d", v, v)
  fls <- outer(dr, file, file.path)
  tt <-   sapply(sub(".csv", "", file),
                 \(x)sprintf("%d0x%d as target",
                  v, v))
  res <- list()
  for (i in seq_along(fls)){
    res[[i]] <- contour_path(fls[i])
    title(tt[i])
  }
  res
}

#f_contours(c("ccd3_43.csv", "oacd3_50.csv"))
#title(outer = TRUE)


# x <- seq(-1, 1, length = 100)
# z <- t(sapply(x, \(i)sapply(x, \(j)fn(c(i,j)))))
# y <- scales::rescale(x, c(0.05, 0.95))
# filled.contour(y,y, z, nlevels = 50)


comparison <- function(n, m, K = 100, result_dir = "data/result", force = FALSE)
{
  name <- sprintf("%s/up%dx%d_%d.txt", result_dir, n, m, K)

  if(!dir.exists(result_dir)) dir.create(result_dir)
  if(!file.exists(name)|force){
    DE1 <- replicate(K, Meta4Design::UniPro(n=n, m=m, NP=100,
                                  itermax=1500, pMut = 0.1, pCR=0.5, pGBest=1, pSelf = 0)$opt)
    DE4 <- replicate(K, Meta4Design::UniPro(n=n, m=m,
                                             NP=100, itermax=1500, pMut = 0.1, pCR=0.5,
                                            pGBest=0.5, pSelf = 0.25)$opt)

#    a <- contour_path(sprintf("data/up%dx%d/oacd3_50.csv", n, m), plot = FALSE)[[1]]
#    opt <- replicate(K, do.call(Meta4Design::UniPro, as.list(c(n=n, m=m, a)))[[1]])
#    res <- list(DE1 = DE1, DE4 = DE4, SO_optim = opt)

#	optimal pMut and pCR based on 100 tries with pCR=1-pMut, 8/24/24
	pMut = 0.15	# for 50X5 and 70x7
	if(n==30 && m==3)	pMut=0.95
	else if(n==50 && m==5)	pMut=0.25
	pCR= 1 - pMut
    opt <- replicate(K, Meta4Design::UniPro(n=n, m=m, NP=100,
                                  itermax=1500, pMut = pMut, pCR=pCR, pGBest=0.95)$opt)
	res <- list(DE1 = DE1, DE4 = DE4, DEoptim = opt)

    write.table(res, file = name)
  }
  boxplot(read.table(name, header = TRUE)[,])#, main = sprintf("%dx%d", n, m))
}


comparisonEGO <- function(n, m, K=50,  name = "oacd3_50", dir = "hetGP1",
                          result_dir = "result", force = FALSE)
{
  load(sprintf("%s/up%dx%d_%d.Rda", result_dir, n, m, K))
  load(sprintf("%s/up%dx%d/%s.Rda", dir, n, m, name))
  result <- get(sprintf("up%dx%d_%s", n, m, name))
  ncores <- parallel::detectCores()
  if(.Platform$OS.type == "windows") ncores <- 1
  name <- sprintf("%s/up%dx%d_%dEGO.Rda", result_dir, n, m, K)
  if(!dir.exists(result_dir)) dir.create(result_dir)
  if(!file.exists(name)){
    EGO <- parallel::mclapply(1:K, \(i)do.call(Meta4Design::UniPro, as.list(result$par))$opt,
                    mc.cores = ncores)%>%unlist()
    save(EGO, file = name)
  }
  load(name)
  boxplot(c(res, EGO = list(EGO)), main = sprintf("%dx%d", n, m))
}

#comparisonEGO(n = 70, m = 7, result_dir = "result", dir="hetGPv10")





