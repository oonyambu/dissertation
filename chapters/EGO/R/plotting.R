suppressMessages(suppressWarnings(library(tidyverse)))

read_data <- function(path, header = TRUE, ..., dir='..'){
  read.table(sub("^/", "", file.path(dir, path)), header = header, ...)
}

persp_plot <- function(fun, x_range, y_range, ...){
  par(mai = c(0.3,0.3,0.1,0.1), mgp = c(0.2,0.2,0))
  x <- do.call(seq, c(x_range, length = list(100)))
  if(missing(y_range)) y <- x
  else y <- do.call(seq, c(y_range, length=list(100)))
  z <- t(sapply(x, \(i)sapply(y, \(j)fun(c(i,j)))))
  persp(x, y, z, ...)
  invisible(list(x=x, y=y, z=z))
}

# box_plot <- function(path, dir = '..'){
#   read_data(path, dir = dir)%>%
#     dplyr::filter(nstep == max(nstep)) %>%
#     ggplot(aes(design, y)) +
#     geom_boxplot()
# }

lineplot <- function(path, transformer = identity, nstep_max = NULL, dir='..'){
  read_data(path, dir = dir) %>%
    dplyr::filter(nstep <= if(is.null(nstep_max)) max(nstep) else nstep_max) %>%
    ggplot(aes(nstep, y, color = design)) +
    stat_summary(geom = 'pointrange',fun.data = ~mean_se(.x, 2)) +
    stat_summary(geom = 'line', fun = mean, linewidth=1)
}


nms <- c('$\\phi_p(D)$', '$\\psi(D)$', '$100CD$', '$100\\phi(D)$', '$\\rho_{ave}$')

comparison <- function(dat){
  cap <- "Comparison of 100 $80\\times 8$ designs using various criteria"
  dat %>%
    select(-iteration)%>%
    summarise(across(everything(), list(mean = mean, median = median)), .by = design) %>%
    pivot_longer(-design, names_to = c('.value', 'summary'), names_sep = "_") %>%
    setNames(c('design', 'summary', nms))%>%
    kableExtra::kable(format = 'latex',digits = 4, escape = FALSE, booktabs = TRUE,
                      caption = cap, label = "comparison")%>%
    kableExtra::row_spec(0, bold = TRUE) %>%
    kableExtra::kable_styling(position = "center",latex_options = "hold_position")%>%
    kableExtra::collapse_rows(1,"top", "full", row_group_label_position = "first")
   # %>%
   #  sub("{l|", "{|l|", x = _, fixed = TRUE) %>%
   #  sub("|r}", "|r|}", x = _, fixed = TRUE) %>%
   #  sub("[!h]", "[h!tb]", x=_, fixed = TRUE)
}


table_results <- function(dat, cap, label){
  dat %>%
    summarise(across(everything(), list(mean = ~mean(.x, na.rm = TRUE),
                                        median = ~median(.x, na.rm = TRUE))), .by = fun) %>%
    pivot_longer(-fun, names_to = c('name', 'name2'), names_sep = "_") %>%
    mutate(value = ifelse(value == min(value, na.rm = TRUE),
                          sprintf("\\textbf{%.4f}", value),
                          sprintf("%.4f", value)), .by = c(fun, name2)) %>%
    pivot_wider() %>%
    kableExtra::kbl(format = 'latex', caption = cap, digits = 4, booktabs = TRUE,
                    label = label, escape = FALSE, tabular = "tabular") %>%
    kableExtra::collapse_rows(1,"top", "full", row_group_label_position = "first") %>%
    # gsub("l|", "l", ., fixed = TRUE) %>%
    # sub("[t]", "", ., fixed = TRUE) %>%
    # sub("{l", "{|l", ., fixed = TRUE) %>%
    # sub("l}", "l|}", ., fixed = TRUE) %>%
 #   gsub("fun|name2", "", .) %>%
    gsub("name2", "", .) %>%
    sub("begin{table}", "begin{table}[h!tb]", ., fixed = TRUE)

}

no_outlier <- function(x){
  y <- 1.5*IQR(x)
  (x+y > quantile(x, .25)) & (x - y < quantile(x, .75))
}

box_plot_rmse <- function(fn, dat, remove_outlier = FALSE){
  filter(dat, fun == fn, if(remove_outlier) no_outlier(value) else TRUE)%>%
    ggplot(aes(design, value)) +
    geom_boxplot() +
    ylab('Normalized RMSE') +
    theme(legend.position="none") +
    theme_bw()
}
