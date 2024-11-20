generate_pdf <- function(file){
  wd <- getwd()
  on.exit(setwd(wd), add = TRUE)
  f <- basename(file)
  all_files <- c(sub("Rnw", "pdf", f), list.files(dirname(file)))
  setwd(dirname(file))
  Sweave(f)
  tools::texi2pdf(sub("Rnw", "tex", f), clean = TRUE)
  file.remove(setdiff(list.files(), all_files))
  system(paste('open', sub("Rnw", "pdf", f)))
  invisible()
}
generate_pdf(commandArgs(TRUE))

