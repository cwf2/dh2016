local.lib <- "/home/vagrant/R/x86_64-pc-linux-gnu-library/3.0"
dir.create(local.lib, recursive=T)

options(repos=c(CRAN="http://cran.rstudio.com"))

install.packages(lib=local.lib, pkgs=c(
  "XML",
  "data.table"
))
