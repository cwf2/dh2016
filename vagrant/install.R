# CRAN repository of choice
repos <- "http://cran.rstudio.com"

# package dependencies for jdmdh scripts
pkgs <- c(
    "XML",
    "data.table",
    "mclust",
    "stringi",
    "tm",
    "topicmodels"
)

# quiet install but echo pkg name
install.no.pbar <- function(pkg) {
	cat(paste("installing R package ", pkg, "\n", sep=""))
	install.packages(pkg, repos=repos, quiet=T)
}

# install required packages
res <- lapply(pkgs, install.no.pbar)