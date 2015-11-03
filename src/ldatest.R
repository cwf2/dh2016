#######

# topics test

one.lda.per.kmeans <- function(ntopics, nclusters, nreps=10, ncores=NA) {

  inner.function <- function(i) {
    t0 <- Sys.time()
    on.exit(
      cat(
        paste(" - [", i, "/", nreps, "]"),
        "...",
        difftime(Sys.time(), t0, units = "min"),
        "minutes\n"
      )
    )
    kmeans(slot(LDA(dtm.tf, k = ntopics), "gamma"), nclusters)$cluster
  }

  cat("Generating", nreps, "reps with", ntopics, "topics and", nclusters, "classes\n")

  if(is.na(ncores)) {
    return(do.call(cbind, lapply(1:nreps, inner.function)))
  } else {
    return(do.call(cbind, mclapply(1:nreps, inner.function, mc.cores = ncores)))
  } 
}

output.dir <- "~/ldatest"
if(! dir.exists(output.dir)) {
	dir.create(output.dir, recursive=T)
}

set.seed(11011981)

lapply(2:20, function(nclusters) {
   output.file <- file.path(output.dir, paste(ntopics, "-", nclusters, ".txt", sep=""))
	write.table(file=output.file, one.lda.per.kmeans(ntopics, nclusters, nreps = 15, ncores))
})

