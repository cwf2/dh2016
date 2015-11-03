stopwords.maxsamples <- 0.5
sample.size <- 50
offset <- 0
max.clusters <- 12
nreps <- 15
output.base.dir <- "data/sampletest"
ncores <- 4

# create output directory
if(! dir.exists(output.base.dir)) {
  dir.create(output.base.dir, recursive = T)
}

sample.test <- function(sample.size = 50, offset = 0) {
  t0 <- Sys.time()
  on.exit(cat(sample.size, "/", offset, ":", difftime(Sys.time(), t0, units = "min"), "minutes\n"))

  output.dir <- file.path(output.base.dir, paste(sample.size, offset, sep = "-"))
  if(dir.exists(output.dir)) {
    unlink(output.dir, recursive=T)
  }
  dir.create(output.dir, recursive = T)

  # Generate samples

  cat("Generating samples: size =", sample.size, "; offset =", offset, "\n")
  samples <- make.samples(sample.size, offset)

  # Generate feature vectors

  cat("Building tf-idf weighted document term matrix\n")

  feat.tfidf <- as.matrix(DocumentTermMatrix(
    x = VCorpus(VectorSource(
      sapply(samples$text, paste, collapse = " ")
    )),
    control = list(
      weighting = weightTfIdf,
      bounds = list(global = c(2, stopwords.maxsamples * nrow(samples)))
    )
  ))

  # Authorship adjustment

  cat ("Trying to remove author signal\n")

  sig.base <- colMeans(feat.tfidf)

  for (name in levels(samples$auth)) {
    cat(" -", name, "\n")
    sig.auth <- colMeans(feat.tfidf[samples$auth == name,]) - sig.base
    feat.tfidf[samples$auth == name,] <- t(t(feat.tfidf[samples$auth == name,]) - sig.auth)
  }

  # Classify

  lapply(2:max.clusters, function(nclusters) {

    output.file <- file.path(output.dir, paste("k", nclusters, sep="-"))

    cat("Generating", nreps, "classifications with k =", nclusters, "\n")

    cluster.one.rep <- function(i) {
      t0 <- Sys.time()
      on.exit(
        cat(
          paste("  [", i, "/", nreps, "]", sep=""),
          "...",
          difftime(Sys.time(), t0, units = "min"),
          "minutes\n"
        )
      )

      sample.class <- kmeans(feat.tfidf, nclusters)$cluster
      verse.class <- rep(NA, nrow(la.verse))
      sample.to.verse <- do.call(
        what = Vectorize(function(start, end) {seq.int(start, end)}),
        args = samples[,.(start, end)]
      )
      verse.class[unlist(sample.to.verse)] <- rep(
        x = sample.class,
        times = lapply(sample.to.verse, length)
      )

      return(verse.class)
    }

    if(is.na(ncores)) {
      invisible(write.table(
        x = do.call(cbind, lapply(1:nreps, cluster.one.rep)),
        file = output.file
      ))
    } else {
      invisible(write.table(
        x = do.call(cbind, mclapply(1:nreps, cluster.one.rep, mc.cores = ncores)),
        file = output.file
      ))
    }
  })
}

for (sample.size in seq(from = 30, to = 70, by = 10)) {
  for (offset in seq(from = 0, to = sample.size - 5, by = 5)) {
    sample.test(sample.size, offset)
  }
}


###

# read results

k <- 6

params <- data.frame(unname(do.call(rbind, sapply(dir("data/sampletest"), strsplit, split="-"))))
colnames(params) <- c("size", "offset")
bar <- do.call(cbind,
  lapply(dir("data/sampletest"), function(param) {
    cat(" ", param, "\n")
    x <- read.table(file.path("data/sampletest", param, paste("k", k, sep="-")))

    apply(combn(ncol(x), 2), 2, function(i) {
      adjustedRandIndex(x[,i[1]], x[,i[2]])
    })
  })
)
