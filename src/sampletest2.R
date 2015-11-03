stopwords.maxsamples <- 0.5
sample.size <- 50
offset <- 0
max.clusters <- 12
nreps <- 15
output.base.dir <- "data/sampletest2/pca"
ncores <- 4

# create output directory
if(! dir.exists(output.base.dir)) {
  dir.create(output.base.dir, recursive = T)
}

sample.test <- function(sample.size = 50, offset = 0, nreps = 15, pca = F, ncores = NA) {
  t0 <- Sys.time()
  on.exit(cat(sample.size, "/", offset, ":", difftime(Sys.time(), t0, units = "min"), "minutes\n"))

  output.dir <- file.path(output.base.dir, 
    paste(sample.size, sprintf("%02i", offset), sep = "-")
  )
  if(dir.exists(output.dir)) {
    unlink(output.dir, recursive = T)
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
  
  # PCA
  
  if (pca) {
    cat ("Performing PCA\n")
    feat.tfidf <- prcomp(feat.tfidf)$x[,1:500]
  }
  
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
      
      # a matrix with one row per sample
      sample.to.verse <- do.call(
        what = Vectorize(function(start, end) {list(seq.int(start, end))}),
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

for (offset in seq(from = 0, to =  45, by = 5)) {
  sample.test(offset = offset, pca = T)
}


###

# read results

k <- 6

params <- data.frame(unname(do.call(rbind, sapply(dir(output.base.dir), strsplit, split="-"))))
colnames(params) <- c("size", "offset")
bar <- do.call(cbind,
  lapply(dir(output.base.dir), function(param) {
    cat(" ", param, "\n")
    x <- read.table(file.path(output.base.dir, param, paste("k", k, sep="-")))

    apply(combn(ncol(x), 2), 2, function(i) {
      adjustedRandIndex(x[,i[1]], x[,i[2]])
    })
  })
)

all.k <- do.call(cbind, 
  mclapply(2:max.clusters, function(k) {
    cat(" - k =", k, "\n")
    x <- do.call(cbind, 
      lapply(dir(output.base.dir), function(param) {
        read.table(file.path(output.base.dir, param, paste("k", k, sep="-")))
      })
    )
    apply(combn(ncol(x), 2), 2, function(i) {
      adjustedRandIndex(x[,i[1]], x[,i[2]])
    })
  }, mc.cores=4)
)

all.class <- do.call(cbind, 
  lapply(dir(output.base.dir), function(param) {
      cat(" ", param, "\n")
      read.table(file.path(output.base.dir, param, paste("k", k, sep="-")))
  })
)
namask <- ! apply(all.class, 1, anyNA)

err.nreps <- 100

roi <- la.verse[, 
  auth == "ovid" & 
  as.numeric(sub(pattern="\\..*", replacement="", x=loc)) > 10
]

# pdf(file = "plot/theb5.unsupervised.pdf", width=10, height=5)
plot.new()
plot.window(xlim=c(1, sum(roi)), ylim=c(-10, 80))

for (i in 1) {
  e <- apply(replicate(err.nreps, sample(ncol(all.class), 2)), 2, function(i) {
    x <- rep(0, sum(namask))
    x[classError(all.class[namask, i[1]], all.class[namask, i[2]])$misclassified] <- 1
    return(x)
  })

  misclass <- rep(NA, nrow(all.class))
  misclass[namask] <- rowSums(e)
  
  lines(misclass[roi], col=i)
}

abline(v = la.verse[roi, which(sub (pattern=".*\\.", replacement="", x=loc) == "1")])

ticsmask <- c(la.verse[roi, which(loc=="1.1")], la.verse[roi, which(as.numeric(sub(pattern=".*\\.", replacement="", x=loc)) %% 100 == 0)])
axis(
  side = 1, 
  at = ticsmask,
  labels = la.verse[roi, loc][ticsmask],
  cex.axis = 0.9
)
axis(side = 2, at = seq(from = 0, to = err.nreps, by = 20))

title(
  main = "Statius, Thebaid 5",
  xlab = "verse",
  ylab = "classification error (%)"
)

points(rep(-10, sum(roi)), col=all.class[roi, 1], pch=15)
mtext(side = 2, text = "class", at = -10, las=1)

dev.off()
