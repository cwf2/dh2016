# Thematic work on Valerius Flaccus

# packages

library(data.table)
library(parallel)
library(tm)
library(topicmodels)
library(mclust)
library(XML)
library(MASS)

source(file.path("src","tesserae.R"))

# functions

load.file <- function(file) {
  # read a single text and produce a data.table

  cat("Reading", file, "... ")

  # parse the Tesserae XML file
  doc <- xmlParse(file)

  # get author and work from Tesserae text id
  tess.id <- xmlGetAttr(getNodeSet(doc, "/TessDocument")[[1]], "id")
  auth <- sub("\\..*", "", tess.id)
  work <- sub(".*\\.", "", tess.id)

  # print a tally when done
  count <- 0
  on.exit(cat(count, "lines\n"))

  # process lines, build a data.table
  rbindlist(
    xpathApply(doc, "//TextUnit", function(node) {
      count <<- count + 1

      data.table(
        auth = auth,
        work = work,
        loc = xmlGetAttr(node, "loc"),
        verse = xmlValue(node)
      )
    })
  )[,
    verse := iconv(verse, from = "ASCII", to = "UTF-8", mark = T)
    ][,
      unitid := .I
      ]
}

load.corpus <- function(index.file) {
  # load a set of Tesserae texts and construct a corpus

  cat("Loading corpus from", index.file, "\n")

  files <- scan(index.file, what="character", sep="\n")
  cat("Reading", length(files), "files\n")

  rbindlist(lapply(files, load.file))
}

make.samples <- function() {
  la.verse[,
    .(auth, work, loc, verse, int.grp = ceiling(unitid / sample.size))
  ][,
    .(
      sampleid = .GRP,
      start = min(.I),
      end = max(.I),
      lstart = loc[which.min(.I)],
      lend = loc[which.max(.I)],
      text = paste(
        unlist(la.stemmer(standardize(unlist(strsplit(verse, " "))))),
        collapse = " "
      )
    ),
    by = .(auth, work, int.grp)
  ][
    nchar(text) > mean(nchar(text)) - 3 * sd(nchar(text))
    & nchar(text) < mean(nchar(text)) + 3 * sd(nchar(text)),
    .(auth = factor(auth), work = factor(work), start, end, sampleid, lstart, lend, text)
  ]
}

cluster.series <- function(x, k = 5, nreps = 10, cores = NA) {
  # generate a series of k-means clusterings

  cat("Generating", nreps, "classifications with k =", k, "\n")

  inner.function <- function(i) {
    t0 <- Sys.time()
    on.exit(
      cat(
        paste(" - [", i, "/", nreps, "]", sep=""),
        "...",
        difftime(Sys.time(), t0, units = "min"),
        "minutes\n")
    )
    kmeans(x, k)$cluster
  }

  if(is.na(cores)) {
    return(do.call(cbind, lapply(1:nreps, inner.function)))
  } else {
    return(do.call(cbind, mclapply(1:nreps, inner.function, mc.cores = cores)))
  }
}

load.scenes <- function(file, author) {
  cat("Loading benchmark themes from", file, "\n")
  
  auth.start <- match(author, la.verse$auth)
  
  getverses <- Vectorize(function(Start, End) {
    paste(
      unlist(
        la.stemmer(
          standardize(
            unlist(
              strsplit(la.verse[Start:End]$verse, " ")
            )
          )
        )
      ), collapse = " "
    )
  })
  
  bench <- data.table(read.table(file, sep="\t", header=T, as.is=T, quote="", fill=T))[,.(Loc, Description, Type)]
  bench[, Loc := sub("\\(?(\\d),(\\d+)\\s*-\\s*(\\d+)\\)?", "\\1.\\2-\\1.\\3", Loc)]
  bench[, Loc := sub("\\(?(\\d),(\\d+)\\s*\\)?$", "\\1.\\2-\\1.\\2", Loc)]
  bench[, Start := match(sub("-.*", "", Loc), la.verse[auth==author, loc]) + auth.start - 1]
  bench[, End := match(sub(".*-", "", Loc), la.verse[auth==author, loc]) + auth.start - 1]
  bench[, Text := getverses(Start, End)]

  return(bench)
}

remapClassToRef <- function(class, ref) {
  cmap <- mapClass(a = ref, b = class)
  
  new.class <- class
  
  for (i in names(cmap$bTOa)) {
    new.class[class == i] <- cmap$bTOa[i]
  }
  
  return(unlist(new.class))
}


#
# main
#

# 0. preliminaries

set.seed(19810111)
sample.size <- 30
stopwords.maxsamples <- 0.5
ntopics <- 50
ncores <- 3

# 1. preprocessing the texts

# load latin corpus
la.verse <- load.corpus("data/texts/la.index.txt")

# load latin stemmer
la.stemmer <- build.stemmer.table(
  stems.file = file.path("data", "tesserae", "la.lexicon.csv"),
  resolve = file.path("data", "tesserae", "la.word.freq")
)

# 2. sampling

# assign samples to verse lines
cat("Sampling\n")

samples <- make.samples()

# 2. Generate feature vectors for samples
#    a) tf-idf weights using tm

cat("Building tf-idf weighted document term matrix\n")

# build samples into a tm-style corpus
tm.corp <- VCorpus(VectorSource(
  sapply(samples$text, paste, collapse = " ")
))

# calculate document-term matrix with tf-idf weights
dtm.tfidf <- DocumentTermMatrix(tm.corp, control=list(
  weighting = weightTfIdf,
  bounds = list(global = c(2, stopwords.maxsamples * nrow(samples)))
))

feat.tfidf <- as.matrix(dtm.tfidf)

cat("Calculating PCA for tf-idf scores\n")
pca.tfidf <- prcomp(feat.tfidf)

# draw a graph showing author distribution
pdf(file.path("plot", "pca.author.tfidf.pdf"), width = 6, height = 6)
plot(pca.tfidf$x,
  col = unclass(samples$auth),
  pch = unclass(samples$auth),
  main = "TF-IDF by author"
)
dev.off()

pdf(file.path("plot", "pca.author.legend.pdf"), width = 4, height = 6)
plot.new()
legend("center",
  legend = levels(samples$auth),
  col = 1:nlevels(samples$auth),
  pch = 1:nlevels(samples$auth)
)
dev.off()

#  2. b) author-adjusted tf-idf weights

sig.base <- colMeans(feat.tfidf)

cat ("Trying to remove author signal\n")
feat.adjusted <- feat.tfidf
for (name in levels(samples$auth)) {
  cat(" -", name, "\n")
  sig.auth <- colMeans(feat.tfidf[samples$auth == name,]) - sig.base
  feat.adjusted[samples$auth == name,] <- t(t(feat.tfidf[samples$auth == name,]) - sig.auth)
}

cat("Calculating PCA for author-adjusted tf-idf scores\n")
proj.pca.adjusted <- predict(pca.tfidf, feat.adjusted)
pca.adjusted <- prcomp(feat.adjusted)

# draw a graph showing author distribution
pdf(file.path("plot", "pca.author.adjusted.pdf"), width = 6, height = 6)
plot(pca.adjusted$x,
  col = unclass(samples$auth),
  pch = unclass(samples$auth),
  main = "TF-IDF: author-adjusted\nby author"
)
dev.off()

# 2. c) LDA

cat("Building tf document-term matrix for LDA\n")
dtm.tf <- DocumentTermMatrix(tm.corp, control=list(
  bounds = list(global = c(2, stopwords.maxsamples * nrow(samples)))
))

cat("Generating topic model with", ntopics, "topics\n")
print(system.time(
  lda <- LDA(dtm.tf, k = ntopics)
))

feat.topics <- slot(lda, "gamma")
pca.topics <- prcomp(feat.topics)

# draw a graph showing author distribution
pdf(file.path("plot", "pca.author.topics.pdf"), width = 6, height = 6)
plot(pca.topics$x,
  col = unclass(samples$auth),
  pch = unclass(samples$auth),
  main = "LDA 50 topics\nby author"
)
dev.off()

#
# 3. Check k-means correlation with authorhsip
#

#  a) tf-idf
print(system.time(
  auth_test.tfidf.cl <- cluster.series(feat.tfidf, k = nlevels(samples$auth), nreps = 10, cores = ncores)
))
auth_test.tfidf <- apply(auth_test.tfidf.cl, 2, adjustedRandIndex, y=samples$auth)

#  b) author-adjusted tf-idf
print(system.time(
  auth_test.adjusted.cl <- cluster.series(feat.adjusted, k = nlevels(samples$auth), nreps = 10, cores = ncores)
))
auth_test.adjusted <- apply(auth_test.adjusted.cl, 2, adjustedRandIndex, y=samples$auth)


#
# 4. Check k-means stability for varying k
#

k.max <- 15
nreps <- 10

#  a) tfidf

cat("Checking stability for tfidf with k = 2 to", k.max, "\n")

kmcl.tfidf.cl <- do.call(cbind, lapply(2:k.max, function(k) {
  cat("k =", k, ": Generating", nreps, "clusters\n")

  t0 <- Sys.time()
  on.exit(cat("k =", k, ":", difftime(Sys.time(), t0, units = "min"),  "minutes\n"))

  cluster.series(feat.tfidf, k, nreps, ncores)
})
)

kmcl.tfidf.k.orig <- rep(2:k.max, each=nreps)

randscores.tfidf.k.orig <- lapply(2:k.max, function(k) {
  combn(which(kmcl.tfidf.k.orig == k), 2, function(i) {
    adjustedRandIndex(kmcl.tfidf.cl[,i[1]], kmcl.tfidf.cl[, i[2]])
  })
})
names(randscores.tfidf.k.orig) <- 2:k.max

randscores.tfidf.k.orig.nobs <- rep(choose(nreps, 2), k.max - 1)

pdf(file.path("plot", "kmeans-authorship.tfidf.pdf"), height = 6, width = 6)
boxplot(randscores.tfidf.k.orig, 
  main = "TF-IDF: k-means stability\n",
  xlab = "number of classes",
  ylab = "Adjusted Rand index",
  cex.axis = 0.9
)
dev.off()


kmcl.tfidf.k.effective <- apply(kmcl.tfidf.cl, 2, function(cl) {
  length(which(table(cl) > 3))
})

randscores.tfidf.k.effective <- lapply(sort(unique(kmcl.tfidf.k.effective)), function(k) {
  combn(which(kmcl.tfidf.k.effective == k), 2, function(i) { 
    adjustedRandIndex(kmcl.tfidf.cl[,i[1]], kmcl.tfidf.cl[, i[2]])
  })
})
names(randscores.tfidf.k.effective) <- sort(unique(kmcl.tfidf.k.effective))

randscores.tfidf.k.effective.nobs <- sapply(table(kmcl.tfidf.k.effective), choose, k = 2)

# boxplot(randscores.tfidf.k.effective, main = "tfidf\neffective k")


#  b) author-adjusted tfidf

cat("Checking stability for adjusted with k = 2 to", k.max, "\n")

kmcl.adjusted.cl <- do.call(cbind, lapply(2:k.max, function(k) {
  cat("k =", k, ": Generating", nreps, "clusters\n")

  t0 <- Sys.time()
  on.exit(cat("k =", k, ":", difftime(Sys.time(), t0, units = "min"),  "minutes\n"))

  cluster.series(feat.adjusted, k, nreps, ncores)
})
)

kmcl.adjusted.k.orig <- rep(2:k.max, each=nreps)

randscores.adjusted.k.orig <- lapply(2:k.max, function(k) {
  combn(which(kmcl.adjusted.k.orig == k), 2, function(i) {
    adjustedRandIndex(kmcl.adjusted.cl[,i[1]], kmcl.adjusted.cl[, i[2]])
  })
})
names(randscores.adjusted.k.orig) <- 2:k.max

randscores.adjusted.k.orig.nobs <- rep(choose(nreps, 2), k.max - 1)

pdf(file.path("plot", "kmeans-authorship.adjusted.pdf"), height = 6, width = 6)
boxplot(randscores.adjusted.k.orig, 
  main = "TF-IDF: author-adjusted\nk-means stability\n",
  xlab = "number of classes",
  ylab = "Adjusted Rand index",
  cex.axis = 0.9
)
dev.off()

kmcl.adjusted.k.effective <- apply(kmcl.adjusted.cl, 2, function(cl) {
  length(which(table(cl) > 9))
})

randscores.adjusted.k.effective <- lapply(sort(unique(kmcl.adjusted.k.effective)), function(k) {
  combn(which(kmcl.adjusted.k.effective == k), 2, function(i) { 
    adjustedRandIndex(kmcl.adjusted.cl[,i[1]], kmcl.adjusted.cl[, i[2]])
  })
})
names(randscores.adjusted.k.effective) <- sort(unique(kmcl.adjusted.k.effective))

randscores.adjusted.k.effective.nobs <- sapply(table(kmcl.adjusted.k.effective), choose, k = 2)

# boxplot(randscores.adjusted.k.effective, main = "tfidf-adjusted\neffective k")


#  c) LDA

cat("Checking stability for adjusted with k = 2 to", k.max, "\n")

kmcl.topics.cl <- do.call(cbind, lapply(2:k.max, function(k) {
  cat("k =", k, ": Generating", nreps, "clusters\n")

  t0 <- Sys.time()
  on.exit(cat("k =", k, ":", difftime(Sys.time(), t0, units = "min"),  "minutes\n"))

  cluster.series(feat.topics, k, nreps, ncores)
})
)

kmcl.topics.k.orig <- rep(2:k.max, each=nreps)

randscores.topics.k.orig <- lapply(2:k.max, function(k) {
  combn(which(kmcl.topics.k.orig == k), 2, function(i) {
    adjustedRandIndex(kmcl.topics.cl[,i[1]], kmcl.topics.cl[, i[2]])
  })
})
names(randscores.topics.k.orig) <- 2:k.max

randscores.topics.k.orig.nobs <- rep(choose(nreps, 2), k.max - 1)

pdf(file.path("plot", "kmeans-authorship.topics.pdf"), height = 6, width = 6)
boxplot(randscores.topics.k.orig, 
  main = "LDA: 50 topics\nk-means stability\n",
  xlab = "number of classes",
  ylab = "Adjusted Rand index",
  cex.axis = 0.9
)
dev.off()

kmcl.topics.k.effective <- apply(kmcl.topics.cl, 2, function(cl) {
  length(which(table(cl) > 9))
})

randscores.topics.k.effective <- lapply(sort(unique(kmcl.topics.k.effective)), function(k) {
  combn(which(kmcl.topics.k.effective == k), 2, function(i) { 
    adjustedRandIndex(kmcl.topics.cl[,i[1]], kmcl.topics.cl[, i[2]])
  })
})
names(randscores.topics.k.effective) <- sort(unique(kmcl.topics.k.effective))

randscores.topics.k.effective.nobs <- sapply(table(kmcl.topics.k.effective), choose, k = 2)



#
# 5. Aeneid 7 graph
#

mask <- samples$auth == "vergil" & (
  substr(samples$lstart, 1, 1) == "7" | substr(samples$lend, 1, 1) == "7")

# create collated graphs

# a. tfidf

for (i in 21:30) {  
  pdf(
    file = file.path("plot", paste("aen7", "tfidf", i, "pdf", sep=".")),
    width = 6,
    height = 4
  )
  
  plot(
    kmcl.tfidf.cl[mask, i],
    col = kmcl.tfidf.cl[mask, i],
    pch = 15,
    cex = 2,
    main = paste("Vergil Aeneid 7\ntfidf [", i, "]", sep=""),
    ylab = "class",
    xlab = "first verse of sample",
    xaxt = "n",
    yaxt = "n",
    ylim = c(1,4)
  )
  mtext(
    text = samples[mask, lstart], 
    at = seq(1:(length(which(mask)))),
    side = 1,
    line = 0.5,
    cex = 0.8,
    las = 2
  )
  dev.off()
}


# b. adjusted

aen.remap <- read.table(
  file = file.path("data", "benchmark", "aen7-30line-remap.txt"),
  head = F,
  sep = "\t"
)

for (i in 21:30) {  
  pdf(
    file = file.path("plot", paste("aen7", "adjusted", i, "pdf", sep=".")),
    width = 6,
    height = 4
  )
  
  oldclass <- kmcl.adjusted.cl[mask, i]
  newclass <- oldclass
  
  for (j in 1:4) {
    newclass[oldclass==j] <- aen.remap[j, i-20]
  }
  
  plot(
    newclass,
    col = newclass,
    pch = 15,
    cex = 2,
    main = paste("Vergil Aeneid 7\ntfidf-adjusted [", i, "]", sep=""),
    ylab = "class",
    xlab = "first verse of sample",
    xaxt = "n",
    yaxt = "n",
    ylim = c(1,4)
  )
  mtext(
    text = samples[mask, lstart], 
    at = seq(1:(length(which(mask)))),
    side = 1,
    line = 0.5,
    cex = 0.8,
    las = 2
  )
  dev.off()
}

# c. LDA

for (i in 21:30) {  
  pdf(
    file = file.path("plot", paste("aen7", "topics", i, "pdf", sep=".")),
    width = 6,
    height = 4
  )
  
  plot(
    kmcl.topics.cl[mask, i],
    col = kmcl.topics.cl[mask, i],
    pch = 15,
    cex = 2,
    main = paste("Vergil Aeneid 7\nlda ", ntopics, " topics [", i, "]", sep=""),
    ylab = "class",
    xlab = "first verse of sample",
    xaxt = "n",
    yaxt = "n",
    ylim = c(1,4)
  )
  mtext(
    text = samples[mask, lstart], 
    at = seq(1:(length(which(mask)))),
    side = 1,
    line = 0.5,
    cex = 0.8,
    las = 2
  )
  dev.off()
}

#
# 6. Mynors' scenes
#

mynors <- load.scenes(
  file = file.path("data", "benchmark", "scenes-aen7-3col.txt"),
  author = "vergil"
)
mynors.tfidf <- as.matrix(DocumentTermMatrix(VCorpus(VectorSource(mynors$Text))))

mynors.pca.tfidf <- predict(pca.tfidf, newdata = mynors.tfidf)


cat ("Trying to remove author signal\n")
feat.adjusted <- feat.tfidf
for (name in levels(samples$auth)) {
  cat(" -", name, "\n")
  sig.auth <- colMeans(feat.tfidf[samples$auth == name,]) - sig.base
  feat.adjusted[samples$auth == name,] <- t(t(feat.tfidf[samples$auth == name,]) - sig.auth)
}
