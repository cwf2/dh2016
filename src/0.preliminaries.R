# Thematic work on Valerius Flaccus

# packages

library(data.table)
library(parallel)
library(tm)
library(topicmodels)
library(mclust)
library(XML)
library(MASS)

source(file.path("src", "tesserae.R"))

# functions

load.file <- function(file) {
  # read a single text and produce a data.table

  cat("Reading", file, "... ")

  # parse the Tesserae XML file
  doc <- xmlParseDoc(file, encoding="UTF-8")

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

      ll = xmlGetAttr(node, "loc")
      vv <- xmlValue(node)
      if (is.na(vv)) { cat("NA at", ll, "\n")}

      data.table(
        auth = auth,
        work = work,
        loc = xmlGetAttr(node, "loc"),
        verse = xmlValue(node)
      )
    })
  )[,
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

make.samples <- function(sample.size = 50, offset = 0) {
  la.verse[,
    .(auth, work, loc, verse, int.grp = ceiling((unitid - offset)/ sample.size))
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


# 1. preprocessing the texts

# load latin corpus
la.verse <- load.corpus(file.path("data", "la.index.txt"))

# load latin stemmer
la.stemmer <- build.stemmer.table(
  stems.file = file.path("data", "tesserae", "la.lexicon.csv"),
  resolve = file.path("data", "tesserae", "la.word.freq")
)
