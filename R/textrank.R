
#' @title Calculate the distance between 2 vectors based on the Jaccard distance
#' @description The jaccard distance computes the percentage of terms in the 2 vectors which are overlapping.
#' @param termsa a character vector of words
#' @param termsb a character vector of words
#' @return The Jaccard distance distance between the 2 vectors
#' @export
#' @examples
#' sentencea <- c("I", "like", "champaign")
#' sentenceb <- c("I", "prefer", "choco")
#' textrank_jaccard(termsa = sentencea, termsb = sentenceb)
textrank_jaccard <- function(termsa, termsb){
  length(intersect(termsa, termsb)) / length(union(termsa, termsb))
}



#' @title Use locality-sensitive hashing to get combinations of sentences which contain words which are in the same minhash bucket
#' @description This functionality is usefull if there are a lot of sentences and most of the sentences have no overlapping
#' words in there. In order not to compute the jaccard distance among all possible combinations of sentences as is
#' done by using \code{\link{textrank_candidates_all}}, we can reduce the combinations of sentences by using the Minhash algorithm.
#' This function sets up the combinations of sentences which are in the same Minhash bucket.
#' @param x a character vector of words or terms
#' @param sentence_id a character vector of identifiers of sentences where the words/terms provided in \code{x} are part of the sentence.
#' The length of \code{sentence_id} should be the same length of \code{x}
#' @param minhashFUN a function which returns a minhash of a character vector. See the examples or look at \code{\link[textreuse]{minhash_generator}}
#' @param bands integer indicating to break down the minhashes in \code{bands} number of bands. Mark that
#' the number of minhash signatures should always be a multiple of the number of local sensitive hashing bands. See the example
#' @return a data.frame with 2 columns textrank_id_1 and textrank_id_2 containing identifiers of sentences \code{sentence_id}
#' which contained terms in the same minhash bucket.
#' This data.frame can be used as input in the \code{\link{textrank}} algorithm.
#' @export
#' @seealso \code{\link{textrank}}
#' @examples
#' library(textreuse)
#' lsh_probability(h = 1000, b = 500, s = 0.1) # A 10 percent Jaccard overlap will be detected well
#'
#' minhash <- minhash_generator(n = 1000, seed = 123456789)
#'
#' data(joboffer)
#' terminology <- subset(joboffer, upos %in% c("NOUN", "ADJ"), select = c("sentence_id", "lemma"))
#' candidates <- textrank_candidates_lsh(x = terminology$lemma, sentence_id = terminology$sentence_id,
#'                                       minhashFUN = minhash, bands = 500)
#' head(candidates)
textrank_candidates_lsh <- function(x, sentence_id, minhashFUN, bands){
  textrank_id_1 <- textrank_id_2 <- sentence_id.left <- sentence_id.right <- hash <- band <- NULL

  stopifnot(length(x) == length(sentence_id))

  ## each hash is put into a bucket
  examplehash <- minhashFUN("detect the n in minhashFUN")
  rows <- length(examplehash) / bands
  if(length(examplehash) %% rows != 0) {
    stop(sprintf("the number of hashes (%s) should be a multiple of bands (%s)", length(examplehash), bands))
  }
  hash_bands <- unlist(lapply(seq_len(bands), FUN=function(i) rep(i, times = rows)))

  sentence_to_bucket <- split(x, sentence_id)
  sentence_to_bucket <- mapply(sentence_id = names(sentence_to_bucket), sentence_to_bucket, FUN=function(sentence_id, words){
    buckets <- data.table(sentence_id = sentence_id,
                          hash = minhashFUN(words),
                          band = hash_bands)
    buckets <- buckets[, list(bucket = digest::digest(object = list(hashes = hash, b = band[1]))), by = list(sentence_id, band)]
    buckets <- buckets[, c("sentence_id", "bucket"), with = FALSE]
    buckets
  }, SIMPLIFY = FALSE)
  sentence_to_bucket <- data.table::rbindlist(sentence_to_bucket)

  candidates <- merge(sentence_to_bucket, sentence_to_bucket, by = "bucket", suffixes = c(".left", ".right"), all.x=TRUE, all.y=FALSE)
  candidates <- candidates[candidates$sentence_id.left != candidates$sentence_id.right, ]
  candidates <- candidates[, textrank_id_1 := ifelse(sentence_id.left < sentence_id.right, sentence_id.left, sentence_id.right)]
  candidates <- candidates[, textrank_id_2 := ifelse(sentence_id.left < sentence_id.right, sentence_id.right, sentence_id.left)]
  candidates <- unique(candidates[, c("textrank_id_1", "textrank_id_2")])
  candidates <- setDF(candidates)
  candidates
}


#' @title Get all combinations of sentences
#' @description Get all combinations of sentences
#' @param x a character vector of sentence identifiers
#' @return a data.frame with 2 columns textrank_id_1 and textrank_id_2 listing up all possible combinations of \code{x}.
#' The columns textrank_id_1 and textrank_id_2 contain identifiers of sentences given in \code{sentence_id}.
#' This data.frame can be used as input in the \code{\link{textrank}} algorithm.
#' @export
#' @seealso \code{\link{textrank}}
#' @examples
#' data(joboffer)
#' candidates <- textrank_candidates_all(unique(joboffer$sentence_id))
#' head(candidates, 50)
textrank_candidates_all <- function(x){
  x <- unique(x)
  x <- setdiff(x, NA)
  candidates <- utils::combn(x = x, m = 2, simplify = FALSE)
  candidates <- lapply(candidates, FUN=function(x){
    list(textrank_id_1 = x[1],
         textrank_id_2 = x[2])
  })
  candidates <- data.table::rbindlist(candidates)
  candidates <- setDF(candidates)
  candidates
}


#' @title Textrank - extract relevant sentences
#' @description The textrank algorithm is a technique to rank sentences in order of importance.\cr
#'
#' In order to find relevant sentences, the textrank algorithm needs 2 inputs:
#' a data.frame (\code{data}) with sentences and a data.frame (\code{terminology})
#' containing tokens which are part of each sentence.\cr
#' Based on these 2 datasets, it calculates the pairwise distance between each sentence by computing
#' how many terms are overlapping (Jaccard distance, implemented in \code{\link{textrank_jaccard}}).
#' These pairwise distances among the sentences are next passed on to Google's pagerank algorithm
#' to identify the most relevant sentences.\cr
#'
#' If \code{data} contains many sentences, it makes sense not to compute all pairwise sentence distances but instead limiting
#' the calculation of the Jaccard distance to only sentence combinations which are limited by the Minhash algorithm.
#' This is implemented in \code{\link{textrank_candidates_lsh}} and an example is show below.
#' @param data a data.frame with 1 row per sentence where the first column
#' is an identifier of a sentence (e.g. textrank_id) and the second column is the raw sentence. See the example.
#' @param terminology a data.frame with with one row per token indicating which token is part of each sentence.
#' The first column in this data.frame is the identifier which corresponds to the first column of \code{data}
#' and the second column indicates the token which is part of the sentence which will be passed on to \code{textrank_dist}.
#' See the example.
#' @param textrank_dist a function which calculates the distance between 2 sentences which are represented by a vectors of tokens.
#' The first 2 arguments of the function are the tokens in sentence1 and sentence2.
#' The function should return a numeric value of length one. The larger the value,
#' the larger the connection between the 2 vectors indicating more strength. Defaults to the jaccard distance (\code{\link{textrank_jaccard}}),
#' indicating the percent of common tokens.
#' @param textrank_candidates a data.frame of candidate sentence to sentence comparisons with columns textrank_id_1 and textrank_id_2
#' indicating for which combination of sentences we want to compute the Jaccard distance or the distance function as provided in \code{textrank_dist}.
#' See for example \code{\link{textrank_candidates_all}} or \code{\link{textrank_candidates_lsh}}.
#' @param max integer indicating to reduce the number of sentence to sentence combinations to compute.
#' In case provided, we take only this max amount of rows from \code{textrank_candidates}
#' @param options_pagerank a list of arguments passed on to \code{\link[igraph]{page_rank}}
#' @param ... arguments passed on to \code{textrank_dist}
#' @seealso \code{\link[igraph]{page_rank}}, \code{\link{textrank_candidates_all}}, \code{\link{textrank_candidates_lsh}}, \code{\link{textrank_jaccard}}
#' @return an object of class textrank
#' which is a list with elements:
#' \itemize{
#' \item sentences: a data.frame with columns textrank_id, sentence and textrank where the textrank is the Google Pagerank importance metric of the sentence
#' \item sentences_dist: a data.frame with columns textrank_id_1, textrank_id_2 (the sentence id) and weight which
#' is the result of the computed distance between the 2 sentences
#' \item pagerank: the result of a call to \code{\link[igraph]{page_rank}}
#' }
#' @export
#' @examples
#' data(joboffer)
#' head(joboffer)
#'
#' sentences <- unique(joboffer[, c("sentence_id", "sentence")])
#' cat(sentences$sentence)
#' terminology <- subset(joboffer, upos %in% c("NOUN", "ADJ"), select = c("sentence_id", "lemma"))
#' head(terminology)
#'
#' ## Textrank for finding the most relevant sentences
#' tr <- textrank(data = sentences, terminology = terminology)
#' summary(tr, n = 2)
#' summary(tr, n = 5, keep.sentence.order = TRUE)
#'
#' ## Using minhash to reduce sentence combinations - relevant if you have a lot of sentences
#' library(textreuse)
#' minhash <- minhash_generator(n = 1000, seed = 123456789)
#' candidates <- textrank_candidates_lsh(x = terminology$lemma, sentence_id = terminology$sentence_id,
#'                                       minhashFUN = minhash, bands = 500)
#' tr <- textrank(data = sentences, terminology = terminology, textrank_candidates = candidates)
#' summary(tr, n = 2)
#'
#' ## You can also reduce the number of sentence combinations by sampling
#' tr <- textrank(data = sentences, terminology = terminology, max = 100)
#' summary(tr, n = 2)
textrank <- function(data, terminology,
                     textrank_dist = textrank_jaccard,
                     textrank_candidates = textrank_candidates_all(data$textrank_id),
                     max = 1000,
                     options_pagerank = list(directed = FALSE),
                     ...){
  textrank_id <- NULL

  stopifnot(sum(duplicated(data[, 1])) == 0)
  data <- as.data.frame(data)
  data.table::setnames(data, old = colnames(data)[1:2], new = c("textrank_id", "sentence"))

  terminology <- as.data.table(terminology)
  data.table::setnames(terminology, old = colnames(terminology)[1:2], new = c("textrank_id", "term"))
  data.table::setkey(terminology, "textrank_id")

  ## Calculate pairwise distances along all sentence combinations
  sentence_dist <- function(id1, id2, distFUN, ...){
    data1 <- terminology[textrank_id == id1, ]
    data2 <- terminology[textrank_id == id2, ]
    if(nrow(data1) == 0 || nrow(data2) == 0){
      w <- 0
    }else{
      w <- distFUN(data1$term, data2$term, ...)
    }
    data.frame(
      textrank_id_1 = id1,
      textrank_id_2 = id2,
      weight = w,
      stringsAsFactors = FALSE)
  }
  sent2sent_distance <- as.data.frame(textrank_candidates)
  if(!missing(max)){
    max <- min(nrow(sent2sent_distance), max)
    sent2sent_distance <- sent2sent_distance[sample.int(n = nrow(sent2sent_distance), size = max), ]
  }
  sent2sent_distance <- mapply(id1 = sent2sent_distance$textrank_id_1,
                               id2 = sent2sent_distance$textrank_id_2, FUN = sentence_dist, MoreArgs = list(distFUN = textrank_dist, ...),
                               SIMPLIFY = FALSE)
  sent2sent_distance <- data.table::rbindlist(sent2sent_distance)

  ## Calculate pagerank
  pr <- igraph::graph_from_data_frame(sent2sent_distance, directed = FALSE)
  options_pagerank$graph <- pr
  pr <- do.call(igraph::page_rank, options_pagerank)

  ## Add pagerank to sentences for having it in the ouput
  data <- merge(data,
                data.frame(textrank_id = names(pr$vector), textrank = as.numeric(pr$vector), stringsAsFactors = FALSE),
                by = "textrank_id", all.x=TRUE, all.y=FALSE, sort = FALSE, suffixes = c("sent.", ""))

  result <- list(
    sentences = data,
    sentences_dist = sent2sent_distance,
    pagerank = pr)
  class(result) <- "textrank"
  result
}

#' @title Extract the most important sentences which were identified with textrank
#' @description Extract the most important sentences which were identified by \code{\link{textrank}}
#' @param object an object of class textrank
#' @param n integer indicating to extract only the top n sentences
#' @param keep.sentence.order logical indicating to keep the sentence order as provided
#' in the original \code{data} argument of the \code{\link{textrank}} function
#' or to order it by the pagerank score. Defaults to FALSE indicating to order by pagerank score.
#' @param ... not used
#' @return a character vector with the top \code{n} most important sentences
#' which were identified by \code{\link{textrank}}
#' @export
#' @seealso \code{\link{textrank}}
summary.textrank <- function(object, n = 3, keep.sentence.order = FALSE, ...){
  pr <- sort(object$pagerank$vector, decreasing = TRUE)
  topsent <- utils::head(names(pr), n)
  out <- object$sentences[object$sentences$textrank_id %in% topsent, ]
  if(!keep.sentence.order){
    out <- out[order(factor(out$textrank_id, levels = topsent), decreasing = FALSE), ]
  }
  out$sentence
}


