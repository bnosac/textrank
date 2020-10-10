
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
#' This data.frame can be used as input in the \code{\link{textrank_sentences}} algorithm.
#' @export
#' @seealso \code{\link{textrank_sentences}}
#' @examples
#' \dontshow{if(require(udpipe) & require(textreuse))
#' \{
#' }
#' library(textreuse)
#' library(udpipe)
#' lsh_probability(h = 1000, b = 500, s = 0.1) # A 10 percent Jaccard overlap will be detected well
#'
#' minhash <- minhash_generator(n = 1000, seed = 123456789)
#'
#' data(joboffer)
#' joboffer$textrank_id <- unique_identifier(joboffer, c("doc_id", "paragraph_id", "sentence_id"))
#' sentences <- unique(joboffer[, c("textrank_id", "sentence")])
#' terminology <- subset(joboffer, upos %in% c("NOUN", "ADJ"), select = c("textrank_id", "lemma"))
#' candidates <- textrank_candidates_lsh(x = terminology$lemma, sentence_id = terminology$textrank_id,
#'                                       minhashFUN = minhash, bands = 500)
#' head(candidates)
#' tr <- textrank_sentences(data = sentences, terminology = terminology,
#'                          textrank_candidates = candidates)
#' summary(tr, n = 2)
#' \dontshow{
#' \}
#' # End of main if statement running only if the required packages are installed
#' }
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
#' This data.frame can be used as input in the \code{\link{textrank_sentences}} algorithm.
#' @export
#' @seealso \code{\link{textrank_sentences}}
#' @examples
#' \dontshow{if(require(udpipe))
#' \{
#' }
#' library(udpipe)
#' data(joboffer)
#' joboffer$textrank_id <- unique_identifier(joboffer, c("doc_id", "paragraph_id", "sentence_id"))
#' candidates <- textrank_candidates_all(unique(joboffer$textrank_id))
#' head(candidates, 50)
#' \dontshow{
#' \}
#' # End of main if statement running only if the required packages are installed
#' }
textrank_candidates_all <- function(x){
  x <- unique(x)
  x <- setdiff(x, NA)
  x_length <- length(x)
  stopifnot(x_length > 1)
  if(x_length < 200){
    candidates <- utils::combn(x = x, m = 2, simplify = FALSE)
    candidates <- lapply(candidates, FUN=function(x){
      list(textrank_id_1 = x[1],
           textrank_id_2 = x[2])
    })
  }else{
    candidates <- lapply(seq(x)[-x_length], function(i){
      data.table::data.table(textrank_id_1 = x[i], textrank_id_2 = x[(i+1L):x_length])
    })
  }
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
#' @return an object of class textrank_sentences
#' which is a list with elements:
#' \itemize{
#' \item sentences: a data.frame with columns textrank_id, sentence and textrank where the textrank is the Google Pagerank importance metric of the sentence
#' \item sentences_dist: a data.frame with columns textrank_id_1, textrank_id_2 (the sentence id) and weight which
#' is the result of the computed distance between the 2 sentences
#' \item pagerank: the result of a call to \code{\link[igraph]{page_rank}}
#' }
#' @export
#' @examples
#' \dontshow{if(require(udpipe))
#' \{
#' }
#' library(udpipe)
#' data(joboffer)
#' head(joboffer)
#' joboffer$textrank_id <- unique_identifier(joboffer, c("doc_id", "paragraph_id", "sentence_id"))
#' sentences <- unique(joboffer[, c("textrank_id", "sentence")])
#' cat(sentences$sentence)
#' terminology <- subset(joboffer, upos %in% c("NOUN", "ADJ"), select = c("textrank_id", "lemma"))
#' head(terminology)
#'
#' ## Textrank for finding the most relevant sentences
#' tr <- textrank_sentences(data = sentences, terminology = terminology)
#' summary(tr, n = 2)
#' summary(tr, n = 5, keep.sentence.order = TRUE)
#' \dontshow{
#' \}
#' # End of main if statement running only if the required packages are installed
#' }
#' \dontrun{
#' ## Using minhash to reduce sentence combinations - relevant if you have a lot of sentences
#' library(textreuse)
#' minhash <- minhash_generator(n = 1000, seed = 123456789)
#' candidates <- textrank_candidates_lsh(x = terminology$lemma, sentence_id = terminology$textrank_id,
#'                                       minhashFUN = minhash, bands = 500)
#' tr <- textrank_sentences(data = sentences, terminology = terminology,
#'                          textrank_candidates = candidates)
#' summary(tr, n = 2)
#' }
#' ## You can also reduce the number of sentence combinations by sampling
#' tr <- textrank_sentences(data = sentences, terminology = terminology, max = 100)
#' tr
#' summary(tr, n = 2)
textrank_sentences <- function(data, terminology,
                     textrank_dist = textrank_jaccard,
                     textrank_candidates = textrank_candidates_all(data$textrank_id),
                     max = 1000,
                     options_pagerank = list(directed = FALSE),
                     ...){
  textrank_id <- NULL

  stopifnot(sum(duplicated(data[, 1])) == 0)
  data <- as.data.frame(data)
  stopifnot(nrow(data) > 1)
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
  class(result) <- "textrank_sentences"
  result
}



#' @title Extract the most important sentences which were identified with textrank_sentences
#' @description Extract the most important sentences which were identified by \code{\link{textrank_sentences}}
#' @param object an object of class textrank_sentences
#' @param n integer indicating to extract only the top n sentences
#' @param keep.sentence.order logical indicating to keep the sentence order as provided
#' in the original \code{data} argument of the \code{\link{textrank_sentences}} function
#' or to order it by the pagerank score. Defaults to FALSE indicating to order by pagerank score.
#' @param ... not used
#' @return a character vector with the top \code{n} most important sentences
#' which were identified by \code{\link{textrank_sentences}}
#' @export
#' @seealso \code{\link{textrank_sentences}}
summary.textrank_sentences <- function(object, n = 3, keep.sentence.order = FALSE, ...){
  pr <- sort(object$pagerank$vector, decreasing = TRUE)
  topsent <- utils::head(names(pr), n)
  out <- object$sentences[object$sentences$textrank_id %in% topsent, ]
  if(!keep.sentence.order){
    out <- out[order(factor(out$textrank_id, levels = topsent), decreasing = FALSE), ]
  }
  out$sentence
}


#' @export
print.textrank_sentences <- function(x, ...){
  cat("Textrank on sentences, showing top 5 most important sentences found:", sep = "\n")
  txt <- summary(x, n = 5)
  txt <- sprintf("  %s. %s", seq_along(txt), txt)
  cat(txt, sep = "\n")
}



#' @title Textrank - extract relevant keywords
#' @description The textrank algorithm allows to find relevant keywords in text.
#' Where keywords are a combination of words following each other. \cr
#'
#' In order to find relevant keywords, the textrank algorithm constructs a word network. This
#' network is constructed by looking which words follow one another.
#' A link is set up between two words if they follow one another, the link gets a higher weight if these 2 words occur
#' more frequenctly next to each other in the text.\cr
#' On top of the resulting network the 'Pagerank' algorithm is applied to get the importance of each word.
#' The top 1/3 of all these words are kept and are considered relevant. After this, a keywords table is constructed
#' by combining the relevant words together if they appear following one another in the text.
#' @param x a character vector of words.
#' @param relevant a logical vector indicating if the word is relevant or not. In the standard textrank
#' algorithm, this is normally done by doing a Parts of Speech tagging and selecting which of the words are
#' nouns and adjectives.
#' @param p percentage (between 0 and 1) of relevant words to keep. Defaults to 1/3.
#' Can also be an integer which than indicates how many words to keep. Specify +Inf if you want to keep all words.
#' @param ngram_max integer indicating to limit keywords which combine \code{ngram_max} combinations of words which follow one another
#' @param sep character string with the separator to \code{paste} the subsequent relevant words together
#' @return an object of class textrank_keywords
#' which is a list with elements:
#' \itemize{
#' \item terms: a character vector of words from the word network with the highest pagerank
#' \item pagerank: the result of a call to \code{\link[igraph]{page_rank}} on the word network
#' \item keywords: the data.frame with keywords containing columns keyword, ngram, freq indicating the keywords found and the frequency of occurrence
#' \item keywords_by_ngram: data.frame with columns keyword, ngram, freq indicating the keywords found and the frequency of occurrence
#' at each level of ngram. The difference with keywords being that if you have a sequence of words e.g. data science consultant, then in the keywords_by_ngram
#' you would still have the keywords data analysis and science consultant, while in the keywords list element you would only have data science consultant
#' }
#' @export
#' @seealso \code{\link[igraph]{page_rank}}
#' @examples
#' data(joboffer)
#' keywords <- textrank_keywords(joboffer$lemma,
#'                               relevant = joboffer$upos %in% c("NOUN", "VERB", "ADJ"))
#' subset(keywords$keywords, ngram > 1 & freq > 1)
#' keywords <- textrank_keywords(joboffer$lemma,
#'                               relevant = joboffer$upos %in% c("NOUN"),
#'                               p = 1/2, sep = " ")
#' subset(keywords$keywords, ngram > 1)
#'
#' ## plotting pagerank to see the relevance of each word
#' barplot(sort(keywords$pagerank$vector), horiz = TRUE,
#'         las = 2, cex.names = 0.5, col = "lightblue", xlab = "Pagerank")
textrank_keywords <- function(x, relevant=rep(TRUE, length(x)), p = 1/3, ngram_max = 5, sep = "-"){
  stopifnot(is.logical(relevant))
  stopifnot(is.character(x))
  stopifnot(length(x) == length(relevant))
  stopifnot(length(p) == 1)
  stopifnot(ngram_max > 1)
  keyword <- freq <- ngram <- NULL

  next_word <- function(x, n){
    data.table::shift(x, n = as.integer(n), type = "lead")
  }
  term_adjacency <- function (x, relevant, order = TRUE, ...) {
    cooc <- weight <- term1 <- term2 <- NULL
    result <- data.table(term1 = x,
                         term2 = next_word(x, n = 1), cooc = 1L)
    result <- result[!is.na(term1) & !is.na(term2) &
                       relevant %in% TRUE & next_word(relevant, n = 1) %in% TRUE, ]
    result <- result[, list(weight = sum(cooc)), by = list(term1, term2)]
    if (order) {
      result <- result[order(weight, decreasing = TRUE), ]
    }
    result
  }

  ## Identify word network by looking at which words are followed by one another
  data <- term_adjacency(x, relevant)

  ## On that network apply Google pagerank
  wordgraph <- igraph::graph_from_data_frame(data, directed = FALSE)
  pr <- igraph::page_rank(graph = wordgraph)

  ## Keep by default 1/3 of the words in the network which have the highest pagerank
  keep_nr <- igraph::vcount(wordgraph)
  if(p == +Inf){
  }else if(p <= 1){
    keep_nr <- keep_nr * p
  }else{
    keep_nr <- min(keep_nr, p)
  }
  keep_nr <- ceiling(keep_nr)
  keywords <- sort(pr$vector, decreasing = TRUE)
  keywords <- head(keywords, keep_nr)
  keywords <- names(keywords)

  ## extract keyword combinations: keywords which are followed by another keyword
  ## this is done for each ngram combination meaning we keep the ngram of the lower n even if the higher n is there
  output_per_ngram <- list()
  keywordcombinations <- data.table(keyword = ifelse(x %in% keywords, x, NA_character_), ngram = 1L)
  stop_already <- is.na(keywordcombinations$keyword)
  output_per_ngram[[1]] <- keywordcombinations
  for(i in 2:ngram_max){
    nextterm <- next_word(x, n = i-1)
    nextterm <- ifelse(nextterm %in% keywords, nextterm, NA_character_)
    stop_already[which(is.na(nextterm))] <- TRUE
    if(sum(stop_already) == length(stop_already)){
      break
    }
    keywordcombinations$keyword <- ifelse(stop_already, keywordcombinations$keyword, paste(keywordcombinations$keyword, nextterm, sep = sep))
    keywordcombinations$ngram <- ifelse(stop_already, keywordcombinations$ngram, keywordcombinations$ngram + 1L)
    output_per_ngram[[i]] <- keywordcombinations[keywordcombinations$ngram == i, ]
  }
  output_per_ngram <- lapply(output_per_ngram, FUN=function(x){
    x <- x[!is.na(keyword), list(freq = .N), by = list(keyword, ngram)]
    x <- x[order(freq, decreasing=TRUE), ]
    x
  })
  output_per_ngram <- data.table::rbindlist(output_per_ngram)
  output_per_ngram <- setDF(output_per_ngram)

  output_keywords <- keywordcombinations[!is.na(keywordcombinations$keyword), ]
  output_keywords <- as.data.table(output_keywords)
  output_keywords <- output_keywords[, list(freq = .N), by = list(keyword, ngram)]
  output_keywords <- output_keywords[order(freq, decreasing = TRUE), ]
  output_keywords <- setDF(output_keywords)

  result <- list(terms = keywords,
                 pagerank = pr,
                 keywords = output_keywords,
                 keywords_by_ngram = output_per_ngram)
  class(result) <- "textrank_keywords"
  result
}

