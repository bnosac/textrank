

textrank_dist <- function(termsa, termsb, type = c("jaccard", "jaccard_bag_of_words")){
  type <- match.arg(type)
  if(type == "jaccard"){
    ## how many % of all terms are overlapping
    length(intersect(termsa, termsb)) / length(union(termsa, termsb))
  }else if(type == "jaccard_bag_of_words"){
    overlap <- intersect(termsa, termsb)
    overlap_minfreq_byterm <- sapply(overlap, function(x) min(sum(x == termsa), sum(x == termsb)), USE.NAMES = FALSE)
    sum(overlap_minfreq_byterm) / log(length(termsb)) + log(length(termsa))
  }
}

#' @title Textrank - extract relevant sentences and keywords
#' @description The textrank algorithm is a technique to rank sentences in order of importance and can also be used to
#' identify the most relevant keywords or key phrases in text.\cr
#'
#' In order to find relevant sentences, the textrank algorithm needs 2 inputs:
#' a data.frame (\code{data}) with sentences and a data.frame (\code{terminology})
#' containing tokens which are part of each sentence.\cr
#' Based on these 2 datasets, it calculates the pairwise distance between each sentence by computing
#' how many terms are overlapping (Jaccard distance). These pairwise distances among the sentences are next passed on to
#' Google's pagerank algorithm to identify the most relevant sentences.\cr\cr
#'
#' For extracting keywords the flow is similar, but is still under construction
#' @param data a data.frame with 1 row per sentence where the first column
#' is an identifier of a sentence (e.g. textrank_id) and the second column is the raw sentence. See the example.
#' @param terminology a data.frame with tokens which are part of the sentence. Where the first column
#' is the identifier which corresponds to the first column of \code{data} and the second column indicates
#' the token which is part of the sentence. See the example.
#' @param method a function which calculates the distance between 2 vectors of tokens. The first 2 arguments of the function
#' are the tokens in sentence1 and sentence2. The function should return a numeric value of length one. The larger the value,
#' the larger the connection between the 2 vectors indicating more strength. Defaults to the jaccard distance,
#' indicating the percent of common tokens. For a large number of sentences, you might be interested in looking at the minhash algorithm
#' which is part of the textreuse package.
#' @param max integer indicating to reduce the number of sentence to sentence combinations to compute.
#' @param options_pagerank a list of arguments passed on to \code{\link[igraph]{page_rank}}
#' @param ... arguments passed on to \code{method}
#' @return an object of class textrank
#' which is a list with elements:
#' \itemize{
#' \item sentences: a data.frame with columns textrank_id, sentence and textrank where the textrank is the Google Pagerank importance metric of the sentence
#' \item sentences_dist: a data.frame with columns textrank_id_1, textrank_id_2 (the sentence id) and weight which
#' is the result of the computed distance between the 2 sentences
#' \item pagerank: the result of a call to \code{\link[igraph]{page_rank}}
#' }
#' @seealso \code{\link[igraph]{page_rank}}
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
textrank <- function(data, terminology, method = textrank_dist, max = 1000, options_pagerank = list(directed = FALSE), ...){
  textrank_id <- NULL

  stopifnot(sum(duplicated(data[, 1])) == 0)
  data <- as.data.frame(data)
  data.table::setnames(data, old = colnames(data)[1:2], new = c("textrank_id", "sentence"))

  terminology <- as.data.table(terminology)
  data.table::setnames(terminology, old = colnames(terminology)[1:2], new = c("textrank_id", "term"))
  data.table::setkey(terminology, "textrank_id")

  ## Calculate pairwise distances along all sentence combinations
  sentence_dist <- function(x, distFUN, ...){
    data1 <- terminology[textrank_id == x[1], ]
    data2 <- terminology[textrank_id == x[2], ]
    if(nrow(data1) == 0 || nrow(data2) == 0){
      w <- 0
    }else{
      w <- distFUN(data1$term, data2$term, ...)
    }
    data.frame(
      textrank_id_1 = x[1],
      textrank_id_2 = x[2],
      weight = w,
      stringsAsFactors = FALSE)
  }
  sent2sent_distance <- utils::combn(x = data$textrank_id, m = 2, simplify = FALSE)
  if(length(sent2sent_distance) > max){
    sent2sent_distance <- sent2sent_distance[sample.int(n = length(sent2sent_distance), size = max)]
  }
  sent2sent_distance <- lapply(sent2sent_distance, FUN = sentence_dist, distFUN = method, ...)
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

#' @title Extract the most important sentences
#' @description Extract the most important sentences which were identified by \code{\link{textrank}}
#' @param object an object of class textrank
#' @param n integer indicating to extract only the top n sentences
#' @param keep.sentence.order logical indicating to keep the sentence order in the original \code{data} element of the \code{\link{textrank}} function
#' or to order it by the pagerank score. Defaults to FALSE - order by pagerank score
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
