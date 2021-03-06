###
# Corpus Readers
###

#' Read in a corpus file.
#' 
#' Converts pre-processed document matrices stored in popular formats to stm
#' format.
#' 
#' This function provides a simple utility for converting other document
#' formats to our own.  Briefly- \code{dtm} takes as input a standard matrix
#' and converts to our format \code{ldac} takes a file path and reads in a
#' document in the sparse format popularized by David Blei's C code
#' implementation of lda.  \code{slam} converts from the
#' \code{simple_triplet_matrix} representation used by the \code{slam} package.
#' This is also the representation of corpora in the popular \code{tm} package
#' and should work in those cases.
#' 
#' \code{dtm} expects a matrix object where each row represents a document and
#' each column represents a word in the dictionary.
#' 
#' \code{ldac} expects a file name or path that contains a file in Blei's LDA-C
#' format. From his ReadMe: "The data is a file where each line is of the form:
#' 
#' [M] [term_1]:[count] [term_2]:[count] ...  [term_N]:[count]
#' 
#' where [M] is the number of unique terms in the document, and the [count]
#' associated with each term is how many times that term appeared in the
#' document.  Note that [term_1] is an integer which indexes the term; it is
#' not a string."
#' 
#' Because R indexes from one, the values of the term indices are incremented
#' by one on import.
#' 
#' \code{slam} expects a \code{\link[slam]{simple_triplet_matrix}} from that
#' package.
#' 
#' \code{Matrix} attempts to coerce the matrix to a
#' \code{\link[slam]{simple_triplet_matrix}} and convert using the
#' functionality built for the \code{slam} package.  This will work for most
#' applicable classes in the \code{Matrix} package such as \code{dgCMatrix}.
#' 
#' Finally the object \code{txtorgvocab} allows the user to easily read in a
#' vocab file generated by the software \code{txtorg}.  When working in English
#' it is straightforward to read in files created by txtorg.  However when
#' working in other languages, particularly Chinese and Arabic, there can often
#' be difficulty reading in the files using \code{\link{read.table}} or
#' \code{\link{read.csv}} This function should work well in those
#' circumstances.
#' 
#' @param corpus An input file or filepath to be processed
#' @param type The type of input file.  We offer several sources, see details.
#' @return \item{documents}{A documents object in our format} \item{vocab}{A
#' vocab object if information is available to construct one}
#' @seealso \code{\link{textProcessor}}, \code{\link{prepDocuments}}
#' @examples
#' 
#' \dontrun{
#' 
#' library(textir)
#' data(congress109)
#' out <- readCorpus(congress109Counts, type="Matrix")
#' documents <- out$documents
#' vocab <- out$vocab
#' }
#' @export
readCorpus <- function(corpus, type=c("dtm", "ldac", "slam", "Matrix","txtorgvocab")) {
  type <- match.arg(type)
  switch(type,
         dtm = read.dtm(corpus),
         ldac = read.ldac(corpus),
         slam = read.slam(corpus),
         txtorgvocab = read.txtorg.vocab(corpus),
         Matrix = read.slam(slam::as.simple_triplet_matrix(corpus)))
}

read.ldac <- function(filename) {
  #Read the .ldac format
  # Based on Jonathan Chang's  reader with addition of zero correction.
  d <- scan(filename, what = "character", sep = "\n")
  d <- chartr(":", " ", d)
  d <- strsplit(d, " ", fixed = TRUE)
  d <- lapply(d, function(x) matrix(as.integer(x[-1]), nrow = 2))
  mapply(function(x) rbind(x[1,]+1, x[2,]), d) #zero correction
}

read.dtm <- function(dtm) {
  #test for and adjust for mispecification
  if("simple_triplet_matrix" %in% class(dtm)) {
    warning("Please use the slam option.  dtm is for dense matrices.")
    read.slam(dtm)
  }
  #convert a standard document-term matrix to list format.
  dtm.mat <- as.matrix(dtm)
  vocab <- colnames(dtm)
  if(any(dtm.mat==0)) {
    #if the dtm is not sparse we have to use a slightly slower method
    #to avoid it coercing back to a matrix
    documents <- lapply(split(dtm.mat, row(dtm.mat)), function(y) {
      rbind(which(y > 0), as.integer(y[y > 0])) }) 
    names(documents) <- NULL #we overwrite the automatically generated labels to match other method
  } else {
    #the more usual sparse matrix case
    documents <- apply(dtm.mat, 1, function(y) {
      rbind(which(y > 0), as.integer(y[y > 0])) })
  }
  return(list(documents=documents, vocab=vocab))
}

read.slam <- function(corpus) {
  #convert a simple triplet matrix to list format.
  if(!inherits(corpus, "simple_triplet_matrix")) stop("corpus is not a simple triplet matrix")
  if ("TermDocumentMatrix" %in% class(corpus)) {
    non_empty_docs <- which(slam::col_sums(corpus) != 0)
    documents <- ijv.to.doc(corpus[,non_empty_docs]$j, corpus[,non_empty_docs]$i, corpus[,non_empty_docs]$v) 
    names(documents) <- corpus[,non_empty_docs]$dimnames$Docs
   } else {
    non_empty_docs <- which(slam::row_sums(corpus) != 0)
    documents <- ijv.to.doc(corpus[non_empty_docs,]$i, corpus[non_empty_docs,]$j, corpus[non_empty_docs,]$v) 
    names(documents) <- corpus[non_empty_docs,]$dimnames$Docs
  }
  vocab <- corpus$dimnames$Terms
  return(list(documents=documents,vocab=vocab))
}

read.txtorg.vocab <- function(filename) {
  return(readLines(filename, encoding="UTF-8",warn=FALSE))
}