#' Make a Word Cloud of enriched pathway terms
#'
#' Use the textmining and wordcloud packages to make a word cloud. Great for using
#' with enriched GO/GSEA term names
#'
#' @param text Character vector of terms to be processed for a word cloud
#' @param cols Colour palette. Can be RColourBrewer palette name, rgb, colour names, or hexadecimal. Default is "Dark2"
#' @param retTerms Boolean if term-document matrix dataframe should be returned. Default=F
#' @param minfreq Minimum frequency for a word to be included in the word cloud. Default=3
#' @param rmwords Character vector of words to be filtered out (in addition to the default english stop words and words shorter than 3 letters) - defaults to c("regulation","process","positive","negative","mediated", "cell","cellular", "protein")
#' @return A Word cloud showing most frequent words in input text, and if indicated, a data frame of word frequencies.
#' @export

PathWC<-function(text, cols="Dark2", retTerms=F, minfreq=3,
                           rmwords=c("regulation","process","positive","negative","mediated", "cell",
                                     "cellular", "protein")){
  text<-sapply(strsplit(text, "%", T),'[[',1)
  docs<-text
  suppressWarnings(docs<-tm::Corpus(tm::VectorSource(text)))
  toSpace<-tm::content_transformer(function(x, pattern) gsub(pattern, " ",x))
  suppressWarnings(docs<-tm::tm_map(docs, toSpace, "&"))#If we need to remove special characters
  suppressWarnings(docs<-tm::tm_map(docs, toSpace, "-"))
  suppressWarnings(docs<-tm::tm_map(docs, toSpace, ";"))
  suppressWarnings(docs<-tm::tm_map(docs, toSpace, ">"))
  suppressWarnings(docs<-tm::tm_map(docs, toSpace, "<"))
  #Remove words shorter than 3 letters
  rmShort<-tm::content_transformer(function(x) gsub('\\b\\w{1,2}\\b','',x))
  suppressWarnings(docs<-tm::tm_map(docs, rmShort))
  suppressWarnings(docs<-tm::tm_map(docs, tm::content_transformer(tolower))) #lowercase
  suppressWarnings(docs<-tm::tm_map(docs, tm::removeNumbers)) #remove numbers
  #Remove common stopwords
  suppressWarnings(docs<-tm::tm_map(docs, tm::removeWords, tm::stopwords("english")))
  if(!is.null(rmwords)){
    suppressWarnings(docs<-tm::tm_map(docs, tm::removeWords, rmwords))
  }
  suppressWarnings(docs<-tm::tm_map(docs, tm::removePunctuation))
  suppressWarnings(docs<-tm::tm_map(docs, tm::stripWhitespace)) #remove whitespace
  #Build term-document matrix
  dtm<-tm::TermDocumentMatrix(docs)
  m<-as.matrix(dtm)
  v<-sort(rowSums(m), decreasing=T)
  d<-data.frame(word=names(v), freq=v)
  #Generate the Word cloud
  set.seed(1234)
  wordcloud::wordcloud(words=d$word, freq=d$freq, min.freq=minfreq,
                       max.words=200, random.order=FALSE, rot.per=0.35,
                       colors=colPal(cols))
  if(isTRUE(retTerms)){
    return(d)
  }
}
