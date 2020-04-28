#' @include regenrichClasses.R
# @importFrom S4Vectors isEmpty
setMethod("isEmpty", signature = "Enrich", 
          definition = function(x){
            nrow(x@allResult) == 0
          })

setMethod("isEmpty", signature = "Score", 
          definition = function(x){
            nrow(x) == 0
          })

setMethod("isEmpty", signature = "TopNetwork", 
          definition = function(x){
            nrow(x@elementset) == 0
          })

#' dimention of `TopNetwork` object
#' @param x a `TopNetwork` object.
#' @return Dimention of regulator-target network edge table.
#' @export
#' @examples 
#' nw = newTopNetwork()
#' dim(nw)
setMethod("dim", signature = "TopNetwork", definition = function(x){
  dim(x@elementset)
})


head.Score = function(x, ...){
  head.default = utils::getFromNamespace("head.default", "utils")
  x1 = head.default(S3Part(x), ...)
  resCall = as.call(c(quote(newScore), x1))
  return(eval(resCall))
}

#' head or tail of Score object
#' 
#' @param x an \code{Score} object.
#' @param ... arguments to be passed to or from other methods.
#' @rdname headTailScore
#' @return Head or tail table of Score object.
#' @export
#' @examples
#' s = newScore(letters, seq(26), seq(26), seq(26), seq(2, 0, len = 26))
#' s1 = head(s)
#' s1
#' 
#' s2 = tail(s)
#' s2
#' 
setMethod("head", signature = signature("Score"), 
          definition = head.Score)


tail.Score = function(x, ...){
  tail.default = utils::getFromNamespace("tail.default", "utils")
  x1 = tail.default(S3Part(x), ...)
  resCall = as.call(c(quote(newScore), x1))
  return(eval(resCall))
}

#' @rdname headTailScore
#' @export
setMethod("tail", signature = signature("Score"), 
          definition = tail.Score)
