#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname hello2
#' @export
hello2 <- function(n) {
  .C("hello", as.integer(n))
  #print(paste0("***",n,"***"))
}
