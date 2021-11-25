suppressVerbose <- function(expr, envir = parent.frame(), suppress = TRUE) {
  if (is.na(suppress)) suppress <- FALSE
  expr <- substitute(expr)
  if (suppress) {
    expr <- bquote({
      capture.output(suppressMessages({
        res <- withVisible(.(expr))
      }), type="output", split=FALSE)
      res
    })
  } else {
    expr <- bquote(withVisible(.(expr)))
  }
  res <- eval(expr, envir=envir, enclos=baseenv())
  if (res$visible) res$value else invisible(res$value)
}


#' @importFrom utils packageVersion
future_version <- local({
  ver <- NULL
  function() {
    if (is.null(ver)) ver <<- packageVersion("future")
    ver
  }
})

assert_future_version <- function() {
  if (future_version() >= "1.22.1") return()
  stop(sprintf("This function requires future (>= 1.22.1). Please update: %s",
               future_version()))
}
