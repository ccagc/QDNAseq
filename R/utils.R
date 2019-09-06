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
