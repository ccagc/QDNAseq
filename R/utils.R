## Until the future package provides a well-defined *apply() API
## we define futurized versions of lapply() and apply() here.
flapply <- function(x, FUN, ...) {
  res <- list()
  for (ii in seq_along(x)) res[[ii]] <- future(FUN(x[[ii]], ...))
  names(res) <- names(x)
  values(res)
}

fapply <- function(X, MARGIN, FUN, ...) {
  fFUN <- function(...) { future(FUN(...)) }
  res <- apply(X, MARGIN=1L, FUN=fFUN, ...)
  res <- values(res)
  sapply(res, FUN=I, simplify=TRUE)
}
