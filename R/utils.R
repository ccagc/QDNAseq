## Until the future package provides a well-defined *apply() API
## we define futurized versions of lapply() and apply() here.
flapply <- function(X, FUN, ...) {
    if (length(X) == 1 || !"future" %in% loadedNamespaces())
        return(lapply(X, FUN, ...))
    res <- list()
    for (ii in seq_along(X)) {
        X_ii <- X[[ii]]
        res[[ii]] <- future::future(FUN(X_ii, ...))
    }
    names(res) <- names(X)
    future::values(res)
}

fapply <- function(X, MARGIN, FUN, ...) {
    if (!"future" %in% loadedNamespaces())
        return(apply(X, MARGIN, FUN, ...))
    fFUN <- function(...) { future::future(FUN(...)) }
    res <- apply(X, MARGIN, FUN=fFUN, ...)
    res <- future::values(res)
    sapply(res, FUN=I, simplify=TRUE)
}

# EOF
