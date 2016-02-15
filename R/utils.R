## Until the future package provides a well-defined *apply() API
## we define futurized versions of lapply() and apply() here.
flapply <- function(X, FUN, ..., seeds=NULL) {
    if (!is.null(seeds)) {
        if (length(seeds) < length(X))
            seeds <- rep_len(seeds, length.out=length(X))
    } else {
        seeds <- rep.int(NA_integer_, times=length(X))
    }
    if (length(X) == 1 || !"future" %in% loadedNamespaces()) {
        res <- list()
        for (ii in seq_along(X)) {
            if (!is.na(seeds[ii]))
                set.seed(seeds[ii])
            res[[ii]] <- FUN(X[[ii]], ...)
        }
        names(res) <- names(X)
        return(res)
    }
    res <- list()
    for (ii in seq_along(X)) {
        X_ii <- X[[ii]]
        seed_ii <- seeds[ii]
        res[[ii]] <- future::future({
            if (!is.na(seed_ii))
                set.seed(seed_ii)
            FUN(X_ii, ...)
        })
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
