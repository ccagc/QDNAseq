sdDiffTrim <- function(x, ..., trim=0.001, scale=TRUE) {
    if (scale)
        x <- x / mean(x, na.rm=TRUE)
    sdDiff(x, ..., trim=trim)
}

expectedVariance <- function(object) {
    expected.variance <- sum(binsToUse(object)) / object$used.reads
    if ("paired.ends" %in% colnames(pData(object))) {
        multiplier <- ifelse(object$paired.ends, 2, 1)
        expected.variance <- expected.variance * multiplier
    }
    expected.variance
}

# EOF
