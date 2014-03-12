sdDiffTrim <- function(x, na.rm=FALSE, diff=1L, trim=0.001, scale=TRUE, ...) {
    if (na.rm) 
        x <- x[!is.na(x)]
    if (scale)
        x <- x / mean(x, na.rm=TRUE)
    if (diff > 0L) 
        x <- diff(x, differences=diff)
    n <- length(x)
    if (trim > 0 && n) {
        lo <- floor(n * trim) + 1
        hi <- n + 1 - lo
        x <- sort.int(x, partial=unique(c(lo, hi)))[lo:hi]
    }
    sd <- sd(x, na.rm=FALSE)
    x <- NULL
    sd / (sqrt(2)^diff)
}

# EOF
