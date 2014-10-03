log2offset <- function(offset=.Machine$double.xmin) {
    # Get offset?
    if (missing(offset)) {
        offset <- getOption("QDNAseq::log2offset", .Machine$double.xmin)
        offset <- as.double(offset);
        stopifnot(is.finite(offset));
        return(offset);
    }
    
    # Reset offset?
    if (is.null(offset)) offset <- eval(formals(log2offset)$offset);
    
    # Set offset
    stopifnot(length(offset) == 1L);
    offset <- as.double(offset);
    stopifnot(is.finite(offset));
    options("QDNAseq::log2offset"=offset);
    
    offset;
}

log2adhoc <- function(x, offset=log2offset(), inv=FALSE) {
    if (!inv) {
        x[x < 0] <- 0
        x <- x + offset
        log2(x)
    } else {
        x <- 2^x
        x - offset
    }
}

unlog2adhoc <- function(y, offset=log2offset()) {
    y <- 2^y
    y - offset
}

sqrtadhoc <- function(x, factor=3/8, offset=0, inv=FALSE) {
    if (!inv) {
        x <- x + offset
        sqrt(x * factor)
    } else {
        x <- x^2 * (factor^-1)
        x - offset
    }  
}

# EOF
