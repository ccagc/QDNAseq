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

sqrtfactor <- function(factor=3/8) {
    # Get factor?
    if (missing(factor)) {
        factor <- getOption("QDNAseq::sqrtfactor", 3/8)
        factor <- as.double(factor);
        stopifnot(is.finite(factor));
        return(factor);
    }

    # Reset factor?
    if (is.null(factor)) factor <- eval(formals(sqrtfactor)$factor);

    # Set factor
    stopifnot(length(factor) == 1L);
    factor <- as.double(factor);
    stopifnot(is.finite(factor));
    options("QDNAseq::sqrtfactor"=factor);

    factor;
}

sqrtoffset <- function(offset=0) {
    # Get offset?
    if (missing(offset)) {
        offset <- getOption("QDNAseq::sqrtoffset", 0)
        offset <- as.double(offset);
        stopifnot(is.finite(offset));
        return(offset);
    }

    # Reset offset?
    if (is.null(offset)) offset <- eval(formals(sqrtoffset)$offset);

    # Set offset
    stopifnot(length(offset) == 1L);
    offset <- as.double(offset);
    stopifnot(is.finite(offset));
    options("QDNAseq::sqrtoffset"=offset);

    offset;
}

sqrtadhoc <- function(x, factor=sqrtfactor(), offset=sqrtoffset(), inv=FALSE) {
    if (!inv) {
        x <- x + offset
        sqrt(x * factor)
    } else {
        x <- x^2 * (factor^-1)
        x - offset
    }
}
