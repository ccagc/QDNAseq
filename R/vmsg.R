vmsg <- function(...) {
    verbose <- getOption("QDNAseq::verbose", TRUE)
    if (is.na(verbose)) {
        args <- list(...)
        appendLF <- TRUE
        if ("appendLF" %in% names(args)) {
            appendLF <- args[["appendLF"]]
            args[["appendLF"]] <- NULL
        }
        cat(unlist(args), sep="")
        if (appendLF)
            cat("\n")
    } else if (verbose) {
        message(...)
    }
}

# EOF
