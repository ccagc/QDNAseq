vmsg <- function(...) {
  verbose <- getOption("QDNAseq::verbose", TRUE)
  if (verbose) message(...)
}

# EOF
