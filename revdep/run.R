#!/usr/bin/env Rscript
library("revdepcheck")
options(warn = 1L)

available_cores <- function() {
  getenv <- function(name) {
    as.integer(Sys.getenv(name, NA_character_))
  }
  getopt <- function(name) {
    as.integer(getOption(name, NA_integer_))
  }
  if (is.finite(n <- getopt("mc.cores") + 1L)) return(n)
  if (is.finite(n <- getopt("Ncpus") + 1L)) return(n)
  if (is.finite(n <- getenv("PBS_NUM_PPN"))) return(n)
  if (is.finite(n <- getenv("SLURM_CPUS_PER_TASK"))) return(n)
  if (is.finite(n <- getenv("NSLOTS"))) return(n)
  1L
}

precheck <- function() {
  ## WORKAROUND: Remove checked pkgs that use file links, which otherwise
  ## produce warnings which are promoted to errors by revdepcheck.
  unlink("revdep/checks/aroma.affymetrix", recursive = TRUE)
}

check <- function() {
  if (file_test("-f", p <- Sys.getenv("R_CHECK_ENVIRON", "~/.R/check.Renviron"))) {
    cat(sprintf("R CMD check will use env vars from %s\n", sQuote(p)))
    cat(sprintf("To disable, set 'R_CHECK_ENVIRON=false' (a fake pathname)\n"))
  }
  
  envs <- grep("^_R_CHECK_", names(Sys.getenv()), value = TRUE)
  if (length(envs) > 0L) {
    cat(sprintf("Detected _R_CHECK_* env vars that will affect R CMD check: %s\n",
                paste(sQuote(envs), collapse = ", ")))
  }

  precheck()
  revdep_check(bioc = TRUE, num_workers = available_cores(),
               timeout = as.difftime(20, units = "mins"), quiet = FALSE)
}


todo <- function() {
  pkgs <- tryCatch(revdep_todo(), error = function(ex) NA)
  if (identical(pkgs, NA)) {
    cat("Revdepcheck has not been initiated\n")
    return()
  }
  pkgs <- subset(pkgs, status == "todo")
  if (nrow(pkgs) == 0) {
    cat("There are no packages on the revdepcheck todo list\n")
  } else {
    cat(sprintf("%d. %s\n", seq_len(nrow(pkgs)), pkgs$package))
  }
}

parse_pkgs <- function(pkgs) {
  pkgs <- unlist(strsplit(pkgs, split = ",", fixed = TRUE))
  pkgs <- gsub("[ \t'\"‘’]", "", pkgs)
  sort(unique(pkgs))
}

revdep_init <- function() {
  if (!revdepcheck:::db_exists(".")) revdepcheck:::db_setup(".")
}

revdep_todo_reset <- function() {
  revdep_init()
  db <- revdepcheck:::db(".")
  df <- data.frame(package = character(0L), stringsAsFactors = FALSE)
  DBI::dbWriteTable(db, "todo", df, overwrite = TRUE, append = FALSE)
}

revdep_children <- local({
  cache <- list()
  function(pkg = NULL) {
    if (is.null(pkg)) pkg <- desc::desc(file = "DESCRIPTION")$get("Package")
    pkgs <- cache[[pkg]]
    if (is.null(pkgs)) {
      pkgs <- revdepcheck:::cran_revdeps(pkg)
      pkgs <- setdiff(pkgs, pkg) ## WORKAROUND
      cache[[pkg]] <- pkgs
    }
    pkgs
  }
})

revdep_pkgs_with_status <- function(status = "error") {
  status <- match.arg(status)
  res <- revdepcheck::revdep_summary()
  field <- switch(status, error = "errors")
  has_status <- vapply(res, FUN = function(x) {
    z <- x[["new"]][[field]]
    is.character(z) && any(nchar(z) > 0)
  }, FUN.VALUE = NA, USE.NAMES = TRUE)
  has_status <- !is.na(has_status) & has_status
  names(has_status)[has_status]
}

revdep_preinstall <- function(pkgs) {
  pkgs <- unique(pkgs)
  lib_paths_org <- lib_paths <- .libPaths()
  on.exit(.libPaths(lib_paths_org))
  lib_paths[1] <- sprintf("%s-revdepcheck", lib_paths[1])
  dir.create(lib_paths[1], recursive = TRUE, showWarnings = FALSE)
  .libPaths(lib_paths)
  message("Triggering crancache builds by pre-installing packages: ",
           paste(sQuote(pkgs), collapse = ", "))
  message(".libPaths():")
  message(paste(paste0(" - ", .libPaths()), collapse = "\n"))
  crancache::install_packages(pkgs)
}

args <- base::commandArgs()
if ("--reset" %in% args) {
  revdep_reset()
} else if ("--todo-reset" %in% args) {
  revdep_todo_reset()
  todo()
} else if ("--todo" %in% args) {
  todo()
} else if ("--add" %in% args) {
  pos <- which("--add" == args)
  pkgs <- parse_pkgs(args[seq(from = pos + 1L, to = length(args))])
  revdep_add(packages = pkgs)
  todo()
} else if ("--rm" %in% args) {
  pos <- which("--rm" == args)
  pkgs <- parse_pkgs(args[seq(from = pos + 1L, to = length(args))])
  revdep_rm(packages = pkgs)
  todo()
} else if ("--add-broken" %in% args) {
  revdep_add_broken()
  todo()
} else if ("--add-all" %in% args) {
  revdep_init()
  pkgs <- revdep_children()
  for (pkg in pkgs) {
    pkgs <- c(pkgs, revdepcheck:::cran_revdeps(pkg))
  }
  pkgs <- unique(pkgs)
  revdep_add(packages = pkgs)
  todo()
} else if ("--add-grandchildren" %in% args) {
  revdep_init()
  pkgs <- NULL
  for (pkg in revdep_children()) {
    pkgs <- c(pkgs, revdepcheck:::cran_revdeps(pkg))
  }
  pkgs <- unique(pkgs)
  revdep_add(packages = pkgs)
  todo()
} else if ("--show-check" %in% args) {
  pos <- which("--show-check" == args)
  pkgs <- parse_pkgs(args[seq(from = pos + 1L, to = length(args))])
  for (pkg in pkgs) {
    for (dir in c("old", "new")) {
      path <- file.path("revdep", "checks", pkg, dir, sprintf("%s.Rcheck", pkg))
      if (!utils::file_test("-d", path)) next
      pathname <- file.path(path, "00check.log")
      cat("-----------------------------------------------\n")
      cat(sprintf("%s (%s):\n", pkg, dir))
      cat("-----------------------------------------------\n")
      bfr <- readLines(pathname, warn = FALSE)
      tail <- tail(bfr, n = 20L)
      writeLines(tail)
    }
  }
} else if ("--list-error" %in% args) {
  cat(paste(revdep_pkgs_with_status("error"), collapse = " "), "\n", sep="")
} else if ("--add-error" %in% args) {
  revdepcheck::revdep_add(packages = revdep_pkgs_with_status("error"))
} else if ("--preinstall-error" %in% args) {
  res <- revdepcheck::revdep_summary()
  revdep_preinstall(revdep_pkgs_with_status("error"))
} else if ("--preinstall" %in% args) {
  pos <- which("--preinstall" == args)
  pkgs <- parse_pkgs(args[seq(from = pos + 1L, to = length(args))])
  revdep_preinstall(pkgs)
} else {
  check()
  revdep_report(all = TRUE)
}
