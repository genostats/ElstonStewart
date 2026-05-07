#' @exportS3Method print es.pedigree
print.es.pedigree <- function(x, ...) {
  cat("An es.pedigree object with", dim(x$st)[1], "individuals\n")
}

