#' Plot pedigrees
#' 
#' @description Plot objects created with \code{es.pedigree}
#' 
#' @param x an object created with \code{es.pedigree}
#' @param ... additionnal arguments including possibly a vector of phenotypes \code{pheno}, and other arguments to be given to \code{plot.pedigree}
#' 
#' 
#' @details This function relies on the function \code{plot.pedigree} provided by the package \code{kinship2}.
#' The phenotype is build assuming that \code{x$pheno} is a list of vectors, the first element
#' of which gives the status. Alternatively, a vector \code{pheno} of appropriate length can be given
#' as an additionnal argument.
#' 
#' 
#' 
#' @seealso \code{\link{es.pedigree}}, \code{\link{plot.pedigree}}
#' 
#' 
#' @examples
#' 
#' ## cf Elston-Stewart vignette for more coments on this example
#' data(conrad2)
#' # creating an es.pedigree object
#' genotypes <- c( rep(list(0:2), 21), 2 )
#' 
#' X <- es.pedigree( id = conrad2$id, father = conrad2$father, mother = conrad2$mother,
#'       sex = conrad2$sex, pheno = rep(0, 22), geno = genotypes )
#' 
#' plot(X)
#' 
#' @exportS3Method plot es.pedigree
plot.es.pedigree <- function(x, ...)
{
  pheno <- if(hasArg("pheno")) list(...)$pheno else x$pheno 
  if(mode(pheno) == "list")
    pheno <- Reduce( function(x,y) c(x,y[1]), pheno, NULL)
  ped.k2 <- pedigree( x$st[,"id"], x$st[,"pere"], x$st[,"mere"], x$st[,"sexe"], pheno )
  plot.pedigree(ped.k2, ...)
}

