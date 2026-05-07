#' Compute a probability function on a pedigree
#' 
#' @description Compute the probability of a pedigree, for a parameter \code{theta} and a given modele
#' 
#' @param ped a pedigree, created by \code{es.pedigree}
#' @param modele a modele
#' @param theta a parameter for the modele
#' @param mem an environment used for memoization
#' 
#' 
#' @details This function runs Elston-Stewart algorithm. If \code{mem} is provided, some intermediate results of previous runs stored in it can be re-used.
#' 
#' 
#' 
#' @return A list with two components, \code{result}: the probability value, and \code{mem}: an object storing intermediate results for future runs
#' @seealso \code{\link{Likelihood}}, \code{\link{es.pedigree}}
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
#' # running Elston-Stewart
#' r <- Elston(X, modele.di, list(p = 0.98))
#' r$result
#' 
#' # using the memoization...
#' r <- Elston(X, modele.di, list(p = 0.99), r$mem)
#' r$result
#' 
#' @export Elston
Elston <- function(ped, modele, theta, mem = new.env())
{
  if(!is(ped,"es.pedigree"))
    stop("Argument ped is not of class es.pedigree")
  Re <- split_ped(ped, mem)
  Elston.on.splitted( Re$result, modele, theta, Re$mem)
}

