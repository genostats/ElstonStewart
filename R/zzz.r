#' The pedigree of Conrad II, Holy Roman Emperor
#' @name conrad2
#' @docType data
#' 
#' @description
#' A pedigree with imbreding used in examples
#' 
#' @details See \code{vignette("Elston-Stewart")} for the name of some other individuals in the pedigree.
#' 
#' 
#' @format
#' A data frame with 22 lines with the following variables.
#' \describe{
#' \item{list("id")}{subject id (22 is Conrad II, 2 is Louis the Stammerer}
#' \item{list("father")}{id of the subject's father}
#' \item{list("mother")}{id of the subject's mother}
#' \item{list("sex")}{1=male, 2=female}
#' }
#' 
#' 
#' 
#' @keywords datasets
#' @examples
#' 
#' data(conrad2)
#' 
NULL

#' Fifty families 
#' @name fams
#' @docType data
#' 
#' @description
#' A set 50 pedigrees
#' @format
#' A data frame with 802 lines with the following variables.
#' \describe{
#' \item{list("fam")}{family id}
#' \item{list("id")}{subject id}
#' \item{list("father")}{id of the subject's father}
#' \item{list("mother")}{id of the subject's mother}
#' \item{list("sex")}{1=male, 2=female}
#' \item{list("genotype")}{coded additively. Many genotypes are unkwnon (\code{NA})}
#' }
#' 
#' @keywords datasets
#' @examples
#' 
#' data(fams)
#' 
NULL

#' @importFrom digest digest
#' @importFrom parallel makePSOCKcluster makeForkCluster clusterExport clusterApply stopCluster 
#' @importFrom methods hasArg is
NULL

