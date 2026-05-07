#' Creates a pedigree
#' @aliases es.pedigree print.es.pedigree
#' 
#' @description Creates a pedigree from vectors giving individuals' id, fathers' id, mothers' is, sex, phenotype and genotype
#' 
#' @param id a numeric vector of unique individuals id. Ids must be different from 0
#' @param father a vector giving the id of the father of each individual (0 for founder)
#' @param mother a vector giving the id of the mother of each individual (0 for founder)
#' @param sex a vector giving the sex of each individual (1 = male, 2 = female)
#' @param pheno a list giving the phenotypes of each individual
#' @param geno a list giving the genotypes of each individual
#' @param famid (facultative) a family id, used only to issue error messages
#' 
#' 
#' @details All vectors must have the same length, including the lists \code{pheno} and \code{geno}. The list \code{geno} is a list
#' of vectors containing the possible genotypes of each individual. The format of \code{pheno} depends on the modele which will
#' be used to compute likelihoods with \code{Elston} or \code{Likelihood}.
#' 
#' 
#' 
#' @return An object of class S3 \code{es.pedigree}
#' @seealso \code{\link{plot.es.pedigree}}
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
#' X  # displays a short information on X
#' 
#' @export es.pedigree
es.pedigree <- function(id, father, mother, sex, pheno, geno, famid)
{
  if(missing(famid)) 
    famid = ""
  ## data check
  n <- length(id)
  if( length(father) != n | length(mother) != n | length(sex) != n | length(pheno) != n | length(geno) != n )
    stop(paste(famid, "id, father, mother, sex, pheno, geno should have same length"))

  if(any(duplicated(id)))
    stop(paste(famid, "duplicated ids"))

  nz.father <- father[father!=0]
  if( length(setdiff(nz.father, id)) != 0 )
     stop(paste(famid, "father ids", setdiff(nz.father, id),"not found in id list"))
  if(any(sex[match(nz.father, id)] != 1))
     stop(paste(famid, "for some individuals indicated as fathers sex is not male (1)"))

  nz.mother <- mother[mother!=0]
  if( length(setdiff(nz.mother, id)) != 0 )
     stop(paste(famid, "mother ids", setdiff(nz.mother, id),"not found in id list"))
  if(any(sex[match(nz.mother, id)] != 2))
     stop(paste(famid, "for some individuals indicated as mothers sex is not female (2)"))

  if( any(xor(mother == 0, father == 0)) )
     stop(paste(famid, "individuals must have both parents in the pedigree, or be founders"))
  ## --

  ped <- list( st = cbind(id = id, pere = father, mere = mother, sexe = sex), pheno = pheno, geno = geno, famid = famid)
  ped$st <- cbind(ped$st, nb.enfants=nb.enfants(ped));
  class(ped) <- "es.pedigree"
  return(ped)
}

