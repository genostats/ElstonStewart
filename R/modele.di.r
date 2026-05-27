#' Base modele for a di-allelic marker
#' 
#' @details \code{modele.di} is a list with components 'name' (value "diallelic"), 'trans' : a function
#' 'trans(of, fa, mo)' giving transmission probabilities for a genotype coded as 0, 1, 2, and 
#' 'p.pheno', a function 'p.pheno(x, g, theta)' which always return 1. 
#' @details See vignette for examples.
#'
#' @export
modele.di <- list(
  name = "diallelic",
  # pr
  proba.g = function(g, theta) {
    if(g == 0)
      return(theta$p**2)
    if(g == 1)
      return(2*theta$p*(1-theta$p))
    if(g == 2)
      return((1-theta$p)**2)
    stop("Unknown Genotype");
  },
  # proba de transmisions
  trans = function(of, fa, mo) {
    if(fa == 0) {
      if(mo == 0 & of == 0) return(1)
      if(mo == 1 & (of == 0 | of == 1)) return(.5)
      if(mo == 2 & of == 1) return(1)
    }
    if(fa == 1) {
      if(mo == 0 & (of == 0 | of == 1)) return(.5)
      if(mo == 1 & (of == 0 | of == 2)) return(.25)
      if(mo == 1 & of == 1) return(.5)
      if(mo == 2 & (of == 1 | of == 2)) return(.5)
    }
    if(fa == 2) {
      if(mo == 0 & of == 1) return(1)
      if(mo == 1 & (of == 1 | of == 2)) return(.5)
      if(mo == 2 & of == 2) return(1)
    }
    return(0)
  },
  # on ignore le paramètre theta 
  p.pheno = function(x, g, theta) 1
)

