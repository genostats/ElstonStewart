
### modele pour un marqueur di-allélique
### (avec un placeholder pour la penetrance)
trans <- function(of, fa, mo)
{
  if(fa == 0)
  {
    if(mo == 0 & of == 0)
      return(1);
    if(mo == 1 & (of == 0 | of == 1))
      return(.5);
    if(mo == 2 & of == 1)
      return(1);
  }
  if(fa == 1)
  {
    if(mo == 0 & (of == 0 | of == 1))
      return(.5);
    if(mo == 1 & (of == 0 | of == 2))
      return(.25)
    if(mo == 1 & of == 1)
      return(.5)
    if(mo == 2 & (of == 1 | of == 2))
      return(.5)
  }
  if(fa == 2)
  {
    if(mo == 0 & of == 1)
      return(1);
    if(mo == 1 & (of == 1 | of == 2))
      return(.5);
    if(mo == 2 & of == 2)
      return(1);
  }
  return(0)
};


modele.di <- list(
   name = "diallelic",
   proba.g = function(g, theta)
    {
      if(g == 0)
        return(theta$p**2)
      if(g == 1)
        return(2*theta$p*(1-theta$p))
      if(g == 2)
        return((1-theta$p)**2)
      stop("Unknown Genotype");
    },
  trans = trans,
  p.pheno = function(x,g,theta) 1
)

