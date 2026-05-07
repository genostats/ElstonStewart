
Likelihood.unrelated <- function(ped, modele, theta)
{
  # cas du pedigree vide
  if(dim(ped$st)[1] == 0) return(1);

  P <- 1;
  for(i in 1:dim(ped$st)[1])
  {
    x <- ped$pheno[[i]];
    g <- ped$geno[[i]];
    S <- 0;
    for(a in g)
      S <- S + modele$proba.g(a,theta)*modele$p.pheno(x,a,theta);
    P <- P*S;
    if(all(P == 0)) return(P)
  }
  P;
}

