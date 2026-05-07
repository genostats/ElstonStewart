
Likelihood.2g <- function(ped, modele, theta, mem = new.env())
{
  # memoization 
  key <- digest(list("Likelihood.2g",ped,modele$name,theta))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  # pas (plus) d'enfants, que des individus non apparentés
  if(sum(ped$st[,"pere"]>0)==0)
    return( list( result=Likelihood.unrelated(ped, modele, theta), mem=mem) );

  eff <- effeuille.sibship(ped, mem);
  r <- eff$result
  mem <- eff$mem

  S <- 0;
 
  for(b in r$gp)
  {
    for(c in r$gm)
    {
      r$ped.x$geno[r$n.parents.x] <- c(b,c);
      lik.2g <- Recall(r$ped.x, modele, theta, mem);
      pro.x  <- lik.2g$result
      mem <- lik.2g$mem
      if(all(pro.x == 0)) next;
      P <- 1;
      for( n.1 in r$n.sib )
      { 
        S1 <- 0;
        x  <- ped$pheno[[n.1]];
        gf <- ped$geno[[n.1]];
        for(a in gf)
          S1 <- S1 + modele$p.pheno(x,a, theta)*modele$trans(a,b,c)
        P <- P*S1;
        if(all(P == 0)) break;
      } 
      S <- S + pro.x*P;
    }
  }
  mem.sv(key, S, mem)
  return( list(result = S, mem = mem) )
}  
