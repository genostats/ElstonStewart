Elston.on.splitted <- function(splitted, modele, theta, mem=new.env())
{
  key <- digest(list("bidouille",splitted, modele, theta))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  if(dim(splitted$n.pivots)[2] == 0) # plus de pivots
  {
    P <- 1;
    for(i in seq_along(splitted$f2G))
    {  
      w <- splitted$f2G[[i]]
      ped <- list( st = splitted$ped$st[w,,drop=FALSE], geno=splitted$ped$geno[w], pheno=splitted$ped$pheno[w] )
      lik.2g <- Likelihood.2g(ped,modele,theta,mem)
      P <- P*lik.2g$result;
      mem <- lik.2g$mem;
    }
    mem.sv(key, P, mem);  
    return( list(result=P, mem=mem) );
  }

  Re <- extract.pivot(splitted, mem);
  extract <- Re$result;
  mem <- Re$mem;

  # enumeration genotypes
  S <- 0;
  for(a in extract$geno.piv)
  {
    numerateur <- modele$proba.g(a,theta)*modele$p.pheno(extract$pheno.piv,a,theta)
    if(all(numerateur == 0)) ## pas la peine de considerer a, combi geno/pheno incompatible
      next;
    numerateur[ numerateur == 0 ] <- 1;
    extract$splitted$ped$geno[extract$n.piv] <- a
    Re <- Recall(extract$splitted, modele, theta, mem)
    P <- Re$result
    mem <- Re$mem
    for(i in seq_along(extract$F2G))
    {
      f2g <- extract$F2G[[i]]$ped
      n.piv.f2g <- extract$F2G[[i]]$n.piv
      f2g$geno[n.piv.f2g] <- a
      lik.2g <- Likelihood.2g(f2g,modele,theta,mem)
      P <- P*lik.2g$result;
      mem <- lik.2g$mem;
    }
    S <- S + P/numerateur
  }
  mem.sv(key, S, mem);
  return( list(result=S, mem=mem) );
}

