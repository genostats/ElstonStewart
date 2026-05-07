extract.pivot <- function(splitted, mem = new.env())
{
  key <- digest(list("extract.pivot",splitted))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  # ----------------------------------------
  # une heuristique pour le choix de pivot
  # prioritaire : pivots ou le nombre de genotypes a envisager est le plus petit
  # L <- sapply(splitted$ped$geno[splitted$n.pivots[1,]], length)
  # if(any(L == 1))
  # {
  #   k <- which(L == 1)[1]
  # }
  # else # sinon : pivots pris dans la famille ou il y a le moins de pivot
  # {
  #   id.piv <- splitted$f2G.pivots[[which.min(sapply(splitted$f2G.pivots, length))]][1]
  #   k <- which( apply((splitted$id.pivots == id.piv),2,any) )
  # }
  ### k <- which.min(L)
  # ----------------------------------------
  k <- 1

  n.piv.1  <- splitted$n.pivots[1,k]
  n.piv.2  <- splitted$n.pivots[2,k]
  id.piv.1 <- splitted$id.pivots[1,k]
  id.piv.2 <- splitted$id.pivots[2,k]

  ph <- splitted$ped$pheno[[n.piv.1]];
  g <- splitted$ped$geno[[ n.piv.1 ]]

  # supression pivot...
  splitted$n.pivots  <- splitted$n.pivots[,-k,drop=FALSE]   
  splitted$id.pivots <- splitted$id.pivots[,-k,drop=FALSE]   

  F2G <- list()
  w.fix <- rep(FALSE, length(splitted$f2G))  # pour reperer dans f2G les familles qui n'ont plus de pivots
  w <- rep(FALSE, dim(splitted$ped$st)[1])      # pour reperer les individus qui appartiennent a ces familles [pour extraction]
  for(i in seq_along(splitted$f2G))
  {
    a <- splitted$f2G.pivots[[i]] 
    splitted$f2G.pivots[[i]] <- a[ a != id.piv.1 & a != id.piv.2 ];
    if( length( splitted$f2G.pivots[[i]] ) == 0 ) # plus de pivot dans cette famille
    {
      w.fix[i] <- TRUE
      w <- w | splitted$f2G[[i]];
      # extraction de la famille 2G
      ped <- list(st = splitted$ped$st[splitted$f2G[[i]],,drop=FALSE],
                  pheno = splitted$ped$pheno[splitted$f2G[[i]]], geno=splitted$ped$geno[splitted$f2G[[i]]] )
      # ou sont les pivots dans cette famille ?
      n.piv <- which( ped$st[,"id"] == id.piv.1 | ped$st[,"id"] == id.piv.2 )
      F2G <- c( F2G, list(list(ped=ped, n.piv=n.piv) ))
    }
  }

  # on efface de splitted les familles f2G qui ont ete extraites ci-dessus
  if(any(w))
  {  
    splitted$ped <- list( st = splitted$ped$st[!w,,drop=FALSE] , geno = splitted$ped$geno[!w], pheno = splitted$ped$pheno[!w] )
    splitted$f2G <- splitted$f2G[!w.fix];
    for(i in seq_along(splitted$f2G)) 
      splitted$f2G[[i]] <- splitted$f2G[[i]][!w]
    splitted$f2G.pivots <- splitted$f2G.pivots[!w.fix];
    splitted$n.pivots <-  matrix(sapply( splitted$id.pivots, function(i) which(splitted$ped$st[,"id"]==i) ), nrow=2)
  }
  n.piv = which( is.element(splitted$ped$st[,"id"], c(id.piv.1, id.piv.2)))
  r <- list(n.piv = n.piv , pheno.piv = ph, geno.piv = g, splitted = splitted, F2G = F2G)
  mem.sv(key, r, mem);
  return(list(result=r, mem=mem))
}

