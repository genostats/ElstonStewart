
effeuille.sibship <- function(ped, mem = new.env())
{
  # memoization 
  key <- digest(list("effeuille.sibship",ped))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  # choix d'une feuille : le premier enfant du pedigree
  # + préparation du pedigree dépouillé de cette feuille
  n   <- which(ped$st[,"pere"]>0)[1]
  p   <- ped$st[n,"pere"];
  m   <- ped$st[n,"mere"];

  n.p <- which(ped$st[,"id"] == p)
  n.m <- which(ped$st[,"id"] == m)
  gp <- ped$geno[[n.p]];
  gm <- ped$geno[[n.m]];

  # on effeuille d'un coup toute la fratrie !!
  n.sib <- which( ped$st[,"pere"] == p & ped$st[,"mere"] == m )

  ped.x <- list( st=ped$st[-n.sib,], geno=ped$geno[-n.sib], pheno=ped$pheno[-n.sib] )
  n.parents.x <- c(which(ped.x$st[,"id"] == p), which(ped.x$st[,"id"] == m))

  r <- list(gp=gp, gm=gm, n.sib=n.sib, ped.x=ped.x, n.parents.x=n.parents.x)
  mem.sv(key, r, mem)
  return( list(result = r, mem = mem) )
}

