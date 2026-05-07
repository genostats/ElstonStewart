# split pedigree for ES

split_ped <- function(ped, mem = new.env())
{
  key <- digest(list("split_ped",ped))
  if(exists(key, envir = mem)) return(list(result= get(key, envir = mem), mem = mem))

  n.vec <- which(nb.enfants(ped) > 0 & ped$st[,"pere"] != 0);
  ped$st <- ped$st[, c('id', 'pere', 'mere', 'sexe') ]

  id.pivots <- matrix(ncol=length(n.vec), nrow=2)
  n.pivots  <- matrix(ncol=length(n.vec), nrow=2)
   
  new.id <- max(ped$st[,"id"])
  new.n  <- dim(ped$st)[1]
  for(i in seq_along(n.vec))
  {
    n <- n.vec[i]

    new.n  <- new.n  + 1;
    new.id <- new.id + 1

    id.pivots[1,i] <- ped$st[n,"id"] # la matrice id.pivot contient les ids des
    id.pivots[2,i] <- new.id         # pivots et de leur copie.
    n.pivots[1,i]  <- n;
    n.pivots[2,i]  <- new.n;

    ped$st <- rbind(ped$st, ped$st[n,]);
    ped$st[new.n,"id"] <- new.id
    ped$st[new.n,"pere"] <- ped$st[n,"pere"]    # le nouveau prend les parents
    ped$st[new.n,"mere"] <- ped$st[n,"mere"]
    ped$st[n,"pere"] <- ped$st[n,"mere"] <- 0   # le pivot garde les enfants mais n'a plus de parents
 
    ped$geno  <- c(ped$geno,  ped$geno[n])
    ped$pheno <- c(ped$pheno, ped$pheno[n])
  }
  
  # on repere les lignes qui correspondent aux familles nucleaires
  # [en fait les familles 2G]
  w <- rep(TRUE, new.n)
  f2G <- list();
  while(any(w))
  {
    b <- min( ped$st[ w , "id"] )  # on attrape un nouvel individu
    L <- 0;
    while(length(b) != L)  # on recupere sa famille nucleaire
    {
      L <- length(b);
      w.b <- is.element(ped$st[,"id"],b)
      # On ajoute les parents de tous les individus de la liste b
      b <- union(b, ped$st[w.b,"pere"])
      b <- union(b, ped$st[w.b,"mere"])
      b <- setdiff(b, 0)
      # on ajoute les enfants de tous les individus de la liste b
      b <- union(b, ped$st[ is.element(ped$st[,"pere"], b) | is.element(ped$st[,"mere"], b) , "id"]);
    }
    f2G <- c(f2G, list(w.b))
    w <- w & (!w.b); 
  } 

  # on prepare le terrain en donnant, pour chaque famille dans f2G, les id des pivots
  x <- as.vector(id.pivots)
  f2G.pivots <- vector('list', length(f2G))
  for(i in seq_along(f2G)) f2G.pivots[[i]] <- x[is.element(x, ped$st[f2G[[i]],"id"])]

  r <- list(ped = ped, id.pivots = id.pivots, n.pivots = n.pivots, f2G = f2G, f2G.pivots=f2G.pivots)
  mem.sv(key, r, mem)
  return( list(result=r, mem=mem));
}

