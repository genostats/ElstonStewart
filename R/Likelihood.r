# On va utiliser l'environnement .mem pour la mémoisation du cluster 
# et celle des pedigree si on tombe sur quelqu'un qui veut pas de
# calcul parallèle
.mem <- new.env()
ClusterEnv <- new.env()


#' Compute the log-likelihood of a parameter
#' 
#' @description Compute the log-likelihood of a parameter \code{theta}, given a list of pedigrees and a modele, using multiple cores and memoization
#' 
#' @param ped.set a list of pedigrees, each created with \code{es.pedigree}
#' @param modele a modele
#' @param theta a parameter for the modele
#' @param n.cores number of cores on which to run the computation
#' @param optim.alloc cf Details
#' @param sum.likelihoods cf Value
#' @param PSOCK Use PSOCK cluster instead of fork cluster (defaults to \code{TRUE} on non-Linux)
#' 
#' 
#' @details
#' The parameter \code{theta} will be given to the functions in \code{modele} to compute the likelihoods.
#' 
#' If \code{n.cores > 1} a cluster is created and left open for futur use with
#' the same parameters \code{ped.set} and \code{modele}. Open clusters can be closed
#' with \code{es.stopCluster()}.
#' 
#' If \code{optim.alloc = TRUE}, the function tries to optimize the distribution of
#' the computation (if \code{n.cores > 1}) between the cluster nodes. This has no effect on the first run
#' for given parameters \code{ped.set} and \code{modele}, but will reduce running time if the algorithm is ran several times
#' with different values of \code{theta} and same parameters \code{ped.set} and \code{modele}.
#' This is typically usefull for likelihood maximization.
#' 
#' 
#' 
#' @return If \code{sum.likelihoods = TRUE}, the function returns a single value, the sum of the
#' log-likelihoods computed for each pedigree of \code{ped.set}. Else, the function
#' returns a vector containing these log-likelihoods.
#' @seealso \code{\link{Elston}}, \code{\link{es.pedigree}}, \code{\link{es.stopCluster}}
#' 
#' 
#' @examples
#' 
#' data(fams)
#' 
#' # this data frame contains various families
#' # getting their famid
#' fam.ids <- unique(fams$fam);
#' 
#' # creating a list of genotypes corresponding to all individuals in fams
#' # NA -> 0, 1 or 2
#' genotypes <- lapply( fams$genotype, function(x) if(is.na(x)) 0:2 else x )
#' 
#' # creating a list of es.pedigree objects
#' X <- vector("list", length(fam.ids))
#' for(i in seq_along(fam.ids))
#' {
#'   w <- which(fams$fam == fam.ids[i])
#'   X[[i]] <- es.pedigree( id = fams$id[w], father = fams$father[w],
#'       mother = fams$mother[w], sex = fams$sex[w], pheno = rep(0, length(w)), 
#'       geno = genotypes[w], famid = fam.ids[i] )
#' }
#' 
#' \dontrun{# computing the likelihood for a single value p
#' Likelihood(X, modele.di, theta = list( p=0.5), n.cores=1 )
#' 
#' # computing the likelihood for a vector p (Elston-Stewart is ran only once!)
#' p <- seq(0,1,length=501)
#' L <- Likelihood(X, modele.di, theta = list( p=p ), n.cores=1 ) 
#' plot( p, exp(L), type="l")
#' 
#' # running an optimization algorithm
#' # Elston-Stewart is ran several times
#' # here we run the algorithm with 2 cores
#' L <- function(p) -Likelihood(X, modele.di, theta = list( p=p ), n.cores=2 ) 
#' optimize(L , c(0.35,0.45) )}
#' 
#' @export Likelihood
Likelihood <- function(ped.set, modele, theta, n.cores = getOption("mc.cores", 2L), optim.alloc = TRUE, sum.likelihoods = TRUE, PSOCK = Sys.info()["sysname"] != "Linux")
{
  if(n.cores == 1)
  {
    lik <- vector("list", length(ped.set))
    for(i in seq_along(ped.set))
    {
      a <- Elston(ped.set[[i]], modele, theta, .mem);
      lik[[i]] <- log(a$result);
      .mem <- a$mem;
    }
    if(sum.likelihoods)
      return(Reduce("+", lik))
    else
      return(Reduce(cbind,lik))
  }
  key <- paste("cluster",digest(list("cluster", ped.set, modele)),sep="")
  if(exists(key, envir = .mem)) 
  {
    # cat("reusing cluster\n")
    x <- get(key, envir = .mem)
    Likelihood.next.runs(theta, x, sum.likelihoods)
  }
  else
  {
    # cat("building cluster\n")
    if(PSOCK) {
      cl <- makePSOCKcluster(n.cores)
      clusterExport(cl, c("likelihood.0", "likelihood.00", "likelihood.0.vec", "likelihood.00.vec", "Elston", "ClusterEnv"), environment(likelihood.0))
    } else
      cl <- makeForkCluster(n.cores)

    x <- Likelihood.first.run(ped.set, modele, theta, cl, optim.alloc, sum.likelihoods, PSOCK)
    mem.sv(key, list(cl=cl, o=x$o), .mem)
    return(x$likelihood)
  }
}






#' Closing clusters
#' 
#' @description This function closes clusters opened by \code{Likelihood}
#' 
#' @param verbose if \code{TRUE}, the function will display a short information message when closing clusters
#' 
#' 
#' 
#' 
#' 
#' @seealso \code{\link{Likelihood}}
#' 
#' 
#' @export es.stopCluster
es.stopCluster <- function(verbose=TRUE)
{
  for(key in ls(pattern="cluster.*", envir = .mem))  
  {
    x <- get(key, envir = .mem)
    if(verbose) cat("stopping one cluster with",length(x$cl),"nodes\n")
    rm(list=key, envir = .mem)
    stopCluster(x$cl)
  }
}

Likelihood.first.run <- function(ped.set, modele, theta, cl, optim.alloc, sum.likelihoods, PSOCK)
{
  n.cores <- length(cl)
  # préparation découpage des données
  n.fam <- length(ped.set)
  X <- vector("list", n.cores)
  Indices <- vector("list", n.cores)
  i1 <- 1
  for(i in seq_along(cl))
  {
    i2 <- round(i*n.fam/n.cores)
    Indices[[i]] <- i1:i2
    i1 <- i2+1
  
    X[[i]] <- ped.set[ Indices[[i]] ]
  }

  # initialise PED et MODELE
  clusterApply(cl, X, function(x) assign( "PED", x, ClusterEnv) )# PED <<- X[[i]] )
  clusterApply(cl, rep(list(modele), n.cores), function(x) assign("MODELE", x, ClusterEnv)) # MODELE <<- x)
  # initialise MEM 
  clusterApply( cl, seq_len(n.cores), function(i) assign("MEMO", replicate( length(get("PED", envir=ClusterEnv)) , new.env()), ClusterEnv)) #  MEMO <<- replicate(length(PED),new.env()))

  if(!optim.alloc)   
  {
    # calcule
    if(sum.likelihoods)
      return( list( likelihood = Reduce("+", clusterApply(cl, rep(list(theta), n.cores), likelihood.0)), o = seq_along(ped.set) ))
    else
      return( list( likelihood = Reduce(cbind, clusterApply(cl, rep(list(theta), n.cores), likelihood.0.vec)),  o = seq_along(ped.set) ))
  }

  # calcule
  if(sum.likelihoods)
    R <- clusterApply(cl, rep(list(theta), n.cores), likelihood.00)
  else
    R <- clusterApply(cl, rep(list(theta), n.cores), likelihood.00.vec)

  time <- Reduce(function(x,y) c(x,y$time), R, NULL)

  # re-réparti
  Indices <- repartition(time, n.cores)
  o <- order( Reduce(c, Indices) )
  Y <- vector("list", n.cores)
  for(i in seq_along(cl))
  {
    X[[i]] <- ped.set[ Indices[[i]] ]
  }
  clusterApply(cl, X, function(x) assign( "PED", x, ClusterEnv) )# PED <<- X[[i]] )
  clusterApply( cl, seq_len(n.cores), function(i) assign("MEMO", replicate( length(get("PED", envir=ClusterEnv)) , new.env()), ClusterEnv))  #MEMO <<- replicate(length(PED),new.env()))
  if(sum.likelihoods)
    return(list( likelihood = Reduce(function(x,y) x + y$likelihood, R, 0), o = o) )


  return(list( likelihood = Reduce(function(x,y) cbind(x, y$likelihood), R,NULL), o = o ));
   
}

Likelihood.next.runs <- function(theta, x, sum.likelihoods)
{
  n.cores <- length(x$cl)
  if(sum.likelihoods)
    return( Reduce("+", clusterApply(x$cl, rep(list(theta), n.cores), likelihood.0)) )
  else
    return( Reduce(cbind, clusterApply(x$cl, rep(list(theta), n.cores), likelihood.0.vec))[,x$o] )
}

repartition <- function(x, n)
{
  Ind <- rep(list(numeric(0)), n)
  S <- rep(0,n)
  I <- order(x, decreasing=TRUE)
  x <- x[I]
  for(i in seq_along(x))
  {
    a <- x[i];
    k <- which.min(S)  # tjs remplir le moins plein...
    S[k] <- S[k] + a
    Ind[[k]] <- c(Ind[[k]], I[i])
  }
  return(Ind)
}


# --------------------------------------------------------------
# Ces fonctions utilisent les variables globales PED MEMO MODELE 
# et renvoient simplement le résultat
# Elles sont destinées à tourner sur les noeuds du cluster seulement.

#MODELE <- modele.pel
likelihood.00 <- function(theta, ped=get("PED", envir = ClusterEnv), modele = get("MODELE", envir = ClusterEnv), memo = get("MEMO", envir = ClusterEnv) ) 
{ 
  T <- numeric(length(ped))
  lik <- vector("list",length(ped))
  for(i in seq_along(ped))
  {
    T[i] <- system.time( {a <- Elston(ped[[i]], modele, theta, memo[[i]]); } )[1]
    gc() # si on ne force pas gc() après chaque calcul ça fausse l'estimation du temps nécessaire...
    lik[[i]] <- log(a$result);
    memo[[i]] <- a$mem; 
    # cat(i," : ", log(a$result), " (", T[i], ")\n", sep='')
  }
  return(list( likelihood=Reduce("+", lik) , time = T ) )
}

likelihood.0 <- function(theta, ped=get("PED", envir = ClusterEnv), modele = get("MODELE", envir = ClusterEnv), memo = get("MEMO", envir = ClusterEnv) ) 
{ 
  lik <- vector("list",length(ped))
  for(i in seq_along(ped))
  {
    a <- Elston(ped[[i]], modele, theta, memo[[i]]); 
    lik[[i]] <- log(a$result);
    memo[[i]] <- a$mem; 
    # cat(i," : ", log(a$result), "\n", sep='')
  }
  return(Reduce("+", lik))
}

likelihood.00.vec <- function(theta, ped=get("PED", envir = ClusterEnv), modele = get("MODELE", envir = ClusterEnv), memo = get("MEMO", envir = ClusterEnv) ) 
{ 
  T <- numeric(length(ped))
  lik <- vector("list",length(ped))
  for(i in seq_along(ped))
  {
    T[i] <- system.time( {a <- Elston(ped[[i]], modele, theta, memo[[i]]); } )[1]
    gc() # si on ne force pas gc() après chaque calcul ça fausse l'estimation du temps nécessaire...
    lik[[i]] <- log(a$result);
    memo[[i]] <- a$mem; 
    # cat(i," : ", log(a$result), " (", T[i], ")\n", sep='')
  }
  return(list( likelihood=Reduce(cbind,lik), time = T ) )
}

likelihood.0.vec <- function(theta, ped=get("PED", envir = ClusterEnv), modele = get("MODELE", envir = ClusterEnv), memo = get("MEMO", envir = ClusterEnv) ) 
{ 
  lik <- vector("list",length(ped))
  for(i in seq_along(ped))
  {
    a <- Elston(ped[[i]], modele, theta, memo[[i]]); 
    lik[[i]] <- log(a$result);
    memo[[i]] <- a$mem; 
    # cat(i," : ", log(a$result), "\n", sep='')
  }
  return(Reduce(cbind,lik))
}


