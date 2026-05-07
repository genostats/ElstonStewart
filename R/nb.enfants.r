# pour pouvoir ajouter une colonne avec le nombre d'enfants !
nb.enfants <- function(ped)
{
  n <- dim(ped$st)[1]
  e <- numeric(n)
  for(i in 1:n)
  {
    if(ped$st[i,"sexe"] == 1) e[i] <- sum( ped$st[,"pere"] == ped$st[i,"id"] )
    if(ped$st[i,"sexe"] == 2) e[i] <- sum( ped$st[,"mere"] == ped$st[i,"id"] )
  }
  return(e);
}

