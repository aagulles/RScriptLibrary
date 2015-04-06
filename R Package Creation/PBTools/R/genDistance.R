#----------------------------------------------------------------
# This function generates genetic distance matrix. This uses the dist.genpop function of adegenet.                      
# PBTools uses this for popgen analysis
#
# ARGUMENTs:
# genindObject - genind object
# method - 1 if Nei's distance; 2 if Angular distance or Edwards' distance; 3 if Coancestrality coefficient or Reynolds' distance
#          4 if Classical Euclidean distance or Rogers' distance; 5 if  Absolute genetics distance or Provesti 's distance
# displayDiag - TRUE if diagonals will be displayed
# displayUpper - TRUE if upper triangle of the matrix will be displayed
#                                                     
# Script Created by: Nellwyn L. Sales
#----------------------------------------------------------------

genDistance<-function(genindObject, method=1, displayDiag=FALSE, displayUpper=FALSE) UseMethod("genDistance")

genDistance.default<-function(genindObject, method=1, displayDiag=FALSE, displayUpper=FALSE) {
  options(width=500)
  library(adegenet)
  
  #convert genindObject to genpop object
  pop <- genind2genpop(genindObject, miss = "0")
  
  dist<-dist.genpop(pop, met=method, diag=displayDiag, upper=displayUpper)
  result<-round(dist, digits=4)
  return(result)
}