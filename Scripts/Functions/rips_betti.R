#!/usr/bin/env Rscript

# Load TDA package
if(require(TDA)==0) install.packages("TDA"); library(TDA)
      
# Load Rips filtration
rips_betti <- function(matriz, mdim){
    # Compute distance matrix
    D=1-matriz
    # Take number of nodes
    l=dim(D)[1]
    # Maximum value of the rips filtration
    mscale <- 1
    # Compute rips diagram
    DiagTri <- ripsDiag(D, mdim, mscale, dist = "arbitrary", printProgress = F)
    # Return results
    return(DiagTri)
}
