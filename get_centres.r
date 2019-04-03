# Originally, the value of the projection on each daughter wavelet is stored at the start of that function's support, meaning that the energy of an individual feature is not concentrated at that feature's location, but mostly far away from it. The issue is resolved by calculating the centre of mass of each daughter's squared values and later shifting each of them back such that they are located at their centre of mass.

rm( list=ls() )
require( "LS2W" )

Nv <- 2**(4:9)

cen2d <- function( mat ){
    gr  <- expand.grid( 1:nrow(mat), 1:ncol(mat) )
    mat <- as.vector( mat ) / sum( mat, na.rm=TRUE )
    cen <- c( sum(gr[,1]*mat, na.rm=TRUE), sum(gr[,2]*mat, na.rm=TRUE ) )
    return( cen )
}

get_centres <- function( N, fam="DaubExPhase", num=2 ){
    data          <- array( dim=c(N,N), data=0 )
    data[N/2,N/2] <- 1
    S <- cddews( data, family=fam, filter.number=num, smooth=FALSE, correct=FALSE )$S

    x <- y <- c()
    for( i in 1:dim(S)[1] ){
        ce <- cen2d( S[i,,] )
        x[i] <- ce[1]
        y[i] <- ce[2] 
    }
    return( list(x=x, y=y) )
}


DBcentres <-  list()
for( N in Nv ){
    NN <- paste0( "N",N )
    DBcentres[[ NN ]] <- list()
    for( nu in 1:10 ){ 
        DBcentres[[ NN ]][["DaubExPhase"]][[nu]] <- get_centres(N, fam="DaubExPhase", num=nu)
        if( nu > 3 ){
             DBcentres[[ NN ]][["DaubLeAsymm"]][[nu]] <- get_centres(N, fam="DaubLeAsymm", num=nu)
        }
    }
}

save( DBcentres, file="DBcentres.rdata" )
