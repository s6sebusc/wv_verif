# Originally, the value of the projection on each daughter wavelet is stored at the start of that function's support, meaning that the energy of an individual feature is not concentrated at that feature's location, but mostly far away from it. The issue is resolved by calculating the centre of mass of each daughter's squared values and later shifting each of them back such that they are located at their centre of mass. Because the boundaries of the wavelet transform are periodic, you can't simply find the centre and shift back in one step. The simplest solution I came up with so far is this iterative procedure - if you have a better idea please let me know.
# Also keep in mind that the i-th Daubechies wavelet at scale j has support length i*( 2**(j+1)-2 ) - (2**j-2), if your domain is smaller than that, this whole thing may not work and you should probably not be using that scale.

rm( list=ls() )
source( "verification_functions.r" )

Nv <- 2**(4:10)

# finds the centre (xC, yC) of a matrix
cen2d <- function( mat ){
    gr  <- expand.grid( 1:nrow(mat), 1:ncol(mat) )
    mat <- as.vector( mat ) / sum( mat, na.rm=TRUE )
    cen <- c( sum(gr[,1]*mat, na.rm=TRUE), sum(gr[,2]*mat, na.rm=TRUE ) )
    return( cen )
}

get_centres <- function( N, fam="DaubExPhase", num=2, maxiter=100 ){
    # get a picture of the basis function
    data          <- array( dim=c(N,N), data=0 )
    data[N/2,N/2] <- 1
    S <- mycddews( data, family=fam, filter.number=num, smooth=FALSE, correct=FALSE )$S
    
    dx <- dy <- rep( 0, dim(S)[1] )
    target   <- c(N/2,N/2) # desired centre of the basis functions
    
    for( i in 1:dim(S)[1] ){ #loop over all scales
        ce    <- c(0,0)
        k     <- 0
        # find necessary shift iteratively 
        # ( one iteration is usually insufficient due to periodic boundaries, this is not elegant but it works )
        while( any( abs( ce - target )>1 ) & k < maxiter ){
            if(k>0){ # update the shifts
                dx[i] <- dx[i] + N/2 - ce[1]
                dy[i] <- dy[i] + N/2 - ce[2]
            }
            # shift the matrix, find its centre
            Si <- shiftmat( S[i,,], dx=dx[i], dy=dy[i] )
            ce <- cen2d( Si )
            k <- k + 1
        }
        # this may not converge if the basis function is much larger than N. In that case, you shouldn't be using it anyways.
        if( any( abs( ce - target )>1 ) ){ 
            warning( paste0("no convergence for ",fam,num,", N=",N,", scale #",i) )
        }
    }
    # convert shift to centre
    x <- N/2 - dx
    y <- N/2 - dy
    return( list(x=x, y=y) )
}

# go over all desired basis functions
DBcentres <-  list()
for( N in Nv ){
    NN <- paste0( "N",N )
    DBcentres[[ NN ]] <- list()
    for( nu in 1:10 ){ 
        print( paste0( "working on D", nu, ", N=", N  ) )
        DBcentres[[ NN ]][["DaubExPhase"]][[nu]] <- get_centres(N, fam="DaubExPhase", num=nu)
        if( nu > 3 ){
             DBcentres[[ NN ]][["DaubLeAsymm"]][[nu]] <- get_centres(N, fam="DaubLeAsymm", num=nu)
        }
    }
}

save( DBcentres, file="DBcentres.rdata" )
