# functions required for the wavelet-based verrification
# two lengthtier functions related to a slightly modified version of the original wavelet transformation routine in LS2W are kept in "wavelet_functions.r" to keep this file reasonably short

require("imager")   # needed only for the shift of the basis functions
require("LS2W")     # contains the basic wavelet transform and bias correction
source( "wavelet_functions.r" ) # modified version of the wavelet transform

# get the positions of the basis function's centres of mass, calculated by "get_centres.r"
load( "DBcentres.rdata" ) 


# shift a matrix by (dx,dy) with periodic boundary conditions
shiftmat <- function( mat, dx, dy) imshift( as.cimg(mat), dx, dy, boundary=2 )[,,1,1]

# shift the output of the RDWT such that each basis function is centred at the centre of mass of its support.
cdd_centre <- function( cd # output of cddews / mycddews
                       ){
    S   <- cd$S
    N   <- dim( S )[2]
    cen <- DBcentres[[ paste0("N",N) ]][[ cd$family ]][[ cd$filter.number ]]
    dx <- N/2 - cen$x
    dy <- N/2 - cen$y
    
    if( length(dx)==0 | length(dy)==0 ){
        stop( paste0( "you forgot to calculate the shifts for N=", N ) )
    }
    
    for( j in 1:dim(S)[1] ) S[j,,] <- shiftmat( S[j,,], dx=dx[j], dy=dy[j] )
    cd$S <- S
    return( cd )
}

# linearly decrease the edges of a field to zero
smooth_borders <- function( field, # input field
                            r      # number of pixels to smooth over 
                           ){
    di <- dim( field )
    mask <- array( dim=di, data=0 )
    mv <- seq( 0,1,,r )
    for( i in 1:r ) mask[ i:(di[1]-i+1), i:(di[2]-i+1) ] <- mv[i]
    return( field*mask )
}

# paste a field of arbitary dimension onto an array of N times N zeroes, smooth at the edges
add_on_square <- function( picture,   # input field 
                           N,         # size of the output square, must be larger than the input
                           rsmooth=0  # number of pixels to smooth over 
                           ){
    
    nx <- nrow(picture)
    ny <- ncol(picture)
    if( nx > N | ny > N ) stop( "give me a bigger square" )
    if( rsmooth>0 ) picture <- smooth_borders( picture, r=rsmooth )
    
    px <- 1:nx + ceiling( (N - nx)/2 ) 
    py <- 1:ny + ceiling( (N - ny)/2 ) 
    
    res <- array( dim=c(N,N), data=0 )
    res[ px, py ] <- res[ px, py ] + picture
    return( list( px=px, py=py, res=res ) )
}

# transforms an intput field into an array of squared wavelet coefficients
# takes care of the boundary conditions, calculation of the bias correction etc.
fld2S <- function( field, # input field
                   N = 2**ceiling( log2( max( dim( field ) )  ) ), # dimension of the square on which the field is pasted for the wavelet transformation
                   rsm = N/20,  # number of pixels at the edge where the field is smoothed
                   fam="DaubExPhase", # wavelet family
                   filnum=2,    # number of vanishing moments
                   BS=TRUE,     # backshift the coefficients such that each basis function is centred at its centre of mass?
                   smn=4,   # number of vanishing moments of the wavelet used for spatial smoothing (see ?cddews)
                   ...  # further arguments passed to mycddews 
                   ){
    # paste onto an NxN square
    fmin  <- min(field)
    field <- field - fmin
    square <- add_on_square( field, N, rsmooth=rsm )
    field  <- square$res + fmin
    
    # get the correction matrix and do the transformation
    AI  <- get_AINV( paste( fam, filnum, sep="_" ), N=N )
    S   <- mycddews( field, 
                     family=fam, 
                     filter.number=filnum, 
                     Ainv=AI, 
                     sm.filter.number=smn, 
                     ... )
    # shift back?
    if(BS) S   <- cdd_centre( S )

    # cut out the original dimensions of the field and return
    S <- S$S[,square$px, square$py]
    return( S )
}

# get the centre of mass of a vector
centreofmass <- function( m, x=1:length(m) ) return( sum(m*x)/sum(m) )


# get the centre of mass at each pixel
fld_centreofmass <- function( fld ){
    fld <- fld * ( 1* fld>0 )
    x   <- 1:dim(fld)[1]
    res <- nor <- array( dim=dim( fld[1,,]), data=0 ) 
    
    for( i in 1:dim(fld)[1] ){
        res <- res + x[i]*fld[i,,]
        nor <- nor + fld[i,,]
    }
    return( res/nor )
}

