# functions required for the wavelet-based verrification
# two lengthtier functions related to a slightly modified version of the original wavelet transformation routine in LS2W are kept in "wavelet_functions.r" to keep this file reasonably short

require("imager")   # needed only for the shift of the basis functions
require("LS2W")     # contains the basic wavelet transform and bias correction
source( "wavelet_functions.r" ) # modified version of the wavelet transform


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
                   smn=4,   # number of vanishing moments of the wavelet used for spatial smoothing (see ?cddews)
                   scales=NULL, # which scales to use: make sure to omit the very large scales
                   smooth=TRUE,
                   ...  # further arguments passed to mycddews 
                   ){
    # paste onto an NxN square
    fmin  <- min(field)
    field <- field - fmin
    square <- add_on_square( field, N, rsmooth=rsm )
    field  <- square$res + fmin
    if( is.null( scales ) ) J <- log2( scales ) else J <- max( scales )
    
    # get the correction matrix and do the transformation
    AI  <- get_AINV( paste( fam, filnum, sep="_" ), N=2**J )
    S   <- mycddews( field, 
                     family=fam, 
                     filter.number=filnum, 
                     Ainv=AI, 
                     sm.filter.number=smn, 
                     scales = scales,
                     smooth = smooth,
                     ... )

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

