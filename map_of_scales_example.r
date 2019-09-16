rm( list=ls() )
source( "verification_functions.r" ) 
require( "SpatialVx" )

# load the data, take the logarithm, transform
data( obs0513 )
N <- 1024
J <- 7
data( ICPg240Locs )
lon  <- ICPg240Locs[,1] 
lat  <- ICPg240Locs[,2]
rain <- log2( obs0513 + 2**-3 )
S    <- fld2S( rain, N=N, scales=1:J, smooth=FALSE, rsm=5)

# average over directions, remove negative values
S   <- ( S[ 1:J,, ] + S[ 1:J+J,, ] + S[ 1:J+2*J,, ] )/3
S   <- S*1*(S>0)

# get map of scales, only show pixels with rain
cen <- fld_centreofmass( S )
cen[obs0513<=0] <- NA

# plot
X11( width=10, height=5 )
par( mfrow=c(1,2) )
image.plot( rain )
par( usr=c( range(lon), range(lat) ) )
map( database="state", add=TRUE )
image.plot( cen )
par( usr=c( range(lon), range(lat) ) )
map( database="state", add=TRUE )

