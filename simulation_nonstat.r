# simulate the non-stationary random fields
rm( list=ls() )

library( "RandomFields" )
library( "foreach" )
library( "doSNOW" )
library( "abind" )
n_cores <- 2


path <- "./"
name <- "02"

sim_nst <- function(r1=.5, 
                    r2=2.5, 
                    theta=0, 
                    dx=0, dy=0, nx=64, # 
                    nf=4 
                ){
    theta <- pi*theta/180
    A <- matrix( c(r1*cos(theta),-r2*sin(theta),r1*sin(theta),r2*cos(theta)),nrow=2)
    x <- seq(-1,1,,nx)
    y <- x + dy
    x <- x + dx
    
    mod <- RMmatern( nu=RMexp( Aniso=A ), scale=1 )
    z   <- RFsimulate( model=RPdirect(mod), x, y, n=nf )
    
    return( z )
}


cl   <- makeCluster(n_cores)
registerDoSNOW(cl)

# select parameters 
N <- 128           # size of the fields
n_cases   <- 200  # number of cases
n_members <- 10     # number of ensemble members
r1    <- c( .5, .5 )
r2    <- c( 1.5, 2.5 )
theta <- 45
th0   <- 80 # fixed threshold


nmod  <- length(r1)
modnames <- paste0( "r1=",r1," r2=",r2 )
gr   <- expand.grid( model=1:nmod )
n    <- dim( gr )[1]
bar  <- txtProgressBar(0,n,style=3)
prog <- function(i) setTxtProgressBar( bar, i )

# simulate without truncation (th=0)
res <- foreach( i=1:n, 
                .options.snow=list(progress=prog),
                .errorhandling="stop",
                .packages=c( "RandomFields", "abind" ),
                .combine=function( x, y ) abind( x, y, along=4 )
               ) %dopar%{
                    RFoptions( cPrintlevel=4, spConform=FALSE, max_variab=128**2, max_chol=128**2 )
                    m <- gr$model[i]
                    sim_nst( r1=r1[m], 
                             r2=r2[m], 
                             theta=theta,
                             nx=N,
                             nf=n_cases*(n_members+1) )
}
stopCluster(cl)

# re-order
rain <- array( dim=c( n_cases, N, N, nmod, n_members + 1 ), data=NA,
               dimnames=list( 1:n_cases, 1:N, 1:N, modnames, 0:n_members ) )
for( j in 1:(n_members+1) ) for( k in 1:nmod){
    rain[ ,,,k,j ] <- aperm( res[,, 1:n_cases + n_cases*( j-1 ) ,k] , c( 3,1,2 ) )
}


# truncate all fields, once randomly and once at a fixed threshold
truncate <- function( x, th ){
    x <- x - quantile( x, th )
    x <- x * ( x>0 )
    return( x )
}
for( i in 1:n_cases ) for( j in 1:(n_members+1) ) for( k in 1:nmod){
    rain[i,,,k,j] <- truncate( rain[i,,,k,j], th0/100 )
}
rain_constant_th <- rain


save( r1, r2, th0, theta, rain_constant_th, name, 
      file=paste0( path, "RRain_const_nst", name, ".rdata" ) )    

