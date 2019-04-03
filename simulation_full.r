# simulate a large number of realizations for our "real" experiment
rm( list=ls() )
source( "stochastic_rain.r" )
library( "foreach" )
library( "doSNOW" )
library( "abind" )
n_cores <- 4

path <- "./"
name <- "03"

cl   <- makeCluster(n_cores)
registerDoSNOW(cl)

# select parameters 
N <- 128           # size of the fields
n_cases   <- 1000  # number of cases
n_members <- 1     # number of ensemble members
nu    <- c( 2.5, 3, 2.5, 3 ) # smoothnesses
sc    <- c( .1, .1, .2, .2 ) # scale parameters
thmin <- 75 # lowest random threshold
th0   <- 80 # fixed threshold
thmax <- 85 # highest random threshold


nmod  <- length(nu)
modnames <- paste0( "nu=",nu," sc=",sc )
gr   <- expand.grid( model=1:nmod, mem=0:n_members )
n    <- dim( gr )[1]
bar  <- txtProgressBar(0,n,style=3)
prog <- function(i) setTxtProgressBar( bar, i )

# simulate without truncation (th=0)
res <- foreach( i=1:n, 
                .options.snow=list(progress=prog),
                .errorhandling="stop",
                .packages=c( "RandomFieldsUtils", "RandomFields", "abind" ),
                .combine=function( x, y ) abind( x, y, along=4 )
               ) %dopar%{
                    m <- gr$model[i]
                    simulate_rain( nu=nu[m], 
                                   scale=sc[m], 
                                   x=1:N, 
                                   y=1:N, 
                                   nf=n_cases, 
                                   th=0 )
}
stopCluster(cl)

# re-order
rain <- array( dim=c( n_cases, N, N, nmod, n_members + 1 ), data=NA,
               dimnames=list( 1:n_cases, 1:N, 1:N, modnames, 0:n_members ) )
for( i in 1:n ) rain[ ,,,gr$model[i],gr$mem[i]+1 ] <- res[,,,i]
rm( res )
gc()

# truncate all fields, once randomly and once at a fixed threshold
truncate <- function( x, th ){
    x <- x - quantile( x, th )
    x <- x * ( x>0 )
    return( x )
}
thv <- c()
rain_constant_th <- rain_random_th <- rain
for( i in 1:n_cases ) for( j in 1:(n_members+1) ) for( k in 1:nmod){
    thr <- runif( 1, thmin, thmax ) 
    rain_constant_th[i,,,k,j] <- truncate( rain_constant_th[i,,,k,j], th0/100 )
    rain_random_th[i,,,k,j]   <- truncate( rain_random_th[i,,,k,j], thr/100 )
    thv <- c( thv, thr )
}

# save all fields and parameters
save( nu, sc, th0, thv, rain, name, 
      file=paste0( path, "RRain_full", name, ".rdata" ) )    


save( nu, sc, th0, thv, rain_constant_th, name, 
      file=paste0( path, "RRain_const", name, ".rdata" ) )    

save( nu, sc, th0, thv, rain_random_th, name, 
      file=paste0( path, "RRain_rand", name, ".rdata" ) )    
      


    
