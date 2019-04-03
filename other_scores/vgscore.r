# A realtively fast fortran-implementation of the stationary Variogram-score from Scheuerer et al 2015: "Variogram-based proper scoring rules for probabilistic forecasts of multivariate quantitie"
# code by Franka Nawrath, used with permission

vgpath <- "./other_scores/"
require( "abind" )


# calculate the p-variogram
vgram <- function(y, R, p=0.5) {
  
  switch(as.character(p),
         '0.5' = { dyn.load( paste0(vgpath, "vg_05.so" ) ) }, 
         '1'   = { dyn.load( paste0(vgpath, "vg_1.so" ) ) },
         '2'   = { dyn.load( paste0(vgpath, "vg_2.so" ) ) },
         stop( 'only p=0.5,1,2 are allowed.' ) 
         )
  
  nx <- nrow(y)
  ny <- ncol(y)
  x <- array(1.5,dim=c(R+1,R+1))
  vg<-.Fortran("vgram_f",nx=as.integer(nx),ny=as.integer(ny), y=y,r=as.integer(R),newr=x, newl=x)
  return(rbind(vg$newl[(R+1):2,],vg$newr)) # vg$newl[R:2,],
}

# average over directions, adapted from vgram.matrix() in the fields package
vgram_sort <- function(holdVG,nx,ny,r){
    
    holdVG[holdVG==0]=NA  
    holdVG[1:r,1]=NA
    d<-outer(-r:r,0:r, function(x,y) sqrt(x**2 + y**2 ))
    holdN<-outer(-r:r,0:r, function(x,y) (nx-abs(x))*(ny-abs(y)))
   
    holdN[holdVG==0]=NA  
    holdN[1:r,1]=NA
   
         
    top<-tapply(holdVG * holdN,d,FUN="mean",na.rm=TRUE)
    bottom<-tapply(holdN,d,FUN="mean",na.rm=TRUE)
    x<-as.numeric(names(bottom))

    vgram<-top/bottom
    vgram<-vgram[x>0 & x<=r]
    
    out <- list( vgram = vgram, holdVG = holdVG, d = x[x>0 & x<=r], d.full=d, holdN=holdN )
    return(out)
}

# caclculate the stationary score
vg_score_stat <- function( y, dat, R, p=0.5, weighted=TRUE  ){
    nx <- nrow(y)
    ny <- ncol(y)
    
    if( length( dim(dat) )==2 ) dat <- abind( dat, along=3 )
    
    nm <-dim(dat)[3]

    vg_y <- vgram_sort( vgram( y, R=R, p=p ), nx, ny, R )
    
    v_y   <- vg_y$vgram
    v_dat <- 0
    for( i in 1:nm ){ 
        v_dat <- v_dat + vgram_sort( vgram( dat[,,i], R=R, p=p ), nx, ny, R )$vgram
    }
    
    if( weighted ) w <- vg_y$d else w <- 1
    
    res <- sum( ( v_y - v_dat/nm )**2/w )
    return( res )
}




