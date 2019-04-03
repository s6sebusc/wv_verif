# Simulation of the stochastic rain fields
# original implementation by Ruediger Hewer, present version for the special case cor(psi, chi, q)=0 by Sebastian Buschow

library('RandomFieldsUtils')
library('RandomFields')


simulate_rain <- function(  nu,          # smoothness paramter 
                            scale,       # size of the features
                            x,y,         # vectors of x and y coordinates
                            nf,          # how many fields do you need?
                            thres=0.8    # fraction of domain without rain
                            ){
    if( nu <= 1) stop( 'nu needs to be greater than one, friend.' )

    # step 1: simulate three independet Matern fields and their first and second derivatives
    RFoptions(cPrintlevel=4,spConform=FALSE,maxGB=6.5)
    # the argument of the modified bessel functions in the Matern covariance is d*( scale/sqrt(2*nu) )
    r     <- scale/sqrt(2*nu)
    
    # the second and third element of RMcurlfree's output are the first derivatives 
    # -> multiplied by r (chain rule)
    # the fourth element contains the second derivatives -> multiply by r**2
    A     <- diag( c( 1, r, r, r**2 ) )
    model <- RMcurlfree( RMmatern(nu=nu), Aniso=A[2:3,2:3] )

    # RMmatrix goes from C to M*C*M^T
    model <- RPcirculant( RMmatrix( model, M=t(A) ), force=TRUE )  
    data  <- RFsimulate( model, x, y, n=3*nf ) 
    data[,,2:3,] <- - data[,,2:3,]
    
    # we interpret those three as humidity, velocity potential and stream function
    # the third dimension of q, chi, psi contains q, dq/dx, dq/dy, div( grad( q ) ) etc.
    q   <- data[ ,,,1:nf ]
    chi <- data[ ,,,1:nf+nf ]
    psi <- data[ ,,,1:nf+2*nf ]

    # step 2: calculate the flux convergence -div( q*( grad(chi) + rot(psi) ) ):
    rain <- -( chi[,,2,] - psi[,,3,] ) * q[,,2,] -  # -( dchi/dx - dpsi/dy ) * dq/dx
             ( chi[,,3,] + psi[,,2,] ) * q[,,3,] -  # -( dchi/dy + dpsi/dx ) * dq/dy
               q[,,1,] * chi[,,4,]                  # - q * div( grad( chi ) )
    
    # step 3: re-order, threshold and return.
    rain <- aperm( rain, c( 3, 1, 2 ) )
    for( i in 1:nf ){
        rain[i,,] <- rain[i,,] - ( quantile( rain[i,,], thres ) )
        rain[i,,] <- rain[i,,] * ( rain[i,,]>0 )
    }
    return(rain)
    
}
