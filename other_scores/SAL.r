# functions to get the original (Wernli et al 2008: "SAL - A novel quality measure for the verification of quantitative precipitation forecasts") and ensmble version (Radanovics et al 2018: "Spatial verification of ensemble precipitation: an ensemble version of SAL") of SAL
# adapted from "SpatialVx"

require( "SpatialVx" )
require( "scoringRules" )

geteSAL <- function( forecast, obs, thf, common_thres=TRUE ){
    if( length(dim(forecast)) == 2 ){
        return( getSAL(forecast, obs, thf, common_thres) )
    }else{
        ne <- dim(forecast)[3]
        feat_list <- list()
        for( i in 1:ne ){
            Vx      <- make.SpatialVx( obs, forecast[,,i])
            th      <- max( obs,na.rm=TRUE ) 
            if( !common_thres ) th <- c( max( forecast[,,i], na.rm=TRUE ), th )
            feat_list[[i]]    <- FeatureFinder( Vx, thresh=th/thf )
        }
        return( esaller( feat_list ) )
    }
}

# modified saller() from SpatialVx to do the ensemble version
esaller <- function (feat_list, d = NULL, distfun = "rdist", ...){
    intRamt <- function(id, x) return(sum(x[id$m], na.rm = TRUE))
    Rmaxer  <- function(id, x) return(max(x[id$m], na.rm = TRUE))
    
    out <- list()
    a <- attributes(feat_list[[1]])
    if (!is.null(a$names)) 
        a$names <- NULL
    attributes(out) <- a
    
    n_ens   <- length(feat_list)
    
    ## get the forecast shit
    Vmod    <- DomRmod <-  0
    cenXhat <- rep(0,2)
    rmod    <- c()
    for( i in 1:n_ens ){
        tmp <- feat_list[[i]]
        y <- tmp$Y.feats
        binY <- im(tmp$Y.labeled)
        binY <- solutionset(binY > 0)
        Y <- tmp$Xhat
        DomRmod <- DomRmod + mean(Y, na.rm = TRUE)/n_ens
        cenXhat <- cenXhat + imomenter(tmp$Xhat)$centroid/n_ens
        RnMod = as.numeric(unlist(lapply(y, intRamt, x = Y)))
        xRmodN = as.numeric(unlist(lapply(y, centdist, y = binY)))
        RmodSum = sum(RnMod, na.rm = TRUE)
        rmod[i] = sum(RnMod * xRmodN, na.rm = TRUE)/RmodSum
        RnMaxMod <- as.numeric(unlist(lapply(y, Rmaxer, x = Y)))
        VmodN <- RnMod/RnMaxMod
        Vmod <- Vmod + sum(RnMod * VmodN, na.rm = TRUE)/(RmodSum*n_ens)
    }
    
    ## get the observation shit
    tmp <- feat_list[[1]]
    x <- tmp$X.feats
    binX <- im(tmp$X.labeled)
    binX <- solutionset(binX > 0)
    X <- tmp$X
    xdim <- dim(X)
    DomRobs <- mean(X, na.rm = TRUE)
    if (is.null(d)) d <- max(xdim, na.rm = TRUE)
    cenX <- imomenter(tmp$X)$centroid
    
    RnObs = as.numeric(unlist(lapply(x, intRamt, x = X)))
    xRobsN = as.numeric(unlist(lapply(x, centdist, y = binX)))
    RobsSum = sum(RnObs, na.rm = TRUE)
    
    numOrig <- sqrt((cenX[1] - cenXhat[1])^2 + (cenX[2] - cenXhat[2])^2)
    RnMaxObs <- as.numeric(unlist(lapply(x, Rmaxer, x = X)))
    VobsN <- RnObs/RnMaxObs
    
    ## calculate scores 
    
    ## S
    Vobs <- sum(RnObs * VobsN, na.rm = TRUE)/RobsSum
    out$S <- 2 * (Vmod - Vobs)/(Vmod + Vobs)
    
    ## A
    A <- 2 * (DomRmod - DomRobs)/(DomRmod + DomRobs)
    out$A <- A
    
    ## L
    L1 = numOrig/d
    
    robs = sum(RnObs * xRobsN, na.rm = TRUE)/RobsSum
    L2 = 2 * abs(rmod - robs)/d
    L2 = 2*crps_sample( y=c(robs/d), dat=rmod/d )
    out$L1 = L1
    out$L2 = L2
    out$L = L1 + L2
    out$L1.alt = NA
    
    class(out) <- "saller"
    return(out)
}

getSAL <- function( forecast, obs, thf, common_thres=TRUE ){
    Vx      <- make.SpatialVx( obs, forecast)
    th      <- max( obs,na.rm=TRUE ) 
    if( !common_thres ) th <- c( max( forecast, na.rm=TRUE ), th )
    feat    <- FeatureFinder( Vx, 
                              thresh=th/thf 
                             )
    return( saller( feat ) )
}
