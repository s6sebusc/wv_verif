# altered version of "cddews()" from LS2W:
# - removed some of the verbosity, made the code more easily readable
# - removed a small bug where the function would always return "DaubExPhase" as its family
# - now accepts a pre-calculated correction matrix Ainv instead of re-calculating it at every step
# - if Ainv isn't supplied as input, it will attempt to load it from harddrive

mycddews <- function ( data,                    # input field, must be 2^N*2^N
                       filter.number = 1,       # number of the wavelet
                       family = "DaubExPhase",  # family, either "DaubExPhase" or "DaubLeAsymm"
                       Ainv=NULL,               # pre-calculated bias correction matrix
                       correct = TRUE,          # do you want to correct the bias or do you like biases?
                       OPLENGTH=35000,          # amount of memory available for the calculation of A^-1
                       
                       # arguments concerning the smoothing procedure, see ?cddews
                       smooth = TRUE, 
                       sm.filter.number = 4, 
                       sm.family = "DaubExPhase", 
                       levels = 3:6, 
                       type = "hard", 
                       policy = "LSuniversal", 
                       by.level = FALSE, 
                       value = 0, 
                       dev = var
                       ) 
{
    nx <- nrow( data )
    nz <- 3*log2( nx )
    if ( nx != ncol( data ) ) stop( "only squares are allowed." )
    if ( abs( round( nz ) - nz ) > 1e-10  ) stop( "only whole powers of 2 please." )
    
    data.wd <- LS2W:::imwd( data, filter.number = filter.number, family = family, type = "station" )
    RawPer  <- getdata( data.wd, switch = "direction" )
    
    # smooth the raw periodogram via wavelet shrinkage?
    if ( smooth ) {
        if( max( levels ) >= log2( nx )  ) stop( "max(levels) must be <= log2( nrow( data ) )" )
        for (i in 1:nz) {
            # transform each level, apply thresholding ("shrinkage"), transform back
            tmp.imwd   <- LS2W:::imwd( RawPer[i,,], filter.number = sm.filter.number, family = sm.family )
            tmp.imwdTH <- LS2W:::threshold.imwd( tmp.imwd, levels = levels, type = type, policy = policy, value = value, by.level = by.level, dev = dev, compression = FALSE, verbose=FALSE )
            tmp.imwr <- imwr( tmp.imwdTH )
            RawPer[i,,] <- tmp.imwr
        }
    }
    # apply the bias correction following Eckley 2010
    if (correct) {
        if( is.null( Ainv ) ) Ainv <- get_AINV( paste( family, filter.number, sep="_" ), 
                                                N=nx, 
                                                OPLENGTH=OPLENGTH )

        tmp  <- matrix( aperm(RawPer), nrow = nz, ncol = nx**2, byrow = TRUE )
        tmp  <- Ainv %*% tmp
        for ( i in 1:nz ) {
            RawPer[i, , ] <- matrix( tmp[i, ], nrow = nx, ncol = nx, byrow = TRUE )
        }
    }
    
    l <- list( S = RawPer, datadim = dim(data), filter.number = filter.number, 
               family = family, structure = "direction", nlevels = data.wd$nlevels, 
               correct = correct, smooth = smooth, sm.filter.number = sm.filter.number, 
               sm.family = sm.family, levels = levels, type = type, 
               policy = policy, date = date() )
    class(l) <- "cddews"
    return(l)
}
    


# calculate the bias correction matrix for a given  wavelet
get_AINV <- function(wv, N, what="Ainv", OPLENGTH=35000 ){
    wv <- strsplit(wv,split="_")[[1]]
    w  <- paste0( wv, collapse="" )
    mat_name <-  paste( w, N, sep="_" ) 
    matfile  <- paste0("Amats/A_",mat_name,".rdata" )
    if( file.exists(matfile) ){
        load( matfile )
    }else{
        cat( "\nlet me just calculate the matrices for",mat_name,"...\n" )
        A    <- D2Amat( J=-log2(N) , 
                        filter.number = as.numeric(wv[2]), 
                        family = wv[1], 
                        switch = "direction", 
                        verbose = FALSE, 
                        OPLENGTH=OPLENGTH)
        AI <- solve(A)
        save( AI, file=matfile )
    }
    return(AI)
}
 
