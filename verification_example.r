rm( list=ls() )
source( "verification_functions.r" ) 
load( "example_fields.rdata" ) # output of simulation_example.r
require( "emdist" )

obs <- obs1 # obs1 has the parameters of M1, obs2 has those of M3

nt <- dim(obs)[1]
J  <- log2( dim(obs)[2] )

# normalize all fields to unit total intensity -> verify structure only
for( i in 1:nt ){
    obs[i,,] <- obs[i,,]/sum( obs[i,,] )
    M1[i,,]  <- M1[i,,]/sum( M1[i,,] )
    M2[i,,]  <- M2[i,,]/sum( M2[i,,] )
    M3[i,,]  <- M3[i,,]/sum( M3[i,,] )
    M4[i,,]  <- M4[i,,]/sum( M4[i,,] )
}

# do the transform, average over space and direction, remove negative energy, normalize
getmeanspec <- function( x ){
    res <- rowMeans( fld2S( x, fam="DaubLeAsymm", filnum=4, smooth=FALSE ) )
    res <- ( res[1:J] + res[1:J+J] + res[1:J+2*J] )/3 
    res <- res*(1*res>0)
    res <- res/sum(res)
    return( res )
}

# transform all the observations and forecasts, calculate the EMD between their spectra
SPemd <- array( dim=c( nt, 4 ) )
for( i in 1:nt ){
    So <- getmeanspec( obs[i,,] )
    for( f in 1:4 ){
        Sf <- getmeanspec( get( paste0("M",f) )[i,,] )
        SPemd[i,f] <- emd( cbind( So, 1:J ), cbind( Sf, 1:J ) )
    }
}

# count how frequently each model "won"
winners   <- apply(SPemd, c(1), which.min)
win_count <- c()
for( i in 1:4 ) win_count[i] <- paste0( "\n M",i,": ", length( which( winners==i ) ), "/", nt )

print( "number of cases where each forecast received the best score:" , quote=FALSE )
cat( win_count, "\n" )
 
