# the full verification experiment from the paper
rm( list=ls() )
source( "verification_functions.r" )
source( "other_scores/vgscore.r" )
source( "other_scores/SAL.r" )

library( "foreach" )
library( "doSNOW" )
library( "abind" )
library( "emdist" )
library( "scoringRules" )
n_cores <- 10
cl   <- makeCluster(n_cores)
registerDoSNOW(cl)

# select parameters
name <- "03"            # name of the input file
fam <- "DaubExPhase"    # wavelet family
num <- 2                # number of vanishing moments
random_th     <- FALSE   # use the randomized thresholds?
correct       <- FALSE   # apply bias correction?
wavelets_only <- TRUE   # do only the wavelet-based scores?
J         <- 5          # largest scale to use
nbr       <- 33         # number of breaks for the scale histogram
RVG       <- 20         # longest lag for the variogram-scores
obs1      <- 1          # first observation
obs2      <- 3          # second observation
path <- "./"
fileend <- name        # end of the output file
if( random_th )           fileend <- paste0( fileend, "_r" )
if( wavelets_only )       fileend <- paste0( fileend, "_wv" )
if( !correct )            fileend <- paste0( fileend, "_raw" )
if( fam!="DaubLeAsymm"  ) fileend <- paste0( fileend, "_EP",num )
if( fam=="DaubLeAsymm" & num!=4  ) fileend <- paste0( fileend, "_LA",num )

# load data
if( random_th ){
    load( file=paste0( path, "RRain_rand", name, ".rdata" ), verbose=TRUE )
    rain <- rain_random_th
    rm( rain_random_th )
    gc()
}else{
    load( file=paste0( path, "RRain_const", name, ".rdata" ), verbose=TRUE )
    rain <- rain_constant_th
    rm( rain_constant_th )
    gc()
}
n_cases <- dim( rain )[1]
nmem    <- dim( rain )[5] - 1
nx      <- dim( rain )[2]
Jpos    <- c( 1:J , 1:J + log2(nx), 1:J + 2*log2(nx) )

# select who gets to be the observation
obs     <- c( rep( obs1, n_cases/2 ), rep( obs2, n_cases/2 ) )

# prepare breaks for the histograms
hbreaks <- seq( 0, J, , nbr )
mids    <- ( hbreaks[1:(nbr-1)] + hbreaks[2:nbr] ) /2

# prepare for parallel calculation
gr      <- expand.grid( case=1:n_cases, mod=1:dim(rain)[4] )
gr$obs  <- obs[gr$case]
n       <- dim(gr)[1]
bar  <- txtProgressBar(0,n,style=3)
prog <- function(i) setTxtProgressBar( bar, i )

res <- foreach( i=1:n, 
                .options.snow=list(progress=prog),
                .errorhandling="stop",
                .packages=c( "LS2W", "emdist", "scoringRules", "SpatialVx", "imager" ),
                .combine=rbind
               ) %dopar%{
    # get observation and forecast
    obs  <- rain[ gr$case[i], , , gr$obs[i], 1]
    forc <- rain[ gr$case[i], , , gr$mod[i], 1:nmem+1]
    if( nmem==1 ) forc <- array( dim=c( dim(forc), 1 ), data=forc )
    
    # normalize
    obs  <- obs/sum(obs)
    for( j in 1:nmem ) forc[,,j] <- forc[,,j]/sum( forc[,,j] )
    
    # transform the observation, get mean spec and scale-histogram
    So   <- fld2S( obs, fam=fam, filnum=num, correct=correct, rsm=0 )[Jpos,,]
    So   <- So*( 1*So>0 )
    So   <- ( So[ 1:J,, ] + So[ 1:J+J,, ] + So[ 1:J+2*J,, ] )/3
    Som  <- rowMeans( So )
    Som <- Som/sum( Som )
    Soc  <- c( fld_centreofmass( So ) )
    ho <- hist( Soc, breaks=hbreaks, plot=FALSE )$density
    
    # transform the forecast, get mean spec and scale-histogram
    Sfm  <- c()
    hf   <- 0
    for( j in 1:nmem ){
        Sf <- fld2S( forc[,,j], fam=fam, filnum=num, correct=correct, rsm=0 )[Jpos,,]
        Sf   <- Sf*( 1*Sf>0 )
        Sf   <- ( Sf[ 1:J,, ] + Sf[ 1:J+J,, ] + Sf[ 1:J+2*J,, ] )/3
        Sfmi <- rowMeans( Sf )
        Sfmi <- Sfmi*( 1*Sfmi>0 )
        Sfmi <- Sfmi/sum( Sfmi )
        Sfm <- rbind( Sfm, Sfmi )
        Sfc <- c( fld_centreofmass( Sf ) )
        hf  <- hf +  hist( Sfc, breaks=hbreaks, plot=FALSE )$density/nmem
    }
    rm( Sf, So )
    gc()
    
    
    # calculate wavelet-based scores
    Hemd  <- emd( cbind( ho, mids ), cbind( hf, mids ) )
    Hcd   <- centreofmass( hf, mids ) - centreofmass( ho, mids )
    if( nmem>1 ){ # ensemble
        SpEn  <- es_sample( y=Som, dat=t(Sfm) )
        scores <- c( SpEn=SpEn, Hemd=Hemd, Hcd=Hcd )
    }else{  # deterministic forecast
        Semd   <- emd( cbind( Som, 1:J ), cbind( c(Sfm), 1:J ) )
        Scd    <- centreofmass( Sfm, 1:J ) - centreofmass( Som, 1:J )
        scores <- c( Semd=Semd, Scd=Scd, Hemd=Hemd, Hcd=Hcd )
    }
    
    if( !wavelets_only ){
        # maybe also get other scores
        VGw_5  <- vg_score_stat( y=obs, dat=forc, R=RVG, p=0.5, weighted=TRUE )
        VG_20  <- vg_score_stat( y=obs, dat=forc , R=RVG, p=2, weighted=FALSE )
        if( nmem>1 ){ # ensemble
            S <- geteSAL( forecast=forc, obs=obs, thf=15, common_thres=TRUE )$S
            scores <- c( scores, Vw5=VGw_5, V20=VG_20, S=S )
        }else{ # deterministic forecast
            S  <- getSAL( forecast=drop(forc), obs=obs, thf=15, common_thres=TRUE )$S
            RMSE <- sqrt( mean( ( drop(forc) - obs )**2 ) )
            scores <- c( scores, Vw5=VGw_5, V20=VG_20, S=S, RMSE=RMSE )
        }
    }
    scores
}
stopCluster(cl)
verif_data  <- data.frame( gr, res )

# check who won, count the number of correct judgements 
winners     <- correct_ans <- array( dim=c( n_cases, dim(res)[2] ), 
                                     data=NA, 
                                     dimnames=list( 1:n_cases, colnames(res) ) )
for( i in 1:n_cases ) for( j in 1:dim(res)[2] ){
    sub              <- res[ which(gr$case==i), j ]
    winners[i,j]     <- which.min( abs(sub) )
    correct_ans[i,j] <- winners[i,j] == obs[i]
}
winners <- data.frame( winners, obs=obs )

correct_ans <- colMeans( correct_ans )
cat( "\n" )
print( correct_ans )

save( verif_data, 
      fam, 
      num ,
      random_th, 
      correct,   
      J,         
      nbr,       
      RVG,  
      correct_ans,
      winners,
      file=paste0( "results/verif",fileend,".rdata" )
      )
