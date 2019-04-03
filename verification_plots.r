# make the plots for the paper
rm( list=ls() )
graphics.off()

require( "gridExtra" )
require( "ggridges" )
require( "ggplot2" )

fileend <- "_nst01"        # end of the input file
thecols <- gray.colors( 4, .7, 0.2 )

load( paste0( "results/verif",fileend,".rdata" ), verbose=TRUE )
scores   <- colnames( verif_data )[-(1:3)]
modnames <- c( "RL", "SmL", "RS", "SmS" )

verif_data <- reshape( verif_data, varying=scores, times=scores, v.names="value",direction="long" )

# figure out, how frequently the correct model was rewarded
winners <- as.data.frame( winners )
sv <- ov <- ms <- hr <- c()
for( o in unique( verif_data$obs ) ) for( s in scores ){
    dat <- subset( verif_data, obs==o&time==s&mod==o )
    ms  <- c( ms, mean(dat$value ) )
    hr  <- c( hr, mean( winners[[s]][ unique( dat$case ) ] == o ) )
    sv  <- c( sv, s )
    ov  <- c( ov, modnames[o] )
}
label_data <- data.frame( ms=ms, obs=ov, time=sv, hr=hr )
winners    <- reshape( winners, varying=scores, times=scores, v.names="value",direction="long" )


# give the models their rightful names
verif_data$obs <- sapply( verif_data$obs, function(x) modnames[x] )
verif_data$mod <- sapply( verif_data$mod, function(x) modnames[x] )
winners$obs   <- sapply( winners$obs, function(x) modnames[x] )
winners$value <- sapply( winners$value, function(x) modnames[x] )

# determine the order of things in the plot
verif_data <- transform( verif_data, 
                         time=factor( time, levels=scores ), 
                         mod=factor( mod, levels=modnames) )
winners <- transform( winners, 
                      time=factor( time, levels=scores ), 
                      value=factor( value, levels=modnames) )
                         
# ridge plot of the scores
p1 <- ggplot( verif_data, aes( y=mod, x=value, col=mod, fill=mod  ) ) + 
      scale_y_discrete( breaks=modnames, labels=modnames ) + 
      geom_density_ridges( show.legend=FALSE, size=.2 ) + 
      facet_grid( paste0( "observations: ",obs ) ~ time, scales="free" ) + 
      scale_color_manual( limits=modnames, values=thecols ) + 
      scale_fill_manual(  limits=modnames, values=thecols ) + 
      labs( x="", y="" )

# bar plot of the winning forecasts
p2 <- ggplot( winners, aes( x=time ) ) + 
      facet_grid( .~paste0( "observations: ",obs ) ) + 
      scale_y_continuous( breaks=seq(1, dim(winners)[1]/(2*length(scores)) ,,5)[-1], 
                          labels=paste0( seq(25,100,25),"%" ) ) + 
      scale_fill_manual(  limits=modnames, values=thecols ) +
      stat_count( geom="bar", aes( fill=value ) ) + 
      coord_flip() +
      labs( x="", y="", fill="" )
      
      
# add labels with "hit-rates"
labtext <- geom_text( size = 3, 
                      data = label_data, 
                      mapping = aes( x=ms, y=obs, label=hr ),
                      hjust = .5,
                      vjust = 0,
                      inherit.aes=FALSE ) 

# add lines at zero (where appropriate)
zeroline <- geom_vline( data=subset( verif_data, value<0 ), 
                        aes( xintercept = 0 ),
                        linetype="dotted" )

# basic look of the plot
style <- theme_minimal() + 
         theme( panel.border = element_rect(colour = "lightgray", 
                                            fill=NA, size=1), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.y = element_text(hjust=0)
                )
# skewed axis
skew_axis <- theme( axis.text.x = element_text(angle=45, vjust=.5, hjust=.5) )
         
# save plots
cairo_pdf( width=9, height=4, file=paste0( "plots/ridges_", fileend, ".pdf" ) )                        
    print( p1 + labtext + zeroline + style + skew_axis )
dev.off()

cairo_pdf( width=7, height=3, file=paste0( "plots/bars_", fileend, ".pdf" ) )                        
    print( p2 +  style )
dev.off()



