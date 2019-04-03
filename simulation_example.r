# simulate a few realizations of the models described in the paper
rm( list=ls() )
source( "stochastic_rain.r" )

x <- y <- 1:128
n <- 64

# observations 1: rough and large-scale
obs1 <- simulate_rain( nu=2.5, scale=0.1, x=x, y=y, nf=n ) 

# observations 2: rough and small-scale
obs2 <- simulate_rain( nu=2.5, scale=0.2, x=x, y=y, nf=n ) 


# "forecast" 1: same as obs1
M1  <- simulate_rain( nu=2.5, scale=0.1, x=x, y=y, nf=n ) 

# "forecast" 2: smooth and large-scale
M2  <- simulate_rain( nu=3,  scale=0.1, x=x, y=y, nf=n ) 

# "forecast" 3: same as obs2
M3  <- simulate_rain( nu=2.5, scale=0.2, x=x, y=y, nf=n ) 

# "forecast" 4: smoother and small-scaled 
M4  <- simulate_rain( nu=3,  scale=0.2, x=x, y=y, nf=n ) 


save( obs1, obs2, M1, M2, M3, M4, file="example_fields.rdata" )
