#Bifurcation diagram
require(rootSolve)
require(shape)

f <- function(x_1,x_2,k_1,k_2,r,alpha)
{
  return(r*x_1*(1-(x_1/k_1)-(alpha*x_2/k_2)) )
}

#This function calculates equilibrium by solving f=0 (using RootSolve package)
#It then calculates stability at equilibrium by first computing f'(x) and then
#evaluating f' at the equilibria.
equilibrium <- function(hi)
{
  Eq <- uniroot.all(f= f,interval=c(0,K),hi=hi)
  
  #Stability of the system.
  eqtype <- vector(length=length(Eq))
  for (i in 1:length(Eq))
  {  
    fprime    <- gradient(f = rate, x=Eq[i], hi=hi) #Find derivative at the root.
    eqtype[i] <- sign(fprime)+2
  }     
  
  return(list(x=Eq, type=eqtype))
}

hseq <- seq(0.0,himax,by=0.01);
h_soln <- numeric(length(hseq))
i <- 1

#Open the plot and plot V=0 root.
plot.new()
plot(0,xlim=range(hseq),ylim=c(0,K),type="n",
     xlab="Harvesting rate (h)",ylab="V*",main="Harvesting Model")

#Loop over grazing rate (h) to calculate equilibrium, stability for each and plot them.
for (hi in hseq) 
{  
  eq  <- equilibrium(hi);
  
  points(rep(hi,length(eq$x)),eq$x,pch='.',
         col=c("black","red","lightgrey")[eq$type],
         bg =c("black","red","lightgrey")[eq$type]) 
  h_soln[i] <- length(eq$type)
  i <- i+1
  
}
