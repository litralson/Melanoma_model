library('fields')

# Growth function
f <- function(x,k,r)
{
  return(r*x*(1-x/k) )
}

 

#To ensure periodic boundary conditions
makePeriodic <- function(xarray)
{
  
  for (j in 2:(n+1)) {
    xarray[1,j] = xarray[n+1,j]
    xarray[n+2,j] = xarray[2,j]
  }
  
  for (i in 2:(n+1)) {
    xarray[i,1] = xarray[i,n+1]
    xarray[i,n+2] = xarray[i,2]
  }
  
  return(xarray)
}



# Initialize the grid
initializeSpace <-function(xarray,k)
{
  for (i in 2:(n+1)) {
    for (j in 2:(n+1)) {
      xarray[i,j] = k; 
    }     
  }
  xarray=makePeriodic(xarray)
}

#Diffusion, based on finite difference method.
diffTerm <- function(xarray)
{
  diffx = xarray[i-1,j]-2*xarray[i,j]+xarray[i+1,j]
  diffy = xarray[i,j-1]-2*xarray[i,j]+xarray[i,j+1]
  return(diffx+diffy) 
}

# Carcinogenic event
carc_event <- function(xarray)
{
  u = sample.int(10,1)
  print(u)
  rn_cells_row = sample.int(50,u)
  rn_cells_col = sample.int(50,u)
  for (i in 1:u) {
      xarray[rn_cells_row[i],rn_cells_col[i]] = 4
    }
  return(xarray)
  }

n = 50
# For normal cells
r1 = 1
k1 =10
dCoeff_1 = 1

# For cancerous cells 
r2 = 1.5 # Look up the values for logistic function
k2 = 15 
dCoeff_2 = 1.5 # make it dependent on neighbouring cells

#Integration related
tmax=200
dt=0.1 #Note that dt < dx^2/2D for the integration scheme to converge! 
dx=1


popCurr_1 = matrix(nrow=n+2,ncol=n+2)
popCurr_1 = initializeSpace(popCurr_1, k1)

popCurr_2 = matrix(nrow=n+2,ncol=n+2)
popCurr_2 = initializeSpace(popCurr_2, 0)  

popPlot=popCurr_1[2:n+1,2:n+1]+popCurr_2[2:n+1,2:n+1]

image.plot(popPlot,zlim=c(0,(k1+k2)))

Sys.sleep(1)

popCurr_2 = carc_event(popCurr_2)
popPlot=popCurr_1[2:n+1,2:n+1]+popCurr_2[2:n+1,2:n+1]
image.plot(popPlot,zlim=c(0,(k1+k2)))

for (t in 2:tmax)
{
  #Numerical integration  
  popPast_1 = popCurr_1
  popPast_2 = popCurr_2
  
  for (i in 2:(n+1) ) {
    for (j in 2:(n+1) ) {
      popCurr_1[i,j] = (popPast_1[i,j] + dt*f(popPast_1[i,j],k1,r1) 
                        + dt*dCoeff_1*diffTerm(popPast_1)/(dx*dx)); 
      popCurr_2[i,j] = (popPast_2[i,j] + dt*f(popPast_2[i,j],k2,r2) 
                        + dt*dCoeff_2*diffTerm(popPast_2)/(dx*dx)); 
    }
  }
  
 
  #Make the array consistent with periodic boundary conditions.
  popCurr_1 = makePeriodic(popCurr_1)
  popCurr_2 = makePeriodic(popCurr_2)
  popPlot = popCurr_1[2:n+1,2:n+1] + popCurr_2[2:n+1,2:n+1];
  
  image.plot(popPlot,zlim=c(0,(k1+k2)))
  
  #Sleep, so that animation is visible.
  Sys.sleep(0.1);
}

# image.plot(popPlot,zlim=c(0,(k1+k2)))