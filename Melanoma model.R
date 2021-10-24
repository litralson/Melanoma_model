library('fields')
library(gsubfn)

# Lotka-Volterra Growth function
f <- function(x_1,x_2,k_1,k_2,r,alpha)
{
  return(r*x_1*(1-(x_1/k_1)-(alpha*x_2/k_2)) )
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

# Make diffusion coefficient of cancerous cells depend on normal cells
dCoeff_C <- function(x,k,d_max)
{
  return(d_max*(1-(x/k)))
}

# Carcinogenic event
carc_event <- function(xarray_C,xarray_N)
{
  u = sample.int(10,1)
  print(u)
  rn_cells_row = sample.int(50,u)
  rn_cells_col = sample.int(50,u)
  mu = 0.2
  for (i in 1:u) {
      xarray_C[rn_cells_row[i],rn_cells_col[i]] = mu*xarray_N[rn_cells_row[i],rn_cells_col[i]]
      xarray_N[rn_cells_row[i],rn_cells_col[i]] = xarray_N[rn_cells_row[i],rn_cells_col[i]] - xarray_C[rn_cells_row[i],rn_cells_col[i]]
    }
  return(list(xarray_C,xarray_N))
  }

n = 50
# For normal cells
r_N = 1
k_N =10
dCoeff_N = 0
alpha_NC = 1
dr = 0.1

# For cancerous cells 
r_C = 1.5 # Look up the values for logistic function
k_C = 15 
dCoeff_C_max = 1.5 # make it dependent on neighbouring cells
alpha_CN = 1
#Integration related
tmax=200
dt=0.1 #Note that dt < dx^2/2D for the integration scheme to converge! 
dx=1


popCurr_N = matrix(nrow=n+2,ncol=n+2)
popCurr_N = initializeSpace(popCurr_N, k_N)

popCurr_C = matrix(nrow=n+2,ncol=n+2)
popCurr_C = initializeSpace(popCurr_C, 0)  



par(2,1)

popPlot_Total = popCurr_N[2:n+1,2:n+1]+popCurr_C[2:n+1,2:n+1]
image.plot(popPlot_Total,zlim=c(0,(k_N+k_C)))

popPlot_Can = popCurr_C[2:n+1,2:n+1]
image.plot(popPlot_Can,zlim=c(0,(k_C)))


Sys.sleep(1)

# Carcinogenic Event
list[popCurr_C, popCurr_N] = carc_event(popCurr_C,popCurr_N )

popPlot_Total=popCurr_N[2:n+1,2:n+1]+popCurr_C[2:n+1,2:n+1]
image.plot(popPlot_Total,zlim=c(0,(k_N+k_C)))

popPlot_Can = popCurr_C[2:n+1,2:n+1]
image.plot(popPlot_Can,zlim=c(0,(k_C)))

for (t in 2:tmax)
{
  #Numerical integration  
  popPast_N = popCurr_N
  popPast_C = popCurr_C
  
  for (i in 2:(n+1) ) {
    for (j in 2:(n+1) ) {
      popCurr_N[i,j] = (popPast_N[i,j] + dt*f(popPast_N[i,j],popPast_C[i,j],k_N,k_C,r1,alpha_NC) 
                        + dt*dCoeff_N*diffTerm(popPast_N)/(dx*dx) - dr*(popPast_N[i,j])); 
      popCurr_C[i,j] = (popPast_C[i,j] + dt*f(popPast_C[i,j],popPast_N[i,j],k_C,k_N,r2,alpha_CN) 
                        + dt*dCoeff_C(popPast_N[i,j],k_N,dCoeff_C_max)*diffTerm(popPast_C)/(dx*dx)); 
    }
  }
  
 
  #Make the array consistent with periodic boundary conditions.
  popCurr_N = makePeriodic(popCurr_N)
  popCurr_C = makePeriodic(popCurr_C)
  
  popPlot_Total = popCurr_N[2:n+1,2:n+1] + popCurr_C[2:n+1,2:n+1];
  image.plot(popPlot_Total,zlim=c(0,(k_N+k_C)))
  
  popPlot_Can = popCurr_C[2:n+1,2:n+1]
  image.plot(popPlot_Can,zlim=c(0,(k_C)))
  
  #Sleep, so that animation is visible.
  Sys.sleep(0.1);
}

# image.plot(popPlot_Total,zlim=c(0,(k_N+k_C)))