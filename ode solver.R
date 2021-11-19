library(deSolve); library(lattice)

f <- function(t,y,p)
{ x_N = y[1]
x_C = y[2]
with(as.list(p),{
  dx_N.dt = r*x_N*(1-(x_N/k_N)-(alpha_CN*x_C/k_C))-dr*x_N 
  dx_C.dt = r*x_C*(1-(x_C/k_C)-(alpha_NC*x_N/k_N))
  return(list(c(dx_N.dt,dx_C.dt)))})
}

# Normal cell parameters
r_N = 0.5
k_N =10
alpha_NC = 1
dr = 0.01

# Cancer cell parameters
r_C = 0.6 
k_C = 15
alpha_CN = 1

p <- c(r_N=r_N, r_C=r_C, k_N=k_N, k_C=k_C,alpha_NC =alpha_NC, alpha_CN =alpha_CN, dr = dr)
times = seq(0,200,0.1)
y0 <- c(x_N =10, x_C=5)
LV.out <- ode(y=y0,times,f,p)

matplot(LV.out[,1],LV.out[,2:3],type = "l",xlab="time", ylab = "Population size")