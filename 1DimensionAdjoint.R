#vector D in this new case has changed, so first thing to do is to update this

#First thing to do is to get vector D and E(D)
#...............................................................................
x <- c(5,7.5,10,12.5,15,17.5,20)

Dpart1 <- function(m){
  len <- length(m)
  Dvals <- rep(0,len)
  for(i in 1:len){
    Dvals[i] = x[i]/3 + cos(x[i])
  }
  return(Dvals)
}

Dpart2 <- function(m){
  len <- length(m)
  Dvals <- rep(0,len)
  for(i in 1:len){
    Dvals[i] = 1/3 - sin(x[i])
  }
  return(Dvals)
}

D <- c(Dpart1(x), Dpart2(x))
expD <- c(4,4,4,4,4,4,4,0,0,0,0,0,0,0)
#...............................................................................
#Now the next thing is to form the big variance matrix
#I know this will be an 8x8 matrix with my setup, and will require different
#subsections, which I need to calculate and then put them together
varu <- 1.2
theta <- 2

varDPart1 <- matrix(nrow = 7, ncol = 7)
for(i in 1:7){
  for(j in 1:7){
    varDPart1[i,j] = varu^2 * exp(-(abs(x[i] - x[j]))^2 / (theta)^2)
  }
}

varDPart2 <- matrix(nrow = 7, ncol = 7)
for(i in 1:7){
  for(j in 1:7){
    varDPart2[i,j] = varu^2 * (2/(theta^2) - 4/(theta^4)*(x[i]-x[j])^2)*exp(-(x[i]-x[j])^2 / (theta^2))
  }
}

varDPart3 <- matrix(nrow = 7, ncol = 7)
for(i in 1:7){
  for(j in 1:7){
    varDPart3[i,j] = varu^2 * 2/(theta)^2 * (x[i] - x[j]) * exp(- (x[i] - x[j])^2 / (theta)^2)
  }
}

varDPart4 <- t(varDPart3)
#Putting this into one big variance matrix
tophalf <- cbind(varDPart1, varDPart3)
bottomhalf <- cbind(varDPart4, varDPart2)
varDmatrix <- rbind(tophalf, bottomhalf)
#.............................................................................

x.vals <- seq(0,23, length.out = 1000)
covmatrixPart1 <- matrix(nrow = 1000, ncol = 7)
for(i in 1:1000){
  for(j in 1:7){
    covmatrixPart1[i,j] = varu^2 * exp(-(x.vals[i] - x[j])^2/(theta)^2) 
  }
}

covmatrixPart2 <- matrix(nrow = 1000, ncol = 7)
for(i in 1:1000){
  for(j in 1:7){
    covmatrixPart2[i,j] = varu^2 * 2/(theta)^2 * (x.vals[i] - x[j]) * exp(- (x.vals[i] - x[j])^2 / (theta)^2)
  }
}

covmatrixPart3 <- matrix(nrow = 1000, ncol = 7)
for(i in 1:1000){
  for(j in 1:7){
    covmatrixPart3[i,j] = - varu^2 * 2/(theta)^2 * (x.vals[i] - x[j]) * exp(- (x.vals[i] - x[j])^2 / (theta)^2)
  }
}

covmatrixPart4 <- matrix(nrow = 1000, ncol = 7)
for(i in 1:1000){
  for(j in 1:7){
    covmatrixPart4[i,j] = varu^2 * (2/(theta^2) - 4/(theta^4)*(x.vals[i]-x[j])^2)*exp(-(x.vals[i]-x[j])^2 / (theta^2))
  }
}

covmatrix <- cbind(covmatrixPart1,covmatrixPart2)

#.............................................................................
#Emulating the derivative of the function, only using the output values
exp.deriv <- rep(0,1000)
var.deriv <- rep(0,1000)
expD <- c(4,4,4,4,4,4,4)

for(i in 1:1000){
  exp.deriv[i] = covmatrixPart3[i,] %*% solve(varDPart1) %*% (Dpart1(x) - expD)
  var.deriv[i] =  0.72 - (covmatrixPart3[i,]%*%solve(varDPart1)%*%covmatrixPart3[i,])
}

plot(x, xlim = c(4,21), ylim = c(-5,5), xlab = "x", ylab = "f'(x)", 
     col = "blueviolet", bg = "blueviolet", pch = 21, cex = 1.5)
lines(x.vals, exp.deriv, col = "blue", lty = 1)
lines(x.vals, exp.deriv + (3*sqrt(var.deriv)), col = "red")
lines(x.vals, exp.deriv - (3*sqrt(var.deriv)), col = "red")
lines(x.vals, 1/3 - sin(x.vals), col = "black")
abline(v = 10, col = 5, lty = 2)
abline(v = 20, col =5, lty = 2)
abline(v=15, col = 5, lty = 2)
abline(v = 5, col = 5, lty = 2)
abline(v = 17.5, col = 5, lty = 2)
abline(v=12.5, col = 5, lty = 2)
abline(v = 7.5, col = 5, lty = 2)
#.................................................................................
#Emulating the derivative of the function, with the full data vector

exp.deriv.both <- rep(0,1000)
var.deriv.both <- rep(0,1000)
expD <- c(4,4,4,4,4,4,4,0,0,0,0,0,0,0)
covmatrix2 <- cbind(covmatrixPart3,covmatrixPart4)

for(i in 1:1000){
  exp.deriv.both[i] = covmatrix2[i,] %*% solve(varDmatrix) %*% (D - expD)
  var.deriv.both[i] = 0.72 - (covmatrix2[i,]%*%solve(varDmatrix)%*%covmatrix2[i,])
}
plot(x, xlim = c(4,21), ylim = c(-10,10), xlab = "x", ylab = "f'(x)", 
     col = "blueviolet", bg = "blueviolet", pch = 21, cex = 1.5)
lines(x.vals, exp.deriv.both, col = "blue")
lines(x.vals, exp.deriv.both + (3*sqrt(var.deriv.both)), col = "red")
lines(x.vals, exp.deriv.both - (3*sqrt(var.deriv.both)), col = "red")
abline(v = 10, col = 5, lty = 2)
abline(v = 20, col =5, lty = 2)
abline(v=15, col = 5, lty = 2)
abline(v = 5, col = 5, lty = 2)
abline(v = 17.5, col = 5, lty = 2)
abline(v=12.5, col = 5, lty = 2)
abline(v = 7.5, col = 5, lty = 2)
lines(x.vals, 1/3 - sin(x.vals), col = "black")
