#function to emulate: x/3 + cos(x)
#derivative function: 1/3 - sin(x)
#vector D in this new case has changed, so first thing to do is to update this

#First thing to do is to get vector D and E(D)
#...............................................................................
x <- c(5,10,15,20)

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
expD <- c(4,4,4,4,0,0,0,0)
#...............................................................................
#Now the next thing is to form the big variance matrix
#I know this will be an 8x8 matrix with my setup, and will require different
#subsections, which I need to calculate and then put them together
varu <- 1
theta <- 3

varDPart1 <- matrix(nrow = 4, ncol = 4)
for(i in 1:4){
  for(j in 1:4){
    varDPart1[i,j] = varu^2 * exp(-(abs(x[i] - x[j]))^2 / (theta)^2)
  }
}

varDPart2 <- matrix(nrow = 4, ncol = 4)
for(i in 1:4){
  for(j in 1:4){
    varDPart2[i,j] = varu^2 * (2/(theta^2) - 4/(theta^4)*(x[i]-x[j])^2)*exp(-(x[i]-x[j])^2 / (theta^2))
  }
}

varDPart3 <- matrix(nrow = 4, ncol = 4)
for(i in 1:4){
  for(j in 1:4){
    varDPart3[i,j] = varu^2 * 2/(theta)^2 * (x[i] - x[j]) * exp(- (x[i] - x[j])^2 / (theta)^2)
  }
}

varDPart4 <- t(varDPart3)
#Putting this into one big variance matrix
tophalf <- cbind(varDPart1, varDPart3)
bottomhalf <- cbind(varDPart4, varDPart2)
varDmatrix <- rbind(tophalf, bottomhalf)
#.............................................................................

#The next thing to do is to form the covariance structure, which is likely to be a little bit trickier
#Like with the D stuff, I need to do a half and half kind of thing and put together
x.vals <- seq(0,23, length.out = 1000)
covmatrixPart1 <- matrix(nrow = 1000, ncol = 4)
for(i in 1:1000){
  for(j in 1:4){
    covmatrixPart1[i,j] = varu^2 * exp(-(x.vals[i] - x[j])^2/(theta)^2) 
  }
}

covmatrixPart2 <- matrix(nrow = 1000, ncol = 4)
for(i in 1:1000){
  for(j in 1:4){
    covmatrixPart2[i,j] = varu^2 * 2/(theta)^2 * (x.vals[i] - x[j]) * exp(- (x.vals[i] - x[j])^2 / (theta)^2)
  }
}

covmatrix <- cbind(covmatrixPart1,covmatrixPart2)
#..............................................................................

exp.x.vals <- rep(0,1000)
var.x.vals <- rep(0,1000)

for(i in 1:1000){
  exp.x.vals[i] = 4 + covmatrix[i,] %*% solve(varDmatrix) %*% (D - expD)
  var.x.vals[i] = varu^2 - (covmatrix[i,]%*%solve(varDmatrix)%*%covmatrix[i,])
}

pdf("1Dwithderiv.pdf", height = 6, width = 7)
plot(x, Dpart1(x), xlim = c(4,21), ylim = c(0,8), xlab = "x", ylab = "f(x)", 
     col = "blueviolet", bg = "blueviolet", pch = 21, cex = 1.5)
lines(x.vals, exp.x.vals, col = "blue")
lines(x.vals, exp.x.vals + (3*sqrt(var.x.vals)), col = "red")
lines(x.vals, exp.x.vals - (3*sqrt(var.x.vals)), col = "red")
lines(x.vals, x.vals/3 + cos(x.vals), col = "black")
dev.off()

#.............................................................................
#Taking just the value of the 4 points, without the derivative information
exp.x4.vals <- rep(0,1000)
var.x4.vals <- rep(0,1000)
expD4 <- c(4,4,4,4)

for(i in 1:1000){
  exp.x4.vals[i] = 4 + covmatrixPart1[i,] %*% solve(varDPart1) %*% (Dpart1(x) - expD4)
  var.x4.vals[i] = varu^2 - (covmatrixPart1[i,]%*%solve(varDPart1)%*%covmatrixPart1[i,])
}

pdf("1Dwithoutderiv.pdf", height = 6, width = 7)
plot(x, Dpart1(x), xlim = c(4,21), ylim = c(0,8), xlab = "x", ylab = "f(x)", 
     col = "blueviolet", bg = "blueviolet", pch = 21, cex = 1.5)
lines(x.vals, exp.x4.vals, col = "blue")
lines(x.vals, exp.x4.vals + (3*sqrt(var.x4.vals)), col = "red")
lines(x.vals, exp.x4.vals - (3*sqrt(var.x4.vals)), col = "red")
lines(x.vals, x.vals/3 + cos(x.vals), col = "black")
dev.off()
