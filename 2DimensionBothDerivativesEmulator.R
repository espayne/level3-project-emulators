#This is similar to the 2D partial derivative emulator
x1 <- seq(-2,2, length = 50)
x2 <- seq(-2,2, length = 50)
dens <- matrix(0,50,50)

for (i in 1:50){
  for (j in 1:50){
    dens[i,j] <- 2*cos(1.2*x1[i]-2) + 3*sin(-0.8*x2[j]+1)
  }
}

require(grDevices) # for colours
#The true contour
filled.contour(x,y, z = dens, color = terrain.colors, asp = NA, xlim = c(-2,2), ylim = c(-2,2), plot.axes={ axis(1, seq(-1,3, by = 1)) ;
  axis(2, seq(-1,4, by = 1))})

#Design points
x1points <- c(0.4919935,-0.434,-0.8231,1.2852674,-1.6543312)
x2points <- c(1.4179598,-0.2974423,0.90911,-0.6044,-1.4505)

#function to emulate: 2cos(1.2x1-2) + 3sin(0.8x2+1)
#x1 partial derivative function: 2.4sin(2-1.2x1)
#x2 partial derivative function -2.4*cos(1-0.8x2)
#vector D in this new case has changed, so first thing to do is to update this

#First thing to do is to get vector D and E(D)
#...............................................................................
#Now we have 3 pieces of data at every design point
Dpart1 <- 2*cos(1.2*x1points-2) + 3*sin(-0.8*x2points + 1)

Dpart2 <- -2.4*sin(1.2*x1points-2)

Dpart3 <- -2.4*cos(1-0.8*x2points) 

D <- c(Dpart1, Dpart2, Dpart3)
#...............................................................................

#Now the next thing is to form the big variance matrix (15 by 15)
#will simulate for each part of the variance matrix first and then put them together
varu <- 0.9
theta <- 1.4

m <- cbind(x1points, x2points)
distmatrix <- as.matrix(dist(m))

varDPart1 <- matrix(nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    varDPart1[i,j] = varu^2 * exp(-(distmatrix[i,j])^2 / (theta)^2)
  }
}

varDPart2 <- matrix(nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    varDPart2[i,j] = varu^2 * (2/(theta^2) - (4/(theta^4)*(x1points[i]-x1points[j])^2))*exp(-(distmatrix[i,j])^2 / (theta^2))
  }
}

varDPart3 <- matrix(nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    varDPart3[i,j] = varu^2 * 2/(theta^2) * (x1points[i]-x1points[j]) * exp(-(distmatrix[i,j])^2 / (theta^2))
  }
}

varDPart5 <-matrix(nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    varDPart5[i,j] = varu^2 * (2/(theta^2) - (4/(theta^4)*(x2points[i]-x2points[j])^2))*exp(-(distmatrix[i,j])^2 / (theta^2))
  }
}

varDPart6 <- matrix(nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    varDPart6[i,j] = varu^2 * 2/(theta^2) * (x2points[i]-x2points[j]) * exp(-(distmatrix[i,j])^2 / (theta^2))
  }
}


varDPart8 <- matrix(nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    varDPart8[i,j] = -1* varu^2 * 4/(theta^4) * (x1points[i]-x1points[j])* (x2points[i]-x2points[j]) *exp(-(distmatrix[i,j])^2 / (theta^2))
  }
}

varDPart4 <- t(varDPart3)
varDPart9 <- t(varDPart8)
varDPart7 <- t(varDPart6)
#Putting this into one big variance matrix
toprow <- cbind(varDPart1, varDPart3, varDPart6)
middlerow <- cbind(varDPart4, varDPart2, varDPart8)
bottomrow <- cbind(varDPart7, varDPart9, varDPart5)
varDmatrix <- rbind(toprow, middlerow, bottomrow)
#...........................................................................

#The next thing to do is to form the covariance structure, again there are three distinct parts
#for this

covmatrixPart1 <- function(est){
  covvec2 <- rep(0,5)
  for(j in 1:5){
    covvec2[j] = varu^2 * exp(-(as.matrix(dist(matrix(c(est, m[j,]), nrow = 2, ncol =2, byrow = TRUE)))[1,2])^2 / (theta)^2)
  }
  return(covvec2)
}

covmatrixPart2 <- function(est){
  covvec2 <- rep(0,5)
  for(j in 1:5){
    covvec2[j] = varu^2 * 2/(theta)^2 *(est[1]-m[j,1])* exp(- (as.matrix(dist(matrix(c(est, m[j,]), nrow = 2, ncol =2, byrow = TRUE)))[1,2])^2 / (theta)^2)
  }
  return(covvec2)
}

covmatrixPart3 <- function(est){
  covvec3 <- rep(0,5)
  for(j in 1:5){
    covvec3[j] = varu^2 * 2/(theta)^2 *(est[2]-m[j,2])* exp(- (as.matrix(dist(matrix(c(est, m[j,]), nrow = 2, ncol =2, byrow = TRUE)))[1,2])^2 / (theta)^2)
  }
  return(covvec3)
}

exp.z.est <- rep(0,2500)
var.z.est <- rep(0,2500)
#Simulation for the Bayes Linear equations
for(i in 1:50){
  for(j in 1:50){
    covvector = c(covmatrixPart1(c(x1[i],x2[j])), covmatrixPart2(c(x1[i],x2[j])), covmatrixPart3(c(x1[i],x2[j])))
    vals = 50 * (i-1) + j
    var.z.est[vals] = varu^2 - (covvector %*% solve(varDmatrix) %*% covvector)
    exp.z.est[vals] = covvector %*% solve(varDmatrix) %*% (D-0)
  }
}


exp.matrix <- matrix(c(exp.z.est), nrow = 50, ncol = 50, byrow = TRUE)
filled.contour(x,y, z = exp.matrix, color = terrain.colors, asp = NA, xlim = c(-2,2), ylim = c(-2,2), plot.axes={ axis(1, seq(-2,2, by = 1)) ;
  axis(2, seq(-2,2, by = 1));
  points(c(xpoints1), c(ypoints1), col = "blueviolet", bg = "blueviolet", pch = 21, cex = 0.7)})


var.matrix <- matrix(c(var.z.est), nrow = 50, ncol = 50, byrow=TRUE)
filled.contour(x,y, z = var.matrix, color = terrain.colors, asp = NA, xlim = c(-2,2), ylim = c(-2,2), plot.axes={ axis(1, seq(-2,2, by = 1)) ;
  axis(2, seq(-2,2, by = 1));
  points(c(xpoints1), c(ypoints1), col = "blueviolet", bg = "blueviolet", pch = 21, cex = 0.7)})
#...................................................................................................................................................
#For the diagnostic plot
diagnostics <- matrix(c(rep(0,2500)), nrow = 50, ncol = 50, byrow = TRUE)

jet.colors <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))

for(i in 1:50){
  for(j in 1:50){
    diagnostics[i,j] = (exp.matrix[i,j] - dens[i,j]) / sqrt(var.matrix[i,j])
  }
}

filled.contour(x, y, diagnostics, color = jet.colors, levels = seq(-3,3,length.out=20),asp = NA, xlim = c(-2,2), ylim = c(-2,2), plot.axes={ axis(1, seq(-2,2, by = 1)) ;
  axis(2, seq(-2,2, by = 1));
  points(c(xpoints1), c(ypoints1), col = "blueviolet", bg = "blueviolet", pch = 21, cex = 0.7)})

