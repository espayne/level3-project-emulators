#2D emulator model, here the point used in the Figures can not be used, since they were randomly simulated
#Instead this code is a 20 point emulator, selecting the points from a uniform distribution on x and y
#However to locate the points via a maximin criteria is not very challenging, and the respective code 
#can be found in the report

x <- seq(-1,3, length = 50)
y <- seq(-1,4, length = 50)
dens <- matrix(0,50,50)

for (i in 1:50){
  for (j in 1:50){
    dens[i,j] <- 2*cos(-0.7*x[i]) + sin(1.5*y[j]) 
  }
}

xpoints <- c(runif(20,-1,3))
ypoints <- c(runif(20,-1,4))

z <- c(rep(0,20))
for(i in 1:20){
  z[i] = 2*cos(-0.7*xpoints[i]) + sin(1.5*ypoints[i])
}

#True function contour
require(grDevices) # for colours, only required once
filled.contour(x,y, z = dens, color = terrain.colors, asp = NA, xlim = c(-1,3), ylim = c(-1,4), plot.axes={ axis(1, seq(-1,3, by = 1)) ;
  axis(2, seq(-1,4, by = 1))})

beta0 <- 1
sigmau <- 0.7
theta <- 1

expz <- rep(beta0,20)

#Need to find the differences between each of the points in 2D, 
#not so trivial as in 1D
m <- cbind(xpoints, ypoints)
distmatrix <- as.matrix(dist(m))

varDmatrix <- matrix(nrow = 20, ncol = 20)
for(i in 1:20){
  for(j in 1:20){
    varDmatrix[i,j] = sigmau^2 * exp(- (distmatrix[i,j])^2 / (theta)^2)
  }
}

covvect <- function(est){
  covvec <- rep(0,20)
  for(j in 1:20){
    covvec[j] = sigmau^2 * exp(-(as.matrix(dist(matrix(c(est, m[j,]), nrow = 2, ncol =2, byrow = TRUE)))[1,2])^2 /(theta)^2)
  }
  return(covvec)
}

var.z.est <- rep(0,2500)
exp.z.est <- rep(0,2500)

for(i in 1:50){
  for(j in 1:50){
    vals = 50 * (i-1) + j
    exp.z.est[vals] = beta0 + covvect(c(x[i],y[j])) %*% solve(varDmatrix) %*% (z - expz)
    var.z.est[vals] = sigmau^2 - covvect(c(x[i],y[j])) %*% solve(varDmatrix) %*% covvect(c(x[i],y[j]))
  }
}

exp.matrix <- matrix(c(exp.z.est), nrow = 50, ncol = 50, byrow = TRUE)
var.matrix <- matrix(c((var.z.est)), nrow = 50, ncol = 50, byrow = TRUE)

#pdf("2Demulatorexp.pdf", height = 6, width = 7)
filled.contour(x,y, z = exp.matrix, color = terrain.colors, asp = NA, xlim = c(-1,3), ylim = c(-1,4), plot.axes={ axis(1, seq(-1,3, by = 1)) ;
  axis(2, seq(-1,4, by = 1));
  points(c(xpoints), c(ypoints), col = "blueviolet", bg = "blueviolet", pch = 21, cex = 0.7)})
#dev.off()

#pdf("2Demulatorvar.pdf", height = 6, width = 7)
filled.contour(x,y, z = var.matrix, color = terrain.colors, asp = NA, xlim = c(-1,3), ylim = c(-1,4), plot.axes={ axis(1, seq(-1,3, by = 1)) ;
  axis(2, seq(-1,4, by = 1));
  points(c(xpoints), c(ypoints), col = "blueviolet", bg = "blueviolet", pch = 21, cex = 0.7)})
#dev.off()