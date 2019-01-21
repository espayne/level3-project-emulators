#This code runs specifically for the emulation shown in Figure 1
#it is clear how to change this code for the parameters and points where 
#the value of the true function is taken.

x <- c(seq(0.1,1,length.out = 10))
D <- rep(0,10)
for(i in 1:10){
  D[i] = 0.5*exp(-3*x[i]^2)*sin(12*x[i])
}

beta0 <- 0
sigmau <- 0.06
theta <- 0.1

expD <- rep(beta0,10)

varDmatrix <- matrix(nrow = 10, ncol = 10)
for(i in 1:10){
  for(j in 1:10){
    varDmatrix[i,j] = sigmau^2 * exp(-(abs(x[i] - x[j]))^2 / (theta)^2)
  }
}

x.vals = c(seq(from = 0, to = 1.1, length.out = 1000))

covmatrix <- matrix(nrow = 1000, ncol = 10)
for(i in 1:1000){
  for(j in 1:10){
    covmatrix[i,j] = sigmau^2 * exp(-(abs(x.vals[i] - x[j]))^2 / (theta)^2)
  }
}

exp.x.vals = c(rep(0,1000))
var.x.vals = c(rep(0,1000))

for(i in 1:1000){
  exp.x.vals[i] = beta0 + covmatrix[i,] %*% solve(varDmatrix) %*% (D - expD)
  var.x.vals[i] = sigmau^2 - (covmatrix[i,]%*%solve(varDmatrix) %*% covmatrix[i,])
}
#pdf("1dimensionemulator.pdf", height = 6, width = 7)
plot(x, D, xlim = c(0.06, 1.04), ylim = c(-0.5,0.5), xlab = "x", 
     ylab = "f(x)", 
     col = "blueviolet", bg = "blueviolet", pch = 21, cex = 1.5)
lines(x.vals, exp.x.vals, col = "blue")
lines(x.vals, exp.x.vals + (3*sqrt(var.x.vals)), col = "red")
lines(x.vals, exp.x.vals - (3*sqrt(var.x.vals)), col = "red")
lines(x.vals, 0.5 * exp(-3* x.vals^2) * sin(12*x.vals), col = "black")
legend(0.6, -0.2, legend=c("True function", "Emulator expectation"), col=c("black", "blue"), lty=1:1, cex=0.8)
#dev.off()
