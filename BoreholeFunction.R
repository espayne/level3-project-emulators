#To do the emulation can just follow the other 2D files, which are uploaded for the toy function
#used in chapter 3. (ie. change the design points and dimensionalities of the structures)

borehole <- function(x){
  rw <- x[1]
  r <- x[2]
  Tu <- x[3]
  Hu <- x[4]
  Tl <- x[5]
  Hl <- x[6]
  L <- x[7]
  Kw <- x[8]
  
  frac1 <- 2*pi*Tu*(Hu-Hl)
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu/Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1/frac2
  return(y)
}

partial_borehole_hu <- function(x){
  rw <- x[1]
  r <- x[2]
  Tu <- x[3]
  Hu <- x[4]
  Tl <- x[5]
  Hl <- x[6]
  L <- x[7]
  Kw <- x[8]
  
  frac1 <- 2*pi*Tu
  frac1a <- 2*Tu/(rw^2 * Kw)
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu/Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- (frac1)/frac2
  return(y)
}

partial_borehole_rw <- function(x){
  rw <- x[1]
  r <- x[2]
  Tu <- x[3]
  Hu <- x[4]
  Tl <- x[5]
  Hl <- x[6]
  L <- x[7]
  Kw <- x[8]
  
  frac1 <- 2*pi*Tu*(Hu-Hl)
  frac1a <- (1/rw) + (4*L*Tu/(Kw*rw^3)) + Tu/(Tl*rw)
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu/Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- (frac1*frac1a)/frac2^2
  return(y)
}

#Colour scheme used
require(grDevices)
blue2 <- colorRampPalette(c("cyan", "deepskyblue", "steelblue","dodgerblue4"))


#Borehole Exploratory 
borehole(c(c(0.1, 7.71, 89335, 1050, 89.55, 760, 1400, 10950)))
par(mfrow = c(2,4))

#Change the values of Tu between its range
TuRange <- seq(from = 63070, to = 115600, by = 1)
yTu <- rep(0,length.out = 52531)
for(i in 1:52531){
  yTu[i] <- borehole(c(c(0.1, 7.71, TuRange[i], 1050, 89.55, 760, 1400, 10950)))
}

plot(TuRange, yTu, type = "l", lty = 1, main = "Tu")


#Changing L
LRange <- seq(from = 1120, to = 1680, by = 1)
yL <- rep(0, length.out = 561)
for(i in 1:561){
  yL[i] <- borehole(c(0.1, 7.71, 89335, 1050, 89.55, 760, LRange[i], 10950))
}

plot(LRange, yL, type = "l", lty = 1, main = "L")

#Changing radius of influence 
rRange <- seq(from = 100, 50000, by = 1)
yr <- rep(0, length.out = 49901)
for(i in 1:49901){
  yr[i] <- borehole(c(0.1, log(rRange[i]), 89335, 1050, 89.55, 760, 1400, 10950))
}

plot(rRange, yr, type = "l", lty = 1,main="r")

#Changing Hu
HuRange <- seq(from = 990, 1100, by = 1)
yHu <- rep(0, length.out = 201)
for(i in 1:201){
  yHu[i] <- borehole(c(0.1, 7.71, 89335, HuRange[i], 89.55, 760, 1400, 10950))
}

plot(HuRange, yHu, type = "l", lty = 1, main = "Hu")

#Changing rw, radius of borehole
r_wRange <- seq(from = 0.05, 0.15, by = 0.001)
yr_w <- rep(0, length.out = 101)
for(i in 1:101){
  yr_w[i] <- borehole(c(r_wRange[i], 7.71, 89335, 1050, 89.55, 760, 1400, 10950))
}

plot(r_wRange, yr_w, type = "l", lty = 1, main="r_w")

#Changing Tl; transmissivity of lower aquifier
TlRange <- seq(from = 63.1, 116, by = 0.1)
yTl <- rep(0, length.out = 530)
for(i in 1:530){
  yTl[i] <- borehole(c(0.1, 7.71, 89335, 1050, TlRange[i], 760, 1400, 10950))
}

plot(TlRange, yTl, type = "l", lty = 1, main = "Tl")

#Changing Hl; potentiometric head of lower aquifier
HlRange <- seq(from = 700, 820, by = 1)
yHl <- rep(0, length.out = 121)
for(i in 1:121){
  yHl[i] <- borehole(c(0.1, 7.71, 89335, 1050, 89.55, HlRange[i], 1400, 10950))
}

plot(HlRange, yHl, type = "l", lty = 1, main = "Hl")

#Changing Kw; hydraulic conductivity of borehole
KwRange <- seq(from = 9855, 12045, by = 1)
yKw <- rep(0, length.out = 2191)
for(i in 1:2191){
  yKw[i] <- borehole(c(0.1, 7.71, 89335, 1050, 89.55, 760, 1400, KwRange[i]))
}

plot(KwRange, yKw, type = "l", lty = 1, main = "Kw")

