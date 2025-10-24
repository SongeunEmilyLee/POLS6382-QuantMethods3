## Tobits, PS 206 Class 10

contrib <- read.table("ps206data4.txt", header=T, sep="\t")

contrib <- na.omit(contrib)

contrib <- as.data.frame(contrib)
attach(contrib)

install.packages("AER")
library(AER)

plot(distance,totcont)

## Exploration of tobit with one independent variable

## tobit

tobit1 <- tobit(totcont ~ distance, left=0, data=contrib)
summary(tobit1)

## linear regression w/ censored dv

regress1 <- lm(totcont ~ distance, data=contrib)
summary(regress1)

## linear regression w/ truncated dv

totcont2 <- totcont[totcont>0]
distance2 <- distance[totcont>0]

regress2 <- lm(totcont2 ~ distance2, data=contrib)
summary(regress2)

plot(distance,totcont)
curve(regress1$coeff[1] + regress1$coeff[2]*x, col="black", lwd="2", add=TRUE)
curve(regress2$coeff[1] + regress2$coeff[2]*x, col="orange", lwd="2", add=TRUE)
curve(tobit1$coefficients[1] + tobit1$coefficients[2]*x, col="red", lwd="2", add=TRUE)
legend(3, 12000, legend=c("OLS", "OLS (y>0)", "y*"), col=c("black","orange","red"), lty=1)


## calculating the expected value for y 

XB <- tobit1$linear.predictors
sigma <- tobit1$scale
IMR <- dnorm(XB/sigma)/pnorm(XB/sigma)

E.y <- pnorm(XB/sigma) * (XB + (sigma * IMR))

## calculating the expected value for y when y>0

E.y0 <- XB + (sigma * IMR)


## Comparison of linear regression, latent variable, E(y), E(y | y>0)
plot(distance,totcont)
curve(regress1$coeff[1] + regress1$coeff[2]*x, col="black", lwd="2", add=TRUE)
curve(regress2$coeff[1] + regress2$coeff[2]*x, col="orange", lwd="2", add=TRUE)
curve(tobit1$coefficients[1] + tobit1$coefficients[2]*x, col="red", lwd="2", add=TRUE)
lines(distance,E.y, type="l", col="green", lwd="2")
lines(distance,E.y0, type="l", col="blue", lwd="2")
legend(3, 12000, legend=c("OLS", "OLS (y>0)", "y*", "E(y)", "E(y | y>0)"), col=c("black","orange","red","green","blue"), lty=1)

## The marginal effect of distance on E(y)

Ey.beta <- pnorm(XB/sigma) * tobit1$coefficients[2]
plot(distance,Ey.beta)

## The marginal effect of distance on E(y | y>0)

Ey0.beta <- tobit1$coefficients[2] * (1 - (IMR * ((XB/sigma) + IMR)))
plot(distance,Ey0.beta)

tobit2 <- tobit(totcont ~ republican+senior+vote90+freshman+distance, left=0, data=contrib)
summary(tobit2)

## Interpretation of the influence of the independent variables on the expected values 
## can be difficult because the effect is nonlinear

## one approach is to calculate the mean of the marginal effects, and see how this mean changes as we change one variable
## we could use the means of our variables, or set up a hypothetical case

hypX1 <- c(1,0,0,65,0,0.3)
hypXB1 <- tobit2$coefficients%*%hypX1

# changing our hypothetical case to a Republican

hypX2 <- c(1,1,0,65,0,0.3)
hypXB2 <- tobit2$coefficients%*%hypX2

## other things we need

sigma <- tobit2$scale
IMR1 <- dnorm(hypXB1/sigma)/pnorm(hypXB1/sigma)
IMR2 <- dnorm(hypXB2/sigma)/pnorm(hypXB2/sigma)

Ey.Dem <- pnorm(hypXB1/sigma) * (hypXB1 + (sigma * IMR1))
Ey.Rep <- pnorm(hypXB2/sigma) * (hypXB2 + (sigma * IMR2))

## Difference in the expected contribution for our hypothetical case, switching from Democrat to Republican

diffEy.party <- Ey.Rep - Ey.Dem

Ey.Dem
Ey.Rep
diffEy.party

## we can do the same calculation for the expected contribution conditional on a positive contribution

Ey0.Dem <- hypXB1 + (sigma * IMR1)
Ey0.Rep <- hypXB2 + (sigma * IMR2)

diffEy0.party <- Ey0.Rep - Ey0.Dem

Ey0.Dem
Ey0.Rep
diffEy0.party

## Using Zelig --- for some reason this does not work.  Bug report has been submitted ...

install.packages("Zelig")
library(Zelig)

tobit3 <- zelig(totcont ~ republican+senior+vote90+freshman+distance, below=0, above=Inf, model="tobit", data=contrib)
summary(tobit3)

x.out <- setx(tobit3, republican = 1, senior = 0, vote90 = 65, freshman = 0, distance = 0.3)

sim3 <- sim(tobit3, x = x.out)
summary(sim3)



## class exercise: complex hypothesis tests 

data("Affairs")
attach(Affairs)
colnames(Affairs)
