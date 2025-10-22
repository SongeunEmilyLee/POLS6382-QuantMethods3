############################################
#POLS6382 Quantitative Research Methods III   
#Lab 7. Tobit and Heckman Selection Model                        
#Ling Zhu, University of Houston           
#Last Update: October 17, 2021  			         
############################################
rm(list = ls())
setwd("/Users/lingzhu/Dropbox/UH Teaching/POLS6382_2021 Fall/2021 Labs/Lab 6 Censored and Truncated Data")
library(ggplot2)
library(GGally)
library(censReg)
library(sampleSelection)
library(VGAM)
library(mvtnorm)

#1. Comparing Tobit and OLS with Simulated Data
N = 10
f = rep(c("s1","s2","s3","s4","s5","s6","s7","s8"),N)
fcoeff = rep(c(-1,-2,-3,-4,-3,-5,-10,-5),N)
set.seed(100) 
x = rnorm(8*N)+1
beta = 5
epsilon = rnorm(8*N,sd = sqrt(1/5))
y.star = x*beta+fcoeff+epsilon ## latent response
y = y.star 
y[y<0] <- 0 ## lefted censored response
simdata<-data.frame(cbind(x,y))

fitols<-lm(y~x)
summary(fitols)
fittobit<-vglm(formula=y~x, family=tobit(Lower=0))
summary(fittobit)
coef(fittobit) # More satisfying estimates

# Compare two regression lines
#pdf(file="olstobit.pdf", width=6, height=6)
plot(x,y)
abline(lm(y~x),col="red",lwd=2,lty=1)
curve(-3.3862557 + 4.5046179 *x, col="orange", lwd="2", add=TRUE)
#dev.off()

# 2.Tobit Model for Left and Right Censored Data
mydata<-read.csv("tobit.csv")
attach(mydata)
summary(mydata$apt)

f <- function(x, var, bw = 15) {
  dnorm(x, mean = mean(var), sd(var)) * length(var)  * bw
}
pdf(file="apt.pdf",height=6, width=8)
p<-ggplot(mydata, aes(x = apt, fill=prog))+theme_bw()
p
p + stat_bin(binwidth=15) +
  stat_function(fun = f, size = 1,
  args = list(var = mydata$apt))
dev.off()

# Check correlations
cor(mydata[, c("read", "math", "apt")])
# Plot matrix
pdf(file="correlation.pdf",height=8, width=8)
ggpairs(mydata[, c("read", "math", "apt")])
dev.off()

# OLS 
olsmodel<-lm(formula=apt~read+math+as.factor(prog))
summary(olsmodel)

# Tobit 
tobitmodel<-vglm(formula=apt~read+math+as.factor(prog), 
                 family=tobit(Upper=800))
tobitmodel<-censReg(apt~read+math+as.factor(prog),right=800,data=mydata)
summary(tobitmodel)

# Consider an alternative specification
tobitmodel2<-vglm(apt ~ read + math, tobit(Upper = 800), 
                  data = mydata)
tobitmodel2<-censReg(apt ~ read + math, right=800, left=200)
# LRT to compare models
(p <- pchisq(2 * (logLik(tobitmodel) - logLik(tobitmodel2)), 
             df = 2, lower.tail = FALSE))


#3. Heckman Sample Selection Models: Simulated Data
# With exclusion restriction
set.seed(0)
eps <- rmvnorm(500, c(0, 0), 
     matrix(c(1, -0.7, -0.7, 1), 2, 2))
# selection
xs <- runif(500)
ys <- xs + eps[, 1] > 0
# outcome
xo <- runif(500)
yoX <- xo + eps[, 2] # Latent variable
yo <- yoX * (ys > 0)
summary(selection(ys ~ xs, yo ~ xo))

# Without exclusion restriction
yoX <- xs + eps[, 2]
yo <- yoX * (ys > 0)
summary(selection(ys ~ xs, yo ~ xs))

# 4. Selection Model
data("Mroz87")
Mroz87$kids <- (Mroz87$kids5 + Mroz87$kids618 > 0)

# Heckman Two-Step Method
selectmod1 <- selection(lfp ~ age + I(age^2) + faminc + kids + educ,
  +   wage ~ exper + I(exper^2) + educ + city, data = Mroz87, 
    method = "2step")
summary(selectmod1)

# ML Estimation
selectmod2 <- selection(lfp ~ age + I(age^2) + faminc + kids + educ,
          +wage ~ exper + I(exper^2) + educ + city, data = Mroz87,
          maxMethod = "BHHH", iterlim = 500) 
summary(selectmod2)

# Lab 7 ends here. 