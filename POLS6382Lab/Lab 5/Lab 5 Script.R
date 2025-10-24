############################################
#POLS6382 Quantitative Methods III   
#Lab 5. Models for Ordinal Data                        
#Ling Zhu, University of Houston           
#Last Update:October 3, 2021			         
############################################
rm(list = ls())
setwd("/Users/lingzhu/Dropbox/UH Teaching/POLS6382_2021 Fall/2021 Labs/Lab 5 Models for Ordinal Data")
library(foreign)
library(ggplot2)
library(gmodels)
library(Hmisc)
library(MASS)
library(ordinal)
library(reshape)
library(stats)

#1. Data and Variables for the Ordered Logit Example
healthcare<-read.dta("gss2012lab5.dta")
attach(healthcare)
# Histogram
#pdf(file="response.pdf",height=6, width=6)
hist(hlthctzen,
     main="Access to public funded health care if one is not a citizen?",
     xlab="Responses"
     )
#dev.off()

# Tab frequencies
data.frame(table(hlthctzen))

# Table proportions
prop.table(table(hlthctzen))

# Histogram with percentages
h = hist(hlthctzen)
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE)

#marginals
prop.table(table(hlthctzen))

#2. Estimate an ordered logit model
model1<-polr(as.factor(hlthctzen)~female+white+age+educ+income+ideology,
             data=healthcare,Hess=TRUE,model=TRUE)
summary(model1)
require(stargazer)
stargazer(model1, type="text")

#store table
(ctable<-coef(summary(model1)))
# get p values
p<-pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
#combined table
(ctable <- cbind(ctable, "p value" = p))

#3. Interpretation
# Odds ratios
exp(coef(model1))

# Mean predicted probabilities, white respondents
predictdata<-cbind(ideology=seq(1,7,length=100),
                   age=mean(age),income=mean(income),educ=mean(educ),
                   female=mean(female),white=1)
opinion.hat<-predict(model1,predictdata,type='prob')
# Mean predicted probabilities, non-white respondents
predictdata2<-cbind(ideology=seq(1,7,length=100),
                   age=mean(age),income=mean(income),educ=mean(educ),
                   female=mean(female),white=0)
opinion.hat2<-predict(model1,predictdata2,type='prob')

# Plot the pps
ideology<-seq(1,7,length=100)
#pdf(file="pp_white.pdf",height=7, width=7)
plot(c(1,7),c(0,1),type='n',
     xlab="Liberal-Conservative Ideology Scale",
     ylab="Predicted Probabilities (y=j)",
     main="Access to public funded health care if one is not a citizen?")
lines(ideology,opinion.hat[1:100,1],lty=1,lwd=3,col="red")
lines(ideology,opinion.hat[1:100,2],lty=2,lwd=3,col="blue")
lines(ideology,opinion.hat[1:100,3],lty=3,lwd=3,col="green")
legend(1,1.1,cex=0.3,c('Disagree','Neither','Agree'),
       lty=1:3,col=c("red","blue","green"))
#dev.off()

#pdf(file="pp_nonwhite.pdf",height=7, width=7)
plot(c(1,7),c(0,1),type='n',
     xlab="Liberal-Conservative Ideology Scale",
     ylab="Predicted Probabilities (y=j)",
     main="Access to public funded health care if one is not a citizen?")
lines(ideology,opinion.hat2[1:100,1],lty=1,lwd=3,col="red")
lines(ideology,opinion.hat2[1:100,2],lty=2,lwd=3,col="blue")
lines(ideology,opinion.hat2[1:100,3],lty=3,lwd=3,col="green")
legend(1,1.1,cex=0.3,c('Disagree','Neither','Agree'),
       lty=1:3,col=c("red","blue","green"))
#dev.off()

# Compare white and non-white respondents' PPs of "disagree"
#pdf(file="compare.pdf",height=7, width=7)
plot(c(1,7),c(0,1),type='n',
     xlab="Liberal-Conservative Ideology Scale",
     ylab="Predicted Probabilities (y=Disagree)",
     main="Access to public funded health care if one is not a citizen?")
lines(ideology,opinion.hat[1:100,1],lty=1,lwd=3,col="red")
lines(ideology,opinion.hat2[1:100,1],lty=2,lwd=3,col="blue")
legend(1,1,cex=0.9,c('White','Non-White'),
       lty=1:3,col=c("red","blue"))
#dev.off()

# Compare white and non-white respondents' PPs of "agree"
#pdf(file="compare2.pdf",height=7, width=7)
plot(c(1,7),c(0,1),type='n',
     xlab="Liberal-Conservative Ideology Scale",
     ylab="Predicted Probabilities (y=Agree)",
     main="Access to public funded health care if R is not a citizen?")
lines(ideology,opinion.hat[1:100,3],lty=1,lwd=3,col="red")
lines(ideology,opinion.hat2[1:100,3],lty=2,lwd=3,col="blue")
legend(1,1,cex=0.9,c('White','Non-White'),
       lty=1:3,col=c("red","blue"))
#dev.off()

