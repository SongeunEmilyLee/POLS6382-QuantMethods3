############################################
#POLS6382 Quantitative Methods III:MLE   
#Lab 4. Binary Choice Models: Extensions                       
#Ling Zhu, University of Houston           
#Last Update: September 26, 2021  			         
############################################
# Laptop
rm(list=ls())
setwd("/Users/lingzhu/Dropbox/UH Teaching/POLS6382_2021 Fall/2021 Labs/Lab 4 Binary Choice Models 2")

library(AER)
library(ggplot2)
library(glmx)
library(heatmapFit)
library(lmtest)
library(maxLik)
library(reshape2)
library(Zelig)
library(MASS)
library(stats)

#1. Hetero Probit Model With Simulated Data
# Simulating DGP
set.seed(48)
n<- 500
x<- rnorm(n)
ystar<-1+x+rnorm(n, sd = exp(x))
y<-ifelse(ystar>0, 1, 0)
simdata<-data.frame(cbind(y,ystar,x))

# visualization
# Visual 1. Scatter plot of simulated x and latent y.
visual1<-ggplot(data=simdata)+
  geom_point(aes(x=x, y=ystar))+
  labs(x="Simulated X (500 Draws from Standard Normal)",
       y="Latent Values of Y")+
  theme_light()
# Visual 2. Scatter plot of simulated x and binary y.
visual2<-ggplot(data=simdata)+
  geom_point(aes(x=x, y=y))+
  labs(x="Simulated X (500 Draws from Standard Normal)",
       y="Simulated Outcomes")+
  theme_light()
# Visual 3. Histogram of the binary dependent variable
visual3<-ggplot(data=simdata)+
  geom_histogram(aes(x=y), binwidth = 0.2, color="black", fill="gray")+
  labs(x="Simulated Binary Outcome Y",
       y="Count")+
  theme_light()
require(ggpubr)
#pdf(file="simdata.pdf", height=5, width=15)
ggarrange(visual1,visual2, visual3, ncol=3, nrow=1)
#dev.off()

# Model1: Probit model 
model1 <- glm(y ~ x, family = binomial(link = "probit"),data=simdata)
# Model2:heteroskedastic model 
model2 <- hetglm(y ~ x | x, family=binomial(link="probit"),data=simdata)
require(stargazer)
stargazer(model1, model2, type="text")

# LR Test
lrtest(model1,model2)

# Model Fit
y<-simdata$y
pred1<-predict(model1,type="response")
pred2<-predict(model2,type="response")

table(true = y,
      predicted = pred1 <= 0.5)

require(heatmapFit)
#pdf(file="heat1",height=6, width=6)
heatmap.fit(y,pred1,reps=1000,legend=FALSE)
#dev.off()
#pdf(file="heat2",height=6, width=6)
heatmap.fit(y,pred2,reps=1000,legend=FALSE)
dev.off()

# 2. Heteroskedastic Probit Model With Real Data
# Labor Force Participation Example from Package AER
load("PSID1976.rda")
#data("PSID1976", package = "AER")
PSID1976$kids <- with(PSID1976, factor((youngkids + oldkids) > 0,
                levels = c(FALSE, TRUE), labels = c("no", "yes")))
PSID1976$fincome <- PSID1976$fincome/10000

#Probit model using glm()
labormodel1<- glm(participation ~ age + I(age^2) + 
              fincome + education + kids,
              data = PSID1976, family = binomial(link = "probit"))

# Heteroskedastic probit model via hetglm() with constant scale
labormodel2 <- hetglm(participation ~ age + I(age^2) + 
                fincome + education + kids | 1,
                data = PSID1976)

#Heteroskedastic Probit model with varying scale by # of kids and income
heterlabor<- hetglm(participation ~ age + I(age^2) + 
              fincome + education + kids | kids + fincome,
              data = PSID1976)

# Likelihood ratio test
lrtest(labormodel2, heterlabor)
stargazer(heterlabor, labormodel2, type="text")
summary(heterlabor)

# Check for predictions
table(true = PSID1976$participation,
      predicted = fitted(labormodel2) <= 0.5)
table(true = PSID1976$participation,
      predicted = fitted(heterlabor) <= 0.5)

# 3. Post-estimation simulation with Zelig
data(turnout,package="Zelig")
# Estimate the model
z.out<-zelig(vote~race+educate+age+I(age^2)+income,
             model="logit",data=turnout)
# Set explanatory variables 
x.low<-setx(z.out, educate=12, age=18:95)
x.high<-setx(z.out, educate=16, age=18:95)
#Simulate quantities of interest:
s.out<-sim(z.out, x=x.low, x1=x.high)
#pdf(file="zelig1.pdf",height=6, width=8)
plot(s.out, xlab="Age in Years",
        ylab="Predited Probability of Voting",
        xlim=c(20,90))
#dev.off()

# 4. Analyzing Rare Events Data with Zelig
data(mid,package="Zelig")
reventmod<-zelig(conflict~major+ contig+power+maxdem+mindem+years, data=mid,
                model="relogit",tau=1042/303772, 
                case.control="prior", bias.correct=TRUE)
# summarize the model output
summary(reventmod)
# Let variable ``power" vary, and hold all other variables constant.
x.out1<-setx(reventmod,power=0:1)
#Simulate quantity of interests
s.out1<-sim(reventmod,x=x.out1)
#pdf(file="rarezelig.pdf",height=6, width=6)
plot(s.out1,ci=95)
#dev.off()

#5. More on simulations: visually weighted logit regression lines
#Simulate some data
n <- 1000
simdata2<- data.frame(X = rnorm(n),
            Female = sample(c(0, 1), n, replace = T))
simdata2$Y <- with(simdata2, X * 2 + Female * 10 + rnorm(n, sd = 5))
simdata2$Y <- (simdata2$Y > 1) * 1
head(simdata2)

# Model:
myModel <- glm(Y ~ X + Female, data = simdata2, family = "binomial")

# Generate predictions
nSims <- 1000
someScenarios <- expand.grid(1,  # Intercept
              seq(min(simdata2$X), max(simdata2$X), len = 100),  # A sequence covering the range of X values
              c(0, 1))  # The minimum and maximum of a binary variable, or some other first difference

simDraws <- mvrnorm(nSims, coef(myModel), vcov(myModel))
simYhats <- plogis(simDraws %*% t(someScenarios))# t():matrix transpose function
dim(simYhats)  # Simulated predictions

# Combine scenario definitions(gender id) and predictions.
predictionFrame <- data.frame(someScenarios, t(simYhats))
# Reshape wide -> long (for ggplot2)
longFrame <- melt(predictionFrame, id.vars = colnames(someScenarios))
head(longFrame)

#pdf(file="weighted.pdf",height=6, width=8)
p1<-ggplot(data = longFrame,
      aes(x = Var2, y = value, group = paste(variable, Var3), 
      colour = factor(Var3)))+
      geom_line(alpha = I(1/sqrt(nSims)))+
      scale_x_continuous("Explanatory Variable X", expand = c(0, 0))+
      scale_y_continuous("Predicted Pr(Y=1)", limits = c(0, 1), expand = c(0, 0))+
      scale_colour_brewer(palette="Set1", labels = c("Male", "Female"))+  # Change colour palette
      guides(colour = guide_legend("Group ID", override.aes = list(alpha = 1)))+  # Avoid an alpha-related legend problem
      ggtitle("The Effect of X on Y by Gender Group")+
      theme_bw()
print(p1)  # This might take a few seconds...
#dev.off()

# Lab 5 ends here

##Appendix: Code to Create a 1:1 Sample
data(mid,package="Zelig")
# create a subset only for cases if conflict=1
subset1<-mid[which(mid$conflict=='1'),]
# create a subset only for cases if conflict=0
subset2<-mid[which(mid$conflict=='0'),]
# creat a sample from subset2, randomly draw 1042 observations
sample1<-subset2[sample(1:nrow(subset2),1042),]
# create a new data frame, with a 1:1 sample design for 1s and 0s
newdata<-rbind(subset1, sample1)
