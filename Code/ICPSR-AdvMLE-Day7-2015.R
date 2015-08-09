#########################################################
# ICPSR "Advanced Maximum Likelihood" 2015 - Survival
# 
# Day seven materials.
#
########################################################

library(RCurl)
library(foreign)
library(gtools)
library(plyr)
library(texreg)
library(survival)
library(logistf)
library(coxphf)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places


# Table

Yeas<-t(c(rep(0,times=212),rep(1,times=219)))
Dems<-t(c(rep(0,times=178),rep(1,times=253)))
table(Yeas,Dems)

# Simulated Logits:

set.seed(7222009)
X<-runif(100,min=-5,max=5)
X<-X[order(X)]
Z<-runif(100,min=-5,max=5)
Y<-ifelse(plogis(X+Z)>0.5,1,0)
Y2<-ifelse(plogis(X+0.5*Z)>0.5,1,0)
Y3<-ifelse(plogis(X+0.1*Z)>0.5,1,0)
Ysep<-ifelse(plogis(X)>0.5,1,0)
Yfit<-glm(Y~X,family="binomial")
Y2fit<-glm(Y2~X,family="binomial")
Y3fit<-glm(Y3~X,family="binomial")
Ysepfit<-glm(Ysep~X,family="binomial")

# Plots:

pdf("Separation.pdf",8,7)
par(mar=c(4,4,2,2))
par(mfrow=c(2,2))
plot(X,Y,pch=19,xlab="X",ylab="Y")
lines(X,plogis(predict(Yfit)),lwd=3)
legend("topleft",inset=0.04,bty="n",cex=1.2,
       legend=c(paste("Beta =", round(Yfit$coefficients[2],digits=2)),
                paste("SE =", round(sqrt(vcov(Yfit))[4],digits=2))))
plot(X,Y2,pch=19,xlab="X",ylab="Y")
lines(X,plogis(predict(Y2fit)),lwd=3)
legend("topleft",inset=0.04,bty="n",cex=1.2,
       legend=c(paste("Beta =", round(Y2fit$coefficients[2],digits=2)),
                paste("SE =", round(sqrt(vcov(Y2fit))[4],digits=2))))
plot(X,Y3,pch=19,xlab="X",ylab="Y")
lines(X,plogis(predict(Y3fit)),lwd=3)
legend("topleft",inset=0.04,bty="n",cex=1.2,
       legend=c(paste("Beta =", round(Y3fit$coefficients[2],digits=2)),
                paste("SE =", round(sqrt(vcov(Y3fit))[4],digits=2))))
plot(X,Ysep,pch=19,xlab="X",ylab="Y")
lines(X,plogis(predict(Ysepfit)),lwd=3)
legend("topleft",inset=0.04,bty="n",cex=1.2,
      legend=c(paste("Beta =", round(Ysepfit$coefficients[2],digits=2)),
            paste("SE =", round(sqrt(vcov(Ysepfit))[4],digits=2))))
dev.off()

# Toy data:

rm(X,Y,Z)
set.seed(7222009)
Z<-rnorm(500)
W<-rnorm(500)
Y<-rbinom(500,size=1,prob=plogis((0.2+0.5*W-0.5*Z)))
X<-rbinom(500,1,0.5)
X<-ifelse(Y==0,0,X)

summary(glm(Y~W+Z+X,family="binomial"))

data<-as.data.frame(cbind(W,X,Y,Z))
write.dta(data,"SepSim.dta") # for the Stata illustration

# Pets data:

PetsURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-git/master/Data/Pets.csv"
temp<-getURL(PetsURL)
Pets<-read.csv(textConnection(temp))
rm(temp)

Pets.1<-glm(petfamily~female+as.factor(married)+as.factor(partyid)
            +as.factor(education),data=Pets,family=binomial)

Pets.2<-glm(petfamily~female+as.factor(married)*female+as.factor(partyid)+
              as.factor(education),data=Pets,family=binomial)

with(Pets, xtabs(~petfamily+as.factor(married)+female))

Pets.Firth<-logistf(petfamily~female+
                    as.factor(married)*female+as.factor(partyid)+
                    as.factor(education),data=Pets)

# Separation in Survival Data...

rm(X)
set.seed(7222009)
X<-rep(0:1,times=100)
T<-abs(rweibull(200,shape=1.2,scale=1/(exp(0+0.2*X))))
C<-rbinom(200,1,0.2)
C<-ifelse(X==0,0,C)

table(C,X)

cox.fit<-coxph(Surv(T,C)~X,method="efron")
summary(cox.fit)

weib.fit<-survreg(Surv(T,C)~X,dist="weibull")
summary(weib.fit)

SIM<-cbind(T,C,X)
SIM<-data.frame(SIM)
firth.fit<-coxphf(SIM,formula=Surv(T,C)~X)
firth.fit

pdf("HeinzePlot.pdf",6,5)
par(mar=c(4,4,4,2))
coxphfplot(SIM,formula=Surv(T,C)~X,profile=~X)
dev.off()

# LHR:

LHRURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-git/master/Data/LHR.csv"
temp<-getURL(LHRURL)
LHR<-read.csv(textConnection(temp))
rm(temp)

keep<-c("id","ccode1","ccode2","year","LHRcluster","X_t0","X_t","X_d",
        "archigosFIRC","archigosFIRClnt","capchange","battletide",
        "thirdpartycfire","index","onedem5","twodem5","twodem5lnt",
        "tie","lndeaths","cfhist","stakes","contiguity","contiguitylnt")
LHR<-LHR[keep]
LHR<-LHR[complete.cases(LHR),]
summary(LHR)

LHR.S<-Surv(LHR$X_t0,LHR$X_t,LHR$X_d)

pdf("LHR-KM.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(LHR.S~1),mark.time=FALSE,lwd=c(3,1,1),
     xlab="Time (in days)",ylab="Survival")
dev.off()

LHR.Cox<-coxph(LHR.S~archigosFIRC+archigosFIRClnt,data=LHR,method="efron",
               iter.max=10000)

table(LHR$archigosFIRC,LHR$X_d)

LHR.CoxF<-coxphf(LHR.S~archigosFIRC+archigosFIRClnt,data=LHR,maxit=1000)

# Rescaling T to years doesn't change anything:

LHR$TstartY<-LHR$X_t0/365.25
LHR$TendY<-LHR$X_t/365.25
LHR.SY<-Surv(LHR$X_t0,LHR$X_t,LHR$X_d)
LHR$archigosFIRClnTY<-LHR$archigosFIRC*(log(LHR$TendY))

LHR.CoxY<-coxph(LHR.SY~archigosFIRC+archigosFIRClnTY,data=LHR,method="efron",
                 iter.max=10000)

LHR.CoxFY<-coxphf(LHR.SY~archigosFIRC+archigosFIRClnTY,data=LHR,maxit=1000)

LHR.CY.out<-extract(LHR.CoxY,include.rsquared=F,include.maxrs=F,
                      include.zph=F,include.nobs=F,
                      include.missings=F)
# LHR.CFY.out<-extract(LHR.CoxFY,include.rsquared=F,include.maxrs=F,
#                     include.zph=F,include.nobs=F,
#                     include.missings=F)
texreg(l=list(LHR.CY.out),
       file="LHR-Years-Table.tex",
       custom.model.names=c("Cox"),
       custom.coef.names=c("FIRC","FIRC x ln(T)"),
       stars=numeric(0),caption="",label="",table=F)

# Then go fix the table by-hand (because there is not yet
# a texreg method for coxphf...)
