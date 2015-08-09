#########################################################
# ICPSR "Advanced Maximum Likelihood" 2015 - Survival
#
# Day eight materials.
#
########################################################

library(RCurl)
library(foreign)
library(gtools)
library(boot)
library(plyr)
library(texreg)
library(survival)
library(pscl)
library(smcure)
library(nltm)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

# Cure Models...
#
# Plot:

t<-seq(0.1,50,by=0.1)
pi<-0.5 # P(cured) = 0.5
lambda<-0.1
S.exp<-exp(-(lambda*t))
S.mix<-pi+((1-pi)*(exp(-lambda*t)))
S.nomix<-exp(log(pi)*(1-(exp(-lambda*t))))

pdf("CureTypeS.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t,S.exp,t="l",lwd=3,xlab="Time",ylab="Survival",
     main=expression(paste("Exponential Hazards with ",
          lambda," = 0.1 and ",pi," = 0.5")))
lines(t,S.mix,lwd=3,lty=2,col="red")
lines(t,S.nomix,lwd=3,lty=3,col="blue")
abline(h=0.5,lty=5,col="grey")
legend("topright",bty="n",lwd=3,lty=c(1,2,3),
       col=c("black","red","blue"),
       c("No Cured Fraction","Mixture Cure","Non-Mixture Cure"))
dev.off()

# Simulate some cured-fraction data:

set.seed(7222009)
X<-rnorm(500)
Z<-rbinom(500,1,0.5)
T<-rweibull(500,shape=1.2,scale=1/(exp(0.5+1*X)))
C<-rbinom(500,1,(0.4-0.3*Z)) # Z increases cure probability
S<-Surv(T,C)

pdf("CureSimKM.pdf",10,5)
par(mar=c(4,4,2,2))
plot(survfit(S~1),mark.time=FALSE,lwd=c(3,1,1),
     xlab="Time",ylab="Survival")
dev.off()

coxph(S~X)
coxph(S~X+Z)

data.cure<-cbind(X,Z,T,C)
data.cure<-data.frame(data.cure)
cure.fit<-smcure(S~X,cureform=~Z,data=data.cure,model="ph")

cure.hat<-predictsmcure(cure.fit,c(rep(mean(X),times=2)),
                        c(0,1),model="ph")
cure.pic<-plotpredictsmcure(cure.hat,type="S",model="ph",
            lwd=c(3,3))

pdf("CureHats.pdf",10,5)
par(mar=c(4,4,2,2))
plotpredictsmcure(cure.hat,type="S",model="ph",
                  lwd=c(3,3),main="Predicted Survival")
legend("topright",inset=0.04,bty="n",lty=c(1,2),cex=1.2,
      c("X at mean, Z = 0","X at mean, Z = 1"))
dev.off()


# Fitting real cure models...sorta.

CFURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-git/master/Data/ceasefiresTC.csv"
temp<-getURL(CFRURL)
CF<-read.csv(textConnection(temp))
rm(temp)

CF<-CF[complete.cases(CF),]

CF.S<-Surv(CF$peace,CF$uncensor)

pdf("CF-KM.pdf",10,5)
par(mar=c(4,4,2,2))
plot(survfit(CF.S~1),mark.time=FALSE,lwd=c(3,1,1),
     xlab="Time (in months)",ylab="Survival")
dev.off()

CF.cox<-coxph(CF.S~tie+imposed+lndeaths+contig+onedem+twodem,
              data=CF,method="efron")
CF.cox

CF.cure1.fit<-smcure(CF.S~tie+lndeaths+imposed,
                    cureform=~contig,data=CF,model="ph",
                    link="logit",emmax=500)
# LNB:

LNB<-read.csv("LNB.csv")
LNB<-LNB[complete.cases(LNB),]
LNB<-LNB[order(LNB$dyad,LNB$year),]

LNB.S<-Surv(LNB$count1-1,LNB$count1,LNB$buofmzmid)
LNB.altS<-Surv(LNB$count1,LNB$buofmzmid)

pdf("LNB-KM.pdf",10,5)
par(mar=c(4,4,2,2))
plot(survfit(LNB.S~1,id=dyad,data=LNB),mark.time=FALSE,lwd=c(3,1,1),
     xlab="Time (in years)",ylab="Survival")
legend("bottomright",inset=0.04,bty="n",cex=1.2,
       c("Kaplan-Meier Survival", "Estimate, Long et al. (2007)"))
dev.off()

LNB.cox<-coxph(LNB.S~relcap+major+jdem+border+wartime+s_wt_glo+
               medarb+noagg+arbcom+organ+milinst+cluster(dyad),
               data=LNB,method="breslow")
LNB.cox

# # Runs for a long time:
LNB.cure<-smcure(LNB.altS~relcap+major+jdem+border+wartime+s_wt_glo+
                   medarb+noagg+arbcom+organ+milinst,
                   cureform=~border,model="ph",data=LNB)

# # Also does not work:
# LNB.cure1<-nltm(LNB.S~relcap+major+jdem+border,
#                 nlt.model="PHPHC",data=LNB)

# Stata code for strsmix:
# 
# stset count1, id(episode) f(buofmzmid==1)
# gen h0=0
# strsmix major jdem border wartime, bhazard(h0) distribution(weibull) \\\
#    link(logistic) k1(relcap major jdem border wartime s_wt_glo medarb \\\
#    noagg arbcom organ milinst)
