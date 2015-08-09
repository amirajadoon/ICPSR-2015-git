#########################################################
# ICPSR "Advanced Maximum Likelihood" 2015 - Survival
#
# Day nine materials.
#
########################################################

library(RCurl)
library(RColorBrewer)
library(colorspace)
library(foreign)
library(gtools)
library(statmod)
library(plyr)
library(survival)
library(coxme)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

# Simulate some frailty-type data...

set.seed(7222009)
G<-1:40        # "groups"
F<-rnorm(40)   # frailties
data<-data.frame(cbind(G,F))
data<-data[rep(1:nrow(data),each=20),]
data$X<-rbinom(nrow(data),1,0.5)
data$T<-rexp(nrow(data),rate=exp(0+1*data$X+(2*data$F)))
data$C<-rbinom(nrow(data),1,0.5)
data<-data[order(data$F),]

S<-Surv(data$T,data$C)

Fcolors<-diverge_hcl(length(F))[rank(F)]

pdf("FrailtyKM.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(S~strata(data$G)),col=Fcolors,mark=20,
     xlab="ln(Time)",ylab="Survival",log="x",xlim=c(0.0001,100))
legend("bottomleft",bty="n",cex=0.9,inset=0,
       c("Reds indicate strata","with larger frailties;",
          "blues with smaller ones"))
dev.off()

cox.noF<-coxph(S~X,data=data)
summary(cox.noF)

weib.noF<-survreg(S~X,data=data,dist="weib")
summary(weib.noF)

cox.F<-coxph(S~X+frailty.gaussian(F),data=data)
summary(cox.F)

weib.F<-survreg(S~X+frailty.gaussian(F),data=data,dist="weib")
summary(weib.F)

# Predicted survival plot:
# 
# plot(survfit(cox.noF),log="x",mark.time=FALSE,lwd=c(3,1,1),
#      xlab="ln(Time)",ylab="Fitted Survival")
# lines(survfit(cox.F),col="red",log="x",mark.time=FALSE,lwd=c(3,1,1))

# Examples using leaders data...

leadURL<-"https://raw.githubusercontent.com/PrisonRodeo/ICPSR-git/master/Data/leaders.csv"
temp<-getURL(leadURL)
lead<-read.csv(textConnection(temp))
rm(temp)

lead<-lead[lead$year<2004,]

lead.S<-Surv(lead$tenstart,lead$tenure,lead$tenureend)

Rs<-as.matrix(lead[,13:17])
lead$region<-factor((Rs %*% 1:ncol(Rs))+1,
                    labels=c("NorthAm",colnames(Rs)))
rm(Rs)

lead.F<-coxph(lead.S~female*region+frailty.gamma(ccode),data=lead)
summary(lead.F)

pdf("leadHats.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(lead.F,se.fit=TRUE),conf.int=TRUE,mark.time=FALSE,
     log="x",lwd=c(2,1,1),col="red",xlab="Time (in days)",
     ylab="Survival")
lines(survfit(lead.S~1),conf.int=FALSE,col="black",
      mark.time=FALSE,lwd=2)
legend("bottomleft",bty="n",inset=0.04,lty=1,lwd=3,col=c("red","black"),
       c("Predicted Survival","Baseline (Univariate) Survival"))
dev.off()

# Mixed-effects

lead.coxME<-coxme(lead.S~female + (1 | ccode/female),data=lead)
lead.coxME

# Stratified vs. Frailty, etc.

pdf("lead-KM.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(lead.S~1,id=leadid,data=lead),mark.time=FALSE,lwd=c(3,1,1),
     xlab="Time (in days)",ylab="Survival",log="x")
dev.off()

# Strata by country

pdf("leadKMcountries.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(lead.S~strata(ccode),id=leadid,data=lead),
     col=brewer.pal(9,"Set1"),log="x",mark.time=FALSE,
     xlab="Time (in days)", ylab="Survival")
dev.off()

# Strata by region:

pdf("leadKMregions.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(lead.S~strata(region),id=leadid,data=lead),
     col=brewer.pal(6,"Set1"),lwd=2,log="x",mark.time=FALSE,
     xlab="Time (in days)", ylab="Survival")
legend("bottomleft",inset=0.02,bty="n",col=brewer.pal(6,"Set1"),
       c("N. America","Latin America","Europe","Africa","Asia",
         "Middle East"),lty=1,lwd=2)
dev.off()

lead.Fstrat<-coxph(lead.S~female*strata(region)+
                     frailty.gamma(ccode),data=lead)
summary(lead.Fstrat)

lead.stratCl<-coxph(lead.S~female*strata(region)+
                       cluster(ccode),data=lead)
summary(lead.stratCl)

lead.FstratCl<-coxph(lead.S~female*strata(region)+frailty.gamma(ccode)+
                     cluster(ccode),data=lead)
# boom
