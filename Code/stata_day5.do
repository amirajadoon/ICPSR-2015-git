* ICPSR Advanced MLE -- 2015 
*
* Day Five Stata code
*
* Note: Results are very slightly different than when using
* R due to differences in algorithms.

insheet using "scotus.csv", clear

summarize

stset service, failure(retire) id(justice)

* Cox model

stcox age pension pagree, nohr efron

* log-log survival plots

xtile agecat = age, n(4)

stphplot, by(age)
stphplot, by(pension)
stphplot, by(pagree)

* Martingale and schoenfeld residuals - refit model

predict mgres, mgale
predict sch*, schoenfeld
predict sca*, scaled

* Test for non-PH

estat phtest, detail

* Plots:

estat phtest, plot(age) yline(0)
estat phtest, plot(pension) yline(0)
estat phtest, plot(pagree) yline(0)

* log-Time interaction

gen lnT = ln(service)
gen ageXlnT = age * lnT

stcox age pension pagree ageXlnT, nohr efron

* Tests:

estat phtest, detail 

*******************************************
* Unobserved heterogeneity

* Simulation:

clear

set obs 500
set seed 7222009
gen W = invnorm(uniform())
gen X = invnorm(uniform())
gen Z = invnorm(uniform())
gen T = (-1/(exp(0+0.5*W+0.5*X-0.6*Z))) * ln(1-uniform()) 
gen C = 1

stset T, failure(C)

* Models:

streg X, nohr dist(weib)
mat p = e(aux_p)
svmat double p, names(p1)
replace p11 = p11[_n-1] if p11[_n-1]!=.
streg X W, nohr dist(weib)
mat p = e(aux_p)
svmat double p, names(p2)
replace p21 = p21[_n-1] if p21[_n-1]!=.
streg X W Z, nohr dist(weib)
mat p = e(aux_p)
svmat double p, names(p3)
replace p31 = p31[_n-1] if p31[_n-1]!=.

* Plot hazards for lambda = 0.02:

range Tsim 1 60
gen Exphat = 0.02*1*((0.02*Tsim)^(1-1))
gen S1hat = 0.02*p11*((0.02*Tsim)^(p11-1))
gen S2hat = 0.02*p21*((0.02*Tsim)^(p21-1))
gen S3hat = 0.02*p31*((0.02*Tsim)^(p31-1))

twoway (line S1hat Tsim, sort lcolor(dkgreen) lwidth(medthick)) /*
 */ (line S2hat Tsim, lcolor(blue) lwidth(medthick) lpattern(dash)) /*
 */ (line S3hat Tsim, lcolor(cranberry) lwidth(medthick) lpattern(shortdash)) /*
 */ (line Exphat Tsim, lcolor(black) lwidth(medthick) /*
 */ lpattern(longdash_dot)), ytitle(Weibull Hazard) xtitle(Time) /*
 */ legend(order(1 "One covariate" 2 "Two covariates" 3 /*
 */ "Correct specification" 4 "True hazard" ) region(lwidth(none)))

* Weibull model on SCOTUS data:

insheet using "scotus.csv", clear
stset service, failure(retire) id(justice)

streg age pension pagree, nohr dist(weibull)

* Add -age- to the shape parameter:

streg age pension pagree, nohr dist(weibull) anc(age)
mat B = e(b)
svmat B, names(Bs)
replace Bs5 = Bs5[_n-1] if Bs5[_n-1]!=.
replace Bs6 = Bs6[_n-1] if Bs6[_n-1]!=.

* Subsequent plot:

range Agesim 32 91
gen Phat = exp(Bs6 + Bs5*Agesim)

twoway (line Phat Agesim, sort lcolor(black) lwidth(medthick)), /*
 */ ytitle(Estimate of p) yline(1, lwidth(medium) lpattern(dash) /*
 */ lcolor(gs8)) xtitle(Age) legend(off)
