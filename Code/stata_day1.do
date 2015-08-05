**Day 1


**In Stata, all survival analysis commands require that we first indicate that the data are 
**“survival-time” data. We do this by “stset-ing” the data.

**The basic command is:

stset durat

**If some observations are censored, then we use:


stset durat, failure(censor)

**If not all the observations in the data entered at the same time, 
**then it is (sometimes) important to record this fact; we do so by indicating, 
**via two variables, when the observation entered and exited the data:

stset durat, failure(eltime2) en(timein) ex(timeout)

**Finally, the id() option is used to indicate that you have time-varying data
**(and possibly that observations experience more than one event):

stset durat, failure(censor) en(timein) ex(timeout) id(id)


****Example using KABL data****

stset durat, failure(ciep12) id(id)

**Cabinet durations (Kaplan-Meier) plot with 95% Greenwood ci's

sts graph, gwood 

**Nelson-Aalen cumulative hazard plot

sts graph, cumhaz ci

**Log-rank test

sts test invest

**Comparing Survival estimates when invest = 0 and invest = 1

sts graph, gwood by(invest)


