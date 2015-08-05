***DAY 4***


**Data is OR.dta on github website**

**Logit examples

gen duration = stop

logit dispute allies contig capratio growth democ trade duration

**fourth order polynomial trend


gen duration2 = duration^2*0.1

gen duration3 = duration^3*0.01

gen duration4 = duration^4*0.001

logit dispute allies contig capratio growth democ trade duration duration2 duration3 duration4

testparm durat*

**Time dummy example

xi i.duration

logit dispute allies contig capratio growth democ trade _Iduration_2 - _Iduration_34 _Iduration_35

testparm _I*

**Cox/Poisson Equivalence

stcox allies contig capratio growth democ trade, nohr

xi: poisson dispute allies contig capratio growth democ trade i.duration

