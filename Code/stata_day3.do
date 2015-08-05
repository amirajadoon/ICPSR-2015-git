***DAY 3***

**Data is OR.dta on github website**

**Cox model, breslow method

stcox allies contig capratio growth democracy trade, nohr breslow bases(basicS) basehc(basich)

**Cox model, hazard ratios
stcox

**Predicted survival curves for contiguous=0 and contiguous=1

stcurve, surv at1(contig=0) at2(contig=1)

**Cox model, robust/sandwich variance estimates

stcox allies contig capratio growth democracy trade, nohr efron robust  

stcox allies contig capratio growth democracy trade, nohr efron robust cluster(dyadid)

