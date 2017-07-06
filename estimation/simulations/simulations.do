
/* Simulation 1: QTTs and Baseline Outcomes */
clear
sim_panel, n(10000) t(2) cov(0)

moments if t==2

* Estimate Average Treatment Effect at Quantiles of Y0
xtile q0 = y0 if t==2, nq(100) 
matrix treatOnQ0 = J(99,1,.)
forv k=1/99 {
 qui su effect if q0==`k'
 matrix treatOnQ0[`k',1] = r(mean)
} 

* Estimate Average Treatment Effect at Quantiles of Baseline Y
xtile q0_pre = y if t==1, nq(100) 
matrix treatOnQ0_pre = J(99,1,.)
forv k=1/99 {
 qui su effect if q0_pre==`k'
 matrix treatOnQ0_pre[`k',1] = r(mean)
} 
keep if t==2
 
* Quantiles of Treatment effects
quantiles effect if t==2
matrix qOft = e(quantiles)
quantiles y treatment if t==2
matrix qtt = e(quantiles)
svmat qtt 
svmat qOft
svmat treatOnQ0
svmat treatOnQ0_pre

gen quantiles = _n if qtt1<.


* Figure 2.
graph twoway (line qtt1 quantiles, lpattern(longdash)) (line qOft1 quantiles, lcolor(black) lwidth(thick)) ///
 (line treatOnQ0_pre1 quantiles, lcolor(black) lpattern(solid)), ///
 graphregion(fcolor(white)) legend(label(1 "QTT") label(2 "Quantiles of F{sub:{&Delta}}") ///
	label(3 "ATE at Quantile of Baseline Y") ) xtitle("Quantile") ///
  ytitle("Treatment Effect") ylabel(,nogrid) note("Note. Based on 10,000 observations drawn from Simulation 1. See appendix 1 for details.")
graph export "${results_dir}\sim_fig2.svg", replace       

    




* Variance of Treatment Effects 
qui su y if treatment==1
local sd1 = r(sd)
qui su y if treatment==0
local sd0 = r(sd)
matrix sdD = J(40,2,.)
local row = 0
forv r=-1(.05)1 {
  local ++row
  matrix sdD[`row',1] = `r'
  matrix sdD[`row',2] = sqrt(`sd1'*`sd1'+`sd0'*`sd0' - 2 * `r' * `sd1' * `sd0')
 }


/* Makarov Bound Simulation */
clear
sim_makarov, n(2000)
graph twoway (kdensity y3 if treat==1, lcolor(black)) ///
 (kdensity y3 if treat==0, lcolor(black) lpattern(longdash)), ///
	graphregion(fcolor(white)) ylabel(,nogrid) xtitle("Outcomes") ytitle("Density") legend(label(1 "F{sub:1},{&sigma}{sub:{&Delta}}=5") label(2 "F{sub:0}"))
	graph export "${results_dir}\sim_fig3_right.svg", replace
	
	
quantiles y1 treat, nq(1000)
matrix qtt1 = e(quantiles)
quantiles y2 treat, nq(1000)
matrix qtt2 = e(quantiles)
quantiles y3 treat, nq(1000)
matrix qtt3 = e(quantiles)

quantiles effect1, nq(1000)
matrix qOft1 = e(quantiles)
quantiles effect2, nq(1000)
matrix qOft2 = e(quantiles)
quantiles effect3, nq(1000)
matrix qOft3 = e(quantiles)


* Makarov Bounds
keep y1 y2 y3 treat
sort treat 
matrix y = J(1000,6,.)
forv k=1/1000 {
 matrix y[`k',1] = y1[`k']
 matrix y[`k',2] = y1[`k'+1000]
 matrix y[`k',3] = y2[`k']
 matrix y[`k',4] = y2[`k'+1000]
 matrix y[`k',5] = y3[`k']
 matrix y[`k',6] = y3[`k'+1000] 
} 
clear
svmat y
matrix makbounds = J(201,7,.)


local r=0
forv d = -10(0.1)10 {
 local ++r
 qui makarov y2 y1 , grid(-10(0.1)10) d(`d')
 matrix makbounds[`r',1] = `d'
 matrix makbounds[`r',2] = r(bound_l)
 matrix makbounds[`r',3] = r(bound_u) 
 qui makarov y4 y3 , grid(-10(0.1)10) d(`d')
 matrix makbounds[`r',4] = r(bound_l)
 matrix makbounds[`r',5] = r(bound_u) 
 qui makarov y6 y5 , grid(-10(0.1)10) d(`d')
 matrix makbounds[`r',6] = r(bound_l)
 matrix makbounds[`r',7] = r(bound_u)  
}

 
svmat makbounds 
matrix inv_makbounds = J(1000,7,.)
forv r=1/1000 {
 local p=`r'/1000
 matrix inv_makbounds[`r',1] = `p'
 qui su makbounds1 if makbounds2>=`p'
 matrix inv_makbounds[`r',2] = r(min)
 qui su makbounds1 if makbounds3>=`p'
 matrix inv_makbounds[`r',3] = r(min)
 
 qui su makbounds1 if makbounds4>=`p'
 matrix inv_makbounds[`r',4] = r(min)
 qui su makbounds1 if makbounds5>=`p'
 matrix inv_makbounds[`r',5] = r(min)
 
 qui su makbounds1 if makbounds6>=`p'
 matrix inv_makbounds[`r',6] = r(min)
 qui su makbounds1 if makbounds7>=`p'
 matrix inv_makbounds[`r',7] = r(min)
 
}

svmat inv_makbounds
ren inv_makbounds1 quantiles
replace quantiles = quantiles*1000
ren inv_makbounds2 makarov_u1
ren inv_makbounds3 makarov_l1
ren inv_makbounds4 makarov_u2
ren inv_makbounds5 makarov_l2
ren inv_makbounds6 makarov_u3
ren inv_makbounds7 makarov_l3

keep quantiles makarov*
keep in 1/1000
replace quantiles = round(quantiles)
svmat qtt1
svmat qtt2
svmat qtt3
svmat qOft1
svmat qOft2
svmat qOft3

graph twoway (line qtt1 qOft1 makarov_u1 makarov_l1 quantiles, lcolor(black black black black) lpattern(solid dot longdash longdash)) ///
			 (line qtt3 qOft3 makarov_u3 makarov_l3 quantiles, lcolor(gs10 gs10 gs10 gs10) lpattern(solid dot longdash longdash)), ///
			  legend(order(1 "QTT, {&sigma}=0" 2 "Quantiles of F{sub:{&Delta}}, {&sigma}=0" 3 "Makarov, {&sigma}=0" 5 "QTT, {&sigma}=5" 6 "Quantiles of F{sub:{&Delta}}, {&sigma}=5" 7 "Makarov, {&sigma}=5") colfirst ) ///
			  graphregion(fcolor(white)) xtitle("Quantiles") ytitle("Treatment Effect") ylabel(,nogrid)
graph export "${results_dir}\sim_fig3_left.svg", replace



/* Table 2. Cross-Sectional Simulations */
local col=0
matrix table = J(9,6,.)
foreach c in 0.5 -1 1 {
	foreach cov in 0.5 0 {
		local ++col
		clear
		sim, c(`c') cov(`cov')
		su delta, d
		matrix table[1,`col']=r(mean)
		matrix table[2,`col']=r(sd)
		matrix table[3,`col']=r(skewness)  
		moments
		matrix table[4,`col']=r(mean)
		matrix table[5,`col']=r(sd)
		matrix table[6,`col']=r(skewness)  

		qui su y if treat==1
		local v1 = r(sd)^2
		qui su y if treat==0
		local v0 = r(sd)^2
		di "Deconvolution Test: " `v1'-`v0'
		gen beta = .
		mata: mallow("y","treat")
		summarize beta, detail
		matrix table[7,`col']=r(mean)
		matrix table[8,`col']=r(sd)
		matrix table[9,`col']=r(skewness)  
	}
}


/* Figure 4.  Simulation 1, b: Panel Data */
foreach t in 4  16  {
	clear
	sim_panel, n(1000) t(`t')	
	local nt = 1000*`t'
	bys id: egen _tmean = mean(y) if treat==1
	bys id: egen _cmean = mean(y) if treat==0
	bys id: egen tmean=max(_tmean)
	bys id: egen cmean=max(_cmean)
	gen betahat = tmean-cmean
	gen y_resid = y-betahat if treat==1
	 replace y_resid = y if treat==0
	bys id: egen alphahat = mean(y_resid)
	gen beta = .
	keep if ever_treat==1
	replace y = y-alpha
	mata: mallow("y","treat")
	bys id: egen seq=seq()
	graph twoway (function y=normalden(x,1.5,1), range(-4 6) lw(thick) lcolor(black)) (kdensity  betahat if seq==1, lcolor(black) lpattern(dot)) (kdensity beta , lcolor(black) lpattern(longdash)) , ///
		legend( label(1 "True Distribution")label(2 "Regression Estimate") label(3 "Deconvolution Estimate")) ///
		graphregion(fcolor(white)) ylabel(,nogrid) xtitle("Treatment Effect") ytitle("Density")
	graph export "${results_dir}\sim_fig4_t`t'.svg", replace	
}




/* Replicate WSDP Analysis  */


* Independent Case: Deconvolution Assumptions Satisfied
clear
sim_gambia, rho(0)
tab year, gen(year)

gen older = floor(t/2)

/* QTTs */
preserve
qui reg y treatment older if t==7 | t==8, robust
predict resid,  resid
replace resid = resid +_b[treatment]*treatment

pctile tr_perc = resid if treatment==1, nq(10)
pctile c_perc = resid if treatment==0, nq(10)
 
gen q = _n 
keep in 1/9
keep q tr_perc c_perc

gen qte1 = tr_perc-c_perc
keep q qte1
tempfile qte1
save `qte1'

restore
keep if ever_treat==1
do "${estimation_dir}\gambia\gambia_simulations.do"


* Makarov Bounds
preserve
keep y treat
sort treat 

/* Need two matrices with treatment and control outcomes */
matrix y1 = J(366,1,.)
matrix y0 = J(122,1,.)
forv r=1/122 {
 matrix y0[`r',1] = y[`r']
}
forv r=123/488 {
 local i=`r'-122
 matrix y1[`i',1] = y[`r']
}


clear
svmat y0
svmat y1
qui makarov y01 y11 , grid(-10(0.1)10) d(0)
di r(bound_l)
di r(bound_u)
restore






/* Deconvolution */

/* 3 times as many treatment as control - so following code triples control observations */
gen n=_n
expand 3
bys n: egen seq=seq()
keep if (seq==1 & treat==1) | treat==0
drop n
/* Run Mallows' Algorithm */
gen beta = .
mata: mallow("y","treat")
ren beta beta_0

gen n=_n
tempfile file1
save `file1'





* ----------------------------------------------------------------- *
* Correlation = 0.1 so deconvolution assumptions "nearly" satisfied *
* ----------------------------------------------------------------- *
clear
sim_gambia, rho(0.1) 

gen year1 = t==1 | t==2
gen year2 = t==3 | t==4
gen year3 = t==5 | t==6
gen year4 = t==7 | t==8

gen older = floor(t/2)


/* QTTs */
preserve
qui reg y treatment older if t==7 | t==8, robust
predict resid,  resid
replace resid = resid +_b[treatment]*treatment

pctile tr_perc = resid if treatment==1, nq(10)
pctile c_perc = resid if treatment==0, nq(10)
 
gen q = _n 
keep in 1/9
keep q tr_perc c_perc

gen qte2 = tr_perc-c_perc
keep q qte2
tempfile qte2
save `qte2'
restore


keep if ever_treat==1
do "${estimation_dir}\gambia\gambia_simulations.do"


* Makarov Bounds
preserve
keep y treat
sort treat 

/* Need two matrices with treatment and control outcomes */
matrix y1 = J(366,1,.)
matrix y0 = J(122,1,.)
forv r=1/122 {
 matrix y0[`r',1] = y[`r']
}
forv r=123/488 {
 local i=`r'-122
 matrix y1[`i',1] = y[`r']
}


clear
svmat y0
svmat y1
qui makarov y01 y11 , grid(-10(0.1)10) d(0)
di r(bound_l)
di r(bound_u)
restore


/* 3 times as many treatment as control - so following code triples control observations */
gen n=_n
expand 3
bys n: egen seq=seq()
keep if (seq==1 & treat==1) | treat==0

/* Run Mallows' Algorithm */
gen beta = .
mata: mallow("y","treat")
ren beta beta_point1
drop n
gen n=_n
tempfile file2
save `file2'




* ------------------------------------------------------------------------- *
* Correlation = 0.8 so treatment effects fairly strongly correlated with Y0 *
* ------------------------------------------------------------------------- *
clear
sim_gambia, rho(0.8) 



gen year1 = t==1 | t==2
gen year2 = t==3 | t==4
gen year3 = t==5 | t==6
gen year4 = t==7 | t==8

gen older = floor(t/2)

/* QTTs */
preserve
qui reg y treatment older if t==7 | t==8, robust
predict resid,  resid
replace resid = resid +_b[treatment]*treatment

pctile tr_perc = resid if treatment==1, nq(10)
pctile c_perc = resid if treatment==0, nq(10)
 
gen q = _n 
keep in 1/9
keep q tr_perc c_perc

gen qte3 = tr_perc-c_perc
keep q qte3
tempfile qte3
save `qte3'
restore

keep if ever_treat==1
do "${estimation_dir}\gambia\gambia_simulations.do"


* Makarov Bounds
preserve
keep y treat
sort treat 

/* Need two matrices with treatment and control outcomes */
matrix y1 = J(366,1,.)
matrix y0 = J(122,1,.)
forv r=1/122 {
 matrix y0[`r',1] = y[`r']
}
forv r=123/488 {
 local i=`r'-122
 matrix y1[`i',1] = y[`r']
}


clear
svmat y0
svmat y1
qui makarov y01 y11 , grid(-10(0.1)10) d(0)
di r(bound_l)
di r(bound_u)
restore



/* 3 times as many treatment as control - so following code triples control observations */
gen n=_n
expand 3
bys n: egen seq=seq()
keep if (seq==1 & treat==1) | treat==0
drop n
/* Run Mallows' Algorithm */
gen beta = .
mata: mallow("y","treat")
ren beta beta_point8

gen n=_n
merge 1:1 n using `file1'
drop _merge
merge 1:1 n using `file2'
drop _merge




* Graph Deconvolution Estimates
graph twoway (function y=normalden(x,1.5,1), range(-4 6) lw(thick) lcolor(black))  (kdensity beta_0 , lcolor(gs10) ) ///	
	(kdensity beta_point1 , lcolor(black) lpattern(longdash)) (kdensity beta_point8 , lcolor(black) lpattern(dot)), ///
		legend( label(1 "True Distribution") label(2 "{&rho}=0") label(3 "{&rho}=0.1") label(4 "{&rho}=0.8")) ///
		graphregion(fcolor(white)) ylabel(,nogrid) xtitle("Treatment Effect") ytitle("Density")
graph export "${results_dir}\sim_fig_10.svg", replace



/* Compare Deciles of Distribution of Impacts to QTTs at each decile*/

 
pctile beta_deciles0 = beta_0, nq(10)
pctile beta_deciles1 = beta_point1, nq(10)
pctile beta_deciles8 = beta_point8, nq(10)


gen beta_deciles_true = 0.2184484 in 1
 replace beta_deciles_true = 0.6583788 in 2
 replace beta_deciles_true = 0.9755995 in 3
 replace beta_deciles_true = 1.2466529 in 4
 replace beta_deciles_true = 1.5000000 in 5
 replace beta_deciles_true = 1.7533471 in 6
 replace beta_deciles_true = 2.0244005 in 7
 replace beta_deciles_true = 2.3416212 in 8
 replace beta_deciles_true = 2.7815516 in 9
  
gen q = _n in 1/9
keep in 1/9
keep q beta_deciles*
merge 1:1 q using `qte1'
drop _merge
merge 1:1 q using `qte2'
drop _merge
merge 1:1 q using `qte3'
drop _merge


graph twoway (connected beta_deciles_true q, lw(thick) lcolor(black) mcolor(black)  lpattern(solid) c(l) m(D) ) ///
    (connected beta_deciles0 q, lcolor(gs10) mcolor(gs10)) ///	
	(connected qte1 q, lcolor(gs10) mcolor(gs10) lpattern(longdash)) , title("{&rho}=0") ///
	legend(off) ///
	graphregion(fcolor(white)) ylabel(,nogrid) xtitle("Treatment Effect") ///
	xtitle("Decile") ytitle("Effect") 	ylabel(0(1)4,nogrid)
graph export "${results_dir}\sim_fig11_rho0.svg", replace	
	
graph twoway (connected beta_deciles_true q, lw(thick) lcolor(black) mcolor(black)  lpattern(solid) c(l) m(D) ) ///
    (connected beta_deciles1 q, lcolor(gs10) mcolor(gs10)) ///	
	(connected qte2 q, lcolor(gs10) mcolor(gs10) lpattern(longdash)) , title("{&rho}=0.1") ///
	legend(off) ///
		graphregion(fcolor(white)) ylabel(,nogrid) xtitle("Treatment Effect") ///
	xtitle("Decile") ytitle("Effect") 	ylabel(0(1)4,nogrid)
graph export "${results_dir}\sim_fig11_rho1.svg", replace	

graph twoway (connected beta_deciles_true q, lw(thick) lcolor(black) mcolor(black)  lpattern(solid) c(l) m(D) ) ///
    (connected beta_deciles8 q, lcolor(gs10) mcolor(gs10)) ///	
	(connected qte3 q, lcolor(gs10) mcolor(gs10) lpattern(longdash)) , title("{&rho}=0.8") ///
	legend( label(1 "True Deciles") label(2 "Deciles of F{sub:{&Delta}}") label(3 "QTT at Deciles")) ///
		graphregion(fcolor(white)) ylabel(,nogrid) xtitle("Treatment Effect") ///
	xtitle("Decile") ytitle("Effect") 	ylabel(0(1)4,nogrid)
graph export "${results_dir}\sim_fig11_rho8.svg", replace	

	

exit
