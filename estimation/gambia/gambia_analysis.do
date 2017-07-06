/* This program produces the analysis in Section 7.2 of distributional analysis toolkit */


use "${data_dir}/gambia/analysis_data.dta", clear

/* Collapse to a single observation per school-year-grade.
 Then can apply Mallow's algorithm without concerns about 
 clustering students within a school. */ 
 
keep if ((grade_student==3 | grade_student==5) & (year==2008 | year==2010)) | ((grade_student==4|grade_student==6) & (year==2009 | year==2011))

gen older = 0
 replace older = 1 if grade_student==5 & (year==2008 | year==2010)
 replace older = 1 if grade_student==6 & (year==2009 | year==2011)
collapse (mean) z_score_total_NL wsd2008 grant2008 control2008 community_support, by(schoolID  older year)


* Restrict sample to balanced panel of schools with no missing observations
gen n = 1
drop if missing(z_score_total_NL)
bys schoolID: egen count=sum(n)
keep if count==8
drop count n

/* Generate Year Indicators */
tab year, gen(y)


/* Create Treatment Indicator */
gen treatment = wsd2008
 replace treatment = 0 if year==2008


/* Limit attention to WSD treatment (drop grant treatment) */
keep if wsd2008==1 | control2008==1


* Save Temporary Version of Data for Later
tempfile data
save `data'


* Stacked Regression *
reg z_score_total_NL treatment i.year older, robust
eststo stacked

estout stacked using "${results_dir}/gambia_reg.txt", replace ///
	keep(treatment) cells(b(star fmt(2)) p( par fmt(2))) ///
	varl(treatment Treatment)  /// 
	collabels(,none) 


* Subgroup Analysis *
gen t_supp = treatment*community_support
gen t_nosupp = treatment*(1-community_support)
reg z_score_total_NL  t_nosupp t_supp community_support i.year older, robust
 eststo bySupp

 estout bySupp using "${results_dir}/gambia_reg.txt", append ///
	keep(t_nosupp t_supp) cells(b(star fmt(2) )   p( par fmt(2)) ) /// 
	stardrop(*) ///
	collabels(,none) 

 
reg z_score_total_NL  treatment t_nosupp community_support i.year older, robust
 eststo bySupp_b

 estout bySupp_b using  ///
	"${results_dir}/gambia_reg.txt", append ///
	keep(t_nosupp) cells(b(star fmt(2))   p(par fmt(2))) ///
	varl(t_nosupp diff) stardrop(*) ///
	collabels(,none) 


* Show treatment group residual variation > control group residual variation (Necessary condition for deconvolution) *
areg z_score_total_NL treatment i.year older, absorb(schoolID)
predict resid,  resid
 replace resid = resid + _b[treatment] if treatment==1
table treatment, c(sd resid)





*************************************************************************************
* Mallows' algorithm code is adapted from the supplementary materials to:                         *
* Manuel Arellano and Stephane Bonhomme:                                            *
* ``Identifying Distributional Characteristics in Random Coefficients Panel Data    *
* Models'' (December 2010)                                                          *
*************************************************************************************


* Define stata output *
gen beta=0
gen densbeta=0
gen alphahat=0
gen betahat=0

gen mean = 0
gen var = 0


* Start mata routine *

mata:

	T=8

	/* controls and other variables */

	ID = st_data(.,"schoolID")

	y = st_data(., "z_score_total_NL")

	z = (st_data(.,"y2"),st_data(.,"y3"),st_data(.,"y4"), st_data(.,"older"))

	x = (st_data(.,"treatment"),J(rows(z),1,1))

	CH = (st_data(.,"year"))

	TR = (st_data(.,"wsd2008"))

	/* Technical: collapse ID so 1,...,# Groups */

	ID2=J(rows(ID),1,1)
	ID2[1,1]=1
	j1=2
	while (j1<=rows(ID)){
		if (ID[j1]==ID[j1-1]){
			ID2[j1]=ID2[j1-1]
			}
		else {
			 ID2[j1]=ID2[j1-1]+1
			}
		j1=j1+1
		}
	ID=ID2
	mata drop ID2


	/* First step: estimate parameters common across schools */

	deltaden=J(cols(z),cols(z),0)
	deltanum=J(cols(z),1,0)
	hz = J(cols(x),cols(z),0)

	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)

		 deltaden = deltaden+z2'q2*z2
		 deltanum = deltanum+z2'q2*y2
		 hz = hz+h2*z2
		 
		 j1=j1+1
		 }

	invdeltaden = invsym(deltaden)
	delta = invdeltaden*deltanum
	delta

	/* Second step: Estimate School Specific Effects Covariates Interacted with Fixed Effects */

	vardeltanum = J(cols(z),cols(z),0)
	meangamma = J(cols(x),1,0)
	varmeangamma1 = J(cols(x),cols(x),0)


	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)
		 
		 /* Standard errors common coefficients */
		 vardeltanum = vardeltanum+z2'*q2*(y2-z2*delta)*(y2-z2*delta)'*q2'*z2
		 
		 /* Overall mean individual effects */
		 meangamma = meangamma+h2*(y2-z2*delta)
		 
		 j1=j1+1
		}

	meangamma=meangamma/maxID
	MG=0
	st_view(MG,.,"mean")
	MG[.] = J(rows(y),1,meangamma[1,1])

	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)
		 
		 /* Standard errors average individual effects */
		 vect = h2*(y2-z2*delta-x2*meangamma)-hz*invdeltaden*z2'*q2*(y2-z2*delta)
		 varmeangamma1 = varmeangamma1+vect*vect'
			 
		 j1=j1+1
		}


	/* Standard Errors for Common coefficients */

	stddelta = sqrt(diagonal(invdeltaden*vardeltanum*invdeltaden'))
	display(" Coefficient -  Std. Dev. ")
	(delta,stddelta)

	/* Standard Error for Average of individual effects */

	stdmeangamma=sqrt(diagonal(varmeangamma1/maxID^2))

	varmeangamma1
	display(" Coefficient -  Std. Dev. ")
	(meangamma,stdmeangamma)

	/* Variances of Effects - Under Different Covariance Assumptions */

	/* Create Matrices */
	invxx=J(cols(x),cols(x),0)
	sigma2u=0
	vargamma=J(cols(x),cols(x),0)
	varsigma2u=0
	varvargammaiid=J(2*cols(x)^2+1,2*cols(x)^2+1,0)
	invxxrobust=J(cols(x),cols(x),0)
	varvargammarobust=J(cols(x)^2,cols(x)^2,0)
	varutvL = J(T^2,1,0)
	varvarutvL=J(T^2,T^2,0)
	varutvgL = J(cols(x)^2,1,0)
	varvargammatvL=J(cols(x)^2,cols(x)^2,0)
	varbeta=J(cols(x)-1,cols(x)-1,0)
	varvarbeta=J((cols(x)-1)^2,(cols(x)-1)^2,0)


	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)
		 q2b = I(T)-J(T,T,1)/T
		 h2b=pinv(q2b*x2[.,1])
		 
		 /* Variance errors, iid */
		 sigma2u = sigma2u+(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))

		 /* Variance individual effects, iid */
		 vargamma = vargamma+h2*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*h2'
		 invxx = invxx+invsym(x2'x2)

		 /* Variance individual effects, robust */
		 invxxrobust = invxxrobust+(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))*invsym(x2'x2)

		 /* Variance beta in first differences, robust */
		 varbeta = varbeta+h2b*q2b*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*q2b'*h2b'/*
		 */ -(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))*h2b*q2b*h2b'

		 j1=j1+1
		}


	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)
		 q2b = I(T)-J(T,T,1)/T
		 h2b=pinv(q2b*x2[.,1])


		 /* std variance errors, iid */
		 varsigma2u = varsigma2u+((y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x)))^2

		 /* std variance individual effects, iid */
		 vect=vec(h2*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*h2')
		 vect=vect\(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))
		 vect=vect\vec(invsym(x2'x2))
		 varvargammaiid = varvargammaiid+vect*vect'

		 /* std variance individual effects, robust */
		 vect=vec(h2*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*h2'/*
				*/-(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))*invsym(x2'x2))
		 varvargammarobust = varvargammarobust+vect*vect'

		 /* std variance beta in first differences, robust */
		 vect = h2b*q2b*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*q2b'*h2b'/*
		 */ -(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))*h2b*q2b*h2b'
		   varvarbeta = varvarbeta+vect*vect'
		 

		 j1=j1+1
		}




	/* Variance errors, iid errors */
	sigma2u=sigma2u/maxID
	stdsigma2u=sqrt(varsigma2u/maxID^2)
	display(" Coefficient -  Std. Dev. ")
	sigma2u,stdsigma2u

	/* Variance individual effects, iid errors */
	vargamma=vargamma/maxID
	vargammaiid=vargamma-sigma2u*invxx/maxID
	M=I(cols(x)^2),(-vec(invxx/maxID)),(-sigma2u*I(cols(x)^2))
	stdvargammaiid = sqrt(diagonal(M*varvargammaiid*M'/maxID^2-vec(vargammaiid)*vec(vargammaiid)'/maxID))
	display(" Coefficient -  Std. Dev. ")
	vec(vargammaiid),stdvargammaiid

	/* Variance individual effects, robust */
	vargammarobust = vargamma-invxxrobust/maxID
	stdvargammarobust = sqrt(diagonal(varvargammarobust/maxID^2-vec(vargammarobust)*vec(vargammarobust)'/maxID))
	display(" Coefficient -  Std. Dev. ")
	vec(vargammarobust),stdvargammarobust

	/* Export Prefered Variance To Stata (out of mata) */
	V=0
	st_view(V,.,"var")
	V[.] = J(rows(y),1,vargammarobust[1,1])


	/* Variance beta in first differences, robust */
	varbeta = varbeta/maxID
	stdvarbeta = sqrt(diagonal(varvarbeta/maxID^2-vec(varbeta)*vec(varbeta)'/maxID))
	display(" Coefficient -  Std. Dev. ")
	vec(varbeta),stdvarbeta



	/* Regression of individual effects */

	gamma=J(maxID,cols(x),.)

	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)
		 
		 /* Fixed effects estimates */
		 gamma[j1,.]=(h2*(y2-z2*delta))'

	j1=j1+1
	}


	j1=1
	maxID=max(ID)
	while (j1<=maxID){

		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)

		 
		 j1=j1+1
		}



	/* Run Mallows' Algorithm (Mallows, 2007; Arellano and Bonhomme, 2012) */

	/* Save Partial Residuals */
	Y=y-z*delta
	Btot=J(0,1,0)

	/* Create A and C Matrices */

	Btot=J(0,1,0)
	Atot=J(0,1,0)
	Ctot=J(0,1,0)

	/*Can Only identify effects among eventually treated schools */
	YY=select(Y,(TR:==1))
	CHCH=select(CH,(TR:==1))

	Y1=select(YY,(CHCH:==2008))
	Y2=select(YY,(CHCH:==2009))
	Y3=select(YY,(CHCH:==2010))
	Y4=select(YY,(CHCH:==2011))

	/* Logic of Next Two Lines:
	   Y4 = Delta+eps4 
	   Y3 = Delta+eps3
	   so Y4-Y3 = eps4-eps3 
	   C terms are of form: epsB-epsA
	   
	   Similarly: 
	   Y1 = eps1 
	   so Y4-Y1 = Delta + eps4-eps1
	   A terms are of form: Delta + epsB-epsA*/
	/* "\" appends vectors */

	C=(Y4-Y3) \ (Y4-Y2) \ (Y3-Y2)
	A=(Y4-Y1)\(Y3-Y1) \ (Y2-Y1)

	Ctot=Ctot\C
	Atot=Atot\A


	C=sort(Ctot,1)
	A=sort(Atot,1)


	/* Set initial condition */
	B=sort(A-C,1)

	/* Iteration */
	j1=1
	while (j1<=2000){
		/* Randomly Sort B */
		indexsort=sort((uniform(rows(B),1),range(1,rows(B),1)),1)
		Bnew=B[indexsort[,2]]
		
		/* Create a new A that matches randomly sorted B with C */
		Atilde=Bnew+C
		
		/* Sort original A based on Atilde */
		D=sort((A,Atilde),2)
		
		/* Save this sorted version of A */
		Anew=D[.,1]
		
		/* Create new draw from treatment effect distribtion */
		B=Anew-C
		
		/* Ignore first 500 iterations to allow the algorithm to "burn in" */
		if (j1>=501){
			Btot=Btot\B
		}
		j1=j1+1
	}

	/* Mean of Distribution of Treatment Effects */
	mean(Btot)
	/* Standard Deviation of Distribution of Treatment Effects */
	sqrt(variance(Btot))


	Btot=Btot[trunc(rows(Btot)*uniform(20000,1))+J(20000,1,1)]
	rows(Btot)

	/* Export the simulated beta to stata */
	stata("count")
	stata("expand 1125")
	stata("keep if _n<=20000")
	BETA=0
	st_view(BETA,.,"beta")
	BETA[.]=Btot


	/* Export fixed effects estimates to stata */
	ALPHA1=0
	st_view(ALPHA1,.,"alphahat")
	ALPHA1[.]=gamma[.,2]\J(20000-rows(gamma),1,.)
	BETA1=0
	st_view(BETA1,.,"betahat")
	BETA1[.]=gamma[.,1]\J(20000-rows(gamma),1,.)

end



******************* END MALLOWS' ALGORITHM CODE ****************************************




* Graphs of Densities*
replace betahat = . if betahat==0 /* Set To Missing for Untreated Schools */

gen xbeta=-2+4*_n/(20000) 
kdensity betahat, at (xbeta) gen(densbetahat) nograph

kdensity beta, at (xbeta) gen(densbetasimul) nograph
graph twoway connected densbetasimul densbetahat xbeta, ms(p p) color(black) c(l l) legend(off) xtitle("beta") ytitle("density")



* Normal Example *
local mean = mean[1]
local sd = sqrt(var[1])
gen density = normalden(xbeta,`mean', `sd')

graph twoway  ( connected densbetasimul xbeta, color(black) lpattern(solid) lcolor(black) ms(i) c(l)) ///
 (connected densbetahat xbeta, lcolor(black) ms(i) c(l) lpattern(dash) color(black)) ///
 (connected density xbeta, color(black) lcolor(black) lpattern(dot) ms(i) c(l)) ,  ///
 legend( label(1 "Deconvolution") label(2 "Estimates") label(3 "Normal")) xtitle("beta") ///
 ytitle("Density") xtitle("Delta-hat") graphregion(fcolor(white)) ylabel(,nogrid)
 graph export "${results_dir}/densities.png", replace

 
* Save deciles for comparison to QTTs at each decile at bottom of program 
pctile effects = beta, nq(11)  
gen decile = _n in 1/10
keep in 1/10
keep decile effects
tempfile deciles
save `deciles'

 
* -------------- *
* Makarov Bounds *
* -------------- *

use `data', clear


keep treatment z_score_total
tempfile reshape
save `reshape'

keep if treatment==0
ren z_score_total_NL z_score_total_NL0
gen i=_n
tempfile control
save `control'

use `reshape', clear
keep if treatment==1
ren z_score_total_NL z_score_total_NL1
gen i=_n
merge 1:1 i using `control'


set obs 1000
gen xbeta=-2+4*_n/(1000) 

gen bound_l = .
gen bound_u = .

forv i=1/1000 {
 local d = xbeta[`i']
 makarov z_score_total_NL1 z_score_total_NL0, grid( -2(.004)4) d(`d')
 qui replace bound_l = r(bound_l) in `i'
 qui replace bound_u = r(bound_u) in `i'
 }
 

graph twoway connected bound_l bound_u xbeta, ms(i i) color(black) c(l l) legend(label(1 "Lower Bound") label(2 "Upper Bound")) xtitle("beta") ///
 ytitle("Density") xtitle("Beta") graphregion(fcolor(white))

 
 
 
 
 
 
/* Compare Deciles of Distribution of Impacts to QTTs at each decile*/

use `data', clear
keep if year==2011
keep z_score_total_NL treatment older 
tempfile qtt_data
save `qtt_data'

local iter 1000
matrix bs = J(`iter',10,0)


/* Bootstrap QTTs - First iteration is real data */
forv z=1/`iter' {
quietly {
* Stacked Regression *
reg z_score_total_NL treatment older , robust

predict resid,  resid
replace resid = resid +_b[treatment]*treatment

pctile tr_perc = resid if treatment==1, nq(11)
pctile c_perc = resid if treatment==0, nq(11)
 
gen q = _n in 1/11
keep in 1/10
keep q tr_perc c_perc

gen qte = tr_perc-c_perc
forv q=1/10 {
 matrix bs[`z',`q'] = qte[`q']
}
}

* Here is where we draw bootstrap sample
clear 
use `qtt_data'
bsample
}

clear
svmat bs

gen true = 0
 replace true = 1 in 1
 
matrix qtt = J(10,6,0)
forv i=1/10 {
 qui su bs`i' if true==1
 matrix qtt[`i',1] = r(mean)
 qui su bs`i',d
 matrix qtt[`i',2] = r(sd)
 sort bs`i'
 matrix qtt[`i',3] = bs`i'[25]
 matrix qtt[`i',4] = bs`i'[975]
 matrix qtt[`i',5] = bs`i'[50]
 matrix qtt[`i',6] = bs`i'[959]
 
} 

clear
svmat qtt

* Bring in deciles of treatment effects
gen decile = _n
merge 1:1 decile using `deciles'

 drop _merge

graph twoway (connected effects decile, lcolor(gs10) mcolor(gs10)  lpattern(solid) c(l) m(D) ) ///
	(connected qtt1 decile, lcolor(black) mcolor(black)  lpattern(solid) c(l) ) ///-
	(connected qtt3 decile, lcolor(black) lpattern(longdash) c(l) m(i)) ///
	(connected qtt4 decile, lcolor(black)  lpattern(longdash) c(l) m(i)), ///
	xtitle("Decile") ytitle("Effect") graphregion(fcolor(white)) ///
	legend(order(1 "Deconvolution Estimates" 2 "QTT" 3 "95% CI for QTT") ) ylabel(,nogrid)
graph save "${results_dir}/qtt.png", replace	
 
