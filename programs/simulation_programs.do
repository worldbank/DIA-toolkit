/******************************************************************************/
/* CIA and Normal Potential Outcomes */
/******************************************************************************/

/* Program for Multivariate Normal Draws */

cap: program drop mkcorr
program define mkcorr
	version 4.0 

	local k = rowsof(P)
	matrix A = cholesky(P)

	local i 1 
	quietly {
		while `i'<=`k' {
			gen c`i' = invnorm(uniform())
			local i=`i'+1
		}
		local i 1
		while `i'<=`k' {
			matrix row = A[`i',.]
			matrix score eps`i' = row
			local i=`i'+1
		}
		local i 1 
		while `i' <= `k' {
			drop c`i'
			local i=`i'+1
		}
	}
end
 
capture program drop moments
program moments, rclass

	qui {
		marksample touse

		reg y treat if `touse'
		local mean = _b[treat]
		su y if `touse' & treat==1, d
		local v1 = r(Var)
		su y if `touse' & treat==0, d
		local v0 = r(Var)
		local sd = sqrt(`v1'+`v0')
		gen y2 = y*y
		gen y3 = y2*y

		su y if `touse' &  treat==1
		local y_1 = r(mean)
		su y if `touse' &  treat==0
		local y_0 = r(mean)
		su y2 if `touse' &  treat==1
		local y2_1 = r(mean)
		su y2 if `touse' &  treat==0
		local y2_0 = r(mean)
		su y3 if `touse' &  treat==1
		local y3_1 = r(mean)
		su y2 if `touse' &  treat==0
		local y3_0 = r(mean)

		local e3 = `y3_1'-3*`y2_1'*`y_0'+3*`y_1'*`y2_0'-`y3_0'
		local skew = (`e3'-3*`mean'*(`sd')^2- `mean'^3)/(`sd')^3
	}
	
	di "Mean: " round(`mean', .01)
	di "SD: "   round(`sd',   .01)
	di "Skew: " round(`skew', .01)
	
	return scalar mean=`mean'
	return scalar sd = `sd'
	return scalar skewness = `skew'
end 

capture program drop sim
program sim , eclass

	syntax [, N(integer 1000) cov(real 0.5) c(real 0.5)]
	set seed 2929

	qui {
		clear
		set obs `n'
		
		matrix P = (1, `cov' \ `cov', 0.5)
		mkcorr
		ren eps2 eps0
		gen x1 = rnormal(1,1)
		gen y1 = 1+`c'*x1+eps1
		gen y0 = 0.5*x1 + eps0
		gen delta = y1-y0
		
		gen rand = uniform()
		 su rand, d
		gen treat = rand>=r(p50)
		gen tx1 = treat*x1
		gen y = y1*treat+y0*(1-treat)
	}
end

capture program drop sim_rankinv
program sim_rankinv , rclass

	syntax [, N(integer 1000) ]
	set seed 2929

	qui {
		clear
		set obs `n'
		
		gen m=rnormal()
		sort m 
		gen rank = _n
		tempfile m
		save `m'
		
		clear
		set obs `n'
		gen rank = _n
		merge 1:1 rank using `m'
		gen y0 = m
		drop _merge m rank
		gen rando = rnormal()
		sort rando
		gen rank = _n
		merge 1:1 rank using `m'
		gen y1 = rank+m
		drop rando _merge m rank
		gen delta = y1-y0
		gen rand = rnormal()
		su rand, d
		gen treat = rand>=r(p50)
		drop rand
		gen y = y1*treat+y0*(1-treat)
	}
end

capture program drop sim_panel
program sim_panel ,eclass 
	syntax [, N(integer 1000) T(integer 4) cov(real 0.5)]
	set seed 2929

	qui {
		clear
		local nt = `n'*`t'
		set obs `nt'
		gen id = ceil(_n/`t')
		gen t = mod(_n,`t')+1
		sort id t
		
		matrix P = (1,`cov' \ `cov', 0.5)
		mkcorr
		gen _delta = eps1 if t==1
		gen _alpha = eps2 if t==1
		bys id: egen delta = max(_delta)
		bys id: egen alpha = max(_alpha)
		 drop _delta _alpha
		gen eps = rnormal()
		gen y1 = alpha+1.5+delta+eps
		gen y0 = alpha+eps
		
		
		gen random = runiform() if t==1
		sort id t
		by id: replace random = random[1]
		
		local treat_after = ceil((`t')/2)
		qui su random if t>=`treat_after', d
		gen treatment =  random>=r(p50) & t>`treat_after'
		bys id: egen ever_treat = max(treatment)
		gen treat = ever_treat*(t>`treat_after')
		gen y = y1*treat+y0*(1-treat)
		gen effect = y1-y0
	
	}
end

capture program drop sim_gambia
program sim_gambia ,eclass 
	syntax [, rho(real 0.5) seed(integer 223417)]
	set seed `seed'

	qui {
		clear
		local n=120
		local t=8
		local nt = `n'*`t'
		set obs `nt'
		gen id = ceil(_n/`t')
		gen t = mod(_n,`t')+1
		sort id t
		gen year = 1 if t==1 | t==2
		 replace year = 2 if t==3 | t==4
		 replace year = 3 if t==5 | t==6
		 replace year = 4 if t==7 | t==8
		
		gen treatment = 0
		 replace treatment = 1 if id<=61 & t>2
		bys id: egen ever_treat = max(treatment)
		

		# delimit ;
		matrix P = 
			(0.5,`rho'*sqrt(0.5) \ 
			`rho'*sqrt(0.5), 1);
		#delimit cr
		mkcorr
		gen _delta = eps2
		bys id: gen delta = _delta[1]
		gen noise=0.1*rnormal()
		gen _eps = eps1
		bys id: gen eps = _eps[1]
		gen epsilon = .5*sqrt(year)*eps+noise
		gen _alpha = rnormal()
		bys id: gen alpha = _alpha[1]
		 drop _delta _alpha
		gen y1 = alpha+epsilon+1.5+delta
		gen y0 = alpha+epsilon
		
	
		gen y = y1*treatment+y0*(1-treatment)
		gen effect = 1.5+delta
	
	}
end

capture program drop sim_makarov
program sim_makarov, eclass
	syntax [, N(integer 1000)]
	set seed 2929
	qui {
		clear
		set obs `n'
		
		gen effect1 = 1
		gen effect2 = rnormal(1,0.5)	
		gen effect3 = rnormal(1,5)
		
		gen eps = rnormal(0,1)
		gen rand = runiform()
		qui su rand, d
		gen treat = (rand>=r(p50))
		gen y1 = treat*effect1+eps
		gen y2 = treat*effect2+eps
		gen y3 = treat*effect3+eps
	}
end

capture program drop quantiles
program quantiles  ,eclass 
	syntax varlist(min=1 max=2 numeric) [if], [nq(real 100)]
	marksample touse
	
	tokenize `varlist'
	local outcome `1'
	local treatment `2'
	
	local nq1 = `nq'-1
	
	// Point estimates for arbitrary quantiles
	
	tempname q b 
	
	if "`treatment'"~="" {
		matrix `q' = J(2, `nq1', .)
	
		forvalues v = 0/1 { 

			_pctile `outcome' if `touse' & `treatment' == `v', nq(`nq')
		
			forvalues k = 1/`nq1' {
				matrix `q'[`v'+1, `k'] = r(r`k')
			}

		}

		matrix `b' = `q'[2,1...] - `q'[1,1...]
	
	}
	if "`treatment'"=="" {
		matrix `q' = J(1, `nq1', .)


			_pctile `outcome' if `touse', nq(`nq')
		
			forvalues k = 1/`nq1' {
				matrix `q'[1, `k'] = r(r`k')
			}

		matrix `b' = `q'[1,1...]
	}
	matrix `b' = `b''
	ereturn matrix quantiles =  `b'
end	

/******************************************************************************/
/* Mallows' Deconvolution Algorithm */
/******************************************************************************/

* This is a crude Mallows' function. 
* It takes variables with partial residuals 
* and treatment status as input. There need to be the same number
* of treatment and control residuals. Also, need to generate a 
* variable beta = . before running.

// Mata Mallows' Deconvolution Function

capture mata: mata drop mallow()
mata: 

	void mallow(string scalar varname1, string scalar varname2){

		// Pull in Residuals

		eps = st_data(.,varname1)
		treat = st_data(.,varname2)

		// Create A and C Matrices

		Btot=J(0,1,0)
		Atot=J(0,1,0)
		Ctot=J(0,1,0)

		A=select(eps,(treat:==1))
		C=select(eps,(treat:==0))

		Ctot=Ctot\C
		Atot=Atot\A

		C=sort(Ctot,1)
		A=sort(Atot,1)


		// Set initial condition
		
		B=sort(A-C,1)

		// Iteration
		
		j1=1
		
		while (j1<=2000){
		
			// Randomly Sort B
			indexsort=sort((uniform(rows(B),1),range(1,rows(B),1)),1)
			Bnew=B[indexsort[,2]]
			
			// Create a new A that matches randomly sorted B with C
			Atilde=Bnew+C
			
			// Sort original A based on Atilde
			D=sort((A,Atilde),2)
			
			// Save this sorted version of A
			Anew=D[.,1]
			
			// Create new draw from treatment effect distribtion
			B=Anew-C
			
			// Ignore first 500 iterations to allow the algorithm to "burn in"
			if (j1>=501){
				Btot=Btot\B
			}
			
			j1=j1+1
		}

		// Mean of Distribution of Treatment Effects
		mean(Btot)
		
		// Standard Deviation of Distribution of Treatment Effects
		sqrt(variance(Btot))

		// Export the simulated beta to stata
		stata("count")
		stata("expand 750")
		BETA=0
		st_view(BETA,.,"beta")
		BETA[.]=Btot
	}
end
