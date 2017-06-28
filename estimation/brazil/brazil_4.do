/**************************************************************************************************/
// CONDITIONAL AVERAGE EFFECTS
/**************************************************************************************************/

run "${estimation_dir}/brazil/brazil_setup.do"

// Prepare data

xtile   initial = vl_proficiencia_bl, nq(2)
replace initial = initial - 1

rename dumm_rp_24_bl repeater
rename dumm_rp_14_bl welfare
rename dumm_rp_49_bl worker
rename dumm_rp_50_bl income

// Regressions

capture program drop myboot
program myboot, eclass

	syntax , stat(name) vce(varname)

	tempname B b V

	matrix `B' = J(1, 18, .)
	local   k  = 0

	qui foreach v of varlist initial repeater worker income welfare female {

		local k  = `k' + 1
		local K1 = 3*(`k' - 1) + 1
		local K2 = 3*(`k' - 1) + 2
		local K3 = 3*(`k' - 1) + 3

		reg outcome c.treatment#i.`v' i.`v', cl(`vce')
		matrix `b' = e(b)
		matrix `V' = e(V)
		
		matrix `B'[1,`K1'] = `b'[1,1]
		matrix `B'[1,`K2'] = `b'[1,2]
		matrix `B'[1,`K3'] = `b'[1,2] - `b'[1,1] - (`stat'[1, `K2'] - `stat'[1, `K1'])
		matrix `B'[1,`K3'] = `B'[1,`K3']/sqrt(`V'[1,1] + `V'[2,2] - 2*`V'[2,1])
		
		local names "`names' `v'_0 `v'_1 `v'_t"	
	}

	matrix colnames `B' = `names'

	ereturn post `B'
end

// Estimation

matrix B = J(1, 18, 0)
myboot, stat(B) vce(cd_escola)
matrix B = e(b)

// Bootstrap

bootstrap, seed(0101) reps(5000) cl(cd_escola) id(bs_id) ///
	saving("${bs_dir}/z_sample.dta", replace) :           ///
	myboot, stat(B) vce(bs_id)

// Confidence region

bscr *_0 *_1 using "${bs_dir}/z_sample.dta", fwer(1)

// Critical values

use "${bs_dir}/z_sample.dta", clear

matrix T0 = J(1, 6, 0)
local  k  = 0
local  l  = 95

qui foreach v in initial repeater worker income welfare female {
	local k = `k' + 1
	gen abs_`v' = abs(_b_`v'_t)
	sum abs_`v', det
	matrix T0[1,`k'] = r(p`l')
}

qui forvalues t = 1/2 {
	bscr *_t, fwer(`t') level(`l') norecenter
	matrix T`t' = e(cv)
}

matrix T = [T0 \ T1 \ T2]'
matrix list T

/**************************************************************************************************/
// CONDITIONAL QUANTILE EFFECTS
/**************************************************************************************************/

run "${estimation_dir}/brazil/brazil_setup.do"

// Prepare data

xtile   initial = vl_proficiencia_bl, nq(2)
replace initial = initial - 1

rename dumm_rp_24_bl repeater
rename dumm_rp_14_bl welfare
rename dumm_rp_49_bl worker
rename dumm_rp_50_bl income

// Conditional quantile effects

capture mata: mata drop cqtt()
mata:
	real rowvector cqtt(string scalar varlist, real vector q)
	{
		real scalar    B, K, N, Q, i00, i01, i10, i11
		real colvector I, R, Y, T, W00, W01, W10, W11
		real rowvector b
		real matrix    S, X
		
		st_view(S = ., ., tokens(varlist))
		
		K = cols(S) ; N = rows(S) ; Q = length(q) ; B = 2*(K - 2)*Q
		Y = S[.,1]
		b = J(1, B, .) ; q = q/100 ; R = range(1, N, 1)
		
		for (i = 3; i <= K; i++) {
		
			i00 = 2*Q*(i - 3) + 1
			i01 = 2*Q*(i - 3) + Q
			i10 = 2*Q*(i - 3) + Q + 1
			i11 = 2*Q*(i - 3) + 2*Q
			
			I = select(R, S[.,i] :!= .)
			
			W00 = (1 :-  S[.,2]):*(1 :- S[.,i])
			W01 = (1 :-  S[.,2]):*S[.,i]
			W10 = S[.,2]:*(1 :- S[.,i])
			W11 = S[.,2]:*S[.,i]
			
			b[(i00..i01)] = mm_quantile(Y[I], W10[I], q) - mm_quantile(Y[I], W00[I], q)
			b[(i10..i11)] = mm_quantile(Y[I], W11[I], q) - mm_quantile(Y[I], W01[I], q)
		}
		
		return(b)
	}
end

// Wrapper for the bootstrap

capture program drop myboot
program myboot, eclass

	tempname b

	local varlist "initial repeater worker income welfare female"
	
	foreach x of local varlist {
		forvalues v = 0/1 {
			foreach q in 10 25 50 75 90 {
				local names "`names' `x'_`v'_`q'"
			}
		}
	}

	local varlist "outcome treatment `varlist'"
	
	mata: st_matrix("`b'", cqtt("`varlist'", (10, 25, 50, 75, 90)))

	matrix colnames `b' = `names' 
	
	ereturn post `b'
end

// Bootstrap

bootstrap, seed(0101) reps(5000) cl(cd_escola) ///
	saving("${bs_dir}/q_sample.dta", replace) : ///
	myboot

// Equality test

use "${bs_dir}/q_sample.dta", clear

foreach v in initial repeater worker income welfare female {
	foreach k in 10 25 50 75 90 {
		gen diff_`v'_`k' = abs(                                         ///
				_b_`v'_1_ `k' - _b_`v'_0_`k'                            ///
				- (`_b_`v'_1_`k'[observed]' - `_b_`v'_0_`k'[observed]') ///
			)
	}
	egen max_`v' = rowmax(diff_`v'_*)
}

keep max_*

matrix obs = [5.612327, 1.406586, 1.558539, 1.799575, 2.370861, 2.317842]
matrix colnames obs = initial repeater worker income welfare female

bscr max_*, fwer(1) level(90) stat(obs) norecenter

