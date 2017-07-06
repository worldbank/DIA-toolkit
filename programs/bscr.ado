/******************************************************************************/
// INTERFACE
/******************************************************************************/

program bscr, eclass
    version 13.1

	local cmdline `"bscr `0'"'
	
    // Load data

    syntax [anything] [using/] [, *]
	
	preserve
    
	if `"`using'"' != "" {
		qui use `anything' using `using', clear
	}
	else if "`anything'" != "" {
		keep `anything'
	}

    // Syntax

    local 0 "`anything', `options'"

    syntax [varlist] [if] [in], ///
       [Stat(name)              ///
        UNIForm                 ///
        FWER(string)            ///
        FDP(string)             ///
		noRECenter 		        ///
		STREAMline 		        ///
        Level(cilevel)          ///
		noTABle]

    // Reserve names

    tempname CR D S

	// Select sample
	
	marksample touse
	
	qui count if `touse'
	
	local N = r(N)
	local K = c(k) - 1

    // Check consistency

    local v = 0
    local v = `v' + ("`uniform'" != "")
	local v = `v' + ("`fwer'" != "")
	local v = `v' + ("`fdp'" != "")

    if `v' > 1 {
        di as err "choose only one type of confidence region"
        error 184
    }
	else if `v' == 0 {
		local fwer "1"
	}
	
	if "`streamline'" != "" & "`fwer'" == "" & "`fdp'" == "" {
		di as error "option streamline requires fwer() or fdp()"
        error 197
	}

	// Check fwer() option
	
	if "`fwer'" != "" {
		if real("`fwer'") == . {
			di as err "option fwer() incorrectly specified"
			exit 125
		}
		if `fwer' != int(`fwer') | `fwer' < 0 {
			di as err "option fwer() incorrectly specified"
			di as err "number of false rejections not a positive integer"
			exit 125
		}
		if `fwer' > `K' {
			di as err "option fwer() incorrectly specified"
			di as err "number of false rejections greater than hypotheses"
			exit 125
		}
	}

	// Check fdp() option
	
	if "`fdp'" != "" {
		if real("`fdp'") == . {
			di as err "option fdp() incorrectly specified"
			exit 125
		}
		if `fdp' < 0 | `fdp' >= 1 {
			di as err "option fdp() incorrectly specified"
			di as err "proportion of false discoveries outside unit interval"
			exit 125
		}
		if `fdp' == 0 {
			local fdp ""
			local fwer = 1
		}
	}

	// Parsing

	local rec = cond("`recenter'"   == "", 1, 0)
	local str = cond("`streamline'" == "", 0, 1)

	// Vector of estimates

	if "`stat'" != "" {
	
		if colsof(`stat') != `K' {
			di as err                                         ///
				"vector of observed values does not conform " ///
				"with bootstrap sample"
			exit 503
		}
		
		local names : colnames `stat'		
	}
	else {
	
		tempname stat
		
		local k = 0
		matrix `stat' = J(1, `K', .)
		
		foreach v of varlist `varlist' {
			local k = `k' + 1
			matrix `stat'[1, `k'] = ``v'[observed]'
			local names "`names' ``v'[colname]'"
		}		
	}

	// Standard errors
	
	mata: st_matrix("`S'", stderr("`varlist'"))
	
	matrix colnames `S' = `names'
	matrix rownames `S' = _se

    // Estimate confidence region

    if "`uniform'" != "" {
		local crtype "uniform"
		local opt `""`varlist'", "`stat'", `rec', `level'"'
        mata: st_matrix("`CR'", cr_unif(`opt'))
    }
	else if "`fwer'" != "" {
		local crtype "`fwer'-FWER"
		local opt `""`varlist'", "`stat'", `rec', `level', `fwer', `str'"'
        mata: st_matrix("`CR'", cr_fwer(`opt'))
    }
	else if "`fdp'" != "" {
		local crtype "`fdp'-FDP"
		local opt `""`varlist'", "`stat'", `rec', `level', `fdp', `str'"'
        mata: st_matrix("`CR'", cr_fdp(`opt'))
    }

	matrix colnames `CR' = `names'

	if   (`rec') matrix rownames `CR' = ll ul
	else         matrix rownames `CR' = cv
	
	// Output matrix

	matrix `D' = [`stat' \ `S' \ `CR']
	matrix colnames `D' = `names'
	
	// Post results
	
	ereturn clear

	ereturn local cmd     "bscr"
	ereturn local cmdline `"`cmdline'"'
	ereturn local crtype  "`crtype'"

	if `rec' {
		if "`table'" == "" {
			estout matrix(`D', transpose fmt(4 4 4 4)), ///
				title("Simultaneous confidence region") ///
				mlabels(, none)                         ///
				collabels("Obs. Coef." "Std. Err." "" "Conf. region")
		}
		ereturn matrix cr `CR'
	}
	else {
	if "`table'" == "" {
			estout matrix(`D', transpose fmt(4 4 4)), ///
				title("Simultaneous critical values") ///
				mlabels(, none)                       ///
				collabels("Obs. Coef." "Std. Err." "Crit. values")
		}
		ereturn matrix cv `CR'
	}

	ereturn matrix se `S'

	ereturn scalar level  = `level'
	ereturn scalar N_reps = `N'

end

/******************************************************************************/
// ESTIMATION
/******************************************************************************/

version 13.1
mata:

	// STANDARD ERRORS
	
	real matrix stderr(string scalar varlist)
	{
		real scalar K
		real rowvector S
        real matrix B

        // Load data and allocate space

        st_view(B = ., ., tokens(varlist))
		K = cols(B)
		S = J(1, K, .)
		
		// Standard errors
		
		for (i=1; i<=K; i++) S[i] = sqrt(quadvariance(B[.,i]))
		
		// Post results
		
        return(S)
    }

    // UNIFORM CONFIDENCE BAND
	
		// Algorithm 3 of Chernozhukov, Fernandez-Val and Melly (2013, ECMA)

    real matrix cr_unif(
        string scalar varlist,
        string scalar beta,
		real   scalar recenter,
        real   scalar level
    ) {
        real scalar    t
        real rowvector b, IQR
        real matrix    B, CR

        // Load data

        st_view(B = ., ., tokens(varlist)) ; b = st_matrix(beta)
		
		// Recenter test statistics
		
		if (recenter) {
			B   = B :- b
			IQR = mm_iqrange(B) / (invnormal(.75) - invnormal(.25))
			B   = abs(B) :/ IQR
		}
		else {
			B = abs(B)
		}

        // Compute critical values and construct confidence region

        t  = mm_quantile(rowmax(B), 1, level/100)
		CR = (recenter) ? (b - t :* IQR \ b + t :* IQR) : J(1, cols(b), t)
		
		// Post results
		
        return(CR)
    }
	
	// BALANCED STEP-DOWN K-FWER CONTROL
	
	real matrix cr_fwer(
		string scalar varlist,
        string scalar beta,
		real   scalar recenter,
        real   scalar level,
		real   scalar k,
		real   scalar streamline
	){
        real rowvector b, cv
        real matrix    B, CR

        // Load data and allocate space

        st_view(B = ., ., tokens(varlist))
		b = st_matrix(beta)
		B = (recenter) ? abs(B :- b) : abs(B)
		
		// Call subroutine
		
		cv = aux_fwer(B, abs(b), level, k, streamline)
		
		// Construct confidence region
		
		CR = (recenter) ? (b - cv \ b + cv) : cv
		
		// Post results

        return(CR)
	}

	// BALANCED STEP-DOWN FDP CONTROL
	
		// Algorithm 8.1 of Romano and Wolf (2010, Ann. of Stat.)
	
	real matrix cr_fdp(
		string scalar varlist,
        string scalar beta,
		real   scalar recenter,
        real   scalar level,
		real   scalar g,
		real   scalar streamline
	){
		real scalar    k, r
        real rowvector b, cv
        real matrix    B, CR

        // Load data and allocate space

        st_view(B = ., ., tokens(varlist))
		b = st_matrix(beta)
		B = (recenter) ? abs(B :- b) : abs(B)

		// Sequence of step-down algorithms
		
		k = 0
		
		do {
			k  = k + 1
			cv = aux_fwer(B, abs(b), level, k, streamline)
			r  = sum(abs(b) :> cv)
		} while (r >= k/g - 1)
		
		// Construct confidence region
		
		CR = (recenter) ? (b - cv \ b + cv) : cv
		
		// Post results

        return(CR)
	}
	
	// BALANCED STEP-DOWN K-FWER CONTROL
	
		// Algorithms 4.1 and 4.2 of Romano and Wolf (2010, Ann. of Stat.)
		// Support function for cr_fwer() and cr_fdp()
	
	real matrix aux_fwer(
		real matrix    B,
        real rowvector b,
        real scalar    level,
		real scalar    k,
		real scalar    streamline
	){
		real scalar    d, in, r0, r1, K
        real rowvector cv0, cv1, i, i0, i1, p, I, R
        real matrix    CR

        // Allocate space
		
		d = cols(b)
		i = range(1, d, 1)'
		K = k-1
		I = J(1, K, 0)

		// Single-step critical values
		
		cv1 = aux_cv(B, level, k)
		R   = abs(b) :> cv1 ; r0 = 0 ; r1 = sum(R)

		// Step-down algorithm
		
			// FWER (k == 1)
			
		while ((r1 >= k) & (r1 > r0) & (r1 < d) & (k == 1)) {

			r0 = r1
			i0 = select(i, 1 :- R)

			cv1[i0] = aux_cv(B[.,i0], level, 1)

			R = abs(b) :> cv1 ; r1 = sum(R)
		}

			// k-FWER (k > 1), all permutations
		
		while ((r1 >= k) & (r1 > r0) & (r1 < d) & (k > 1) & (streamline == 0)) {
		
			r0 = r1 ; a = d - r0
			in = comb(r1,K)
			i0 = select(i, 1 :- R)
			i1 = select(i, R)

			for (x=1; x<=in; x++) {

				_aux_comb(I,r1,K,x)

				cv0 = aux_cv(B[.,(i0, i1[I])], level, k)

				if   (x == 1) cv1[i0] = cv0[(1..a)]
				else          cv1[i0] = colmax((cv1[i0] \ cv0[(1..a)]))

			}

			R = abs(b) :> cv1 ; r1 = sum(R)
		}
	
			// k-FWER (k > 1), streamlined
		
		if (streamline != 0) {
			p = 1 :- mm_relrank(B, 1, abs(b))
			p = order(p',1)'
		}
		
		while ((r1 >= k) & (r1 > r0) & (r1 < d) & (k > 1) & (streamline != 0)) {
		
			r0 = r1 ; a = d - r0
			i0 = select(i, 1 :- R)
			i1 = p[(1..r1)][(1..K)]
			
			cv0     = aux_cv(B[.,(i0,i1)], level, k)
			cv1[i0] = cv0[(1..a)]

			R = abs(b) :> cv1 ; r1 = sum(R)
		}

		// Post results

        return(cv1)
	}

    // SINGLE-STEP CRITICAL VALUES FOR BALANCED K-FWER CONTROL
	
		// Remark 3.1 of Romano and Wolf (2010, Ann. of Stat.)
		// Single-step balanced critical values to control the k-FWER
		// Support function for aux_fwer()
	
	real rowvector aux_cv(real matrix B, real scalar level, real scalar k)
	{
		real colvector L
		real matrix    H
		
		H = mm_ecdf(B)
		L = (k == 1) ? rowmax(H) : aux_kmax(H, k)

		return(mm_quantile(B, 1, mm_quantile(L, 1, level/100)))
	}

	// COMBINATION ACCORDING TO LEXICOGRAPHIC INDEX
	
		// Algorithm 515 of ACM TOMS (Buckles and Lybanon, 1977) 
		// Translated from C code by user nlucaroni,
			// available at Stack Overflow (questions 561 and 127704)
		// We build the x-th combination of k out of n elements,
			// filling an array c of size k
		// Support function for aux_cv()

	void _aux_comb(
		real rowvector c,
		real scalar    n,
		real scalar    k,
		real scalar    x
	) {
		real scalar j, r
		
		if (cols(c) != k) _error(3200)
		
		if (k == 1) {
			c[.] = x
		}
		else {
		
			j = 0 ; r = 0
			
			for (i=1; i<k; i++) {
			
				c[i] = (i != 1) ? c[i-1] : 0
				
				do {
					c[i] = c[i] + 1
					r = comb(n-c[i],k-i)
					j = j + r
				} while (j < x)
				
				j = j - r
			}
			
			c[k] = c[k-1] + x - j
		}
	}
	
	// K-TH LARGEST ELEMENT OF EACH ROW OF A MATRIX
	
		// Support function for aux_cv()

	real colvector aux_kmax(real matrix X, real scalar k) {
		real scalar    R
		real rowvector x
		real colvector K
		
		R = rows(X) ; K = J(R, 1, .) ; x = J(cols(X) , 1, .)
		
		for (i=1; i<=R; i++) {
			x[.] = sort(X[i,.]', -1)
			K[i] = x[k]
		}
		
		return(K)
	}
	
end
