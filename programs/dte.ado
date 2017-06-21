/******************************************************************************/
// INTERFACE, PARSING AND BOOTSTRAP
/******************************************************************************/

program dte, eclass
    version 13.1
	
	if replay() {
		if ("`e(cmd)'" != "dte") error 301
		syntax , Level(cilevel)
		ereturn display, level(`level')
		exit
	}

    syntax anything [fw aw pw iw/] [if] [in] [, ///
           Statistics(namelist)                 ///
           Percentiles(numlist)                 ///
           CDF(numlist)                         ///
		   PROB(varname)                        ///
           EFFECT(string)                       ///
           MODEL(string)                        ///
           TRIM(real 0.0)                       ///
           REPS(integer 1000)                   ///
		   TITLE(string)                        ///
		   *]

    // Reserve names

    tempvar w

    // Check requirements (Moremata)

    capture findfile lmoremata.mlib
    if _rc {
        di as err "-moremata- is required; type {stata ssc install moremata}"
        error 499
    }

    // Parse input variables (outcome, treatment, instrument and controls)

    _iv_parse `anything'

	local outcome    `s(lhs)'
	local treatment  `s(endog)'
	local control    `s(exog)'
	local instrument `s(inst)'
	
	if `: list sizeof treatment' > 1 {
		di as err `"syntax is "y t [X]" or "y (t = z [X])""'
		exit 198
	}
	if "`treatment'" != "" & "`control'" != "" {
		di as err `"syntax is "y t [X]" or "y (t = z [X])""'
		exit 198
	}
	
	if ("`control'" != "")    gettoken treatment  control : control
	if ("`instrument'" != "") gettoken instrument control : control
	
	// Check input variables
	
	foreach v in treatment instrument {
		if "``v''" != "" {
			qui count if ``v'' != 0 & ``v'' != 1
			if `r(N)' > 1 {
				di as err "`v' not binary"
				exit 149
			}
		}
	}

	// Switches
	
	local p = ("`prob'" != "")
	local x = ("`control'" != "")
	local z = ("`instrument'" != "")
	local r = `p' + `x' + `z'
	
	// Mark sample

    marksample touse

    markout `touse' `outcome'
    markout `touse' `treatment'

    if (`p' == 1) markout `touse' `prob'
	if (`x' == 1) markout `touse' `control'
	if (`z' == 1) markout `touse' `instrument'

    // Check consistency

    if `p' == 1 & `x' == 1 {
        di as err "option prob() may not be combined with control variables"
        exit 198
    }
	if `r' == 0 & "`effect'" != "" {
        di as err "option effect() only available for IPW; see help file"
        exit 198
    }
    if `x' == 0 & "`model'" != "" {
        di as err "option model() requires control variables"
        exit 197
    }
    if `r' == 0 & `trim' != 0.0 {
        di as err "option trim() only available for IPW and IV; see help file"
        exit 198
    }
	
	// Parse weight
	
	if ("`exp'" != "") local weight "[`weight' = `exp']"

	// Parse objects

    if "`statistics'" == "" & "`percentiles'" == "" & "`cdf'" == "" {
        local statistics = "mean"
    }

    if "`statistics'" == "" {
        local obj "J(0, 0, NULL)"
    }
    else {
        local names "`names' `statistics'"
        local statistics : subinstr local statistics " " "(), &the_", all
        local obj "`obj' (&the_`statistics'())"
    }

    if "`percentiles'" == "" {
        local obj "`obj', J(0, 0, .)"
    }
    else {
        local temp : subinstr local percentiles "." "_", all
        local temp : subinstr local temp " " " p_", all
        local names "`names' p_`temp'"
        local percentiles : subinstr local percentiles " " ", ", all
        local obj "`obj', (`percentiles')"
    }

    if "`cdf'" == "" {
        local obj "`obj', J(0, 0, .)"
    }
    else {
        local temp : subinstr local cdf "." "_", all
        local temp : subinstr local temp " " " d_", all
        local names "`names' d_`temp'"
        local cdf : subinstr local cdf " " ", ", all
        local obj "`obj', (`cdf')"
    }
	
	// Parse effect

    if `r' == 1 {
	
        if ("`effect'" == "") local effect "treated"
		
		local ze "complier treated"
		local xe "current population treated"
		
		if `: list sizeof effect' > 1 {
            di as err "option effect() incorrectly specified"
            exit 198
        }		
        if `z' == 1 & `: list posof `"`effect'"' in ze' == 0 {
            di as err "option effect() incorrectly specified"
            exit 198
        }
		if `x' == 1 & `: list posof `"`effect'"' in xe' == 0 {
            di as err "option effect() incorrectly specified"
            exit 198
        }
    }


	// Parse bootstrap options

	if `reps' > 1 {
		if ("`title'" == "") local title "Distributional treatment effects"
		local boot "bootstrap, `options' reps(`reps') title(`title') :"
	}
	
	// Set up estimation
	
	local e_opts "outcome(`outcome')"
	local e_opts "`e_opts' treatment(`treatment')"
	local e_opts "`e_opts' object(`obj')"
	local e_opts "`e_opts' names(`names')"

    // Set up reweighting

    if `p' == 1 {

        local w_opts "gen(`w')"
        local w_opts "`w_opts' treatment(`treatment')"
        local w_opts "`w_opts' prob(`prob')"
        local w_opts "`w_opts' effect(`effect')"
		local w_opts "`w_opts' trim(`trim')"
		
		if (`z' == 1) local w_opts "`w_opts' instrument(`instrument')"

        _s_weight `weight' if `touse', `w_opts'

		local w_opts "reweight(`w')"
		local weight ""
    }
    else if `x' == 1 {

		if ("`model'" == "") local model = "logit"

        gettoken model w_opts : model, parse(",")

		local w_opts "options(`w_opts')"
		local w_opts "`w_opts' model(`model')"
        local w_opts "`w_opts' control(`control')"
        local w_opts "`w_opts' effect(`effect')"
		local w_opts "`w_opts' trim(`trim')"
		
		if (`z' == 1) local w_opts "`w_opts' instrument(`instrument')"		
    }
	else if `z' == 1 {
		local w_opts "instrument(`instrument')"
		local w_opts "`w_opts' effect(`effect')"
	}

	// Estimation

	`boot' _s_`z'_`x' `weight' if `touse', `e_opts' `w_opts'
	
	// Post results

	ereturn local cmdline `"`0'"'
	ereturn local command `"`0'"'
	ereturn local cmd     "dte"
	ereturn local cmdname "dte"

	ereturn repost

end

/******************************************************************************/
// REWEIGHTING FOR EXOGENOUS CONDITIONAL PROBABILITY
/******************************************************************************/

program _s_weight
    version 13.1

    syntax [fw aw pw iw/] [if] [in],    ///
            Generate(name)              ///
            TREATment(varname numeric)  ///
			PROB(varname numeric)       ///
           [INSTrument(varname numeric) ///
            TRIM(real 0.0)              ///
            EFFECT(name)]

    // Mark sample

    marksample touse

    // Construct appropriate weight

    qui if "`instrument'" != "" & "`effect'" == "complier" {
		gen `generate' = (2*`treatment' - 1) * ///
			(`instrument' -`prob')/(`prob'*(1 - `prob'))
    }
	qui else if "`instrument'" != "" & "`effect'" == "treated" {
        gen `generate' = `treatment' + ///
			(1 - `treatment')*(`prob' - `instrument')/(1 - `prob')
    }
    qui else if "`effect'" == "population" {
        gen `generate' = `treatment'/`prob' + (1 - `treatment')/(1 - `prob')
    }
    qui else if "`effect'" == "treated" {
        gen `generate' = `treatment' + (1 - `treatment')*`prob'/(1 - `prob')
    }
	
	// Trim selection probability

    qui replace `generate' = . if `prob' < `trim' | `prob' > 1 - `trim'
	qui replace `generate' = . if `touse' == 0

    // Adjust for existing weights

    if ("`exp'" != "") qui replace `generate' = `generate' * `exp'

end

/******************************************************************************/
// ESTIMATION SUBROUTINES
/******************************************************************************/

// NO INSTRUMENT
// NO OR EXOGENOUS REWEIGHTING

program _s_0_0, eclass
    version 13.1

    syntax [fw aw pw iw/] [if] [in],   ///
            OUTCOME(varname numeric)   ///
            TREATment(varname numeric) ///
            OBJect(string)             ///
            NAMes(string) 			   ///
		   [REWEIGHT(varname numeric)]

    // Reserve names

    tempvar T0 T1
    tempname b

    // Mark sample

    marksample touse

    // Split treated and untreated observations

    qui gen `T0' = `touse'*(1 - `treatment')
    qui gen `T1' = `touse'*`treatment'

    // Check number of observations

    qui count if `T0'
    if (`r(N)' == 0) error 2000
    local N0 = `r(N)'

    qui count if `T1'
    if (`r(N)' == 0) error 2000
    local N1 = `r(N)'

    local N = `N0' + `N1'

    // Estimation
	
	if      ("`reweight'" != "") local object `"`object', "`reweight'""'
	else if ("`exp'" != "")      local object `"`object', "`exp'""'
    
    mata: st_matrix("`b'", m_dte("`outcome'", "`T0'", "`T1'", `object'))
    matrix colnames `b' = `names'
	
	// Return results
	
    ereturn post `b', esample("`touse'") depname("`outcome'") obs(`N')

end

// INSTRUMENT
// NO REWEIGHTING

program _s_1_0, eclass
    version 13.1

    syntax [fw aw pw iw/] [if] [in],    ///
            OUTCOME(varname numeric)    ///
            TREATment(varname numeric)  ///
            INSTrument(varname numeric) ///
            OBJect(string)              ///
            NAMes(string)

    // Reserve names

    tempvar T0 T1 w
    tempname b

    // Mark sample

    marksample touse
	
	// Split treated and untreated observations

    qui gen `T0' = `touse'*(1 - `treatment')
    qui gen `T1' = `touse'*`treatment'

    // Check number of observations

    qui count if `T0'
    if (`r(N)' == 0) error 2000
    local N0 = `r(N)'

    qui count if `T1'
    if (`r(N)' == 0) error 2000
    local N1 = `r(N)'

	local N = `N0' + `N1'
	
    // Parse weight

    if ("`exp'" != "") local weight "[`weight' = `exp']"
	
	// Reweight

    summarize `instrument' `weight', meanonly
	
	// Construct appropriate weight
	
	qui if "`effect'" == "complier" {
		gen `w' = (2*`T1' - 1)*(`instrument' - r(mean))/(r(mean)*(1 - r(mean)))
    }
	qui else {
        gen `w' = `T1' + `T0'*(r(mean) - `instrument')/(1 - r(mean))
    }

    // Adjust for existing weights

    if ("`exp'" != "") replace `w' = `w' * `exp'

    // Estimation

    mata: st_matrix("`b'", m_dte("`outcome'", "`T0'", "`T1'", `object', "`w'"))
    matrix colnames `b' = `names'
	
	// Return results
	
    ereturn post `b', esample("`touse'") depname("`outcome'") obs(`N')

end

// NO INSTRUMENT
// REWEIGHTING

program _s_0_1, eclass
    version 13.1

    syntax [fw aw pw iw/] [if] [in],    ///
            OUTCOME(varname numeric)    ///
            TREATment(varname numeric)  ///
            OBJect(string)              ///
            NAMes(string)               ///
            EFFECT(name)                ///
            CONTROL(varlist numeric fv) ///
            MODEL(string)               ///
           [OPTIONS(string)             ///
            TRIM(real 0.0)]

    // Reserve names

    tempvar T0 T1 w
    tempname b

    // Mark sample

    marksample touse

    // Parse weight

    if ("`exp'" != "") local weight "[`weight' = `exp']"

    // Reweight

    qui `model' `treatment' `control' `weight' if `touse' `options'
    qui predict `w'
    qui replace `w' = . if `w' < `trim' | `w' > 1 - `trim'

	// Split treated and untreated observations

	markout `touse' `w'
	
    gen `T0' = `touse'*(1 - `treatment')
    gen `T1' = `touse'*`treatment'
	
	// Check number of observations

    qui count if `T0'
    if (`r(N)' == 0) error 2000
    local N0 = `r(N)'

    qui count if `T1'
    if (`r(N)' == 0) error 2000
    local N1 = `r(N)'

    local N = `N0' + `N1'
	
	// Construct appropriate weight

    qui if ("`effect'" == "population") replace `w' = `T1'/`w' + `T0'/(1 - `w')
    qui else replace `w' = `T1' + `T0'*`w'/(1 - `w')

     // Adjust for existing weights

    if ("`exp'" != "") replace `w' = `w' * `exp'

    // Estimation

    mata: st_matrix("`b'", m_dte("`outcome'", "`T0'", "`T1'", `object', "`w'"))
    matrix colnames `b' = `names'
	
	// Return results
	
    ereturn post `b', esample("`touse'") depname("`outcome'") obs(`N')

end

// INSTRUMENT
// REWEIGHTING

program _s_1_1, eclass
    version 13.1

    syntax [fw aw pw iw/] [if] [in],    ///
            OUTCOME(varname numeric)    ///
            TREATment(varname numeric)  ///
            INSTrument(varname numeric) ///
            OBJect(string)              ///
            NAMes(string)               ///
            EFFECT(name)                ///
            CONTROL(varlist numeric fv) ///
            MODEL(string)               ///
           [OPTIONS(string)             ///
            TRIM(real 0.0)]

    // Reserve names

    tempvar T0 T1 w
    tempname b

    // Mark sample

    marksample touse

    // Parse weight

    if ("`exp'" != "") local weight "[`weight' = `exp']"

    // Reweight

    qui `model' `instrument' `control' `weight' if `touse' `options'
    qui predict `w'
    qui replace `w' = . if `w' < `trim' | `w' > 1 - `trim'
	
	// Split treated and untreated observations

	markout `touse' `w'
	
    gen `T0' = `touse'*(1 - `treatment')
    gen `T1' = `touse'*`treatment'

    // Check number of observations

    qui count if `T0'
    if (`r(N)' == 0) error 2000
    local N0 = `r(N)'

    qui count if `T1'
    if (`r(N)' == 0) error 2000
    local N1 = `r(N)'

    local N = `N0' + `N1'
	
	// Construct appropriate weight
	
	qui if "`effect'" == "complier" {
		gen `w' = (2*`T1' - 1)*(`instrument' - `w')/(`w'*(1 - `w'))
    }
	qui else {
        gen `w' = `T1' + `T0'*(`w' - `instrument')/(1 - `w')
    }

     // Adjust for existing weights

    if ("`exp'" != "") replace `w' = `w' * `exp'

    // Estimation

    mata: st_matrix("`b'", m_dte("`outcome'", "`T0'", "`T1'", `object', "`w'"))
    matrix colnames `b' = `names'
	
	// Return results
	
    ereturn post `b', esample("`touse'") depname("`outcome'") obs(`N')

end

/******************************************************************************/
// MATA SUBROUTINES
/******************************************************************************/

version 13.1
mata:

    // ESTIMATION

    real rowvector m_dte(
        string scalar outcome,
        string scalar untreated,
        string scalar treated,
        pointer(real scalar function) matrix stat,
        real matrix qte,
        real matrix cdf,
        | string scalar weight
    ) {
        real scalar N_s, N_q, N_c, d0, d1
        real colvector Y0, Y1, W0, W1, D

        // Number of objects

        N_s = length(stat) ; N_q = length(qte) ; N_c = length(cdf)

        // Outcome data

        st_view(Y0 = ., ., tokens(outcome), untreated)
        st_view(Y1 = ., ., tokens(outcome), treated)

        // Weights

        if (args() == 6) {
            W0 = 1
            W1 = 1
        }
        else {
            st_view(W0 = ., ., tokens(weight), untreated)
            st_view(W1 = ., ., tokens(weight), treated)
        }

        // Estimation

        d1 = 0
        D  = J(1, N_s + N_q + N_c, .)

        if (N_s > 0) {
            for (i = 1; i <= N_s; i++) {
                d1   = d1 + 1
                D[i] = (*stat[i])(Y1, W1) - (*stat[i])(Y0, W0)
            }
        }
        if (N_q > 0) {
            qte = qte/100
            d0 = d1+1
            d1 = d1 + N_q
            D[(d0 .. d1)] = mm_quantile(Y1, W1, qte) - mm_quantile(Y0, W0, qte)
        }
        if (N_c > 0) {
            d0 = d1+1
            d1 = d1 + N_c
            D[(d0 .. d1)] = mm_relrank(Y1, W1, cdf) - mm_relrank(Y0, W0, cdf)
        }

        return(D)
    }

    // SUMMARY STATISTICS

    real scalar the_mean(real colvector x, | real colvector w) {
        if (args()==1) w = 1
        return(mean(x,w))
    }
    real scalar the_var(real colvector x, | real colvector w) {
        if (args()==1) w = 1
        return(quadvariance(x,w))
    }
    real scalar the_sd(real colvector x, | real colvector w) {
        if (args()==1) w = 1
        return(sqrt(quadvariance(x,w)))
    }
	real scalar the_skewness(real colvector x, | real colvector w) {
		real vector xx
        if (args()==1) w = 1
        xx = x :- mean(x)
		return(mean(xx:^3) / (mean(xx:^2)^(3/2)))
    }
	real scalar the_kurtosis(real colvector x, | real colvector w) {
        real vector xx
        if (args()==1) w = 1
        xx = x :- mean(x)
		return(mean(xx:^3) / (mean(xx:^2)^2))
    }
    real scalar the_gini(real colvector x, | real colvector w) {
        if (args()==1) w = 1
        return(mm_gini(x,w))
    }
    real scalar the_cv(real colvector x, | real colvector w) {
        real colvector Y
        if (args()==1) w = 1
        Y = quadmeanvariance(x, w)
        return(sqrt(Y[1,1])/Y[2,1])
    }
    real scalar the_iqrange(real colvector x, | real colvector w) {
        if (args()==1) w = 1
        return(mm_iqr(x,w))
    }

    // QUANTILE RATIOS

    real scalar qratio(
        real scalar n, real scalar d, real colvector x, | real colvector w
    ) {
        real matrix Q
        if (args()==1) w = 1
        Q = mm_quantile(x, w, (n, d))
        return(Q[1]/Q[2])
    }

    real scalar the_pr_90_10(real colvector x, | real colvector w) {
        if (args()==1) w = 1
        return(qratio(0.9, 0.1, x, w))
    }
    real scalar the_pr_50_10(real colvector x, | real colvector w) {
        if (args()==1) w = 1
        return(qratio(0.5, 0.1, x, w))
    }
    real scalar the_pr_90_50(real colvector x, | real colvector w) {
        if (args()==1) w = 1
        return(qratio(0.9, 0.5, x, w))
    }
    real scalar the_pr_75_25(real colvector x, | real colvector w) {
        if (args()==1) w = 1
        return(qratio(0.75, 0.25, x, w))
    }

end
