/**************************************************************************************************/
// POWER CALCULATIONS
/**************************************************************************************************/

run "${estimation_dir}/brazil/brazil_setup.do"

// Parameters

local size   = 0.05
local power  = 0.8
local effect "+ 1"

// Prepare data

keep if treatment == 0

gen y0 = outcome
gen y1 = outcome `effect'

// Power calculations

forvalues k = 0/1 {
    pctile q`k' = y`k', nq(100) genp(t`k')
    kdens y`k', adaptive bw(sjpi) at(q`k') gen(f`k') nograph
    gen s`k' = t`k'*(100 - t`k')/(100*f`k')^2
}

gen p = s1/(s0 + s1)
gen n = round((((invnormal(1-`size'/2) + invnormal(`power'))/(q1-q0))^2)*(s0/(1-p) + s1/p))

keep t0 p n
keep if inlist(t0, 10, 20, 30, 40, 50, 60, 70, 80, 90)

ed

/**************************************************************************************************/
// MHT ADJUSTMENT
/**************************************************************************************************/

use "${bs_dir}/b_sample.dta", clear

global vars = ""

foreach v in mean p_10 p_20 p_30 p_40 p_50 p_60 p_70 p_80 p_90 {
    global vars = "$vars _b_`v'"
}

foreach v of varlist * {
    replace `v' = abs(`v' - ``v'[observed]')
}

foreach v of varlist $vars {
    qui sum `v', det
    di `r(p95)'
}

bscr,       norecenter
bscr $vars, norecenter
