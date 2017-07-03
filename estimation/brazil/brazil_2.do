/**************************************************************************************************/
// REPLICATE COLUMN 1 OF TABLE 3 OF BRUHN ET AL. (2016)
/**************************************************************************************************/

run "${estimation_dir}/brazil/brazil_setup.do"

reg  outcome treatment                                   , cl(cd_escola)
areg outcome treatment                                   , cl(cd_escola) a(pair_fe)
areg outcome treatment baseline missing_bl i.female_coded, cl(cd_escola) a(pair_fe)

/**************************************************************************************************/
// DISTRIBUTIONAL EFFECTS
/**************************************************************************************************/

run "${estimation_dir}/brazil/brazil_setup.do"

// Unweighted

dte outcome treatment,                        ///
    stat(mean sd pr_75_25) p(5/95) cdf(33/81) ///
    cluster(cd_escola) reps(5000) seed(0101)  ///
    saving("${bs_dir}/b_sample.dta", replace)

estat bootstrap, percentile
bscr using "${bs_dir}/b_sample.dta", fwer(1)

// Weighted for gender and baseline outcome (appendix)

dte outcome treatment i.female_coded i.missing_bl c.baseline##c.baseline, ///
    stat(mean sd pr_75_25) p(5/95) cdf(33/81)                             ///
    cluster(cd_escola) reps(5000) seed(0101)                              ///
    saving("${bs_dir}/w_sample.dta", replace)

estat bootstrap, percentile
bscr using "${bs_dir}/w_sample.dta", fwer(1)
