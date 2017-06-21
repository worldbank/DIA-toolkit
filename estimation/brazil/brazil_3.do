/**************************************************************************************************/
// BOUNDS ON FEATURES OF THE DISTRIBUTION OF TREATMENT EFFECTS
/**************************************************************************************************/

run "${estimation_dir}/brazil/brazil_setup.do"

// Standard deviation of individual treatment effects and P(Y1 - Y0 > 0)

gen nontreatment = 1 - treatment

mata:
	st_view(Y0 = ., ., tokens("outcome"), "nontreatment")
	st_view(Y1 = ., ., tokens("outcome"), "treatment")
	
	y  = range(20, 105, 0.01)	
	F0 = mm_relrank(Y0, 1, y)
	F1 = mm_relrank(Y1, 1, y)
	D  = F1 - F0
	
	st_matrix("F0", F0)
	st_matrix("F1", F1)
	
	st_numscalar("mlb", - min((0, min(F1 - F0))))
	st_numscalar("mub", 1 - max((0, max(F1 - F0))))
	
	st_numscalar("V0", quadvariance(Y0))
	st_numscalar("V1", quadvariance(Y1))
end

di _n 																			///
	"Bounds on features of the distribution of treatment effects" _n            ///
	"Share of winners: [" round(100*mlb, 0.01) ", " round(100*mub, 0.01) "]" _n ///
	"Standard deviation: [" round(sqrt(V0 + V1 - 2*sqrt(V0*V1)), 0.01) ", " 	///
		round(sqrt(V0 + V1 + 2*sqrt(V0*V1)), 0.01) "]"  _n  					///
	"Std. dev. (pos. corr.): [" round(sqrt(V0 + V1 - 2*sqrt(V0*V1)), 0.01) ", " ///
		round(sqrt(V0 + V1), 0.01) "]"
