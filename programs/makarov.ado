cap program drop makarov
program define makarov , rclass
	version 12
	
	syntax varlist(min=2 max=2 numeric) , Grid(numlist ) [d(real 0.0)]

	gettoken var1 var0 : varlist
	
	tempvar gridmat 

	numlist "`grid'"
	local nlist `r(numlist)'
	
	foreach el of local nlist {
		matrix `gridmat' = nullmat(`gridmat'),`el'
	}
	
	matrix `gridmat' = `gridmat''
	matrix list `gridmat'

	qui count if `var1'<.
	local n1 = r(N)
	qui count if `var0'<.
	local n0 = r(N)
	 
	mata: m_makarov("`var1'", "`var0'", "`gridmat'", `d', `n1', `n0') 

	return scalar bound_l = r(bound_l)
	return scalar bound_u = r(bound_u)
end 
 
capture mata: mata drop m_makarov()
mata: 

	void m_makarov(
		string scalar var1,
		string scalar var0,
		string scalar grid,
		real   scalar d,
		real   scalar n1,
		real   scalar n0
	){
		G = st_matrix(grid)
		d = J(rows(G),1,d)
		Gd = G-d
		Y1 = (st_data(.,var1))[1::n1,1]
		Y0 = (st_data(.,var0))[1::n0,1]

		F1=J(rows(G),1,0)
		F0=J(rows(G),1,0)
		one = J(rows(G),1,1)

		for (i=1; i<=rows(G); i++) {

		F1[i] = mean(Y1 :<=G[i,1])
		F0[i] = mean(Y0 :<=Gd[i,1])
		}


		diff_l = F1-F0
		diff_u = F1-F0+one

		bound_l = max((diff_l \ 0))
		st_numscalar("r(bound_l)", bound_l)
		bound_u = min((diff_u \ 1))
		st_numscalar("r(bound_u)", bound_u)
	}
end



