/******************************************************************************/
// SET UP DATA (FIRST ROUND)

use "${data_dir}/brazil/school_intervention_panel_final.dta", clear

keep if round == 0 & !missing(vl_proficiencia_fup) & !missing(treatment)

bys pair_all : egen flag = mean(treatment) if !missing(vl_proficiencia_fup)

gen pair_fe    = cond(flag == 0 | flag == 1, 0, pair_all)
gen outcome    = vl_proficiencia_fup
gen baseline   = cond(missing(vl_proficiencia_bl), 0, vl_proficiencia_bl)
gen missing_bl = missing(vl_proficiencia_bl)

