/* This program cleans and combines the data files from the WSDP evaluation */


* Covariates *

use "${data_dir}/gambia/gambia_headteacher_cleaned2008.dta", clear

gen building_good = BUILDING_HT==1
gen twoShifts = Q202_HT==1

gen students = Q203AM_HT if twoShifts==0
 replace students = 0.5*(Q203AM_HT+Q203PM_HT) if twoShifts==1
gen students_ms = missing(students)
 replace students = 0 if students_ms==1
gen teachers = Q204AM_HT if twoShifts==0
 replace teachers = 0.5*(Q204AM_HT+Q204PM_HT) if twoShifts==1
ren Q207_HT classrooms

gen trackExpenses = Q304_HT == 1 | Q304_HT==2
gen staffConduct = Q404_HT==1 | Q404_HT==2
gen teacherConduc = Q420_HT==1 | Q420_HT==2
gen pta_scheduled = Q508_HT==8888 
gen community_support = Q500_HT==1
gen dev_plan = Q509_HT==1 | Q509_HT==2
keep schoolID building_good-dev_plan

replace schoolID = (schoolID +111)/11 

tempfile Xs
save `Xs'




* Save Treatment Status from Head Teacher Baseline Survey *
* schoolID is unique within year *


use "${data_dir}/gambia/gambia_headteacher_cleaned2008.dta", clear
gen year = 2008
ren CL_CODE_HT subcluster

keep schoolID treatment year subcluster

replace schoolID = (schoolID +111)/11 
tempfile id2008
save `id2008'

clear
use "${data_dir}/gambia/gambia_headteacher_followup2009.dta"
 gen year = 2009
 ren T_RECEIVE_HT treatment
  replace treatment = 3 if treatment==2
 bys schoolID: egen seq=seq()
 keep if seq==1
 replace schoolID = (schoolID +111)/11 
 keep schoolID treatment year
 tempfile id2009
 save `id2009'
 
clear
use "${data_dir}/gambia/gambia_HeadTeacher2010_PUF.dta"
 gen year = 2010
 keep schoolID treatment year
 tempfile id2010
 save `id2010'

clear
use "${data_dir}/gambia/gambia_HeadTeacher2011_PUF.dta"
 ren t_receiv treatment
 ren cde_sch schoolID
 ren sub_clus subcluster
 gen year=2011
 keep schoolID treatment year subcluster
tempfile id2011
save `id2011'
 


use `id2008'
 append using `id2009'
 append using `id2010'
 append using `id2011'
 

gen wsd = treatment==1
gen grant = treatment==2
gen wsd08 = wsd==1 & year==2008
gen grant08 = grant==1 & year==2008
bys schoolID: egen wsd2008 = max(wsd08)
bys schoolID: egen grant2008 = max(grant08)
gen control2008 = wsd2008==0 & grant2008==0
drop wsd08 grant08

 
save "${data_dir}/gambia/treatment.dta", replace




* ------------------------------------------------------- *
* Make Numeracy/Literacy Test Files Uniform, then Combine *
* ------------------------------------------------------- * 

use "${data_dir}/gambia/gambia_numlit_cleaned2008.dta", clear

ren BOY_GIRL_NL genderStudent_nl 
ren GRADE_NL grade_student_nl

forv i=1/32 {
 egen c`i' = mode(Q`i'_NL)
 gen n`i' = Q`i'_NL==c`i'
 drop c`i' 
}
gen n=1
local i=1
quietly {
foreach var of varlist WM1_NL MW_55_NL {
  bys `var': egen count=sum(n)
  egen mode = max(count)
  gsort -count
  gen correct = `var'[1]
  gen l`i' = `var'==correct
  drop count mode correct
  local i=`i'+1
 }
}
drop n
egen score_num_NL = rowtotal(n*)
egen score_lit_NL = rowtotal(l*)
egen score_total_NL = rowtotal(n* l*)
 replace score_num_NL = 100*score_num_NL / 32
 replace score_lit_NL = 100*score_lit_NL / 55
 replace score_total_NL = 100*score_total_NL / 87
 
foreach var of varlist score* {
 qui su `var'
 gen z_`var' = (`var'-r(mean))/r(sd)
}

gen year = 2008
* Moussa Blimpo told us this fix to school IDs across years on April 3rd, 2016
replace schoolID = (schoolID +111)/11 

save "${data_dir}/gambia/gambia_numeracyliteracy_cleaned2008.dta", replace


use  "${data_dir}/gambia/gambia_numeracyliteracy_followup2009.dta", clear
ren studentNumber_NL qnumb 
ren BOY_GIRL_NL genderStudent_nl 
ren GRADE_NL grade_student_nl
ren tract trackedStudents_nl

replace schoolID = (schoolID +111)/11 
forv i=1/32 {
 egen c`i' = mode(Q`i'_NL)
 gen n`i' = Q`i'_NL==c`i'
 drop c`i' 
}
gen n=1
local i=1
quietly {
foreach var of varlist AT_NL - MW_55 {
  bys `var': egen count=sum(n)
  egen mode = max(count)
  gsort -count
  gen correct = `var'[1]
  gen l`i' = `var'==correct
  drop count mode correct
  local i=`i'+1
 }
}
drop n
egen score_num_NL = rowtotal(n*)
egen score_lit_NL = rowtotal(l*)
egen score_total_NL = rowtotal(n* l*)
 replace score_num_NL = 100*score_num_NL / 32
 replace score_lit_NL = 100*score_lit_NL / 55
 replace score_total_NL = 100*score_total_NL / 87
 
foreach var of varlist score* {
 su `var'
 gen z_`var' = (`var'-r(mean))/r(sd)
}

 
gen year = 2009
append using "${data_dir}/gambia/gambia_numeracyliteracy_cleaned2008.dta"
append using "${data_dir}/gambia/gambia_StudentLiteracyNumeracy2010_PUF.dta"
 replace year = 2010 if missing(year)
 replace genderStudent_nl = genderStudent_nl10
 drop genderStudent_nl10
append using "${data_dir}/gambia/gambia_StudentLiteracyNumeracy2011_PUF.dta"
 replace year = 2011 if missing(year)
 replace schoolID = IDschool if year==2011
 replace genderStudent_nl = genderStudent_nl11
 drop genderStudent_nl11

 
* Merge with Treatment Status * 
 

merge m:1 schoolID year using "${data_dir}/gambia/treatment.dta", gen(mergeTreat)
 *keep if mergeTreat==3
merge m:1 schoolID  using `Xs', gen(mergeX)


gen grade4 = grade==4

drop if missing(treatment)


save "${data_dir}/gambia/analysis_data.dta", replace



exit
