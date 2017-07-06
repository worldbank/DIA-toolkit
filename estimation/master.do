/******************************************************************************/
// NOTES
/******************************************************************************/

* This program organizes the empirical exercises in the DIA Toolkit.
* Last updated: June 17, 2017

/******************************************************************************/
// MISE EN PLACE
/******************************************************************************/

// Main directory

global main_dir = "`c(pwd)'\.."

// Directory macros

global estimation_dir "${main_dir}/estimation"
global program_dir    "${main_dir}/programs"
global results_dir    "${main_dir}/results"
global bs_dir         "${main_dir}/bootstrap"
global data_dir       "${main_dir}/data"

// Other preliminaries

clear all
set more off
set matsize 11000

//Packages required

if 1 { //change to 0 if you have all these packages installed
	
	ssc install kdens
	ssc install moremata
}

/******************************************************************************/
// SIMULATIONS
/******************************************************************************/

* Last updated: June 15, 2017
* Author: Jonathan M.V. Davis
* Email: jonmvdavis@gmail.com

// Load user-written commands

run "${program_dir}/makarov.ado"

// Load user-written commands

run "${program_dir}/simulation_programs.do"

// Analysis

do "${estimation_dir}/simulations/simulations.do"

/******************************************************************************/
// SECTION 7.1
/******************************************************************************/

* Financial literacy in Brazil

* Last updated: June 15, 2017
* Author: Luca Bittarello
* Email: luca.bittarello@gmail.com

* The data and original study are available from the AEJ AE:
* http://dx.doi.org/10.1257/app.20150149

* Original study:
* M. Bruhn, L. de Souza Leao, A. Legovini, R. Marchetti and B. Zia (2016):
* American Economic Journal: Applied Economics 8(4), 256-295.

// Load user-written commands

run "${program_dir}/dte.ado"
run "${program_dir}/bscr.ado"

// Analysis

* One must adjust line 11 and rerun file brazil_5.do
* to fully replicate our power calculations

do "${estimation_dir}/brazil/brazil_2.do"
do "${estimation_dir}/brazil/brazil_3.do"
do "${estimation_dir}/brazil/brazil_4.do"
do "${estimation_dir}/brazil/brazil_5.do"

/******************************************************************************/
// SECTION 7.2
/******************************************************************************/

* Whole School Development Program in The Gambia

* Last updated: June 15, 2017
* Author: Jonathan M.V. Davis
* Email: jonmvdavis@gmail.com

* The data are available from the World Bank:
* http://microdata.worldbank.org/index.php/catalog/2523

* Original study:
* M. Blimpo, D.K. Evans and N. Lahire (2015): WB Policy Research WP 7238

// Analysis

do "${estimation_dir}/gambia/gambia_setup.do"
do "${estimation_dir}/gambia/gambia_analysis.do"

/******************************************************************************/
// EXIT
/******************************************************************************/

exit
