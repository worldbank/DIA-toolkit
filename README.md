# DIA-toolkit
This repository contains all the program codes developed in the "Distributional Impact Analysis: Toolkit and Illustrations of Impacts Beyond the Average Treatment Effect" by Guadalupe Bedoya (World Bank), Luca Bittarello (Northwestern University), Jonathan Davis (University of Chicago), and Nikolas Mittag (CERGE-EI) available at [the World Bank Policy Research Working Paper Series]( http://documents.worldbank.org/curated/en/292901499351272899/pdf/WPS8139.pdf)


## How to Use the DIA-Toolkit

To replicate the simulations in the toolkit and the analysis in Part 7, start on `/estimation/master.do`. This master file sets global paths, sets the current directory  and loads auxiliary programs.

### Replicating Simulations throughout Toolkit

All simulated data are created by user-written commands in `/programs/simulation_programs.do`. 
The simulation results in the toolkit are run by `/estimation/simulations/simulations`.

### Replicating Section 7.1

Five files are available under `/estimation/brazil/`:

- `brazil_setup.do` sets up the data (the baseline dataset is available from [AEJ-EP](https://www.aeaweb.org/articles?id=10.1257/app.20150149)).
- `brazil_2.do` replicates Subsection 7.1.2 (treatment effects on summary statistics, weighted and unweighted).
- `brazil_3.do` replicates Subsection 7.1.3 (bounds on features of the distribution of treatment effects).
- `brazil_4.do` replicates Subsection 7.1.4 (conditional effects).
- `brazil_5.do` replicates Subsection 7.1.5 (power calculations). Adjust line 11 to replicate the different columns of Table 6.

These files involve bootstrapping, which may be lengthy. The folder `/results/bootstrap/` contains our bootstrap samples, which you can use to replicate our tables without running the algorithms.

### Replicating Section 7.2

Three files are avaiilable under `/estimation/gambia/`:

- `gambia_setup.do` sets up the data (which was provided by M. Blimpo, D.K. Evans and N. Lahire (2015): WB Policy Research WP 7238).
- `gambia_analysis.do` replicates all of Section 7.2
- `gambia_simulations.do` replicates the deconvolution estimates using simulated data designed to mimic our sample from The Gambia. (This file is called by `/estimation/simulations/simulations.do`.

### Programs

These programs are available in `/programs/`.

- `dte.do` computes changes the distribution of marginal outcomes and its features. `dte.sthlp` documents the syntax (from Stata, run `help dte`).
- `bscr.do` computes simultaneous confidence regions and critical values from a bootstrap sample of estimators according to Romano and Wolf (2010). `bscr.sthlp` documents the syntax (from Stata, run `help bscr`).

## How to Contribute 
If you are familiar with GitHub pelase feel free to fork the repository and submit a *pull request* with your suggested edits.

An easy but still very efficient way to provide any feedback on the code in this repository is to create an *issue* in GitHub. You can read *issues* submitted by other users or create a new *issue* in the top menu below [**worldbank**/**DIA-toolkit**](https://github.com/worldbank/DIA-toolkit) at [https://github.com/worldbank/DIA-toolkit](https://github.com/worldbank/DIA-toolkit). While the word *issue* has a negative connotation outside GitHub, it can be used for any kind of feedback. If you have identified a bug or an enhancement please create a new *issue* to let us know. Your feedback is very important to us. Please read already existing *issues* to check whether someone else has made the same suggestion or reported the same error before creating a new *issue*.

## License

This software is released under the Creative Commons Attribution 4.0 International license. You can view a license summary here:
https://creativecommons.org/licenses/by/4.0/


## How to Reference the DIA-Toolkit

Use the following reference for attribution:

Bedoya, Guadalupe, Luca Bittarello, Jonathan Davis, and Nikolas Mittag. "Distributional Impact Analysis: Toolkit and Illustrations of Impacts Beyond the Average Treatment Effect." [World Bank Policy Research Working Paper Series]( http://documents.worldbank.org/curated/en/292901499351272899/pdf/WPS8139.pdf), No. WPS 8139. Washington, D.C. (July 2017)

