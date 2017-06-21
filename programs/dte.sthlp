{smcl}
{* *! version 1.0.0  10jan2017}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{viewerjumpto "Syntax" "dte##syntax"}{...}
{viewerjumpto "Description" "dte##description"}{...}
{viewerjumpto "Options" "dte##options"}{...}
{viewerjumpto "Remarks" "dte##remarks"}{...}
{viewerjumpto "Examples" "dte##examples"}{...}
{viewerjumpto "Stored results" "dte##results"}{...}
{viewerjumpto "References" "dte##references"}{...}
{title:Title}

{phang}
{bf:dte} {hline 2} Unconditional distributional treatment effects for binary
	treatments


{marker syntax}{...}
{title:Syntax}

{p 4}Baseline

{p 8 17 2}
{cmdab:dte}
{it:{help varlist:outcome}}
{it:{help varlist:treatment}}
{ifin}
{weight}
[{cmd:,} {it:options}]

{p 4}Inverse probability weighting

{p 8 17 2}
{cmdab:dte}
{it:{help varlist:outcome}}
{it:{help varlist:treatment}}
{it:{help varlist:controls}}
{ifin}
{weight}
[{cmd:,} {it:options}]

{p 4}Instrumental variable

{p 8 17 2}
{cmdab:dte}
{it:{help varlist:outcome}}
{cmd:(}{it:{help varlist:treatment}}
{cmd:=} 
{it:{help varlist:instrument}}
[{it:{help varlist:controls}}]{cmd:)}
{ifin}
{weight}
[{cmd:,} {it:options}]

{synoptset 27 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{cmdab:s:tatistics:(}{it:{help dte##dte:statname}} [{it:...}]{cmd:)}}estimate
	effects on specified statistics{p_end}
{synopt:{opth p:ercentiles(numlist)}}estimate effects on specified
	percentiles{p_end}
{synopt:{opth cdf:(numlist)}}estimate effects on the CDF at specified outcome
	values{p_end}

{syntab:IPW and IV}
{synopt:{cmdab:effect(}{it:{help dte##dte:effecttype}}{cmd:)}}effect of
	interest; default is {cmd:effect(treated)}{p_end}
{synopt:{opth prob(varname)}}declare exogenous selection probability{p_end}
{synopt:{opt m:odel(model [, options])}}control the selection model; default is
	{cmd:model(logit)}{p_end}
{synopt:{opt trim(#)}}set trimming threshold; default is {cmd:trim(0)}{p_end}

{syntab :Bootstrap}
{synopt:{opt reps(#)}}perform # bootstrap replications; default is
	{cmd:reps(1000)}{p_end}
{synopt :{it:{help bootstrap:bootstrap_options}}}bootstrap
	options{p_end}{synoptline}
{p2colreset}{...}
{p 4 6 2}
{cmd:weight}s require option {cmd:force}; see {help weight} and
	{help bootstrap}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:dte} estimates unconditional effects of a binary treatment on features
of the distribution of outcome. Standard errors are computed via the bootstrap.

{pstd}
It offers two methods for applications in which the probability of treatment
take-up differs across observations: inverse probability weighting (Firpo, 2007;
Donald and Hsu, 2014; Firpo and Pinto, 2015) and IV (Frölich and Melly,
2013a,b). Only binary instruments are allowed.

{pstd}
Both IPW and IV require reweighting the sample according to some selection
probability. Write T for the treatment, Z for the instrument and X for
controls. For IPW, this probability is P(T = 1 | X). For IV, it is
P(Z = 1 | X). {cmd:dte} estimates a selection model by default. It is
also possible to declare an exogenous selection probability.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{cmd:statistics(}{it:statname} [{it:...}]{cmd:)}
specifies summary statistics for which treatment effects should be computed. It
is possible to specify multiple statistics, separated by spaces. Available
statistics are:

{marker statname}{...}
{synoptset 17}{...}
{synopt:{space 4}{it:statname}}Definition{p_end}
{space 4}{synoptline}
{synopt:{space 4}{opt mean}}mean{p_end}
{synopt:{space 4}{opt sd}}standard deviation{p_end}
{synopt:{space 4}{opt var}}variance{p_end}
{synopt:{space 4}{opt cv}}coefficient of variation{p_end}
{synopt:{space 4}{opt skewness}}skewness{p_end}
{synopt:{space 4}{opt kurtosis}}kurtosis{p_end}
{synopt:{space 4}{opt gini}}Gini index{p_end}
{synopt:{space 4}{opt iqr}}interquartile range{p_end}
{synopt:{space 4}{opt pr_90_10}}ratio between 90th and 10th percentiles{p_end}
{synopt:{space 4}{opt pr_50_10}}ratio between 50th and 10th percentiles{p_end}
{synopt:{space 4}{opt pr_90_50}}ratio between 90th and 50th percentiles{p_end}
{synopt:{space 4}{opt pr_75_25}}ratio between 75th and 25th percentiles{p_end}
{space 4}{synoptline}
{p2colreset}{...}

{phang}
{opth percentiles(numlist)}
specifies percentiles for which treatment effects should be computed.

{phang}
{opth cdf(numlist)}
specifies outcome values for which treatment effects on the CDF should be
	computed.

{dlgtab:IPW and IV}

{phang}
{cmd:effect(}{it:effecttype}{cmd:)}
specifies the effect of interest. See Imbens and Angrist (1994) and Firpo and
Pinto (2015) for definitions. Available effects are:

{marker effect}{...}
{synoptset 17}{...}
{synopt:{space 4}{it:effecttype}}Definition{p_end}
{space 4}{synoptline}
{synopt:{space 4}{opt compliers}}effect on compliers (IV only){p_end}
{synopt:{space 4}{opt current}}current effect (IPW only){p_end}
{synopt:{space 4}{opt population}}overall effect on the population
	(IPW only){p_end}
{synopt:{space 4}{opt treated}}effect on the treated (IV and IPW);
	default{p_end}
{space 4}{synoptline}
{p2colreset}{...}

{phang}
{space 4}For IV, the effect on the treated pertains to treated compliers, unless
the instrument satisfies one-sided noncompliance (Frölich and Melly, 2013a).

{phang}
{opth prob(varname)}
specifies that the selection probability is exogenous and contained in
{it:varname}.

{phang}
{opt model(model [, options])}
specifies the selection model. The default is {cmd:model(logit)}. Any model is
allowed, as long as the post-estimation command {cmd:predict} is available
afterwards.

{phang}
{opt trim(#)}
controls trimming. Observations are discarded if the appropriate conditional
probability is below # or above 1 - #. For IPW, this probability is the
propensity score, P(T = 1 | X). For IV, it is P(Z = 1 | X). Trimming helps
enforce common support (Crump et al., 2009), but it alters the composition of
the sample. The default is 0 (no trimming).

{dlgtab:Bootstrap}

{phang}
{opt reps(#)}
sets the number of bootstrap replicates. {cmd:reps(0)} suppresses the bootstrap
altogether; only point estimates are reported. Default is {cmd:reps(0)}.

{marker remarks}{...}
{title:Remarks}

{pstd}
Options {cmd:statistics()}, {cmd:percentiles()} and {cmd:cdf()} may be
combined. If neither is specified, {cmd:dte} defaults to {cmd:statistics(mean)}.

{pstd}
The option {cmd:prob()} declares an exogenous selection probability. If the
instrument is treatment assignment, for example, it is a known feature of the
experiment and does not need to be estimated. If the selection probability is
instead estimated, the bootstrap will yield incorrect inference.

{marker examples}{...}
{title:Examples}

{phang}{cmd:. dte y t, stat(mean sd) perc(10(10)90)}{p_end}

{phang}{cmd:. dte y t x1, model(hetprobit, het(x2)) cdf(0)}{p_end}

{phang}{cmd:. dte y (t = z x1##x2)}{p_end}


{marker results}{...}
{title:Stored results}

{phang} See {help bootstrap}.


{marker references}{...}
{title:References}

{phang}
Crump, R.K., V.J. Hotz, G.W. Imbens and O.A. Mitnik (2009): "Dealing with
limited overlap in estimation of average treatment effects",
{it:Biometrika} 96(1), 187–199.

{phang}
Donald, S.G., and Y.C. Hsu (2014): "Estimation and inference for distribution
functions and quantile functions in treatment effect models", {it:Journal of
Econometrics} 178(3), 383-397.

{phang}
Firpo, S., and C. Pinto (2015): "Identification and estimation of distributional
impacts of interventions using changes in inequality measures",
{it:Journal of Applied Econometrics} 31(3), 457–486.

{phang}
Frölich, M., and B. Melly (2013a): "Identification of treatment effects on the
treated with one-sided non-compliance," {it:Econometric Reviews} 32(3), 384–414.

{phang}
Frölich, M., and B. Melly (2013b): "Unconditional quantile treatment effects
under endogeneity", {it:Journal of Business & Economic Statistics} 31(3),
346–357.

{phang}
Imbens, G.W., and J.D. Angrist (1994): "Identification and estimation of local
average treatment effects," {it:Econometrica} 62(2), 467-475.

