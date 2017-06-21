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
{bf:bscr} {hline 2} Bootstrap simultaneous confidence regions


{marker syntax}{...}
{title:Syntax}

{p 4}Baseline

{p 8 17 2}
{cmdab:bscr}
[{varlist}]
[{cmd:using} {it:{help filename}}]
{ifin}
[{cmd:,} {it:options}]


{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opth s:tat(vector)}}observed values for each statistic{p_end}
{synopt:{opt unif:orm}}uniform confidence band{p_end}
{synopt:{opt fwer(#)}}control the #-FWER{p_end}
{synopt:{opt fdp(#)}}control the #-FDP{p_end}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt stream:line}}streamlined algorithm for FWER or FDP control{p_end}
{synopt:{opt norec:enter}}do not recenter test statistics; seldom used{p_end}
{synopt:{opt notab:le}}suppress display{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{cmd:weight}s are not allowed.


{marker description}{...}
{title:Description}

{pstd}
{cmd:bscr} constructs simultaneous confidence regions from a bootstrap sample of
statistics. Three choices of confidence set are available. Given the option
{cmd:uniform}, {cmd:bscr} constructs a uniform confidence band, which
accommodates functional hypotheses. The algorithm is due to Chernozhukov,
Fernández-Val and Melly (2013). The options {cmd:fwer()} and {cmd:fdp()} yield
joint confidence intervals, which adjust critical values for multiple hypothesis
testing according to the selected criterion. The k-FWER is the probability of
k false rejections. The k-FDP is the probability of a share k of false
rejections. The intervals are balanced, in the sense that each has the same
asymptotic marginal coverage probability. {cmd:bscr} implements the step-down
algorithm of Romano and Wolf (2010).

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{cmd:stat(vector)} specifies the observed value of each statistic (that is, the
value of the statistic using the original dataset).

{phang}
{opt uniform} displays a uniform confidence band.

{phang}
{opt fwer(#)} displays simultaneous confidence regions that controls the #-FWER.

{phang}
{opt fdp(#)} displays simultaneous confidence regions that controls the #-FDP.

{phang}
{opt level(#)} sets the confidence level.

{phang}
{opt streamline} calls a streamlined algorithm for FWER or FDP control. The
streamlined algorithm is faster, but the resulting confidence intervals are only
asymptotically valid.

{phang}
{opt norecenter} suppresses the recentering of boostrap statistics around their
observed value. This option should only be used when the null hypotheses have
been imposed on the bootstrap sample (e.g., it contains {it:t}-statistics). Instead
of confidence regions, {cmd:bscr} will return critical values.

{phang}
{opt notable} suppresses the display of results.


{marker remarks}{...}
{title:Remarks}

{pstd}
The options {cmd:uniform}, {cmd:fwer()} and {cmd:fdp()} may not be combined. Default
is {cmd:fwer(1)}. Options {cmd:fwer(1)} and {cmd:fdp(0)} are equivalent.

{pstd}
Note that {it:varlist} should cointain variable names as they appear in the
bootstrap sample. 


{marker examples}{...}
{title:Examples}

{phang}{cmd:. bscr using z_sample, fwer(1) norecenter}{p_end}

{phang}{cmd:. bscr using z_sample, fdp(0.1) level(99)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:regress} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N_reps)}}number of complete bootstrap replications{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:bscr}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(crtype)}}type of confidence region{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(cr)}}confidence region{p_end}
{synopt:{cmd:e(cv)}}critical values{p_end}
{synopt:{cmd:e(se)}}standard errors{p_end}

{marker references}{...}
{title:References}

{phang}
Chernozhukov, V., I. Fernández-Val and B. Melly (2013): "Inference on
counterfactual distributions", {it:Econometrica} 81(6), 2205-268.

{phang}
Romano, J.P. and M. Wolf (2013): "Balanced control of generalized error rates",
{it:Annals of Statistics} 38(1), 598-633.

