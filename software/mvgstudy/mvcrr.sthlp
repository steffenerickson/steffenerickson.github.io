{smcl}
{* *! version 1.3.2  15aug2026}{...}
{viewerjumpto "Syntax" "mvcrr##syntax"}{...}
{viewerjumpto "Description" "mvcrr##description"}{...}
{viewerjumpto "Options" "mvcrr##options"}{...}
{viewerjumpto "Remarks" "mvcrr##remarks"}{...}
{viewerjumpto "Examples" "mvcrr##examples"}{...}
{viewerjumpto "Stored results" "mvcrr##results"}{...}
{title:Title}

{phang}
{bf:mvcrr} {hline 2} Projected construct-relevant reliability (CRR) after
{helpb mvgstudy}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:mvcrr} {cmd:,} {opt o:bject(facetname)} {opt nl(#)}
[{opt a:var(varname)} {opt j:var(varname)} {opt p:var(varname)}
{opt nu0(#)} {opt nu1(#)} {opt nutt(#)} {opt sigmae(#)}
{opt n0:var(varname)} {opt n1:var(varname)} {opt rho_lp(#)}
{opt pw:means} {opt b:ootstrap} {opt ci:_level(#)}]

{p 4 6 2}
{cmd:mvcrr} is a post-estimation command.  It requires a prior {cmd:mvgstudy}
run with {it:exactly three} outcome variables and either the one-facet nested
two-effect design {it:(object lesson|object)} — the standard mode — or the
one-effect design {it:(object)} alone, which triggers
{it:single-replication mode} (see the remarks below), for example{p_end}

{p 8 12 2}{cmd:. mvgstudy (Ahat Jhat pihat = p l|p)}{p_end}
{p 8 12 2}{cmd:. mvgstudy (Ahat Jhat pihat = p)}{space 6}{it:(single-replication mode)}{p_end}

{p 4 6 2}
where, for each person-lesson, {cmd:Ahat} is the estimated negative-class
(false-positive) intercept, {cmd:Jhat} the estimated class separation, and
{cmd:pihat} the estimated prevalence, computed upstream from gold-labeled
utterances.  Like {helpb mvdstudy}, {cmd:mvcrr} reads the persistent Mata
object {cmd:c}; run {cmd:mata drop c} when finished.


{marker description}{...}
{title:Description}

{pstd}
{cmd:mvcrr} computes the {it:projected construct-relevant reliability},

{p 8 8 2}
CRR = lambda x Erho2_DCF:PL,

{pstd}
the share of observed person-score variance that is {it:both} reproducible
across lessons {it:and} driven by true prevalence, evaluated at the planned
number of lessons {opt nl(#)}.  It uses the generalizability coefficient that
assumes {it:both} person-level and lesson-level differential classifier
functioning (DCF).  After the universe-score variance cancels between the two
factors, the computed quantity is

{p 8 8 2}
CRR = Jbar^2*s_pi /
[ Jbar^2*s_pi + s_aP + s_jP*(mupi^2 + s_pi)
+ ( s_aL + (Jbar^2 + s_jP)*s_piL + s_jL*(mupi^2 + s_pi + s_piL)
+ sigma2_eps ) / n_L ]

{pstd}
where s_aP, s_jP, s_pi are the person-level (P) variance components of the
intercept, separation, and prevalence (diagonal of the {cmd:emcp} matrix for
the object effect), s_aL, s_jL, s_piL are the corresponding lesson-within-person
(L:P) components, and Abar, Jbar, mupi are the sample means of the three
outcome variables — by default the observation-weighted grand means over all
person-lesson rows of the stored G-study data, or person-weighted means under
{opt pwmeans} (see Options).

{pstd}
Because the intercept Abar cancels from every variance ratio yet biases
absolute prevalence estimates directly, CRR is always reported together with
Abar; treat a low intercept as a companion requirement whenever absolute
prevalence, rather than relative standing, is the reported quantity.

{pstd}
Negative variance-component estimates are truncated at zero (with a note).  If
truncation removes both person-level DCF components, lambda = 1 by construction
and the reported CRR is an {it:upper bound}, not an estimate.


{marker options}{...}
{title:Options}

{phang}
{opt object(facetname)} specifies the object of measurement (e.g., {cmd:p});
it must equal the first effect of the {cmd:mvgstudy} design.  Required.

{phang}
{opt nl(#)} specifies n_L, the planned number of lessons at which CRR is
projected.  Required; must be a positive integer.

{phang}
{opt avar(varname)}, {opt jvar(varname)}, {opt pvar(varname)} identify which
outcome variable is the intercept (A), the separation (J), and the prevalence.
Defaults: the first, second, and third variable of the {cmd:mvgstudy} varlist.

{phang}
{opt sigmae(#)} supplies the mean realization variance sigma2_eps directly.
May not be combined with {opt nu0()}/{opt nu1()}/{opt nutt()}.

{phang}
{opt nu0(#)}, {opt nu1(#)}, {opt nutt(#)} compute sigma2_eps by the closed
form

{p 12 12 2}
sigma2_eps = [ mupi*nu1 + (1-mupi)*nu0 + mupi*(1-mupi)*Jbar^2 ] / nutt

{pmore}
with {opt nu0(#)} and {opt nu1(#)} the pooled within-class score variances
(nu_0^2, nu_1^2) and {opt nutt(#)} the planned utterances per lesson.  All
three must be given together.

{pmore}
{it:Default when neither is specified:} sigma2_eps = 0.  See the remarks below.

{phang}
{opt n0var(varname)}, {opt n1var(varname)} ({it:single-replication mode only;}
{it:required there}) name the per-object counts of gold-negative and
gold-positive utterances, read from the dataset in memory (which must still
hold the estimation sample).  They drive the sampling-error disattenuation of
the person-level components.

{phang}
{opt rho_lp(#)} ({it:single-replication mode only}) sets the
lesson-within-person prevalence variance as a proportion of the person-level
prevalence variance, sigma2_pi_LP = rho_lp x sigma2_pi.  Default is
{cmd:rho_lp(1)}, the dissertation's value.

{phang}
{opt pwmeans} computes Abar, Jbar, and mu_pi as {it:person-weighted} means —
the mean over objects of measurement of the within-object means — instead of
the default observation-weighted grand mean over all person-lesson rows.  The
two are identical when the design is balanced (every object contributes the
same number of lessons).  Under unbalance the default grand mean weights
objects by the number of lessons they contribute; specify {opt pwmeans} to
weight every object equally, matching the convention of validation samples in
which each object contributes one pooled utterance set.  The option propagates
through {opt bootstrap}: per-replicate person-weighted means are stored by
{cmd:mvgstudy}'s bootstrap and used for the CIs.  The means definition in use
is always shown in the output and stored in {cmd:r(pwmeans)}.

{phang}
{opt bootstrap} computes bootstrap BCa confidence intervals for CRR, Abar,
lambda, and Erho2_DCF:PL.  Requires that {cmd:mvgstudy} was run with its
{opt bootstrap} option: {cmd:mvcrr} propagates the stored G-study bootstrap and
jackknife replicates — both the covariance components {it:and} the per-replicate
outcome means (Abar, Jbar, mupi) — through the CRR formula, so no additional
resampling is performed and the option adds essentially no runtime.

{phang}
{opt ci_level(#)} sets the confidence level for the bootstrap intervals;
default is {cmd:ci_level(95)}.


{marker srmode}{...}
{title:Remarks: single-replication mode}

{pstd}
When the G study has the one-effect design {it:(A J prev = object)} — one row
per object of measurement, no lesson replication (e.g., a validation subsample
with all of a teacher's coded utterances pooled into a single set) — the
lesson-within-person components cannot be estimated, and {cmd:mvcrr} switches
to single-replication mode.  It then computes the dissertation's person-only
criterion, CRR = lambda x {it:Erho2_DCF:P}: the across-object variances of
Ahat, Jhat, and the prevalence (from {cmd:emcp1}) are {it:disattenuated} by
subtracting mean per-object sampling variances — nu0^2*mean(1/n0) for the
intercept, nu1^2*mean(1/n1) + nu0^2*mean(1/n0) for the separation, and the
binomial mean(pi(1-pi)/n) for the prevalence — with negatives truncated at
zero; the lesson side enters as D-study parameters
(sigma2_pi_LP = {opt rho_lp(#)} x sigma2_pi; lesson-level DCF components 0)
and sigma2_eps by the {opt nu0()}/{opt nu1()}/{opt nutt()} closed form, which
is {it:required} in this mode ({opt sigmae()} is not allowed, and there is no
double-counting concern because no lesson-level components are estimated).
{opt bootstrap} is not yet supported in this mode; {opt pwmeans} is a no-op
(one row per object) and is ignored with a note.  The coefficient is returned
as {cmd:r(erho2_dcfp)}, and the sampling-variance corrections as
{cmd:r(sampvar)}.

{marker remarks}{...}
{title:Remarks: double-counting of utterance-sampling noise}

{pstd}
When the lesson-level (L:P) components are estimated from per-lesson
Ahat/Jhat/pihat values, those inputs already carry within-lesson
(utterance-sampling) error, so the estimated L:P components {it:absorb} the
noise that the formula's separate sigma2_eps term represents.  Supplying a
nonzero sigma2_eps on top of raw L:P components counts that noise twice.

{pstd}
{cmd:mvcrr} therefore defaults to sigma2_eps = 0: the noise travels inside the
estimated L:P components instead.  Both live in the same /n_L error block, so
the total error budget is approximately preserved; the residual distortion is
slightly conservative (CRR-lowering), which is the safe direction for a
selection criterion.

{pstd}
Users who have {it:disattenuated} their per-lesson estimates upstream should
supply {opt sigmae(#)} or the {opt nu0()}/{opt nu1()}/{opt nutt()} closed form
to restore the textbook decomposition.  Specifying either signals that the
L:P components are noise-free, and a note to that effect is printed.  The
resolved sigma2_eps and its source are always displayed and stored, so every
reported CRR is reproducible from its stated inputs.


{marker examples}{...}
{title:Examples}

{phang}{cmd:. use mvcrrexampledata.dta, clear}{p_end}
{phang}{cmd:. mvgstudy (Ahat Jhat pihat = p l|p)}{p_end}
{phang}{cmd:. mvcrr, object(p) nl(4)}{p_end}

{pstd}Closed-form realization variance for disattenuated inputs
(nu_0^2 = .01, nu_1^2 = .02, 150 utterances per lesson):{p_end}

{phang}{cmd:. mvcrr, object(p) nl(4) nu0(.01) nu1(.02) nutt(150)}{p_end}

{pstd}Bootstrap BCa confidence intervals (G study must be bootstrapped
first):{p_end}

{phang}{cmd:. mvgstudy (Ahat Jhat pihat = p l|p), bootstrap reps(1000) seed(90210)}{p_end}
{phang}{cmd:. mvcrr, object(p) nl(4) bootstrap}{p_end}

{pstd}Unbalanced corpus (unequal lessons per teacher), weighting every teacher
equally in the means:{p_end}

{phang}{cmd:. mvcrr, object(p) nl(4) pwmeans}{p_end}
{phang}{cmd:. mvcrr, object(p) nl(4) pwmeans bootstrap}{p_end}

{pstd}Clean up when finished:{p_end}

{phang}{cmd:. mata drop c}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:mvcrr} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(crr)}}projected construct-relevant reliability (CRR){p_end}
{synopt:{cmd:r(Abar)}}mean false-positive intercept Abar{p_end}
{synopt:{cmd:r(lambda)}}prevalence-driven share of universe-score variance{p_end}
{synopt:{cmd:r(erho2_dcfpl)}}Erho2_DCF:PL, the generalizability coefficient
under person- and lesson-level DCF (standard mode){p_end}
{synopt:{cmd:r(erho2_dcfp)}}Erho2_DCF:P (single-replication mode; replaces
{cmd:r(erho2_dcfpl)}){p_end}
{synopt:{cmd:r(rho_lp)}}(single-replication mode) rho_lp used{p_end}
{synopt:{cmd:r(Jbar)}}mean class separation{p_end}
{synopt:{cmd:r(mupi)}}mean prevalence{p_end}
{synopt:{cmd:r(sigmae)}}resolved sigma2_eps{p_end}
{synopt:{cmd:r(nl)}}planned number of lessons n_L{p_end}
{synopt:{cmd:r(pwmeans)}}1 if person-weighted means were used, 0 otherwise{p_end}
{synopt:{cmd:r(trunc)}}number of components truncated at zero{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(components)}}2 x 3 variance components as used
(rows P, L_P; columns A, J, prevalence){p_end}
{synopt:{cmd:r(crr_table)}}({opt bootstrap} only) 4 x 4 bootstrap summary:
rows crr, Abar, lambda, erho2_dcfpl; columns estimate, se, ci_lo, ci_hi{p_end}
{synopt:{cmd:r(sampvar)}}(single-replication mode) 1 x 3 mean per-object
sampling variances subtracted in the disattenuation{p_end}


{title:Reference}

{pstd}
Erickson, S.  Dissertation, chapter 2: reliability of LLM-classifier-based
observational measures under differential classifier functioning (DCF).
CRR and the DCF-extended generalizability coefficients are defined there.


{title:Also see}

{psee}
{helpb mvgstudy}, {helpb mvdstudy}
