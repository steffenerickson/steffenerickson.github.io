{smcl}
{* *! version 1.3.2  15aug2026}{...}
{viewerdialog mvdstudy "dialog _mvdstudy"}{...}
{viewerjumpto "Syntax" "mvdstudy##syntax"}{...}
{viewerjumpto "Description" "mvdstudy##description"}{...}
{viewerjumpto "Options" "mvdstudy##options"}{...}
{viewerjumpto "Bootstrap" "mvdstudy##bootstrap"}{...}
{viewerjumpto "Fixed-facet designs" "mvdstudy##fixedfacet"}{...}
{viewerjumpto "Examples" "mvdstudy##examples"}{...}
{viewerjumpto "Stored results" "mvdstudy##results"}{...}
{viewerjumpto "Reference" "mvdstudy##reference"}{...}
{vieweralsosee "mvgstudy" "help mvgstudy"}{...}

{p2col:{bf:mvdstudy}} Decision study for univariate and multivariate multifaceted designs  {p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:mvdstudy} {cmd:,}
{opt o:bject}{cmd:(}{it:name}{cmd:)}
{opt e:rrortype}{cmd:(}{it:string}{cmd:)}
{cmd:(} {opt facet:num}{cmd:(}{it:numlist}{cmd:)} {cmd:|} {opt cur:rent} {cmd:)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt o:bject(name)}}object of measurement (required){p_end}
{synopt:{opt e:rrortype(string)}}{cmd:relative} or {cmd:absolute} error model (required){p_end}
{synopt:{opt facet:num(numlist)}}number of levels per non-object facet for the D-study sweep{p_end}
{synopt:{opt cur:rent}}use the observed G-study sample sizes (single-point D-study){p_end}
{synopt:{opt comp:ositeweights(matname)}}column vector of weights for composite score reliability{p_end}
{synopt:{opt fix(facetname n ...)}}treat one or more non-object facets as fixed at {it:n} levels{p_end}
{syntab:Bootstrap}
{synopt:{opt boot:strap}}propagate G-study bootstrap distribution through D-study{p_end}
{synopt:{opt ci_level(#)}}confidence level as a percentage; default is {cmd:ci_level(95)}{p_end}
{synoptline}

{p 8 10 2}
{opt facetnum()} and {opt current} are mutually exclusive; exactly one must be specified.


{marker description}{...}
{title:Description}

{pstd}{opt mvdstudy} is a post-estimation command for {helpb mvgstudy}. It uses the variance and covariance component estimates from a G-study to project reliability coefficients (generalizability coefficient Eρ² or phi coefficient Φ) and the associated error and true-score variances under a proposed measurement design.{p_end}

{pstd}Two modes are available:{p_end}

{phang2}
{bf:Current mode} ({opt current}): projects reliability at the sample sizes observed in the G-study data. Returns a single row per outcome variable.{p_end}

{phang2}
{bf:Sweep mode} ({opt facetnum()}): projects reliability over a range of facet-level counts specified by the user. Returns one row per combination of facet counts.{p_end}

{pstd}{bf:Note on unbalanced nested designs.} In a nested design where the number of nested levels varies (e.g., a different number of items per task), {opt current} mode uses {it:max(levels per nesting group)} for the nested facet — the value returned by {cmd:levelsof} after the within-group renorming in {cmd:init_inputs_direct}. For example, with tasks containing 4, 5, and 7 items, {opt current} reports Eρ² as if every task had 7 items. Users wanting reliability at the {it:mean} observed number of items should switch to {opt facetnum()} sweep mode and pass the mean explicitly.{p_end}

{pstd}{cmd:mvdstudy} requires the internal Mata object {cmd:c} created by {cmd:mvgstudy}. It persists across calls, so multiple D-studies (e.g., relative then absolute, or several sweep ranges) can be run without repeating the G-study. After all D-studies are complete, release the object with:{p_end}

{p 8 10 2}
{stata mata drop c}

{pstd}{bf:Design specification note.} The {it:termlist} passed to {cmd:mvgstudy} must explicitly include the interaction or residual term for {cmd:mvdstudy} to compute error variances correctly. For a p×i design, use {cmd:(x1 = p i p#i)}, not {cmd:(x1 = p i)}.{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt o:bject(name)} specifies the object of measurement — the facet whose variance is the numerator of the reliability coefficient. Must match a facet name used in {cmd:mvgstudy}. Required.

{phang}
{opt e:rrortype(string)} specifies the error model. Must be {cmd:relative} or {cmd:absolute}.
Relative errors (interactions involving the object of measurement) yield the generalizability coefficient Eρ².
Absolute errors (all non-object variances) yield the phi coefficient Φ.
See Brennan (2001, Ch. 3) for a full discussion.

{phang}
{opt facet:num(numlist)} specifies the D-study sample sizes for each non-object facet, in the order they appear in the {cmd:mvgstudy} equation. Each value {it:n} generates projected reliability at n_1 = 1, 2, ..., {it:n} levels for that facet. For a design with two non-object facets, provide two values. Mutually exclusive with {opt current}.

{phang}
{opt cur:rent} uses the observed G-study sample sizes as the single D-study design point. Produces one row of output per outcome variable. Mutually exclusive with {opt facetnum()}.

{phang}
{opt comp:ositeweights(matname)} specifies the name of a k×1 Stata matrix of weights for computing composite score reliability across outcome variables. If omitted in a multivariate design, reliability is computed separately for each outcome variable. Not applicable for univariate designs.

{phang}
{opt fix(facetname n ...)} treats one or more non-object facets as fixed at specific numbers of levels, implementing the mixed-design D-study formulas of Brennan (2001, Ch. 4). The option takes alternating facet names and positive integers, e.g., {cmd:fix(i 5)} or {cmd:fix(i 5 h 3)}.

{pmore}
A fixed facet is not randomly sampled from a universe; its variance components augment the remaining effects rather than contributing to measurement error. For example, fixing the item facet at {it:n_i} levels removes item sampling variance from the error formula: σ²_pi/n_i is added to the universe score variance, and σ²_pih/n_i is added to the person×rater interaction. The item variance σ²_i/n_i is discarded as a fixed constant (the same for all persons) and does not appear in either error formula.

{pmore}
{opt facetnum()} then sweeps only the remaining non-fixed, non-object facets. After the D-study completes, the original design is restored so that subsequent {cmd:mvdstudy} calls are unaffected.

{pmore}
{it:Validation:} exits with error code 198 if any named facet is not present in the design, or if the named facet is the object of measurement.

{dlgtab:Bootstrap}

{phang}
{opt bootstrap} propagates the G-study bootstrap distribution (stored by {cmd:mvgstudy, bootstrap}) through the D-study equations to produce bootstrap standard errors and BCa confidence intervals for all D-study quantities (error variance, true-score variance, and the reliability coefficient). Requires that {cmd:mvgstudy} was previously called with the {cmd:bootstrap} option; otherwise exits with error code 198.

{phang}
{opt ci_level(#)} specifies the confidence level for BCa intervals, as a percentage. The default is {cmd:ci_level(95)}. This is independent of the {cmd:ci_level()} specified on {cmd:mvgstudy} — D-study and G-study CIs can differ.


{marker bootstrap}{...}
{title:Bootstrap D-study Details}

{pstd}When {cmd:bootstrap} is specified, {cmd:mvdstudy} replays each stored G-study bootstrap replicate through the D-study projection and computes bias-corrected accelerated (BCa) confidence intervals for three D-study quantities:{p_end}

{p 8 10 2}
1. {bf:var_delta} (relative) or {bf:var_Delta} (absolute) — error variance{p_end}
{p 8 10 2}
2. {bf:var_tau} — true-score variance{p_end}
{p 8 10 2}
3. {bf:erho2} (relative) or {bf:phi} (absolute) — reliability coefficient{p_end}

{pstd}The BCa acceleration is estimated from the jackknife samples stored during the G-study bootstrap. No additional resampling is performed at the D-study stage.{p_end}

{pstd}{bf:Current mode:} the output matrix has 3 rows (one per quantity) and 4 columns ({cmd:estimate}, {cmd:se}, {cmd:ci_lo}, {cmd:ci_hi}).{p_end}

{pstd}{bf:Sweep mode:} for {it:n} projected sample-size combinations, the output matrix has 3×{it:n} rows: the three quantities repeat for each projected design, labeled {cmd:r}{it:i}{cmd:_var_delta}, {cmd:r}{it:i}{cmd:_var_tau}, {cmd:r}{it:i}{cmd:_erho2}, etc.{p_end}


{marker fixedfacet}{...}
{title:Fixed-Facet Designs}

{pstd}When a researcher commits to using a specific set of items (e.g., a fixed test form), those items are not randomly sampled from a universe of items. The {cmd:fix()} option implements the fixed-facet D-study described in Brennan (2001, Ch. 4) by augmenting the remaining variance components and removing the fixed facet from the active design.{p_end}

{pstd}{bf:Augmentation rule.} For each variance component associated with an effect that contains a fixed facet:{p_end}
{p 8 10 2}1. Compute the divisor = product of the fixed {it:n} values for each fixed facet present in that effect.{p_end}
{p 8 10 2}2. Identify the target effect = the original effect with the fixed-facet tokens removed.{p_end}
{p 8 10 2}3. If the target is empty (the effect consists entirely of fixed facets): discard the component — it is a constant offset, not measurement error.{p_end}
{p 8 10 2}4. Otherwise: add the component / divisor to the target effect's component; remove the original effect from the active design.{p_end}

{pstd}{bf:Example.} For a p×i×h design with object p, fixing item at n_i levels:{p_end}

{col 5}{hline 65}
{col 5}Absorbed effect{col 25}Contributes to{col 45}Amount
{col 5}{hline 65}
{col 5}i{col 25}discarded{col 45}fixed constant
{col 5}p#i{col 25}p (universe score var){col 45}σ²_pi / n_i
{col 5}i#h{col 25}h{col 45}σ²_ih / n_i
{col 5}p#i#h{col 25}p#h{col 45}σ²_pih / n_i
{col 5}{hline 65}

{pstd}The active design becomes p×h. The resulting reliability formula is:{p_end}

{p 8 10 2}σ²(τ) = σ²_p + σ²_pi / n_i{p_end}
{p 8 10 2}σ²(δ) = (σ²_ph + σ²_pih / n_i) / n_h{p_end}
{p 8 10 2}Eρ² = σ²(τ) / (σ²(τ) + σ²(δ)){p_end}

{pstd}Fixing items always produces Eρ² ≥ the fully-random Eρ² at the same n_i and n_h, because item sampling variance is removed from the error denominator.{p_end}

{pstd}{bf:State restoration.} After each fixed-facet D-study, {cmd:mvdstudy} automatically restores the original design. Subsequent calls (e.g., to change {it:n_i} or to run a fully-random sweep) always operate on the original components.{p_end}


{marker examples}{...}
{title:Examples}

{pstd}{cmd:mvgstudyexampledata.dta} is a simulated dataset with 50 persons (p), 6 items (i), and 4 raters (h) measured on two outcomes (x1 and x2) in a fully-crossed p×i×h design (1,200 observations). It was generated from a multivariate random effects model using seed 90210; see {cmd:gen_mvgstudyexampledata.do} for the data-generating code and population parameters.{p_end}

{hline}
{pstd}D-study at observed sample sizes (n_p=50, n_i=6, n_h=4){p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) current}{p_end}

{hline}
{pstd}D-study sweep: items 1..6, raters 1..4 (24 design points){p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) facetnum(6 4)}{p_end}
{p 8 10 2}{stata matrix rel = r(x1)}{p_end}

{hline}
{pstd}Composite score D-study (60% x1, 40% x2){p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h)}{p_end}
{p 8 10 2}{stata matrix weights = (0.6 \ 0.4)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) facetnum(6 4) compositeweights(weights)}{p_end}
{p 8 10 2}{stata matrix rel = r(composite)}{p_end}

{hline}
{pstd}Bootstrap D-study: BCa CIs for Eρ² at observed sample sizes{p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h), bootstrap reps(200) seed(42)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) current bootstrap ci_level(95)}{p_end}
{p 8 10 2}{stata matrix boot_d = r(dboot_x1)}{p_end}
{p 8 10 2}{stata matlist boot_d, format(%8.4f)}{p_end}

{hline}
{pstd}Bootstrap D-study: sweep mode with BCa CIs{p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h), bootstrap reps(200) seed(42)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) facetnum(6 4) bootstrap ci_level(95)}{p_end}

{hline}
{pstd}Phi coefficient (absolute error) with bootstrap{p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h), bootstrap reps(200) seed(42)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(absolute) current bootstrap ci_level(95)}{p_end}

{hline}
{pstd}Fixed-facet D-study: fix items at observed count (n_i=6), sweep raters 1..8{p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) facetnum(8) fix(i 6)}{p_end}

{pstd}Returns 8 rows, one per rater count. Item sampling variance is removed from the error formula because items are treated as fixed.{p_end}

{hline}
{pstd}Compare fixed vs. random items at matched sample sizes{p_end}

{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) facetnum(8) fix(i 6)}{p_end}
{p 8 10 2}{stata matrix fixed = r(x1)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) facetnum(6 8)}{p_end}
{p 8 10 2}{stata matrix random = r(x1)}{p_end}
{pstd}Row 4 of {cmd:fixed} (n_h=4) versus row 44 of {cmd:random} (n_i=6, n_h=4) gives Eρ² with fixed versus random item sampling at the same observed sample sizes.{p_end}

{hline}
{pstd}Fixed-facet D-study at observed sample sizes{p_end}

{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) current fix(i 6)}{p_end}


{marker results}{...}
{title:Stored Results}

{pstd}
{cmd:mvdstudy} stores the following in {cmd:r()}:

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Matrices}{p_end}
{synopt:{cmd:r(}{it:varname}{cmd:)}}D-study projection matrix for outcome {it:varname}. Columns: one per non-object facet (sample sizes) plus {cmd:error_var}, {cmd:true_var}, and the reliability coefficient. Rows: one per projected design point.{p_end}
{synopt:{cmd:r(composite)}}D-study projection matrix for the composite score (when {opt compositeweights()} is specified).{p_end}

{pstd}When {cmd:bootstrap} is specified, {cmd:mvdstudy} additionally stores:{p_end}

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Matrices (bootstrap)}{p_end}
{synopt:{cmd:r(dboot_}{it:varname}{cmd:)}}bootstrap summary for outcome {it:varname}: (3×n_rows)×4 matrix with rows labeled {cmd:var_delta}/{cmd:var_Delta}, {cmd:var_tau}, {cmd:erho2}/{cmd:phi} (current mode) or {cmd:r}{it:i}{cmd:_}{it:qty} (sweep mode), and columns {cmd:estimate}, {cmd:se}, {cmd:ci_lo}, {cmd:ci_hi}.{p_end}


{marker reference}{...}
{title:References}

{phang}
Brennan, R. L. (2001). {it:Generalizability theory}. Springer. https://doi.org/10.1007/978-1-4757-3456-0
{p_end}

{phang}
Efron, B., & Tibshirani, R. J. (1993). {it:An introduction to the bootstrap}. CRC Press.
{p_end}

{phang}
Li, G., Michaelides, M. P., & Haertel, E. (2023). Bootstrap confidence intervals for generalizability theory variance components. {it:PLOS ONE}, 18(7), e0288069.
{p_end}
