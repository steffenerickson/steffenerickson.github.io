{smcl}
{* *! version 1.3.2  15aug2026}{...}
{viewerdialog mvgstudy "dialog _mvgstudy"}{...}
{viewerjumpto "Syntax" "mvgstudy##syntax"}{...}
{viewerjumpto "Description" "mvgstudy##description"}{...}
{viewerjumpto "Options" "mvgstudy##options"}{...}
{viewerjumpto "Bootstrap options" "mvgstudy##bootstrap"}{...}
{viewerjumpto "Termlist expansion" "mvgstudy##expansion"}{...}
{viewerjumpto "Unbalanced designs" "mvgstudy##unbalanced"}{...}
{viewerjumpto "Known limitations" "mvgstudy##limitations"}{...}
{viewerjumpto "Post estimation" "mvgstudy##postestimation"}{...}
{viewerjumpto "Examples" "mvgstudy##examples"}{...}
{viewerjumpto "Stored results" "mvgstudy##results"}{...}
{viewerjumpto "Reference" "mvgstudy##reference"}{...}
{vieweralsosee "mvdstudy" "help mvdstudy"}{...}
{vieweralsosee "[MV] manova" "help manova"}{...}

{p2col:{bf:mvgstudy}} Variance and covariance component estimation for multifaceted designs  {p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:mvgstudy} {cmd:(}{depvarlist} {cmd:=} {it:termlist}{cmd:)} {ifin}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Bootstrap}
{synopt:{opt boot:strap}}activate bootstrap standard errors and BCa confidence intervals{p_end}
{synopt:{opt reps(#)}}number of bootstrap replications; default is {cmd:reps(1000)}{p_end}
{synopt:{opt ci_level(#)}}confidence level as a percentage; default is {cmd:ci_level(95)}{p_end}
{synopt:{opt seed(#)}}random-number seed for reproducibility{p_end}
{synopt:{opt boot:unit(namelist)}}facet(s) to resample; default is the object of measurement{p_end}
{synoptline}

{p 8 10 2}
where {it:termlist} is a factor variable list of nested and crossed facets:

{p 8 10 2}
o Variables are assumed to be categorical

{p 8 10 2}
o The {cmd:#} symbol indicates crossing

{p 8 10 2}
o The {cmd:|} symbol indicates nesting: the term to the left is nested within the term to the right

{p 8 10 2}
o For fully-crossed or purely-nested designs, only the residual term need be given; {cmd:mvgstudy} expands it to the full termlist automatically. See {help mvgstudy##expansion:Termlist expansion}.
{p_end}


{marker description}{...}
{title:Description}

{pstd}{opt mvgstudy} estimates variance and covariance components for multifaceted generalizability designs, following Brennan (2001). It supports any combination of crossed and nested facets, univariate and multivariate outcome variables, balanced and unbalanced designs, and data with missing outcome values. {p_end}

{pstd}For {bf:balanced designs}, the command uses the SSCP-based method-of-moments estimator derived from a MANOVA decomposition (Brennan, 2001, Ch. 9). {p_end}

{pstd}For {bf:unbalanced designs} (unequal cell sizes), the command uses the exact CP-terms method described in Brennan (2001, Sec. 11.1.3). This is a method-of-moments estimator that constructs a linear system from the observed cross-product matrices of cell means (CP matrices) and a coefficient matrix (K) derived from cell counts. The system is solved exactly — no harmonic-mean approximation is made. The estimates are unbiased and consistent, though they carry higher sampling variance than balanced-design estimates. {p_end}

{pstd}When {depvarlist} contains multiple variables, the command returns covariance component matrices (EMCP matrices) for each effect in {it:termlist}. When {depvarlist} contains one variable, the returned matrices are 1×1 and contain the corresponding variance components. {p_end}

{pstd}{bf:Missing outcome values:} Rows with any missing value in {depvarlist} are automatically excluded before estimation. The reported degrees of freedom and cell counts reflect the filtered data. {p_end}


{marker options}{...}
{title:Options}

{dlgtab:Bootstrap}

{phang}
{opt bootstrap} activates bootstrap standard errors and bias-corrected accelerated (BCa) confidence intervals for all EMCP covariance component matrices. When specified, {cmd:mvgstudy} runs {cmd:reps()} bootstrap replications by cluster-level resampling and stores the distribution for use by {helpb mvdstudy}. Bootstrap is strictly opt-in; all existing output is unchanged when this option is omitted.

{phang}
{opt reps(#)} specifies the number of bootstrap replications. The default is {cmd:reps(1000)}. Values of 200–500 are sufficient for SE estimation; 1000+ recommended for stable CI endpoints.

{phang}
{opt ci_level(#)} specifies the confidence level as a percentage. The default is {cmd:ci_level(95)}. Must be strictly between 1 and 99.

{phang}
{opt seed(#)} sets the random-number seed before drawing bootstrap samples, ensuring reproducibility. If omitted, results will vary across runs.

{phang}
{opt bootunit(namelist)} specifies which facet or facets to use as the resampling unit. Each named facet must appear in {it:termlist}. The default is the first (leftmost) single-facet term in {it:termlist}, which by convention is the object of measurement. Specifying multiple facets (e.g., {cmd:bootunit(p i)}) resamples all listed facets simultaneously via Cartesian product (Li et al., 2023).

{p 8 10 2}
{it:Bootstrap unit options (Li et al., 2023 taxonomy):}

{col 12}{hline 55}
{col 12}Specification{col 32}Resamples{col 50}Fixes
{col 12}{hline 55}
{col 12}{cmd:bootunit(p)}{col 32}persons{col 50}items, residual
{col 12}{cmd:bootunit(i)}{col 32}items{col 50}persons, residual
{col 12}{cmd:bootunit(p i)}{col 32}persons and items{col 50}residual
{col 12}{cmd:bootunit(p i h)}{col 32}all three facets{col 50}—
{col 12}{hline 55}


{marker bootstrap}{...}
{title:Bootstrap Details}

{pstd}{bf:BCa intervals.} The bootstrap uses bias-corrected accelerated (BCa) confidence intervals (Efron & Tibshirani, 1993, Ch. 14). The bias-correction term z0 is estimated from the proportion of bootstrap replicates below the point estimate. The acceleration term a is estimated from a leave-one-cluster-out jackknife at the {cmd:bootunit()} level. BCa intervals are preferred over percentile intervals because they correct for bias and skewness, which are common for near-zero variance component estimates.{p_end}

{pstd}{bf:Cluster-level resampling.} For each bootstrap replicate, unique values of each bootunit facet are sampled with replacement. All rows matching the sampled facet-level values are selected and their facet-level labels are re-indexed sequentially. For balanced designs this preserves cell structure, so the fast SSCP path is used for all replicates. For unbalanced designs the CP-terms path is used; replicates with rank-deficient K matrices are automatically skipped.{p_end}

{pstd}{bf:Nested facets.} When a single bootunit facet is nested within another (e.g., {cmd:bootunit(h)} in a design where items are nested in raters), the command uses within-parent-group resampling: child-facet levels are resampled independently within each combination of parent-facet values. The jackknife for BCa acceleration is performed at the first-parent level.{p_end}

{pstd}{bf:Output.} One {cmd:r(emcp_table_}{it:eff}{cmd:)} matrix per effect, with rows named {it:v_vp} for all k² outcome pairs and columns {cmd:estimate}, {cmd:se}, {cmd:ci_lo}, {cmd:ci_hi}. Effect names containing {cmd:#} or {cmd:|} are sanitized (replaced with {cmd:_}) to satisfy Stata matrix name rules.{p_end}

{pstd}{bf:D-study propagation.} The stored bootstrap distribution is automatically available to {helpb mvdstudy} via the {cmd:bootstrap} option on that command. See {helpb mvdstudy} for details.{p_end}


{marker expansion}{...}
{title:Termlist Expansion}

{pstd}For fully-crossed or purely-nested designs, specifying only the residual term is sufficient. When {it:termlist} contains a single token, {cmd:mvgstudy} automatically derives all higher-level effects and constructs the full termlist before estimation. When a full multi-token termlist is already provided, the expansion code is a no-op.{p_end}

{pstd}{bf:Fully crossed designs} (residual contains {cmd:#} only): all 2^n - 1 non-empty subsets of the n facets are generated.{p_end}

{col 12}{hline 46}
{col 12}Shorthand{col 30}Expands to
{col 12}{hline 46}
{col 12}{cmd:p#i}{col 30}{cmd:p i p#i}
{col 12}{cmd:p#i#h}{col 30}{cmd:p i h p#i p#h i#h p#i#h}
{col 12}{hline 46}

{pstd}{bf:Purely nested designs} (single facet before {cmd:|}): the nesting suffix is expanded recursively, then the nested-facet term is appended.{p_end}

{col 12}{hline 46}
{col 12}Shorthand{col 30}Expands to
{col 12}{hline 46}
{col 12}{cmd:p|i}{col 30}{cmd:i p|i}
{col 12}{cmd:p|i|h}{col 30}{cmd:h i|h p|i|h}
{col 12}{cmd:p|i#h}{col 30}{cmd:i h i#h p|i#h}
{col 12}{hline 46}

{pstd}{bf:Ambiguous mixed designs}: when {cmd:#} appears before the first {cmd:|} (e.g., {cmd:p#i|h}), the shorthand is ambiguous — the same residual string could represent at least two distinct designs. {cmd:mvgstudy} exits with error code 198 and instructs the user to write the full termlist explicitly.{p_end}


{marker unbalanced}{...}
{title:Unbalanced Designs}

{pstd}When cell sizes vary across conditions, {cmd:mvgstudy} automatically detects the imbalance and applies the exact CP-terms estimator from Brennan (2001, Sec. 11.1.3) instead of the balanced SSCP method. No user action is required. {p_end}

{pstd}{bf:Method.} Let α denote a score effect (e.g., person×item interaction) and β range over the effects coarser than α. For each β, the K coefficient K[β,α] equals the sum over β-cells b of the squared sum of cell counts within b for the γ-cells of α, divided by the squared total count in b. The CP matrix for α is the unweighted sum of outer products of the α-cell means. The linear system K × EMCP_stack = CP_stack is solved directly (LU decomposition). {p_end}

{pstd}{bf:Stability.} If K is rank-deficient, the command falls back to the harmonic-mean approximation and issues a warning. In practice, well-posed designs with reasonable cell sizes solve without difficulty. {p_end}

{pstd}{bf:Validated designs.} The following eight designs have been verified against hand-calculated estimates across five random seeds (balanced) and 12% random-deletion scenarios (unbalanced). The last listed component for each design is the residual, denoted (e).{p_end}

{col 5}{hline 74}
{col 5}Design{col 14}Termlist{col 44}Variance components
{col 5}{hline 74}
{col 5}p×i{col 14}{cmd:p i p#i}{col 44}p, i, p#i (e)
{col 5}i:p{col 14}{cmd:p i|p}{col 44}p, i|p (e)
{col 5}p×i×h{col 14}{cmd:p i h p#i p#h i#h p#i#h}{col 44}p, i, h, p#i, p#h, i#h, p#i#h (e)
{col 5}p×(i:h){col 14}{cmd:p h i|h p#h p#i|h}{col 44}p, h, i|h, p#h, p#i|h (e)
{col 5}(i:p)×h{col 14}{cmd:p h i|p p#h h#i|p}{col 44}p, h, i|p, p#h, h#i|p (e)
{col 5}i:(p×h){col 14}{cmd:p h p#h i|p#h}{col 44}p, h, p#h, i|p#h (e)
{col 5}(i×h):p{col 14}{cmd:p i|p h|p h#i|p}{col 44}p, i|p, h|p, h#i|p (e)
{col 5}i:h:p{col 14}{cmd:p h|p i|h|p}{col 44}p, h|p, i|h|p (e)
{col 5}{hline 74}

{pstd}{bf:Validity under random missingness.} Each design was tested with approximately 12% of observations removed at random. Deletion is at the observation level (one row = one person-condition cell), with each row independently flagged for removal with probability 0.12. A minimum of three observations per level of each facet is enforced: any flagged row that would leave a facet level with fewer than three observations is restored. The resulting datasets have approximately 88% of the balanced observations, with unequal cell sizes across all facets. The CP-terms estimator produced non-missing, finite EMCP matrices in all cases across five random seeds. {p_end}

{pstd}{bf:Sensitivity to minimal imbalance.} A separate test removed exactly one observation from a balanced p×i dataset (50 persons × 6 items = 300 observations), leaving one person with five rather than six item scores. This 0.3% reduction in data triggers the CP-terms path and verifies that the estimator gives sensible results even for near-balanced data. Maximum relative deviation from the fully balanced EMCP estimates was less than 5% in all cases. {p_end}


{marker limitations}{...}
{title:Known Limitations}

{pstd}{bf:Facet name prefix collision.} The internal routine that detects which facets belong to an effect uses substring matching. If one facet name is a prefix of another (e.g., facets named {bf:p} and {bf:pr} in the same design), false matches will occur and EMCP estimates will be incorrect. Use facet names that are not prefixes of one another. Single-letter names are safe provided no name is a prefix of another name in the same termlist.{p_end}

{pstd}{bf:Rank-deficient K matrix.} With extreme imbalance or very sparse cell structures the K coefficient matrix may be rank-deficient. The command falls back to the harmonic-mean approximation automatically and issues a warning; those results may be biased relative to the exact CP-terms solution.{p_end}


{marker postestimation}{...}
{title:Post Estimation Commands}

{pstd}The variance and covariance components estimated by {cmd:mvgstudy} can be used in many subsequent analyses. The {helpb mvdstudy} command projects reliability coefficients over varying facet sample sizes (decision study). The {helpb mvcrr} command computes the projected construct-relevant reliability (CRR) for LLM-classifier measurement designs, from a G-study of per-lesson intercept, separation, and prevalence estimates. The internal Mata object {cmd:c} persists across calls so multiple post-estimation commands can be run without repeating the G-study. To release it after all post-estimation analyses are complete, run:{p_end}

{p 8 10 2}
{stata mata drop c}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Commands}{p_end}
{synopt:{helpb mvdstudy}} Decision study for univariate and multivariate multifaceted designs{p_end}
{synopt:{helpb mvcrr}} Projected construct-relevant reliability (CRR) with optional bootstrap BCa CIs{p_end}


{marker results}{...}
{title:Stored Results}

{pstd}
{cmd:mvgstudy} stores the following in {cmd:r()}:

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Matrices}{p_end}
{synopt:{cmd:r(df)}}degrees of freedom for each effect in {it:termlist}{p_end}
{synopt:{cmd:r(P)}}upper triangular matrix of EMCP equation coefficients (SSCP path){p_end}
{synopt:{cmd:r(emcp}{it:i}{cmd:)}}covariance component matrix for the {it:i}th effect in {it:termlist}{p_end}

{pstd}When {cmd:bootstrap} is specified, {cmd:mvgstudy} additionally stores:{p_end}

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Matrices (bootstrap)}{p_end}
{synopt:{cmd:r(emcp_table_}{it:eff}{cmd:)}}bootstrap summary for effect {it:eff}: k²×4 matrix with rows {it:v_vp} for all outcome pairs and columns {cmd:estimate}, {cmd:se}, {cmd:ci_lo}, {cmd:ci_hi}. Effect names with {cmd:#} or {cmd:|} are sanitized to {cmd:_}.{p_end}


{marker examples}{...}
{title:Examples}

{pstd}{cmd:mvgstudyexampledata.dta} is a simulated dataset with 50 persons (p), 6 items (i), and 4 raters (h) measured on two outcomes (x1 and x2) in a fully-crossed p×i×h design (1,200 observations). It was generated from a multivariate random effects model using seed 90210; see {cmd:gen_mvgstudyexampledata.do} in the same directory for the data-generating code and population parameters.{p_end}

{hline}
{pstd}G-study: persons by items by raters, two outcomes{p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h)}{p_end}

{hline}
{pstd}Bootstrap BCa confidence intervals for all components (boot-p, B=200){p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h), bootstrap reps(200) seed(42) ci_level(95)}{p_end}
{p 8 10 2}{stata matrix boot_p = r(emcp_table_p)}{p_end}
{p 8 10 2}{stata matlist boot_p, format(%8.4f)}{p_end}

{hline}
{pstd}Bootstrap with rater resampling (boot-h){p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h), bootstrap reps(200) seed(42) bootunit(h)}{p_end}

{hline}
{pstd}G-study then D-study with bootstrap CIs{p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 x2 = p i h p#i p#h i#h p#i#h), bootstrap reps(200) seed(42)}{p_end}
{p 8 10 2}{stata mvdstudy, object(p) errortype(relative) current bootstrap ci_level(95)}{p_end}

{hline}
{pstd}Univariate G-study (single outcome){p_end}

{p 8 10 2}{stata use mvgstudyexampledata.dta, clear}{p_end}
{p 8 10 2}{stata mvgstudy (x1 = p i h p#i p#h i#h p#i#h)}{p_end}


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
