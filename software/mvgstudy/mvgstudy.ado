
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
*! mvgstudy / mvdstudy / mvcrr  version 1.3.3  15aug2026
*! (Phase 12 adds mvcrr — projected construct-relevant reliability;
*!  Phase 12b adds bootstrap BCa CIs for CRR; Phase 12c adds the pwmeans
*!  option — person-weighted means for unbalanced designs; Phase 12d adds
*!  single-replication mode — lambda x Erho2_DCF:P with disattenuation, for
*!  designs with no lesson replication.
*!  The prior release, v1.2.1, is preserved unmodified as mvgstudy_old.ado.)
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Main Stata Routine
//----------------------------------------------------------------------------//
program mvgstudy, rclass
	version 19
	syntax anything [if] [in] [, Bootstrap REPs(integer 1000) CI_level(integer 95) SEED(string) BOOTunit(string)]

	marksample touse
	tempname facetlevels df_mat P_mat _cp_slvd

	//---//
	// Validate input
	//---//
	if !regexm(`"`anything'"', "\(.*=.*\)") {
		di as error "mvgstudy: invalid syntax -- expected (varlist = termlist)"
		exit 198
	}

	//---//
	// Set up
	// c persists after this program exits so that mvdstudy can access it.
	// Run "mata drop c" manually when all d-studies are complete.
	//---//
	mata c = mvgstudy()
	parse_equation `anything'
	local tempvarlist `r(vars)'
	local effects `r(effects)'
	mata expand_if_residual()
	mata st_local("effects",invtokens(c.sortbylength(tokens(st_local("effects")),0)))
	foreach var of varlist `tempvarlist' {
		local varlist : list varlist | var
	}

	mata facets = uniqrows(tokens(subinstr(subinstr(invtokens(st_local("effects")),"|"," "), "#"," "))')
	mata st_local("facets",invtokens(facets'))
	foreach facet of local facets {
		qui levelsof `facet' if `touse' == 1
		matrix `facetlevels' = (nullmat(`facetlevels') \ `r(r)')
	}
	//---//
	// Direct SSCP computation
	//---//
	mata Y = st_data(., tokens(st_local("varlist")), "`touse'")
	mata Z = st_data(., tokens(st_local("facets")), "`touse'")
	mata c.init_inputs_direct(Y, Z, tokens(st_local("effects")), tokens(st_local("facets")), tokens(st_local("varlist")), st_matrix("`facetlevels'"))
	mata mata drop Y Z facets

	//---//
	// Mata mvgstudy main routine
	//---//
	mata c.mvgstudy_main_routine()
	mata for (loc=asarray_first(c.covcomps); loc!=NULL; loc=asarray_next(c.covcomps, loc)) st_matrix(asarray_key(c.covcomps, loc),asarray_contents(c.covcomps, loc))
	mata st_matrix("`df_mat'",c.df)
	mata st_matrix("`P_mat'",c.P)

	//---//
	// Results
	//---//
	local lengtheffectlist : list sizeof local(effects)
	foreach x of local varlist {
		local name "`x'"
		local names `" `names' "`name'" "'
	}
	local _neg_warn 0
	forvalues i = 1/`lengtheffectlist' {
		matrix rownames emcp`i' = `names'
		matrix colnames emcp`i' = `names'
		matlist emcp`i' , twidth(20) title("`:word `i' of `effects'' Component")
		// Check for negatives before return matrix clears the local copy
		mata: st_numscalar("_min_diag_`i'", min(diagonal(st_matrix("emcp`i'"))))
		if scalar(_min_diag_`i') < 0 local _neg_warn 1
		return matrix emcp`i' = emcp`i'
	}
	if `_neg_warn' {
		di as text "(note: one or more variance component estimates are negative;" ///
		           " this is expected for near-zero components in unbalanced designs;" ///
		           " see Brennan 2001, Ch. 3)"
	}

	return matrix P = `P_mat'
	return matrix df = `df_mat'
	return local varlist `varlist'
	return local effects `effects'
	mata: st_numscalar("`_cp_slvd'", c.cp_solved)
	return scalar cp_solved = `_cp_slvd'

	//---//
	// Phase 10a: Bootstrap CIs (opt-in)
	//---//
	if ("`bootstrap'" != "") {
		// Validate scalar options
		if (`reps' < 1) {
			di as error "mvgstudy: reps() must be >= 1"
			exit 198
		}
		if (`ci_level' <= 0 | `ci_level' >= 100) {
			di as error "mvgstudy: ci_level() must be between 1 and 99"
			exit 198
		}
		if ("`seed'" != "") {
			confirm integer number `seed'
			set seed `seed'
		}

		// Resolve bootunit: default to first single-facet effect (object of measurement)
		if ("`bootunit'" == "") {
			mata st_local("bootunit", tokens(subinstr(subinstr(c.effects[1], "|", " "), "#", " "))[1])
		}
		// Validate: every bootunit facet must be a design facet
		foreach bu_fac of local bootunit {
			if !`:list bu_fac in facets' {
				di as error "mvgstudy: bootunit facet '`bu_fac'' is not in the design facets (`facets')"
				exit 198
			}
		}

		local ci_alpha = (100 - `ci_level') / 100
		di as text _newline "Running bootstrap (B=`reps', boot-`bootunit', `ci_level'% BCa CIs)..."
		mata c.run_bootstrap(`reps', `ci_alpha', tokens("`bootunit'"))
		di as text "Bootstrap complete."

		// Display and return bootstrap tables
		forvalues i = 1/`lengtheffectlist' {
			local effname   : word `i' of `effects'
			// Sanitize effect name for use as Stata matrix name (# and | not allowed)
			local effname_c = subinstr(subinstr("`effname'", "#", "_", .), "|", "_", .)
			tempname _btab`i'
			mata c.push_boot_r_matrix("emcp`i'", "`_btab`i''", `ci_alpha')
			di as text _newline "Bootstrap Summary: `effname' Component (`ci_level'% BCa CI)"
			matlist `_btab`i'', twidth(20) format(%9.4f)
			return matrix emcp_table_`effname_c' = `_btab`i''
		}
	}

end

program mvdstudy, rclass
	version 19
	syntax, Object(string) Errortype(string) ///
	        [FACETnum(numlist integer >=1) CURRent COMPositeweights(string) ///
	         Bootstrap CI_level(real 95) FIX(string)]

	//---//
	// Validate inputs
	//---//
	if ("`errortype'" == "relative") local etype = 0
	else if ("`errortype'" == "absolute") local etype = 1
	else {
		di as error "mvdstudy: errortype must be 'relative' or 'absolute'"
		exit 198
	}

	// facetnum and current are mutually exclusive; at least one required
	if ("`facetnum'" != "" & "`current'" != "") {
		di as error "mvdstudy: facetnum() and current are mutually exclusive"
		exit 198
	}
	if ("`facetnum'" == "" & "`current'" == "") {
		di as error "mvdstudy: specify either facetnum() or current"
		exit 198
	}

	// compositeweights validation
	if ("`compositeweights'" != "") {
		capture confirm matrix `compositeweights'
		if _rc {
			di as error "mvdstudy: compositeweights matrix '`compositeweights'' not found"
			exit 198
		}
		mata st_local("_nvarlist", strofreal(length(c.varlist)))
		if (rowsof(`compositeweights') != `_nvarlist') | (colsof(`compositeweights') != 1) {
			di as error "mvdstudy: compositeweights must be a `_nvarlist'x1 column vector"
			exit 198
		}
	}

	//---//
	// Phase 11: parse fix() and apply fixed-facet augmentation
	//---//
	local n_fixed 0
	if "`fix'" != "" {
		local n_fix_toks : word count `fix'
		if mod(`n_fix_toks', 2) != 0 {
			di as error "mvdstudy: fix() requires alternating facet names and integers, e.g., fix(i 5)"
			exit 198
		}
		local n_fixed = `n_fix_toks' / 2
		tempname _fixed_ns
		matrix `_fixed_ns' = J(1, `n_fixed', .)
		local fixed_names ""
		forvalues fk = 1/`n_fixed' {
			local _fn : word `=2*`fk'-1' of `fix'
			local _fv : word `=2*`fk'' of `fix'
			capture confirm integer number `_fv'
			if _rc | `_fv' <= 0 {
				di as error "mvdstudy: fix() values must be positive integers (got '`_fv'' for facet '`_fn'')"
				exit 198
			}
			local fixed_names "`fixed_names' `_fn'"
			matrix `_fixed_ns'[1, `fk'] = `_fv'
		}
		local fixed_names = strtrim("`fixed_names'")
		mata c.apply_fixed_facets("`object'", tokens("`fixed_names'"), st_matrix("`_fixed_ns'"))
	}

	//---//
	// Build projectionnum vector and observe_only flag
	//---//
	if ("`current'" != "") {
		// Single-point at observed G-study sample sizes
		tempname _fnum
		matrix `_fnum' = J(1, 1, 0)   // dummy; init_dstudyinputs ignores when observe_only=1
		if ("`compositeweights'" == "") {
			mata c.init_dstudyinputs("`object'", `etype', st_matrix("`_fnum'"), 1)
		}
		else {
			mata c.init_dstudyinputs("`object'", `etype', st_matrix("`_fnum'"), 1, ///
			                         st_matrix("`compositeweights'"))
		}
	}
	else {
		// Sweep mode: convert numlist to a Stata matrix row vector
		local nfnum : word count `facetnum'
		tempname _fnum
		matrix `_fnum' = J(1, `nfnum', .)
		forvalues fk = 1/`nfnum' {
			matrix `_fnum'[1, `fk'] = `:word `fk' of `facetnum''
		}
		if ("`compositeweights'" == "") {
			mata c.init_dstudyinputs("`object'", `etype', st_matrix("`_fnum'"), 0)
		}
		else {
			mata c.init_dstudyinputs("`object'", `etype', st_matrix("`_fnum'"), 0, ///
			                         st_matrix("`compositeweights'"))
		}
	}

	//---//
	// D-study routine
	// c persists across calls so mvdstudy can be called multiple times
	// (e.g., relative then absolute) without re-running mvgstudy.
	// Run "mata drop c" manually after all d-studies are complete.
	// Wrap body in capture so fixed-facet state is restored on error.
	//---//
	capture noisily {
		mata c.mvdstudy_main_routine()
		mata c.export_projections()

		//---//
		// Results
		//---//
		foreach v of local vars {
			local names
			foreach x of local colnames {
				local name "`x'"
				local names `" `names' "`name'" "'
			}
			mat colnames `v' = `names'
			matlist `v' , names(c) title("`v'")
			return matrix `v' = `v'
		}
	}
	local _body_rc = _rc
	if `_body_rc' != 0 {
		if `n_fixed' > 0 {
			capture mata c.restore_fixed_facets()
		}
		exit `_body_rc'
	}

	//---//
	// Phase 10d: Bootstrap D-study CIs (opt-in)
	//---//
	if ("`bootstrap'" != "") {
		capture noisily {
			if (`ci_level' <= 0 | `ci_level' >= 100) {
				di as error "mvdstudy: ci_level must be between 1 and 99"
				exit 198
			}
			mata: st_numscalar("_boot_B_ds", c.boot_B)
			if (missing(scalar(_boot_B_ds)) | scalar(_boot_B_ds) <= 0) {
				di as error "mvdstudy: bootstrap requires mvgstudy to have been run with the bootstrap option first"
				exit 198
			}
			if ("`errortype'" == "relative") local coeff_name "erho2"
			else local coeff_name "phi"
			local ci_alpha_ds = (100 - `ci_level') / 100
			di as text _newline "Running D-study bootstrap (`ci_level'% BCa CIs, B=`=scalar(_boot_B_ds)')..."
			mata c.run_dstudy_bootstrap(`ci_alpha_ds')
			di as text "D-study bootstrap complete."
			local _dbidx 0
			foreach v of local vars {
				local ++_dbidx
				tempname _dbtab`_dbidx'
				mata c.push_dstudy_r_matrix("`v'", "`_dbtab`_dbidx''", "`coeff_name'")
				di as text _newline "Bootstrap D-study: `v' (`ci_level'% BCa CI)"
				matlist `_dbtab`_dbidx'', names(c) format(%9.4f)
				return matrix dboot_`v' = `_dbtab`_dbidx''
			}
		}
		local _boot_rc = _rc
		if `_boot_rc' != 0 {
			if `n_fixed' > 0 {
				capture mata c.restore_fixed_facets()
			}
			exit `_boot_rc'
		}
	}

	//---//
	// Phase 11: restore original design after fixed-facet D-study
	//---//
	if `n_fixed' > 0 {
		mata c.restore_fixed_facets()
	}

end

//----------------------------------------------------------------------------//
// Phase 12: mvcrr — Projected Construct-Relevant Reliability
//
// CRR = lambda x Erho2_DCF:PL, the share of observed person-score variance that
// is both reproducible across lessons and driven by true prevalence, under
// person- AND lesson-level differential classifier functioning (DCF).
// Reference: Erickson dissertation ch. 2 (chapter_02_redraft_v5.md), eqs. at
// lines 241 (lambda), 267 (Erho2_DCF:PL), 387 (CRR), 393 (sigma2_eps).
//
// Requires a prior mvgstudy run with exactly three outcome variables
// (A-hat, J-hat, prevalence-hat) and the two-effect design
// (object lesson|object), e.g.  mvgstudy (Ahat Jhat pihat = p l|p).
// Reads the persistent mata object c; does not modify it.
//----------------------------------------------------------------------------//
program mvcrr, rclass
	version 19
	syntax, Object(string) ///
	        [NL(integer 0) Avar(string) Jvar(string) Pvar(string) ///
	         nu0(string) nu1(string) NUTT(string) SIGMAE(string) ///
	         N0var(string) N1var(string) RHO_lp(string) ///
	         PWmeans Bootstrap CI_level(real 95)]

	tempname _nv C_crr C_ab C_lam C_er C_jb C_mp C_se C_tr C_ub CMP

	//---//
	// Validate that mvgstudy results exist and match the required design
	//---//
	capture mata: st_numscalar("`_nv'", length(c.varlist))
	if _rc {
		di as error "mvcrr: no mvgstudy results in memory; run mvgstudy first"
		exit 301
	}
	if scalar(`_nv') != 3 {
		di as error "mvcrr: requires an mvgstudy run with exactly 3 outcome variables (A J prevalence); found `=scalar(`_nv')'"
		exit 198
	}
	mata: st_local("_crr_neff", strofreal(length(c.effects)))
	if `_crr_neff' != 2 & `_crr_neff' != 1 {
		di as error "mvcrr: requires the design (object lesson|object), or (object) alone for single-replication mode; found `_crr_neff' effects"
		exit 198
	}
	// Phase 12d: one-effect design => single-replication (SR) mode
	local srmode = (`_crr_neff' == 1)
	mata: st_local("_crr_eff1", c.effects[1])
	if "`_crr_eff1'" != "`object'" {
		di as error "mvcrr: object(`object') does not match the design's object effect '`_crr_eff1''"
		exit 198
	}
	if !`srmode' {
		mata: st_local("_crr_eff2", c.effects[2])
		if !regexm("`_crr_eff2'", "^[^|#]+\|`object'$") {
			di as error "mvcrr: second effect '`_crr_eff2'' is not of the form lesson|`object'"
			exit 198
		}
	}
	if `nl' < 1 {
		di as error "mvcrr: nl() is required and must be a positive integer"
		exit 198
	}

	//---//
	// Resolve variable roles (default order: A, J, prevalence)
	//---//
	mata: st_local("_crr_vl", invtokens(c.varlist))
	if ("`avar'" == "") local avar : word 1 of `_crr_vl'
	if ("`jvar'" == "") local jvar : word 2 of `_crr_vl'
	if ("`pvar'" == "") local pvar : word 3 of `_crr_vl'
	foreach rv in avar jvar pvar {
		if !`:list `rv' in _crr_vl' {
			di as error "mvcrr: variable '``rv''' is not among the mvgstudy outcome variables (`_crr_vl')"
			exit 198
		}
	}
	if ("`avar'" == "`jvar'" | "`avar'" == "`pvar'" | "`jvar'" == "`pvar'") {
		di as error "mvcrr: avar(), jvar(), and pvar() must name three distinct variables"
		exit 198
	}

	//---//
	// Resolve sigma2_eps (realization variance)
	// Default sigmae = 0: when the L:P components are estimated from per-lesson
	// data they already absorb utterance-sampling noise, so adding a separate
	// sigma2_eps would count that noise twice (see build_mvcrr_plan.md, sec. 5).
	//---//
	local sig_mode 0
	local s1 0
	local s2 0
	local s3 0
	if ("`sigmae'" != "" & ("`nu0'" != "" | "`nu1'" != "" | "`nutt'" != "")) {
		di as error "mvcrr: sigmae() may not be combined with nu0()/nu1()/nutt()"
		exit 198
	}
	if ("`nu0'" != "" | "`nu1'" != "" | "`nutt'" != "") {
		if ("`nu0'" == "" | "`nu1'" == "" | "`nutt'" == "") {
			di as error "mvcrr: nu0(), nu1(), and nutt() must be specified together"
			exit 198
		}
		foreach o in nu0 nu1 nutt {
			capture confirm number ``o''
			if _rc {
				di as error "mvcrr: `o'() must be a number"
				exit 198
			}
		}
		if `nu0' < 0 | `nu1' < 0 {
			di as error "mvcrr: nu0() and nu1() must be nonnegative"
			exit 198
		}
		if `nutt' <= 0 {
			di as error "mvcrr: nutt() must be positive"
			exit 198
		}
		local sig_mode 2
		local s1 `nu0'
		local s2 `nu1'
		local s3 `nutt'
	}
	else if ("`sigmae'" != "") {
		capture confirm number `sigmae'
		if _rc {
			di as error "mvcrr: sigmae() must be a nonnegative number"
			exit 198
		}
		if `sigmae' < 0 {
			di as error "mvcrr: sigmae() must be a nonnegative number"
			exit 198
		}
		local sig_mode 1
		local s1 `sigmae'
	}

	//---//
	// Phase 12d: single-replication mode option validation.
	// SR mode reproduces the paper's person-only CRR (lambda x Erho2_DCF:P):
	// person components disattenuated for sampling error; lesson side entered
	// as D-study parameters (rho_lp, closed-form sigma2_eps).
	//---//
	if `srmode' {
		if ("`bootstrap'" != "") {
			di as error "mvcrr: bootstrap is not supported in single-replication mode"
			exit 198
		}
		if ("`sigmae'" != "") {
			di as error "mvcrr: single-replication mode requires nu0()/nu1()/nutt(); sigmae() is not allowed"
			exit 198
		}
		if `sig_mode' != 2 {
			di as error "mvcrr: single-replication mode requires nu0(), nu1(), and nutt()"
			exit 198
		}
		if ("`n0var'" == "" | "`n1var'" == "") {
			di as error "mvcrr: single-replication mode requires n0var() and n1var() (per-object gold-negative/-positive counts)"
			exit 198
		}
		confirm numeric variable `n0var'
		confirm numeric variable `n1var'
		if ("`rho_lp'" == "") local rho_lp 1
		capture confirm number `rho_lp'
		if _rc {
			di as error "mvcrr: rho_lp() must be a nonnegative number"
			exit 198
		}
		if `rho_lp' < 0 {
			di as error "mvcrr: rho_lp() must be a nonnegative number"
			exit 198
		}
		if ("`pwmeans'" != "") {
			di as text "(note: pwmeans ignored in single-replication mode — one row per object)"
		}
	}
	else if ("`n0var'`n1var'`rho_lp'" != "") {
		di as error "mvcrr: n0var()/n1var()/rho_lp() apply only to single-replication mode (one-effect design)"
		exit 198
	}

	//---//
	// Compute
	// pw = 1: means of within-object means (person-weighted); default is the
	// observation-weighted grand mean.  Identical under balanced designs.
	//---//
	local pw = cond("`pwmeans'" == "" | `srmode', 0, 1)
	if `srmode' {
		tempname SVM
		mata c.compute_crr_sr("`avar'", "`jvar'", "`pvar'", `nl', ///
		                      `nu0', `nu1', `nutt', `rho_lp', ///
		                      "`n0var'", "`n1var'", "`CMP'", "`SVM'", ///
		                      "`C_crr' `C_ab' `C_lam' `C_er' `C_jb' `C_mp' `C_se' `C_tr' `C_ub'")
	}
	else {
		mata c.compute_crr("`avar'", "`jvar'", "`pvar'", `nl', `sig_mode', ///
		                   `s1', `s2', `s3', `pw', "`CMP'", ///
		                   "`C_crr' `C_ab' `C_lam' `C_er' `C_jb' `C_mp' `C_se' `C_tr' `C_ub'")
	}

	//---//
	// Display: CRR reported together with A-bar (companion requirement)
	//---//
	local sig_src "0 (default; utterance-sampling noise absorbed in L:P components)"
	if `sig_mode' == 1 local sig_src "user-supplied via sigmae()"
	if `sig_mode' == 2 local sig_src "closed form from nu0(), nu1(), nutt()"
	local mean_src "observation-weighted grand means"
	if `pw' local mean_src "person-weighted means (pwmeans; mean of within-`object' means)"

	if `srmode' {
		local coeflab "Erho2_DCF:P"
		local modelab "single-replication mode: person-level DCF, disattenuated"
	}
	else {
		local coeflab "Erho2_DCF:PL"
		local modelab "person- and lesson-level DCF"
	}
	di as text _newline "Projected Construct-Relevant Reliability (`modelab')"
	di as text "{hline 72}"
	di as text "  CRR   (lambda x `coeflab')" _column(39) "= " as result %9.4f scalar(`C_crr')
	di as text "  A-bar (mean false-pos. intercept)   = " as result %9.4f scalar(`C_ab')
	di as text "{hline 72}"
	di as text "  lambda                              = " as result %9.4f scalar(`C_lam')
	di as text "  `coeflab'" _column(39) "= " as result %9.4f scalar(`C_er')
	di as text "  J-bar (mean class separation)       = " as result %9.4f scalar(`C_jb')
	di as text "  mu_pi (mean prevalence)             = " as result %9.4f scalar(`C_mp')
	di as text "  sigma2_eps (realization variance)   = " as result %9.6f scalar(`C_se')
	di as text "        source: `sig_src'"
	di as text "  n_L (planned lessons)               = " as result %9.0g `nl'
	if `srmode' {
		di as text "  rho_LP (sigma2_pi_LP / sigma2_pi)   = " as result %9.4f `rho_lp'
	}
	di as text "  A-bar/J-bar/mu_pi from: `mean_src'"
	di as text ""
	matlist `CMP', twidth(8) format(%12.6f) ///
		title("Variance components as used (post-truncation)")
	if `srmode' {
		di as text ""
		matlist `SVM', twidth(8) format(%12.6f) ///
			title("Mean per-object sampling variances subtracted (disattenuation)")
	}

	if scalar(`C_tr') > 0 {
		di as text "(note: `=scalar(`C_tr')' negative variance component(s) truncated at zero)"
	}
	if scalar(`C_ub') == 1 {
		di as text "(note: person-level DCF components truncated to zero; lambda = 1 by" ///
		           " construction and CRR is an upper bound, not an estimate)"
	}
	if !`srmode' & `sig_mode' > 0 & scalar(`C_se') > 0 {
		di as text "(note: sigmae > 0 assumes the lesson-level (L:P) components are" ///
		           " disattenuated for utterance-sampling error; otherwise this noise" ///
		           " is counted twice)"
	}

	//---//
	// Phase 12b: bootstrap BCa CIs for CRR (opt-in)
	// Propagates the G-study bootstrap reps (components AND per-rep means)
	// through the CRR formula; requires mvgstudy , bootstrap first.
	//---//
	if ("`bootstrap'" != "") {
		if (`ci_level' <= 0 | `ci_level' >= 100) {
			di as error "mvcrr: ci_level() must be between 1 and 99"
			exit 198
		}
		tempname _bB CRT
		capture mata: st_numscalar("`_bB'", c.boot_B)
		if _rc | missing(scalar(`_bB')) | scalar(`_bB') <= 0 {
			di as error "mvcrr: bootstrap requires mvgstudy to have been run with the bootstrap option first"
			exit 198
		}
		local ci_alpha = (100 - `ci_level') / 100
		di as text _newline "Running CRR bootstrap (`ci_level'% BCa CIs, B=`=scalar(`_bB')')..."
		mata c.run_crr_bootstrap("`avar'", "`jvar'", "`pvar'", `nl', `sig_mode', ///
		                         `s1', `s2', `s3', `pw', `ci_alpha', "`CRT'")
		di as text "CRR bootstrap complete."
		di as text ""
		matlist `CRT', format(%9.4f) ///
			title("Bootstrap Summary: CRR (`ci_level'% BCa CI)")
		return matrix crr_table = `CRT'
	}

	//---//
	// Returned results
	//---//
	return scalar crr         = scalar(`C_crr')
	return scalar Abar        = scalar(`C_ab')
	return scalar lambda      = scalar(`C_lam')
	if `srmode' {
		return scalar erho2_dcfp = scalar(`C_er')
		return scalar rho_lp     = `rho_lp'
		return matrix sampvar    = `SVM'
	}
	else {
		return scalar erho2_dcfpl = scalar(`C_er')
	}
	return scalar Jbar        = scalar(`C_jb')
	return scalar mupi        = scalar(`C_mp')
	return scalar sigmae      = scalar(`C_se')
	return scalar nl          = `nl'
	return scalar pwmeans     = `pw'
	return scalar trunc       = scalar(`C_tr')
	return matrix components  = `CMP'

end

//----------------------------------------------------------------------------//
// Stata Sub Routines
//----------------------------------------------------------------------------//
program parse_equation, rclass
	version 19
	syntax anything

	mata: parse_equation_mata(st_local("anything"))
	return local vars `vars'
	return local effects `effects'
end

//----------------------------------------------------------------------------//
// Mata Code
//----------------------------------------------------------------------------//
mata

// Standalone functions

void parse_equation_mata(string scalar equation)
{
	string rowvector	toks
	real rowvector		eq_pos
	real scalar			i

	toks   = tokens(subinstr(subinstr(equation,"(",""),")",""))
	eq_pos = strpos(toks,"=")
	i = 1
	while (eq_pos[i] < 1) ++i
	st_local("vars",    invtokens(toks[1, 1..i-1]))
	st_local("effects", invtokens(toks[1, i+1..cols(toks)]))
}

void expand_if_residual()
{
	string scalar		effects_str, residual
	string rowvector	pipe_toks, canon_facets, result

	effects_str = st_local("effects")
	if (length(tokens(effects_str)) != 1) return

	residual = tokens(effects_str)[1]

	if (strpos(residual, "|") > 0) {
		pipe_toks = tokens(residual, "|")
		if (strpos(pipe_toks[1], "#") > 0) {
			errprintf("mvgstudy: ambiguous termlist '%s'\n", residual)
			errprintf("  Multiple facets appear before '|'; this could be a\n")
			errprintf("  fully-nested design or a partially-crossed design.\n")
			errprintf("  Specify the full termlist explicitly.\n")
			exit(198)
		}
	}

	canon_facets = sort(tokens(subinstr(subinstr(residual,"|"," "),"#"," "))', 1)'
	result = _expand_residual(residual, canon_facets)
	st_local("effects", invtokens(result))
}

string rowvector _expand_residual(string scalar residual, string vector canon)
{
	string rowvector	pipe_toks, facs, nesting_parts, result, beta
	string scalar		nesting_suffix
	real scalar			n, b, j, k

	if (strpos(residual, "|") == 0) {
		facs   = tokens(subinstr(residual, "#", " "))
		n      = length(facs)
		result = J(1, 0, "")
		for (b=1; b<2^n; b++) {
			beta = J(1, 0, "")
			for (j=1; j<=n; j++) {
				if (mod(floor(b/(2^(j-1))),2)==1) beta = (beta, facs[j])
			}
			result = (result, _canon_join(beta, canon))
		}
		return(result)
	}
	else {
		pipe_toks     = tokens(residual, "|")
		nesting_parts = J(1, 0, "")
		for (k=3; k<=length(pipe_toks); k=k+2) nesting_parts = (nesting_parts, pipe_toks[k])
		nesting_suffix = invtokens(nesting_parts, "|")
		result = (_expand_residual(nesting_suffix, canon), pipe_toks[1] + "|" + nesting_suffix)
		return(result)
	}
}

string scalar _canon_join(string rowvector subset, string vector canon)
{
	string rowvector	out
	real scalar			fi, j, found

	out = J(1, 0, "")
	for (fi=1; fi<=length(canon); fi++) {
		found = 0
		for (j=1; j<=length(subset); j++) {
			if (canon[fi]==subset[j]) found = 1
		}
		if (found) out = (out, canon[fi])
	}
	return(invtokens(out, "#"))
}

// True iff facet 'fac' appears as a whole token in effect string 'eff'.
// Splits eff on "#" and "|" so it works for crossed/nested effect forms.
// Why: strpos(eff, fac) gives false positives when one facet name is a
// non-delimited substring of another (e.g. facets p and pr; "p" is found
// in "pr#i" by strpos but is not actually a facet of that effect).
real scalar _facet_in_effect(string scalar eff, string scalar fac)
{
	string rowvector toks
	real scalar      t
	toks = tokens(subinstr(subinstr(eff, "|", " "), "#", " "))
	for (t = 1; t <= length(toks); t++) {
		if (toks[t] == fac) return(1)
	}
	return(0)
}

// True iff every facet in eff_A is also a facet in eff_B (subset of facets).
// Splits both on "#" and "|".  Used by get_p() to decide which P[row,col] cells
// to retain.  Replaces the old concatenated-string strpos check, which mis-
// handled cases like p,r subset of p,t,r ("pr" is not a substring of "ptr").
real scalar _effect_subset(string scalar eff_A, string scalar eff_B)
{
	string rowvector toks_A, toks_B
	real scalar      a, b, found
	toks_A = tokens(subinstr(subinstr(eff_A, "|", " "), "#", " "))
	toks_B = tokens(subinstr(subinstr(eff_B, "|", " "), "#", " "))
	for (a = 1; a <= length(toks_A); a++) {
		found = 0
		for (b = 1; b <= length(toks_B); b++) {
			if (toks_A[a] == toks_B[b]) found = 1
		}
		if (!found) return(0)
	}
	return(1)
}

// Class Definition

class mvgstudy
{

	//----------//
	public:
	//---------//

	// Gstudy routines
	void init_inputs()
	void init_inputs_direct()
	void mvgstudy_main_routine()

	// D Study routines
	void init_dstudyinputs()
	void mvdstudy_main_routine()

	// Push projection matrices to Stata and set colnames/vars locals for display
	void export_projections()

	// Phase 10: Bootstrap CIs (boot-p/i/pi/etc.; Li 2023)
	void run_bootstrap()
	void push_boot_r_matrix()
	string rowvector _nesting_parents()   // returns parent facets if bootunit is nested

	// Phase 10d: D-study bootstrap propagation
	void run_dstudy_bootstrap()
	void push_dstudy_r_matrix()

	// Phase 11: fixed-facet D-study augmentation
	void apply_fixed_facets()
	void restore_fixed_facets()

	// Phase 12: projected construct-relevant reliability (CRR)
	void compute_crr()
	void run_crr_bootstrap()
	void compute_crr_sr()   // Phase 12d: single-replication mode

	// Mata functions used in stata wrapper program
	string matrix sortbylength()
	real vector flip_select()

	// Gstudy objects
	string vector   effects
	string vector 	facets
	string vector 	varlist
	real vector		facetlevels
	real vector 	flproducts
	real vector 	df
	real matrix     P
	real matrix     compsstacked
	transmorphic  	sscplist
	transmorphic  	mscplist
	transmorphic  	covcomps

	// Phase 8: CP terms fields
	real matrix     K_mat
	transmorphic    cpcache
	real scalar     cp_solved

	// Data (stored for direct SSCP computation and Phase 4 residuals)
	real matrix     Y_data
	real matrix     Z_data

	// Phase 10: Bootstrap fields
	real scalar     boot_B
	real scalar     boot_ci_alpha
	real scalar     boot_cluster_col
	transmorphic    boot_store
	transmorphic    jack_store
	real scalar     quiet_          // 1 = suppress per-rep messages (set in bootstrap loop)

	// Phase 12b: per-rep outcome means (B x k and n_jack x k), stored by
	// run_bootstrap so run_crr_bootstrap can propagate Abar/Jbar/mupi
	real matrix     boot_means
	real matrix     jack_means

	// Phase 12c: per-rep person-weighted outcome means (mean of within-object
	// means, grouped by the first facet of effects[1]); differ from
	// boot_means/jack_means only under unbalance
	real matrix     boot_pwmeans
	real matrix     jack_pwmeans

	// Phase 10d: D-study bootstrap summary
	transmorphic    dstudy_boot_summary   // asarray: varname → (n_rows*3)×4 matrix

	// Phase 11: saved state for fixed-facet restoration
	string vector   orig_effects
	transmorphic    orig_covcomps
	string vector   orig_facets
	real vector     orig_facetlevels

	// Dstudy objects
	string scalar objmeasurement
	real scalar errortype
	real vector projectionnum
	real scalar observe_only
	real colvector weights
	string matrix projectionmat
	string rowvector uniqueffects
	transmorphic  	projections

	//----------//
	private:
	//----------//

	// Gstudy internal routines
	void get_facetlevelproducts()
	void get_p()
	void get_mscplist()
	void get_compsstacked()
	void emcpmatrixprocedure()

	// Direct SSCP engine (Phase 1+3+4)
	void             get_df_balanced()
	void             get_df_unbalanced()
	void             get_sscp_all()
	real matrix      get_sscp_for_effect()
	real matrix      get_sscp_nested()
	void             _cache_cellmeans()
	void             get_eff_facetlevels()
	real scalar      is_balanced()
	string scalar    _canonical_key()
	transmorphic     _build_idx_map()
	transmorphic     cmcache
	transmorphic     szcache
	real scalar      is_bal
	real matrix      eff_fl_mat

	// Phase 8: CP terms engine
	void             filter_missing_values()
	void             get_cp_all()
	void             get_cp_coefficients()
	void             solve_cp_system()
	string scalar    _flat_key()
	string rowvector _flat_facets()
	string rowvector _setdiff()
	real rowvector   _find_facet_cols_in_key()
	real scalar      _kcoef_mu_alpha()
	real scalar      _kcoef_beta_alpha()
	real scalar      _ncells()
	real matrix      _get_means_only()

	real scalar relativevariance()
	real scalar profilevar()
	real scalar errorprofilevar()

	// Phase 10: BCa helpers
	real scalar    bca_acceleration()
	real rowvector bca_ci()

	// Phase 12: CRR internals
	real rowvector _crr_resolve_idx()
	real rowvector _crr_core()
	real rowvector _pw_means()
	real scalar    _crr_obj_col()

	// Dstudy internal routines
	void absolute_errors_effects()
	void relative_errors_effects()
	real matrix errorprojections()
	void combine_matrices()
	real matrix rowwise_nonzero_product()
	real matrix compositevariance()
	real vector extract_bracket_numbers()

	real vector errorvariances
	real scalar  objmeasurevariance

}

// Public Routines
void mvgstudy::init_inputs(transmorphic sscplist,
                           string vector effects,
                           string vector facets,
                           string vector varlist,
                           real vector facetlevels,
                           real df)
{
	this.sscplist    = sscplist
	this.effects     = effects
	this.facets      = facets
	this.varlist     = varlist
	this.facetlevels = facetlevels
	this.df          = df
}

void mvgstudy::mvgstudy_main_routine()
{
	get_facetlevelproducts()
	get_p()
	if (cp_solved) return
	get_mscplist()
	get_compsstacked()
	emcpmatrixprocedure()
}

void mvgstudy::init_dstudyinputs(string scalar objmeasurement,
                                 real scalar errortype,
                                 real vector projectionnum_vec,
                                 real scalar observe_only_flag,
                               | real vector weights_in)
{
	real scalar  fi
	real vector  design_nonobj_n

	this.objmeasurement = objmeasurement
	this.errortype      = errortype
	this.observe_only   = observe_only_flag
	if (args()>4) this.weights = weights_in
	else          this.weights = J(0,0,.)

	// Build observed-n vector for non-object facets (canonical order)
	design_nonobj_n = J(1, 0, .)
	for (fi=1; fi<=length(facets); fi++) {
		if (facets[fi] != objmeasurement) design_nonobj_n = (design_nonobj_n, facetlevels[fi])
	}

	if (observe_only_flag) {
		// Use observed G-study sample sizes directly
		this.projectionnum = design_nonobj_n
	}
	else {
		if (length(projectionnum_vec) == 1) {
			// Replicate single value to all non-object facets
			this.projectionnum = J(1, length(design_nonobj_n), projectionnum_vec[1,1])
		}
		else if (length(projectionnum_vec) == length(design_nonobj_n)) {
			this.projectionnum = projectionnum_vec
		}
		else {
			errprintf("mvdstudy: facetnum expects %g value(s) for %g non-object facet(s)\n",
			          length(design_nonobj_n), length(design_nonobj_n))
			exit(198)
		}
	}
}

void mvgstudy::mvdstudy_main_routine()
{
	real scalar		i, j
	string scalar	name
	real vector 	variances
	real matrix 	result

	projections = asarray_create()
	if (weights == J(0,0,.)) {
		for (i=1;i<=length(varlist);i++) {
			variances = J(length(effects),1,.)
			for (j=1;j<=length(effects);j++) {
				name = "emcp" + strofreal(j)
				variances[j] = asarray(covcomps, name)[i,i]
			}
			if (errortype == 0) relative_errors_effects(variances)
			else if (errortype == 1) absolute_errors_effects(variances)
			result = errorprojections()
			asarray(projections,varlist[i],result)
		}
	}
	else {
		variances = compositevariance()
		if (errortype == 0) relative_errors_effects(variances)
		else if (errortype == 1) absolute_errors_effects(variances)
		result = errorprojections()
		asarray(projections,"composite",result)
	}
}

void mvgstudy::export_projections()
{
	string scalar	stringlist, colnames
	transmorphic	loc

	stringlist = ""
	for (loc=asarray_first(projections); loc!=NULL; loc=asarray_next(projections, loc)) {
		stringlist = stringlist + " " + asarray_key(projections, loc)
		st_matrix(asarray_key(projections, loc), asarray_contents(projections, loc))
	}
	colnames = invtokens((uniqueffects , "Error" , "True" , "Rel."))
	st_local("vars",     stringlist)
	st_local("colnames", colnames)
}

string matrix mvgstudy::sortbylength(string vector M, real scalar descending)
{
	real matrix		idx

	idx = sort((strlen(M)', (1::length(M))), descending ? (-1, 2) : (1, 2))
	return(M[idx[,2]'])
}

real vector mvgstudy::flip_select(real vector mask) return(mm_cond(mask:>0,mask:-mask,mask:+1))


// Private Routines
// (1) Covariance component calculation

void mvgstudy::get_facetlevelproducts()
{
	real scalar 		flp, x, i, j

	flproducts = J(length(effects),1,.)
	for(i=1;i<=length(effects);i++) {
		flp = 1
		for(j=1;j<=length(facets);j++) {
			if (!_facet_in_effect(effects[i], facets[j])) {
				x = is_bal ? facetlevels[j] : eff_fl_mat[i,j]
			}
			else x = 1
			flp = flp * x
		}
		flproducts[i] = flp
	}
}

void mvgstudy::get_p()
{
	real scalar 		row, col, n_eff

	// Work directly on this.effects (with "#" and "|" preserved); use the
	// token-set helper _effect_subset() to test α ⊆ β.  The earlier
	// concatenated-stripped form ("p#r" → "pr") collided with standalone
	// facet names (e.g. a facet literally called "pr") and also failed the
	// {p,r} ⊆ {p,t,r} test because "pr" is not a consecutive substring of "ptr".
	n_eff = length(effects)
	P     = J(n_eff, n_eff, 0)
	for (col = 1; col <= n_eff; col++) {
		row = col
		while (row >= 1) {
			P[row, col] = flproducts[col]
			--row
		}
	}
	// Residual (last effect) contains every facet, so every row is trivially
	// a subset of it; skipping col = n_eff leaves its column as all-ones,
	// which matches the closed-form EMS coefficient for the residual.
	for (row = 1; row <= n_eff; row++) {
		for (col = 1; col <= n_eff - 1; col++) {
			if (!_effect_subset(effects[row], effects[col])) P[row, col] = 0
		}
	}
}

void mvgstudy::get_mscplist()
{
	real scalar		i
	real matrix 	temp

	mscplist = asarray_create()
	for(i=1;i<=length(effects);i++) {
		temp = asarray(sscplist, strofreal(i)):/df[i]
		asarray(mscplist, strofreal(i), temp)
	}
}

void mvgstudy::get_compsstacked()
{
	real scalar 	i, ncols

	ncols = cols(asarray(mscplist,"1"))
	compsstacked = J(0,ncols,.)
	for(i=1;i<=length(effects);i++) {
		compsstacked = compsstacked \ asarray(mscplist, strofreal(i))
	}
}

void mvgstudy::emcpmatrixprocedure()
{
	string scalar	name
	real scalar 	skip, i, j
	real matrix 	rangemat, starting, ending, inv_p, temp, res

	covcomps = asarray_create()
	skip     = cols(compsstacked)
	rangemat = range(0, rows(compsstacked), skip)
	starting = rangemat[1..(rows(rangemat)-1)]:+1
	ending   = rangemat[2..rows(rangemat)]
	inv_p    = luinv(P)
	for(i=1;i<=cols(P);i++){
		res = J(cols(compsstacked),cols(compsstacked),0)
		for(j=1;j<=cols(P);j++) {
			temp = compsstacked[starting[j]..ending[j], 1..cols(compsstacked)]
			res  = res + (temp * inv_p[i,j])
		}
		name = "emcp" + strofreal(i)
		asarray(covcomps, name, res)
	}
}


// (2) Direct SSCP Engine

void mvgstudy::init_inputs_direct(real matrix Y,
                                   real matrix Z,
                                   string vector effects_in,
                                   string vector facets_in,
                                   string vector varlist_in,
                                   real vector   facetlevels_in)
{
	real matrix     Z_norm
	real scalar     jj_n, kk_n
	real colvector  uniq_n, sel
	real scalar     warn_exp, warn_obs, warn_pct, warn_small_n

	Y_data = Y
	// Normalize each facet column to sequential 1-indexed ranks so the
	// packed-integer row formula works even when variable values have gaps.
	Z_norm = Z
	for (jj_n=1; jj_n<=cols(Z); jj_n++) {
		uniq_n = uniqrows(Z[.,jj_n])
		for (kk_n=1; kk_n<=rows(uniq_n); kk_n++) {
			sel = selectindex(Z[.,jj_n] :== uniq_n[kk_n])
			Z_norm[sel, jj_n] = J(rows(sel), 1, kk_n)
		}
	}
	Z_data           = Z_norm
	this.effects      = effects_in
	this.facets       = facets_in
	this.varlist      = varlist_in
	this.facetlevels  = facetlevels_in
	cmcache           = asarray_create()
	szcache           = asarray_create()

	cp_solved = 0
	is_bal = is_balanced()

	if (is_bal) {
		get_df_balanced()
		get_sscp_all()
	}
	else {
		warn_exp = 1
		for (jj_n=1; jj_n<=length(facetlevels_in); jj_n++) warn_exp = warn_exp * facetlevels_in[jj_n]
		warn_obs = rows(asarray(szcache, invtokens(facets_in, "#")))
		warn_pct = round(100 * (warn_exp - warn_obs) / warn_exp)
		if (!quiet_) {
			printf("{text}Note: unbalanced design detected ({res}%g%%{text} of fully-saturated cells missing).{smcl}\n", warn_pct)
			printf("{text}      Using exact CP-terms method (Brennan, 2001, Sec. 11.1.3).{smcl}\n")
		}
		filter_missing_values()
		get_df_unbalanced()
		cp_solved = 0
		get_cp_all()
		get_cp_coefficients()
		solve_cp_system()
		get_eff_facetlevels()
		if (!cp_solved) {
			// Rank-deficient fallback: revert to harmonic-mean approximation
			get_sscp_all()
		}
	}

	init_inputs(sscplist, effects_in, facets_in, varlist_in, facetlevels_in, df)
}

// Populate this.df.
// Crossed effect "f1#f2":         df = product of (n_f - 1)
// Nested effect "f1|f2|f3|...":   df = product(n_nesting) * product(n_nested - 1)
//   where nesting = all levels after the first |
void mvgstudy::get_df_balanced()
{
	string rowvector	active_facets, nested_facs, nesting_facs, toks
	real scalar			i, j, fi, d, k

	df = J(length(effects), 1, .)
	for (i=1; i<=length(effects); i++) {
		if (strpos(effects[i], "|") > 0) {
			toks         = tokens(effects[i], "|")
			nested_facs  = tokens(subinstr(toks[1], "#", " "))
			// Collect every nesting level: toks[3], toks[5], ... (handles R|T|P, S|R|T|P, etc.)
			nesting_facs = J(1, 0, "")
			for (k=3; k<=length(toks); k=k+2) {
				nesting_facs = (nesting_facs, tokens(subinstr(toks[k], "#", " ")))
			}
			d = 1
			for (j=1; j<=length(nesting_facs); j++) {
				for (fi=1; fi<=length(facets); fi++) {
					if (facets[fi] == nesting_facs[j]) {
						d = d * facetlevels[fi]
						break
					}
				}
			}
			for (j=1; j<=length(nested_facs); j++) {
				for (fi=1; fi<=length(facets); fi++) {
					if (facets[fi] == nested_facs[j]) {
						d = d * (facetlevels[fi] - 1)
						break
					}
				}
			}
		}
		else {
			active_facets = tokens(subinstr(effects[i], "#", " "))
			d = 1
			for (j=1; j<=length(active_facets); j++) {
				for (fi=1; fi<=length(facets); fi++) {
					if (facets[fi] == active_facets[j]) {
						d = d * (facetlevels[fi] - 1)
						break
					}
				}
			}
		}
		df[i] = d
	}
}

// Unbalanced df using observed cell counts.
// Crossed effect "A#B...":  Möbius formula on observed cell counts.
// Nested effect "A#B|C...": Möbius over nested facets conditioned on nesting facets.
//   For |nested_facs|=1 reduces to n_all - n_nesting (the prior two-term formula).
void mvgstudy::get_df_unbalanced()
{
	string rowvector	toks, nested_facs, nesting_facs, active_facets, beta_facets, beta_d
	string scalar		beta_str, all_str, nesting_str, key_d
	real scalar			i, b, j, nbeta, nfacets_alpha, sign, df_val, k
	real scalar			b_d, j_d, nbeta_d, sign_d, df_val_d

	df = J(length(effects), 1, .)

	for (i=1; i<=length(effects); i++) {
		if (strpos(effects[i], "|") > 0) {
			// Nested: Möbius over nested_facs conditioned on nesting_facs
			// df = Σ_{S⊆nested_facs} (-1)^{|nested_facs|-|S|} × n_cells(S∪nesting_facs)
			// Reduces to n_all - n_nesting when |nested_facs|=1 (unchanged for D2,D6,D8)
			toks         = tokens(effects[i], "|")
			nested_facs  = tokens(subinstr(toks[1], "#", " "))
			nesting_facs = J(1, 0, "")
			for (k=3; k<=length(toks); k=k+2) {
				nesting_facs = (nesting_facs, tokens(subinstr(toks[k], "#", " ")))
			}
			df_val_d = 0
			for (b_d=0; b_d<2^length(nested_facs); b_d++) {
				beta_d  = J(1, 0, "")
				nbeta_d = 0
				for (j_d=1; j_d<=length(nested_facs); j_d++) {
					if (mod(floor(b_d / (2^(j_d-1))), 2) == 1) {
						beta_d  = (beta_d, nested_facs[j_d])
						nbeta_d++
					}
				}
				sign_d = (-1)^(length(nested_facs) - nbeta_d)
				key_d  = _canonical_key((beta_d, nesting_facs))
				_cache_cellmeans(key_d)
				if (key_d == "") {
					df_val_d = df_val_d + sign_d * 1
				}
				else {
					df_val_d = df_val_d + sign_d * rows(asarray(cmcache, key_d))
				}
			}
			df[i] = df_val_d
		}
		else {
			// Crossed: Möbius on observed cell counts
			active_facets = tokens(subinstr(effects[i], "#", " "))
			nfacets_alpha = length(active_facets)
			df_val = 0
			for (b=0; b<2^nfacets_alpha; b++) {
				beta_facets = J(1, 0, "")
				nbeta = 0
				for (j=1; j<=nfacets_alpha; j++) {
					if (mod(floor(b / (2^(j-1))), 2) == 1) {
						beta_facets = (beta_facets, active_facets[j])
						nbeta++
					}
				}
				sign = (-1)^(nfacets_alpha - nbeta)
				if (nbeta == 0) {
					df_val = df_val + sign * 1
				}
				else {
					beta_str = _canonical_key(beta_facets)
					_cache_cellmeans(beta_str)
					df_val = df_val + sign * rows(asarray(cmcache, beta_str))
				}
			}
			df[i] = df_val
		}
	}
}

// Populate this.sscplist.
// Balanced: compute n-1 effects directly, residual by subtraction (numerically stable).
// Unbalanced: compute all effects directly (subtraction unstable with unequal cells).
void mvgstudy::get_sscp_all()
{
	real scalar		i, k
	real matrix		Yc, sscp_total, sscp_sum

	sscplist = asarray_create()
	k = cols(Y_data)

	if (is_bal) {
		for (i=1; i<=length(effects)-1; i++) {
			asarray(sscplist, strofreal(i), get_sscp_for_effect(effects[i]))
		}
		Yc = Y_data :- mean(Y_data)
		sscp_total = Yc' * Yc
		sscp_sum = J(k, k, 0)
		for (i=1; i<=length(effects)-1; i++) {
			sscp_sum = sscp_sum + asarray(sscplist, strofreal(i))
		}
		asarray(sscplist, strofreal(length(effects)), sscp_total - sscp_sum)
	}
	else {
		// Unbalanced: compute every effect directly including the last
		for (i=1; i<=length(effects); i++) {
			asarray(sscplist, strofreal(i), get_sscp_for_effect(effects[i]))
		}
	}
}

// Dispatch: nested effects use within-group formula; crossed use Möbius.
real matrix mvgstudy::get_sscp_for_effect(string scalar effect)
{
	string rowvector	alpha_active, beta_facets
	real scalar			nfacets_alpha, ncells_alpha, n_outside, k_out
	real scalar			b, nsubs, j, jj, nbeta, sign, fi, beta_row, mult
	real matrix			entry_alpha, entry_beta, cell_idx_alpha
	real matrix			means_beta, adj, grand
	real rowvector		beta_in_alpha_cols, beta_level_counts, beta_vals
	string scalar		beta_str
	transmorphic		idx_map

	if (strpos(effect, "|") > 0) return(get_sscp_nested(effect))

	alpha_active  = tokens(subinstr(effect, "#", " "))
	nfacets_alpha = length(alpha_active)
	k_out         = cols(Y_data)

	_cache_cellmeans(effect)
	entry_alpha    = asarray(cmcache, effect)
	cell_idx_alpha = entry_alpha[., 1..nfacets_alpha]
	ncells_alpha   = rows(entry_alpha)

	adj   = J(ncells_alpha, k_out, 0)
	nsubs = 2^nfacets_alpha

	for (b=0; b<nsubs; b++) {
		// Build beta = subset of alpha_active indicated by bitmask b
		beta_facets = J(1, 0, "")
		nbeta = 0
		for (j=1; j<=nfacets_alpha; j++) {
			if (mod(floor(b / (2^(j-1))), 2) == 1) {
				beta_facets = (beta_facets, alpha_active[j])
				nbeta++
			}
		}
		sign = (-1)^(nfacets_alpha - nbeta)

		if (nbeta == 0) {
			// Grand mean: same correction for every alpha-cell
			_cache_cellmeans("")
			grand = asarray(cmcache, "")  // 1 × k
			for (j=1; j<=ncells_alpha; j++) {
				adj[j,.] = adj[j,.] + sign * grand[1,.]
			}
		}
		else {
			beta_str = invtokens(beta_facets, "#")
			_cache_cellmeans(beta_str)
			entry_beta  = asarray(cmcache, beta_str)
			means_beta  = entry_beta[., nbeta+1..cols(entry_beta)]

			// Column of cell_idx_alpha that corresponds to each beta facet
			beta_in_alpha_cols = J(1, nbeta, .)
			for (jj=1; jj<=nbeta; jj++) {
				for (j=1; j<=nfacets_alpha; j++) {
					if (beta_facets[jj] == alpha_active[j]) {
						beta_in_alpha_cols[jj] = j
						break
					}
				}
			}

			// Level counts for beta facets (for cellrow formula)
			beta_level_counts = J(1, nbeta, .)
			for (jj=1; jj<=nbeta; jj++) {
				for (fi=1; fi<=length(facets); fi++) {
					if (beta_facets[jj] == facets[fi]) {
						beta_level_counts[jj] = facetlevels[fi]
						break
					}
				}
			}

			// For unbalanced designs, build a packed-key → row lookup so that
			// missing cells don't cause out-of-range subscripts.
			if (!is_bal) {
				idx_map = _build_idx_map(entry_beta[., 1..nbeta], beta_level_counts)
			}

			// Accumulate Möbius term for each alpha-cell
			for (j=1; j<=ncells_alpha; j++) {
				beta_vals = cell_idx_alpha[j, beta_in_alpha_cols]
				beta_row = 1
				mult = 1
				for (jj=nbeta; jj>=1; jj--) {
					beta_row = beta_row + (beta_vals[jj]-1) * mult
					mult     = mult * beta_level_counts[jj]
				}
				if (!is_bal) beta_row = asarray(idx_map, beta_row)
				adj[j,.] = adj[j,.] + sign * means_beta[beta_row,.]
			}
		}
	}

	if (is_bal) {
		n_outside = rows(Y_data) / ncells_alpha
		return(n_outside * (adj' * adj))
	}
	else {
		// Unbalanced: weight each cell's outer product by its own cell size
		real matrix		sscp_u
		real colvector	sz_u
		real scalar		ju
		sscp_u = J(k_out, k_out, 0)
		sz_u   = asarray(szcache, effect)
		for (ju=1; ju<=ncells_alpha; ju++) {
			sscp_u = sscp_u + sz_u[ju] * adj[ju,.]' * adj[ju,.]
		}
		return(sscp_u)
	}
}

// Compute and cache cell means (and cell indices) for one effect level.
// Cache entry format: (cell_idx_cols | means_cols).
// For the grand mean (""), cache entry is just means (1 × k).
// szcache[effect] stores the ncells × 1 vector of cell observation counts.
void mvgstudy::_cache_cellmeans(string scalar effect)
{
	string rowvector	active_facets
	real rowvector		facet_cols
	real matrix			Z_sub, sorted_Z, sorted_Y, cell_idx, means
	real colvector		ord, sizes
	real scalar			i, j, g, start_row, r_row, ncells

	if (asarray_contains(cmcache, effect)) return

	if (effect == "") {
		asarray(cmcache, "", mean(Y_data))
		asarray(szcache, "", rows(Y_data))
		return
	}

	active_facets = tokens(subinstr(effect, "#", " "))
	facet_cols    = J(1, length(active_facets), .)
	for (i=1; i<=length(active_facets); i++) {
		for (j=1; j<=length(facets); j++) {
			if (active_facets[i] == facets[j]) {
				facet_cols[i] = j
				break
			}
		}
	}

	Z_sub    = Z_data[., facet_cols]
	ord      = order(Z_sub, (1..cols(Z_sub)))
	sorted_Z = Z_sub[ord, .]
	sorted_Y = Y_data[ord, .]

	// Count unique cells
	ncells = 1
	for (r_row=2; r_row<=rows(sorted_Z); r_row++) {
		if (rowsum(sorted_Z[r_row,.] :!= sorted_Z[r_row-1,.]) > 0) ncells++
	}

	cell_idx  = J(ncells, cols(Z_sub), .)
	means     = J(ncells, cols(Y_data), .)
	sizes     = J(ncells, 1, .)
	g         = 1
	start_row = 1
	for (r_row=2; r_row<=rows(sorted_Z); r_row++) {
		if (rowsum(sorted_Z[r_row,.] :!= sorted_Z[r_row-1,.]) > 0) {
			cell_idx[g,.] = sorted_Z[start_row,.]
			means[g,.]    = mean(sorted_Y[start_row..r_row-1,.])
			sizes[g]      = r_row - start_row
			g++
			start_row = r_row
		}
	}
	cell_idx[g,.] = sorted_Z[start_row,.]
	means[g,.]    = mean(sorted_Y[start_row..rows(sorted_Y),.])
	sizes[g]      = rows(sorted_Y) - start_row + 1

	asarray(cmcache, effect, (cell_idx, means))
	asarray(szcache, effect, sizes)
}

// Build packed-integer → actual-row lookup for a sparse cell-means matrix.
// cell_idx: ncells × nfacs matrix of integer cell indices (already 1-indexed).
// level_counts: 1 × nfacs vector of facet level counts (for the key formula).
// Returns an asarray: packed_key → row_in_cell_idx.
transmorphic mvgstudy::_build_idx_map(real matrix cell_idx,
                                       real rowvector level_counts)
{
	transmorphic	idx_map
	real scalar		r, key, mult, jj

	idx_map = asarray_create("real")
	for (r=1; r<=rows(cell_idx); r++) {
		key  = 1
		mult = 1
		for (jj=cols(cell_idx); jj>=1; jj--) {
			key  = key  + (cell_idx[r,jj] - 1) * mult
			mult = mult * level_counts[jj]
		}
		asarray(idx_map, key, r)
	}
	return(idx_map)
}

// Return canonical (facets-order) "#"-joined key for a set of facet names.
string scalar mvgstudy::_canonical_key(string rowvector fac_names)
{
	string rowvector	ordered
	real scalar			fi, j

	ordered = J(1, 0, "")
	for (fi=1; fi<=length(facets); fi++) {
		for (j=1; j<=length(fac_names); j++) {
			if (facets[fi] == fac_names[j]) {
				ordered = (ordered, facets[fi])
				break
			}
		}
	}
	return(invtokens(ordered, "#"))
}

// Return 1 if all fully-saturated cells have equal size and no cells are missing.
real scalar mvgstudy::is_balanced()
{
	string scalar	all_key
	real colvector	sizes
	real scalar		fi, expected

	all_key = invtokens(facets, "#")
	_cache_cellmeans(all_key)
	sizes = asarray(szcache, all_key)
	expected = 1
	for (fi=1; fi<=length(facetlevels); fi++) expected = expected * facetlevels[fi]
	return(min(sizes) == max(sizes) & rows(sizes) == expected)
}

// Within-group SSCP for a nested effect "A|B", "A|B|C", "A|B|C|D", etc.
// adj[i] = mean_{all facets at cell i} - mean_{all nesting facets at cell i}
// SSCP = n_outside * adj' * adj
real matrix mvgstudy::get_sscp_nested(string scalar effect)
{
	string rowvector	toks, nested_facs, nesting_facs, all_in_eff, all_eff_facs
	string scalar		all_effect, nesting_effect
	real scalar			k_out, i, j, fi, jj, k, nesting_row, mult, n_all_cells, n_outside
	real matrix			entry_all, entry_nesting, cell_idx_all, means_all, means_nesting, adj
	real rowvector		nesting_in_all_pos, nesting_level_counts
	transmorphic		nest_idx_map
	// Möbius loop variables (used when length(nested_facs) > 1)
	real scalar			b_s, j_s, nbeta_s, sign_s, nkey_s, key_row_s, mult_s, cnt_s
	string rowvector	beta_s, combined_s
	string scalar		key_s
	real matrix			entry_s, means_s
	real rowvector		key_in_all_cols_s, key_level_counts_s
	transmorphic		idx_map_s

	toks         = tokens(effect, "|")
	nested_facs  = tokens(subinstr(toks[1], "#", " "))
	// Collect every nesting level: toks[3], toks[5], ... (handles R|T|P, S|R|T|P, etc.)
	nesting_facs = J(1, 0, "")
	for (k=3; k<=length(toks); k=k+2) {
		nesting_facs = (nesting_facs, tokens(subinstr(toks[k], "#", " ")))
	}
	k_out        = cols(Y_data)

	// Build combined effect key in canonical (this.facets) order
	all_in_eff   = (nested_facs, nesting_facs)
	all_eff_facs = J(1, 0, "")
	for (fi=1; fi<=length(facets); fi++) {
		for (j=1; j<=length(all_in_eff); j++) {
			if (facets[fi] == all_in_eff[j]) {
				all_eff_facs = (all_eff_facs, facets[fi])
				break
			}
		}
	}
	all_effect     = invtokens(all_eff_facs, "#")
	nesting_effect = invtokens(nesting_facs, "#")

	_cache_cellmeans(all_effect)
	_cache_cellmeans(nesting_effect)

	entry_all     = asarray(cmcache, all_effect)
	entry_nesting = asarray(cmcache, nesting_effect)
	n_all_cells   = rows(entry_all)
	cell_idx_all  = entry_all[., 1..length(all_eff_facs)]
	means_all     = entry_all[., length(all_eff_facs)+1..cols(entry_all)]
	means_nesting = entry_nesting[., length(nesting_facs)+1..cols(entry_nesting)]

	if (length(nested_facs) == 1) {
		// Simple nesting (|nested_facs|=1): existing single-subtraction path (unchanged)
		// Column positions in cell_idx_all for each nesting facet (for row-lookup)
		nesting_in_all_pos   = J(1, length(nesting_facs), .)
		nesting_level_counts = J(1, length(nesting_facs), .)
		for (jj=1; jj<=length(nesting_facs); jj++) {
			for (i=1; i<=length(all_eff_facs); i++) {
				if (nesting_facs[jj] == all_eff_facs[i]) {
					nesting_in_all_pos[jj] = i
					break
				}
			}
			for (fi=1; fi<=length(facets); fi++) {
				if (facets[fi] == nesting_facs[jj]) {
					nesting_level_counts[jj] = facetlevels[fi]
					break
				}
			}
		}

		// For unbalanced designs build a lookup from packed-key → row in means_nesting
		if (!is_bal) {
			nest_idx_map = _build_idx_map(
				entry_nesting[., 1..length(nesting_facs)], nesting_level_counts)
		}

		// Compute deviations from nesting-group mean
		adj = J(n_all_cells, k_out, 0)
		for (i=1; i<=n_all_cells; i++) {
			nesting_row = 1
			mult = 1
			for (jj=length(nesting_facs); jj>=1; jj--) {
				nesting_row = nesting_row + (cell_idx_all[i, nesting_in_all_pos[jj]] - 1) * mult
				mult        = mult * nesting_level_counts[jj]
			}
			if (!is_bal) nesting_row = asarray(nest_idx_map, nesting_row)
			adj[i,.] = means_all[i,.] - means_nesting[nesting_row,.]
		}
	}
	else {
		// Crossed-nested (|nested_facs|>1): Möbius inclusion-exclusion over nested_facs
		// adj[abc] = Σ_{S⊆nested_facs} (-1)^{|nested_facs|-|S|} × mean_{S∪nesting_facs}[abc]
		// Fixes the inflated single-subtraction that causes sign flips under unbalance (D4,D5,D7)
		adj = J(n_all_cells, k_out, 0)
		for (b_s=0; b_s<2^length(nested_facs); b_s++) {
			// Enumerate subset beta_s of nested_facs via bitmask b_s
			beta_s  = J(1, 0, "")
			nbeta_s = 0
			for (j_s=1; j_s<=length(nested_facs); j_s++) {
				if (mod(floor(b_s / (2^(j_s-1))), 2) == 1) {
					beta_s  = (beta_s, nested_facs[j_s])
					nbeta_s++
				}
			}
			sign_s = (-1)^(length(nested_facs) - nbeta_s)

			// Combined effect: (beta_s ∪ nesting_facs) in canonical order
			combined_s = (beta_s, nesting_facs)
			key_s      = _canonical_key(combined_s)
			_cache_cellmeans(key_s)
			entry_s = asarray(cmcache, key_s)
			nkey_s  = nbeta_s + length(nesting_facs)
			means_s = entry_s[., nkey_s+1..cols(entry_s)]

			// Column positions in cell_idx_all and level counts, in canonical key order
			key_in_all_cols_s  = J(1, nkey_s, .)
			key_level_counts_s = J(1, nkey_s, .)
			cnt_s = 0
			for (fi=1; fi<=length(all_eff_facs); fi++) {
				for (j=1; j<=length(combined_s); j++) {
					if (all_eff_facs[fi] == combined_s[j]) {
						cnt_s++
						key_in_all_cols_s[cnt_s] = fi
						for (jj=1; jj<=length(facets); jj++) {
							if (facets[jj] == all_eff_facs[fi]) {
								key_level_counts_s[cnt_s] = facetlevels[jj]
								break
							}
						}
						break
					}
				}
			}

			// For unbalanced designs build packed-key → row lookup
			if (!is_bal) {
				idx_map_s = _build_idx_map(entry_s[., 1..nkey_s], key_level_counts_s)
			}

			// Accumulate Möbius term for each all-cell
			for (i=1; i<=n_all_cells; i++) {
				key_row_s = 1
				mult_s    = 1
				for (jj=nkey_s; jj>=1; jj--) {
					key_row_s = key_row_s + (cell_idx_all[i, key_in_all_cols_s[jj]] - 1) * mult_s
					mult_s    = mult_s * key_level_counts_s[jj]
				}
				if (!is_bal) key_row_s = asarray(idx_map_s, key_row_s)
				adj[i,.] = adj[i,.] + sign_s * means_s[key_row_s,.]
			}
		}
	}

	if (is_bal) {
		n_outside = rows(Y_data) / n_all_cells
		return(n_outside * (adj' * adj))
	}
	else {
		real matrix	sscp_n
		real colvector	sz_n
		real scalar		in_n
		sscp_n = J(k_out, k_out, 0)
		sz_n   = asarray(szcache, all_effect)
		for (in_n=1; in_n<=n_all_cells; in_n++) {
			sscp_n = sscp_n + sz_n[in_n] * adj[in_n,.]' * adj[in_n,.]
		}
		return(sscp_n)
	}
}


// Compute effective facet-level counts (harmonic mean of distinct absent-facet
// values per alpha-cell) for each (effect, facet) pair.  Stored in eff_fl_mat
// (n_effects × n_facets); entries for facets present in the effect are set to 1.
void mvgstudy::get_eff_facetlevels()
{
	real scalar		i, j, fi, n_obs, na, ncells_alpha, r, n_dist
	real scalar		hm_sum, alpha_changed, fac_changed
	string rowvector	alpha_facs
	string scalar		key
	real matrix		entry, Z_sorted
	real rowvector	sort_cols, curr_alpha, alpha_col_idx
	real colvector	ord

	eff_fl_mat = J(length(effects), length(facets), 1)
	n_obs = rows(Z_data)

	for (i=1; i<=length(effects); i++) {
		alpha_facs = tokens(subinstr(subinstr(effects[i], "|", " "), "#", " "))
		key = _canonical_key(alpha_facs)
		_cache_cellmeans(key)
		entry = asarray(cmcache, key)
		ncells_alpha = rows(entry)
		na = length(alpha_facs)

		// Column indices in Z_data for each alpha facet (in canonical order)
		alpha_col_idx = J(1, na, .)
		for (fi=1; fi<=na; fi++) {
			for (j=1; j<=length(facets); j++) {
				if (alpha_facs[fi] == facets[j]) {
					alpha_col_idx[fi] = j
					break
				}
			}
		}

		for (j=1; j<=length(facets); j++) {
			if (_facet_in_effect(effects[i], facets[j])) continue  // facet present in effect

			// Sort Z_data by (alpha_cols, absent_fac_col j)
			sort_cols = (alpha_col_idx, j)
			ord      = order(Z_data[., sort_cols], (1..length(sort_cols)))
			Z_sorted = Z_data[ord, sort_cols]

			// Scan: count distinct absent-facet values per alpha-cell
			hm_sum    = 0
			n_dist    = 1
			curr_alpha = Z_sorted[1, 1..na]

			for (r=2; r<=n_obs; r++) {
				alpha_changed = (rowsum(Z_sorted[r, 1..na] :!= curr_alpha) > 0)
				if (alpha_changed) {
					hm_sum     = hm_sum + 1/n_dist
					n_dist     = 1
					curr_alpha = Z_sorted[r, 1..na]
				}
				else {
					fac_changed = (Z_sorted[r, na+1] != Z_sorted[r-1, na+1])
					if (fac_changed) n_dist++
				}
			}
			hm_sum = hm_sum + 1/n_dist  // last alpha-cell

			eff_fl_mat[i, j] = ncells_alpha / hm_sum
		}
	}
}

// (3) Reliability Projections

real matrix mvgstudy::compositevariance()
{
	real matrix 	M
	real scalar 	i
	string scalar 	name
	real vector 	compositevariances

	compositevariances = J(length(effects),1,.)
	for (i=1;i<=length(effects);i++) {
		name = "emcp" + strofreal(i)
		M = asarray(covcomps, name)
		compositevariances[i] = weights'*M*weights
	}
	return(compositevariances)
}

void mvgstudy::absolute_errors_effects(real vector variances)
{
	real matrix 	s1, s2
	string matrix 	erroreffects
	real scalar     fi, ri, found_in_proj
	string rowvector canonical_facets

	s1 = effects':== objmeasurement
	s2 = flip_select(s1)
	objmeasurevariance = select(variances,s1)
	erroreffects       = select(effects',s2)
	errorvariances     = select(variances,s2)
	// All non-object effects absorbed (e.g., every non-object facet was fixed
	// via mvdstudy fix() option).  No error sources remain.
	if (rows(erroreffects) == 0) {
		projectionmat = J(0, 1, "")
		uniqueffects  = J(1, 0, "")
		return
	}
	projectionmat = strtrim(subinstr(subinstr(subinstr(erroreffects,"|"," "), "#"," "),objmeasurement," "))

	// Build canonical uniqueffects: non-object facets that appear in projectionmat,
	// in the canonical order established by c.facets (alphabetical from uniqrows).
	canonical_facets = J(1, 0, "")
	for (fi=1; fi<=length(facets); fi++) {
		if (facets[fi] == objmeasurement) continue
		found_in_proj = 0
		for (ri=1; ri<=rows(projectionmat); ri++) {
			// projectionmat[ri] is space-delimited facet names; compare tokens
			// (avoids strpos false-positives when one facet name is a substring of another)
			if (sum(tokens(projectionmat[ri]) :== facets[fi]) > 0) {
				found_in_proj = 1
				break
			}
		}
		if (found_in_proj) canonical_facets = (canonical_facets, facets[fi])
	}
	uniqueffects = canonical_facets
}

void mvgstudy::relative_errors_effects(real vector variances)
{
	real matrix 	s1, s2
	real colvector  s3
	string matrix 	erroreffects
	real scalar     fi, ri, found_in_proj, ei
	string rowvector canonical_facets

	s1 = effects':==objmeasurement
	s2 = flip_select(s1)
	objmeasurevariance = select(variances,s1)
	erroreffects       = select(effects',s2)
	errorvariances     = select(variances,s2)
	// Token-based membership: keep only error effects that actually contain
	// the object of measurement as a facet.  Avoids substring-prefix false
	// positives (e.g. objmeasurement "p" wrongly matching effect "pr#i").
	s3 = J(rows(erroreffects), 1, 0)
	for (ei = 1; ei <= rows(erroreffects); ei++) {
		s3[ei] = _facet_in_effect(erroreffects[ei], objmeasurement)
	}
	erroreffects       = select(erroreffects,s3)
	errorvariances     = select(errorvariances,s3)
	// All non-object effects absorbed (e.g., every non-object facet was fixed
	// via mvdstudy fix() option, or no error effect contained the object).
	if (rows(erroreffects) == 0) {
		projectionmat = J(0, 1, "")
		uniqueffects  = J(1, 0, "")
		return
	}
	projectionmat = strtrim(subinstr(subinstr(subinstr(erroreffects,"|"," "), "#"," "),objmeasurement," "))

	// Build canonical uniqueffects: non-object facets that appear in projectionmat,
	// in the canonical order established by c.facets (alphabetical from uniqrows).
	canonical_facets = J(1, 0, "")
	for (fi=1; fi<=length(facets); fi++) {
		if (facets[fi] == objmeasurement) continue
		found_in_proj = 0
		for (ri=1; ri<=rows(projectionmat); ri++) {
			// projectionmat[ri] is space-delimited facet names; compare tokens
			// (avoids strpos false-positives when one facet name is a substring of another)
			if (sum(tokens(projectionmat[ri]) :== facets[fi]) > 0) {
				found_in_proj = 1
				break
			}
		}
		if (found_in_proj) canonical_facets = (canonical_facets, facets[fi])
	}
	uniqueffects = canonical_facets
}

real matrix mvgstudy::errorprojections()
{
	string matrix 	res, temp
	real matrix 	first, zero, result, M, numres
	real scalar 	colspan, difference, i, j, fi, fj
	transmorphic 	inner, mats, outerarray, loc
	real vector 	nums, div, temp2

	// uniqueffects is set canonically by relative/absolute_errors_effects
	colspan = length(uniqueffects)

	// --- All non-object facets absorbed (e.g., every non-object facet fixed) ---
	// No error sources remain; reliability is 1 by construction.
	if (colspan == 0) {
		// Single row: (error_var=0, true_var, reliability=1).  Column count is
		// 3 (uniqueffects contributes 0) — matches "Error True Rel." stripe.
		return((0, objmeasurevariance, 1))
	}

	// --- Single-point mode (observe_only=1) ---
	if (observe_only) {
		// Build nums by looking up observed facetlevels for each facet in uniqueffects.
		// This avoids any indexing mismatch between projectionnum and uniqueffects.
		nums = J(1, length(uniqueffects), .)
		for (fi=1; fi<=length(uniqueffects); fi++) {
			for (fj=1; fj<=length(facets); fj++) {
				if (facets[fj] == uniqueffects[fi]) {
					nums[fi] = facetlevels[fj]
					break
				}
			}
		}
		div  = J(rows(projectionmat), 1, 1)
		for (i=1; i<=rows(projectionmat); i++) {
			for (fi=1; fi<=length(uniqueffects); fi++) {
				// Token compare on the space-delimited projectionmat string
				if (sum(tokens(projectionmat[i]) :== uniqueffects[fi]) > 0) {
					div[i] = div[i] * nums[fi]
				}
			}
		}
		result = nums, colsum(errorvariances:/div), objmeasurevariance, ///
		         (objmeasurevariance / (objmeasurevariance + colsum(errorvariances:/div)))
		return(result)
	}

	// --- Sweep mode (observe_only=0) ---
	res = J(0,colspan,"")
	for (i=1;i<=rows(projectionmat);i++) {
		difference = colspan - cols(tokens(projectionmat[i]))
		if (difference !=0) temp = tokens(projectionmat[i]), J(1,difference,"")
		else temp = tokens(projectionmat[i])
		res = res \ temp
	}
	outerarray = asarray_create()
	for (i=1;i<=colspan;i++) {
		inner = asarray_create()
		asarray(outerarray,uniqueffects[i],inner)
	}
	for (i=1;i<=colspan;i++) {
		numres = (res:==uniqueffects[i])
		for (j=1;j<=projectionnum[i];j++) {
			M = numres * j
			asarray(asarray(outerarray,uniqueffects[i]),strofreal(j),M)
		}
	}
	// Use uniqueffects directly as group_keys to guarantee canonical column ordering
	mats  = asarray_create()
	first = asarray(asarray(outerarray, uniqueffects[1]), "1")
	zero  = first * 0
	combine_matrices(outerarray,"", zero, 1, uniqueffects, projectionnum, mats)
	result = J(0,cols(first)+3,.)

	for (loc=asarray_first(mats); loc!=NULL; loc=asarray_next(mats, loc)) {
		nums   = extract_bracket_numbers(asarray_key(mats, loc))
		div    = rowwise_nonzero_product(asarray_contents(mats, loc))
		temp2  = nums, colsum(errorvariances:/div), objmeasurevariance, (objmeasurevariance / (objmeasurevariance + colsum(errorvariances:/div)))
		result = result \ temp2
	}

	return(sort(result,(1..(cols(result)-3))))
}

void mvgstudy::combine_matrices(transmorphic outerarray,
                                string scalar path,
                                real matrix acc,
                                real scalar depth,
                                string vector group_keys,
                                real vector projectionnum,
                                transmorphic mats)
{
	real scalar 	num_groups, i
	string scalar	group, new_path
	real matrix		M

	num_groups = length(group_keys)
	if (depth > num_groups) {
		// Base case: all groups processed
		asarray(mats,path,acc)
		return
	}
	group = group_keys[depth]
	for (i = 1; i <= projectionnum[depth]; i++) {
		M        = asarray(asarray(outerarray, group), strofreal(i))
		new_path = path + group + "[" + strofreal(i) + "] "
		combine_matrices(outerarray, new_path, acc + M, depth + 1, group_keys, projectionnum, mats)
	}
}

real matrix mvgstudy::rowwise_nonzero_product(real matrix X)
{
	real matrix 	result
	real scalar 	i, j

	result = J(rows(X), 1, 1)
	for (i = 1; i <= rows(X); i++) {
		for (j = 1; j <= cols(X); j++) {
			if (X[i, j] != 0) {
				result[i] = result[i] * X[i, j]
			}
		}
	}
	return(result)
}

real vector mvgstudy::extract_bracket_numbers(string scalar s)
{
	string scalar 	pattern
	real vector 	results
	real scalar 	pos, match

	pattern = "\[([0-9]+)\]"
	results = J(1, 0, .)
	pos = 1

	while (regexm(substr(s, pos, .), pattern)) {
		match   = strtoreal(regexs(1))
		results = results , match
		pos     = pos + strpos(substr(s, pos, .), "]")
	}

	return(results)
}

// (3) Profile Variances

// (4) Phase 8: CP Terms Engine

// Remove rows where any outcome variable is missing (listwise deletion).
void mvgstudy::filter_missing_values()
{
	real colvector  keep
	real scalar     r, n_dropped

	keep = J(rows(Y_data), 1, 1)
	for (r = 1; r <= rows(Y_data); r++) {
		if (missing(Y_data[r, .]) > 0) keep[r] = 0
	}
	n_dropped = rows(Y_data) - sum(keep)
	if (n_dropped == 0) return
	Y_data = select(Y_data, keep)
	Z_data = select(Z_data, keep)
	if (!quiet_) printf("{text}Note: %g observations excluded due to missing values in outcome variables.{smcl}\n", n_dropped)
	// Clear caches so they are recomputed on clean data
	cmcache = asarray_create()
	szcache = asarray_create()
}

// Compute CP(α) = M_α' * M_α for each effect α and store in cpcache.
void mvgstudy::get_cp_all()
{
	string scalar    eff, fkey
	real matrix      means_mat
	real rowvector   grand_mean
	real scalar      i

	cpcache = asarray_create()

	_cache_cellmeans("")
	grand_mean = asarray(cmcache, "")
	asarray(cpcache, "", grand_mean' * grand_mean)

	for (i = 1; i <= length(effects); i++) {
		eff   = effects[i]
		fkey  = _flat_key(eff)
		_cache_cellmeans(fkey)
		means_mat = _get_means_only(fkey)
		asarray(cpcache, eff, means_mat' * means_mat)
	}
}

// Build the (m+1) × (m+1) K coefficient matrix from observed cell counts.
void mvgstudy::get_cp_coefficients()
{
	real scalar   m, i, j, n_plus
	string scalar eff_beta, eff_alpha

	m     = length(effects)
	K_mat = J(m+1, m+1, 0)
	n_plus = rows(Y_data)

	// Row 1: grand mean equation  ECP(μ) = μμ' + Σ_α K[μ,α] σ(α)
	K_mat[1, 1] = 1
	for (j = 1; j <= m; j++) {
		K_mat[1, j+1] = _kcoef_mu_alpha(effects[j], n_plus)
	}

	// Rows 2..m+1: effect equations  ECP(β) = n_β μμ' + Σ_α K[β,α] σ(α)
	for (i = 1; i <= m; i++) {
		eff_beta       = effects[i]
		K_mat[i+1, 1]  = _ncells(eff_beta)
		for (j = 1; j <= m; j++) {
			eff_alpha       = effects[j]
			K_mat[i+1, j+1] = _kcoef_beta_alpha(eff_beta, eff_alpha)
		}
	}
}

// K[μ, α] = Σ_{α-cells} n²_a / n_+²
real scalar mvgstudy::_kcoef_mu_alpha(string scalar alpha, real scalar n_plus)
{
	string scalar   fkey
	real colvector  sz

	fkey = _flat_key(alpha)
	_cache_cellmeans(fkey)
	sz = asarray(szcache, fkey)
	return(sum(sz :^ 2) / n_plus^2)
}

// K[β, α] using Brennan Eq. 11.37 (C-terms):
//   K = Σ_{β-cells b} Σ_{γ-cells within b} n²_bγ / n²_b
// where γ = flat_facets(α) not in flat_facets(β).
real scalar mvgstudy::_kcoef_beta_alpha(string scalar beta, string scalar alpha)
{
	string rowvector  flat_a, flat_b, gamma
	string scalar     b_key, bg_key
	real matrix       entry_bg, entry_b
	real colvector    sz_bg, sz_b
	real rowvector    b_col_in_bg, b_vals
	real scalar       K_val, b, bg, n_b, match, fi

	flat_a = _flat_facets(alpha)
	flat_b = _flat_facets(beta)
	gamma  = _setdiff(flat_a, flat_b)

	b_key = _canonical_key(flat_b)

	// Special case: α ⊆ β (γ empty) → K = n_cells_β
	if (length(gamma) == 0) {
		_cache_cellmeans(b_key)
		return(rows(asarray(cmcache, b_key)))
	}

	// General case
	bg_key = _canonical_key((flat_b, gamma))
	_cache_cellmeans(b_key)
	_cache_cellmeans(bg_key)

	entry_b  = asarray(cmcache, b_key)
	entry_bg = asarray(cmcache, bg_key)
	sz_b     = asarray(szcache, b_key)
	sz_bg    = asarray(szcache, bg_key)

	b_col_in_bg = _find_facet_cols_in_key(flat_b, bg_key)

	K_val = 0
	for (b = 1; b <= rows(entry_b); b++) {
		n_b    = sz_b[b]
		b_vals = entry_b[b, 1..length(flat_b)]
		for (bg = 1; bg <= rows(entry_bg); bg++) {
			match = 1
			for (fi = 1; fi <= length(flat_b); fi++) {
				if (entry_bg[bg, b_col_in_bg[fi]] != b_vals[fi]) match = 0
			}
			if (match) K_val = K_val + sz_bg[bg]^2 / n_b^2
		}
	}
	return(K_val)
}

// Solve K * X = CP_stack; populate covcomps; set cp_solved flag.
void mvgstudy::solve_cp_system()
{
	real scalar    m, k, v, vp, idx, i
	real matrix    CP_stack, X_all, comp_mat

	m = length(effects)
	k = cols(Y_data)

	CP_stack = J(m+1, k*k, .)
	for (v = 1; v <= k; v++) {
		for (vp = 1; vp <= k; vp++) {
			idx = (v-1)*k + vp
			CP_stack[1, idx] = asarray(cpcache, "")[v, vp]
			for (i = 1; i <= m; i++) {
				CP_stack[i+1, idx] = asarray(cpcache, effects[i])[v, vp]
			}
		}
	}

	if (rank(K_mat) < rows(K_mat)) {
		if (!quiet_) {
			printf("{err}Warning: CP-terms coefficient matrix is rank-deficient; ")
			printf("falling back to harmonic-mean approximation.{smcl}\n")
		}
		cp_solved = 0
		return
	}

	X_all = luinv(K_mat) * CP_stack

	covcomps = asarray_create()
	for (i = 1; i <= m; i++) {
		comp_mat = J(k, k, .)
		for (v = 1; v <= k; v++) {
			for (vp = 1; vp <= k; vp++) {
				idx = (v-1)*k + vp
				comp_mat[v, vp] = X_all[i+1, idx]
			}
		}
		asarray(covcomps, "emcp" + strofreal(i), comp_mat)
	}

	cp_solved = 1
}

// Return canonical (#-joined, this.facets order) key for all flat facets of effect.
string scalar mvgstudy::_flat_key(string scalar effect)
{
	string rowvector facs
	facs = tokens(subinstr(subinstr(effect, "|", " "), "#", " "))
	return(_canonical_key(facs))
}

// Return sorted unique facet names from an effect string (alphabetical rowvector).
string rowvector mvgstudy::_flat_facets(string scalar effect)
{
	return(sort(uniqrows(tokens(subinstr(subinstr(effect, "|", " "), "#", " "))'), 1)')
}

// Return elements of A not in B (string rowvectors).
string rowvector mvgstudy::_setdiff(string rowvector A, string rowvector B)
{
	string rowvector result
	real scalar      i, j, inB

	result = J(1, 0, "")
	for (i = 1; i <= length(A); i++) {
		inB = 0
		for (j = 1; j <= length(B); j++) {
			if (A[i] == B[j]) inB = 1
		}
		if (!inB) result = (result, A[i])
	}
	return(result)
}

// For each facet in flat_b, return its column index in the facet-index portion
// of the cmcache entry keyed by bg_key (which is in this.facets canonical order).
real rowvector mvgstudy::_find_facet_cols_in_key(string rowvector flat_b,
                                                  string scalar   bg_key)
{
	string rowvector bg_facets
	real rowvector   cidx
	real scalar      i, j

	bg_facets = tokens(subinstr(bg_key, "#", " "))
	cidx = J(1, length(flat_b), .)
	for (i = 1; i <= length(flat_b); i++) {
		for (j = 1; j <= length(bg_facets); j++) {
			if (flat_b[i] == bg_facets[j]) cidx[i] = j
		}
	}
	return(cidx)
}

// Return number of observed cells for an effect (using canonical flat key).
real scalar mvgstudy::_ncells(string scalar effect)
{
	string scalar fkey
	fkey = _flat_key(effect)
	_cache_cellmeans(fkey)
	return(rows(asarray(cmcache, fkey)))
}

// Return only the means columns from a cmcache entry (strips facet index columns).
real matrix mvgstudy::_get_means_only(string scalar key)
{
	real matrix  entry
	real scalar  n_fac_cols

	if (key == "") return(asarray(cmcache, key))
	entry      = asarray(cmcache, key)
	n_fac_cols = length(tokens(subinstr(key, "#", " ")))
	return(entry[., n_fac_cols+1..cols(entry)])
}

// (5) Phase 10: Bootstrap CIs

// General multi-facet bootstrap (Li 2023: boot-p, boot-i, boot-pi, etc.).
// Returns the parent facets of fac if it appears on the LHS of any nested effect ("|").
// E.g., effects={"p","h|p","i|h|p"}: _nesting_parents("h",...) = {"p"}, "i"→{"h","p"}, "p"→{}.
string rowvector mvgstudy::_nesting_parents(string scalar    fac,
                                             string rowvector effects_in)
{
	real scalar      ei, pipe_pos
	string scalar    eff, lhs_str, rhs_str
	string rowvector lhs_facs, rhs_facs, parent_facs

	parent_facs = J(1, 0, "")

	for (ei = 1; ei <= length(effects_in); ei++) {
		eff      = effects_in[ei]
		pipe_pos = strpos(eff, "|")
		if (pipe_pos == 0) continue

		lhs_str  = substr(eff, 1, pipe_pos - 1)
		rhs_str  = substr(eff, pipe_pos + 1, .)
		lhs_facs = tokens(subinstr(lhs_str, "#", " "))
		rhs_facs = tokens(subinstr(subinstr(rhs_str, "|", " "), "#", " "))

		if (sum(lhs_facs :== fac) > 0) parent_facs = (parent_facs, rhs_facs)
	}

	if (length(parent_facs) > 0) parent_facs = uniqrows(parent_facs')'

	return(parent_facs)
}

// bootunit_facs: string rowvector of facet names to resample simultaneously.
// For top-level (non-nested) bootunit: global Cartesian-product resampling (Phase 10a/b).
// For nested bootunit (single facet): within-parent-group resampling (Phase 10c).
// BCa jackknife: at bootunit level (non-nested) or first-parent level (nested).
void mvgstudy::run_bootstrap(real scalar       B,
                              real scalar       ci_alpha,
                              string rowvector  bootunit_facs)
{
	real scalar     nbu, k_out, m, ei, b, fc, fi, pair_idx, j_p, v_i, vp_i
	real scalar     combo, total_combos, remainder
	real scalar     uses_cp, n_boot_failed, n_jack_failed
	real colvector  sel, match_vec
	real matrix     Y_boot, Z_boot, Z_rows, emcp_b, boot_mat, jack_mat
	real vector     bu_cols, bu_n, bu_strides, jack_fl, b_idx
	class mvgstudy  scalar tmp
	string scalar   eff_name
	transmorphic    uniq_store, sv_store
	// Phase 10c: nested-aware resampling
	string rowvector   bu_parents
	real scalar        is_nested_bu, n_jack_reps, jack_col_idx, pc, xi, pi_n, n_X_group
	real colvector     parent_cols_nested, jack_uniq_c, sel_pc, X_vals_in_group, sampled_X_group
	real matrix        parent_combos
	// Phase 12c: column of Z_data for person-weighted mean grouping
	real scalar        pw_col_

	k_out         = cols(Y_data)
	m             = length(effects)
	uses_cp       = !is_bal          // original data uses CP path (unbalanced)
	n_boot_failed = 0
	n_jack_failed = 0

	// If bootunit_facs is empty, default to first single-facet effect
	if (length(bootunit_facs) == 0) {
		bootunit_facs = tokens(subinstr(subinstr(effects[1], "|", " "), "#", " "))[1..1]
	}
	nbu = length(bootunit_facs)

	// Resolve column indices in Z_data and unique value counts per bootunit facet
	bu_cols = J(1, nbu, 0)
	bu_n    = J(1, nbu, 0)
	for (fi = 1; fi <= nbu; fi++) {
		for (fc = 1; fc <= length(facets); fc++) {
			if (facets[fc] == bootunit_facs[fi]) {
				bu_cols[fi] = fc
				bu_n[fi]    = rows(uniqrows(Z_data[., fc]))
				break
			}
		}
	}
	boot_cluster_col = bu_cols[1]  // default; may differ from jackknife col for nested bootunit

	// Pre-store unique values for each bootunit facet (string keys: "1", "2", ...)
	uniq_store = asarray_create()
	for (fi = 1; fi <= nbu; fi++) {
		asarray(uniq_store, strofreal(fi), uniqrows(Z_data[., bu_cols[fi]]))
	}

	// Phase 10c: detect if bootunit is a nested facet (single bootunit only)
	bu_parents   = J(1, 0, "")
	is_nested_bu = 0
	if (nbu == 1) bu_parents = _nesting_parents(bootunit_facs[1], effects)
	if (length(bu_parents) > 0) {
		is_nested_bu = 1
		printf("{text}Note: %s is nested within %s; using within-group cluster bootstrap.\n",
		       bootunit_facs[1], invtokens(bu_parents, " "))
		parent_cols_nested = J(length(bu_parents), 1, 0)
		for (pi_n = 1; pi_n <= length(bu_parents); pi_n++) {
			for (fc = 1; fc <= length(facets); fc++) {
				if (facets[fc] == bu_parents[pi_n]) {
					parent_cols_nested[pi_n] = fc
					break
				}
			}
		}
		{
			real rowvector _pc_rv
			_pc_rv        = parent_cols_nested'
			parent_combos = uniqrows(Z_data[., _pc_rv])
		}
		// Jackknife at first parent level (immediate nesting context)
		jack_col_idx = parent_cols_nested[1]
		jack_uniq_c  = uniqrows(Z_data[., jack_col_idx])
		n_jack_reps  = rows(jack_uniq_c)
	}
	else {
		// Non-nested: jackknife at bootunit level (existing behavior)
		jack_col_idx = bu_cols[1]
		jack_uniq_c  = asarray(uniq_store, "1")
		n_jack_reps  = bu_n[1]
	}

	// Strides for Cartesian-product enumeration (used by non-nested path only)
	bu_strides    = J(1, nbu, 1)
	for (fi = nbu-1; fi >= 1; fi--) bu_strides[fi] = bu_strides[fi+1] * bu_n[fi+1]
	total_combos  = bu_strides[1] * bu_n[1]

	// Initialize stores (B bootstrap reps, N_jack jackknife reps at first-facet level)
	boot_store = asarray_create()
	jack_store = asarray_create()
	for (ei = 1; ei <= m; ei++) {
		eff_name = "emcp" + strofreal(ei)
		asarray(boot_store, eff_name, J(B,            k_out * k_out, .))
		asarray(jack_store, eff_name, J(n_jack_reps,  k_out * k_out, .))
	}
	// Phase 12b: per-rep outcome means (for CRR bootstrap propagation)
	boot_means = J(B,           k_out, .)
	jack_means = J(n_jack_reps, k_out, .)
	// Phase 12c: person-weighted per-rep means, grouped by the first facet of
	// effects[1] (the canonical object of measurement)
	pw_col_      = _crr_obj_col()
	boot_pwmeans = J(B,           k_out, .)
	jack_pwmeans = J(n_jack_reps, k_out, .)

	b_idx = J(1, nbu, 0)

	//----//
	// Bootstrap replications
	//----//
	for (b = 1; b <= B; b++) {
		if (is_nested_bu) {
			// Phase 10c: within-parent-group resampling for nested bootunit facet.
			// For each parent combination, independently resample the child facet within that group.
			Y_boot = J(0, k_out, .)
			Z_boot = J(0, cols(Z_data), .)
			for (pc = 1; pc <= rows(parent_combos); pc++) {
				match_vec = J(rows(Z_data), 1, 1)
				for (pi_n = 1; pi_n <= length(bu_parents); pi_n++) {
					match_vec = match_vec :& (Z_data[., parent_cols_nested[pi_n]] :== parent_combos[pc, pi_n])
				}
				sel_pc = selectindex(match_vec)
				if (rows(sel_pc) == 0) continue
				X_vals_in_group = uniqrows(Z_data[sel_pc, bu_cols[1]])
				n_X_group       = rows(X_vals_in_group)
				sampled_X_group = X_vals_in_group[ceil(uniform(n_X_group, 1) * n_X_group), .]
				for (xi = 1; xi <= n_X_group; xi++) {
					sel    = selectindex(Z_data[sel_pc, bu_cols[1]] :== sampled_X_group[xi])
					Z_rows = Z_data[sel_pc[sel], .]
					Z_rows[., bu_cols[1]] = J(rows(Z_rows), 1, xi)
					Y_boot = Y_boot \ Y_data[sel_pc[sel], .]
					Z_boot = Z_boot \ Z_rows
				}
			}
		}
		else {
			// Global Cartesian-product resampling (existing path for top-level bootunit facets)
			sv_store = asarray_create()
			for (fi = 1; fi <= nbu; fi++) {
				asarray(sv_store, strofreal(fi),
					asarray(uniq_store, strofreal(fi))[ceil(uniform(bu_n[fi], 1) * bu_n[fi]), .])
			}

			Y_boot = J(0, k_out, .)
			Z_boot = J(0, cols(Z_data), .)

			for (combo = 1; combo <= total_combos; combo++) {
				remainder = combo - 1
				for (fi = 1; fi <= nbu; fi++) {
					b_idx[fi] = floor(remainder / bu_strides[fi]) + 1
					remainder  = mod(remainder, bu_strides[fi])
				}

				match_vec = J(rows(Z_data), 1, 1)
				for (fi = 1; fi <= nbu; fi++) {
					match_vec = match_vec :& (Z_data[., bu_cols[fi]] :== asarray(sv_store, strofreal(fi))[b_idx[fi]])
				}
				sel = selectindex(match_vec)
				if (rows(sel) == 0) continue

				Z_rows = Z_data[sel, .]
				for (fi = 1; fi <= nbu; fi++) {
					Z_rows[., bu_cols[fi]] = J(rows(Z_rows), 1, b_idx[fi])
				}

				Y_boot = Y_boot \ Y_data[sel, .]
				Z_boot = Z_boot \ Z_rows
			}
		}

		tmp.quiet_ = 1
		tmp.init_inputs_direct(Y_boot, Z_boot, effects, facets, varlist, facetlevels)
		tmp.mvgstudy_main_routine()

		// For unbalanced designs: skip rep if K was rank-deficient (avoids mixing estimators)
		if (uses_cp & !tmp.cp_solved) {
			n_boot_failed++
			continue
		}

		// Phase 12b/12c: store rep means (post-filter data, aligned with components)
		boot_means[b, .]   = mean(tmp.Y_data)
		boot_pwmeans[b, .] = _pw_means(tmp.Y_data, tmp.Z_data[., pw_col_])

		for (ei = 1; ei <= m; ei++) {
			eff_name = "emcp" + strofreal(ei)
			emcp_b   = asarray(tmp.covcomps, eff_name)
			boot_mat = asarray(boot_store, eff_name)
			pair_idx = 0
			for (v_i = 1; v_i <= k_out; v_i++) {
				for (vp_i = 1; vp_i <= k_out; vp_i++) {
					pair_idx++
					boot_mat[b, pair_idx] = emcp_b[v_i, vp_i]
				}
			}
			asarray(boot_store, eff_name, boot_mat)
		}
	}

	//----//
	// Jackknife for BCa acceleration: leave-one-level-out at jackknife level.
	// Non-nested bootunit: jackknife at the bootunit facet level (existing behavior).
	// Nested bootunit: jackknife at the first parent level (captures cluster influence).
	//----//
	jack_fl               = facetlevels
	jack_fl[jack_col_idx] = n_jack_reps - 1

	for (j_p = 1; j_p <= n_jack_reps; j_p++) {
		sel    = selectindex(Z_data[., jack_col_idx] :!= jack_uniq_c[j_p])
		Y_boot = Y_data[sel, .]
		Z_boot = Z_data[sel, .]

		tmp.quiet_ = 1
		tmp.init_inputs_direct(Y_boot, Z_boot, effects, facets, varlist, jack_fl)
		tmp.mvgstudy_main_routine()

		// Skip jackknife rep if CP failed — jack_mat[j_p,.] stays missing; bca_acceleration handles it
		if (uses_cp & !tmp.cp_solved) {
			n_jack_failed++
		}
		else {
			// Phase 12b/12c: store jackknife rep means
			jack_means[j_p, .]   = mean(tmp.Y_data)
			jack_pwmeans[j_p, .] = _pw_means(tmp.Y_data, tmp.Z_data[., pw_col_])
			for (ei = 1; ei <= m; ei++) {
				eff_name = "emcp" + strofreal(ei)
				emcp_b   = asarray(tmp.covcomps, eff_name)
				jack_mat = asarray(jack_store, eff_name)
				pair_idx = 0
				for (v_i = 1; v_i <= k_out; v_i++) {
					for (vp_i = 1; vp_i <= k_out; vp_i++) {
						pair_idx++
						jack_mat[j_p, pair_idx] = emcp_b[v_i, vp_i]
					}
				}
				asarray(jack_store, eff_name, jack_mat)
			}
		}
	}

	// Warn if any reps had rank-deficient K (unbalanced designs only)
	if (n_boot_failed > 0) {
		printf("{text}Note: %g of %g bootstrap reps had a rank-deficient CP coefficient matrix and were skipped.\n", n_boot_failed, B)
		if (n_boot_failed > B / 2) {
			printf("{err}Warning: more than half of bootstrap reps failed — SEs may be unreliable.{smcl}\n")
		}
	}

	boot_B        = B
	boot_ci_alpha = ci_alpha
}

// BCa acceleration from jackknife leave-one-cluster-out estimates.
real scalar mvgstudy::bca_acceleration(real colvector jack_vals)
{
	real colvector diffs
	real scalar    sum2, sum3

	if (rows(jack_vals) < 2) return(0)
	diffs = mean(jack_vals) :- jack_vals
	sum2  = sum(diffs :^ 2)
	sum3  = sum(diffs :^ 3)
	if (sum2 == 0) return(0)
	return(sum3 / (6 * sum2^1.5))
}

// BCa confidence interval for one parameter element.
real rowvector mvgstudy::bca_ci(real colvector boot_vals,
                                 real scalar    theta_hat,
                                 real scalar    a,
                                 real scalar    ci_alpha)
{
	real scalar    B, prop_less, z0, zalpha_lo, zalpha_hi, alpha1, alpha2, k_lo, k_hi
	real colvector sorted_boot

	B = rows(boot_vals)
	if (B == 0) return((., .))

	sorted_boot = sort(boot_vals, 1)

	prop_less = mean(boot_vals :< theta_hat)
	if (prop_less <= 0) prop_less = 1 / (2 * B)
	if (prop_less >= 1) prop_less = 1 - 1 / (2 * B)
	z0 = invnormal(prop_less)

	zalpha_lo = invnormal(ci_alpha / 2)
	zalpha_hi = invnormal(1 - ci_alpha / 2)
	alpha1    = normal(z0 + (z0 + zalpha_lo) / (1 - a * (z0 + zalpha_lo)))
	alpha2    = normal(z0 + (z0 + zalpha_hi) / (1 - a * (z0 + zalpha_hi)))

	if (missing(alpha1) | alpha1 <= 0) alpha1 = 0.001
	if (alpha1 >= 1)                   alpha1 = 0.999
	if (missing(alpha2) | alpha2 <= 0) alpha2 = 0.001
	if (alpha2 >= 1)                   alpha2 = 0.999

	k_lo = max((1, floor(alpha1 * B)))
	k_hi = min((B, ceil(alpha2 * B)))

	return((sorted_boot[k_lo], sorted_boot[k_hi]))
}

// Phase 10d: Propagate G-study bootstrap reps through D-study; store BCa summaries
// per projection variable in dstudy_boot_summary (asarray: varname → (n_rows*3)×4).
// Must be called after init_dstudyinputs(), mvdstudy_main_routine(), and
// run_bootstrap() have completed.  Uses existing projections asarray for n_rows
// (avoids re-calling mvdstudy_main_routine with original covcomps unnecessarily).
void mvgstudy::run_dstudy_bootstrap(real scalar ci_alpha_ds)
{
	real scalar      m, ei, b, k_out, j_p, q_n, n_jack_total, n_rows, n_quant
	real scalar      rr_q, qq_q, rep_valid, src_c1, src_c2, dst_c1, dst_c2
	real scalar      theta_hat, a_bca, se_bca, n_proj_cols
	real matrix      boot_mat, jack_mat, proj_result, bs_mat_v, jk_mat_v, summ
	real matrix      boot_check, jack_check
	real rowvector   flat_row, ci_bca
	real colvector   boot_vals, jack_vals
	string scalar    eff_name, varname
	transmorphic     orig_covcomps, dstudy_bs_raw, dstudy_jk_raw, loc

	k_out        = cols(Y_data)
	m            = length(effects)
	n_quant      = 3
	n_jack_total = rows(asarray(jack_store, "emcp1"))

	// Save original covcomps
	orig_covcomps = asarray_create()
	for (ei = 1; ei <= m; ei++) {
		eff_name = "emcp" + strofreal(ei)
		asarray(orig_covcomps, eff_name, asarray(covcomps, eff_name))
	}

	// Determine projection structure from the existing projections asarray
	// (already computed by mvdstudy_main_routine() in the Stata wrapper).
	loc         = asarray_first(projections)
	proj_result = asarray_contents(projections, loc)
	n_rows      = rows(proj_result)
	n_proj_cols = cols(proj_result)
	flat_row    = J(1, n_rows * n_quant, .)

	// Initialize raw bootstrap and jackknife stores
	dstudy_bs_raw = asarray_create()
	dstudy_jk_raw = asarray_create()
	for (loc = asarray_first(projections); loc != NULL; loc = asarray_next(projections, loc)) {
		varname = asarray_key(projections, loc)
		asarray(dstudy_bs_raw, varname, J(boot_B,        n_rows * n_quant, .))
		asarray(dstudy_jk_raw, varname, J(n_jack_total,  n_rows * n_quant, .))
	}

	// Bootstrap loop: reconstruct covcomps from boot_store each rep
	boot_check = asarray(boot_store, "emcp1")
	for (b = 1; b <= boot_B; b++) {
		rep_valid = !missing(boot_check[b, 1])
		if (!rep_valid) continue

		for (ei = 1; ei <= m; ei++) {
			eff_name = "emcp" + strofreal(ei)
			boot_mat = asarray(boot_store, eff_name)
			asarray(covcomps, eff_name, rowshape(boot_mat[b, .], k_out))
		}

		mvdstudy_main_routine()

		for (loc = asarray_first(projections); loc != NULL; loc = asarray_next(projections, loc)) {
			varname     = asarray_key(projections, loc)
			proj_result = asarray_contents(projections, loc)
			src_c1      = cols(proj_result) - n_quant + 1
			src_c2      = cols(proj_result)
			for (rr_q = 1; rr_q <= n_rows; rr_q++) {
				dst_c1 = (rr_q - 1) * n_quant + 1
				dst_c2 = rr_q * n_quant
				flat_row[1, dst_c1..dst_c2] = proj_result[rr_q, src_c1..src_c2]
			}
			bs_mat_v = asarray(dstudy_bs_raw, varname)
			bs_mat_v[b, .] = flat_row
			asarray(dstudy_bs_raw, varname, bs_mat_v)
		}
	}

	// Jackknife loop: reconstruct covcomps from jack_store each rep
	jack_check = asarray(jack_store, "emcp1")
	for (j_p = 1; j_p <= n_jack_total; j_p++) {
		rep_valid = !missing(jack_check[j_p, 1])
		if (!rep_valid) continue

		for (ei = 1; ei <= m; ei++) {
			eff_name = "emcp" + strofreal(ei)
			jack_mat = asarray(jack_store, eff_name)
			asarray(covcomps, eff_name, rowshape(jack_mat[j_p, .], k_out))
		}

		mvdstudy_main_routine()

		for (loc = asarray_first(projections); loc != NULL; loc = asarray_next(projections, loc)) {
			varname     = asarray_key(projections, loc)
			proj_result = asarray_contents(projections, loc)
			src_c1      = cols(proj_result) - n_quant + 1
			src_c2      = cols(proj_result)
			for (rr_q = 1; rr_q <= n_rows; rr_q++) {
				dst_c1 = (rr_q - 1) * n_quant + 1
				dst_c2 = rr_q * n_quant
				flat_row[1, dst_c1..dst_c2] = proj_result[rr_q, src_c1..src_c2]
			}
			jk_mat_v = asarray(dstudy_jk_raw, varname)
			jk_mat_v[j_p, .] = flat_row
			asarray(dstudy_jk_raw, varname, jk_mat_v)
		}
	}

	// Restore original covcomps and re-run D-study to restore projections
	for (ei = 1; ei <= m; ei++) {
		eff_name = "emcp" + strofreal(ei)
		asarray(covcomps, eff_name, asarray(orig_covcomps, eff_name))
	}
	mvdstudy_main_routine()

	// Compute BCa summary for each projection variable
	dstudy_boot_summary = asarray_create()
	for (loc = asarray_first(projections); loc != NULL; loc = asarray_next(projections, loc)) {
		varname     = asarray_key(projections, loc)
		proj_result = asarray_contents(projections, loc)
		bs_mat_v    = asarray(dstudy_bs_raw, varname)
		jk_mat_v    = asarray(dstudy_jk_raw, varname)
		summ        = J(n_rows * n_quant, 4, .)

		for (q_n = 1; q_n <= n_rows * n_quant; q_n++) {
			rr_q      = ceil(q_n / n_quant)
			qq_q      = mod(q_n - 1, n_quant) + 1
			theta_hat = proj_result[rr_q, n_proj_cols - n_quant + qq_q]

			boot_vals = select(bs_mat_v[., q_n], !missing(bs_mat_v[., q_n]))
			jack_vals = select(jk_mat_v[., q_n], !missing(jk_mat_v[., q_n]))

			se_bca = (rows(boot_vals) > 1 ? sqrt(variance(boot_vals)) : .)
			a_bca  = (rows(jack_vals) >= 2 ? bca_acceleration(jack_vals) : 0)
			if (rows(boot_vals) > 0) {
				ci_bca = bca_ci(boot_vals, theta_hat, a_bca, ci_alpha_ds)
			}
			else {
				ci_bca = J(1, 2, .)
			}
			summ[q_n, .] = (theta_hat, se_bca, ci_bca)
		}
		asarray(dstudy_boot_summary, varname, summ)
	}
}

// Push D-study bootstrap summary to Stata as stata_matname with row/col stripes.
// coeff_name: "erho2" (relative) or "phi" (absolute).
void mvgstudy::push_dstudy_r_matrix(string scalar varname,
                                     string scalar stata_matname,
                                     string scalar coeff_name)
{
	real scalar      n_rows, n_quant, q_n, rr_q, qq_q
	real matrix      summ
	string matrix    row_stripe, col_stripe
	string rowvector quant_names
	string scalar    err_name

	summ    = asarray(dstudy_boot_summary, varname)
	n_quant = 3
	n_rows  = rows(summ) / n_quant

	if (coeff_name == "erho2") err_name = "var_delta"
	else                       err_name = "var_Delta"

	quant_names = (err_name, "var_tau", coeff_name)

	row_stripe = J(rows(summ), 2, "")
	for (rr_q = 1; rr_q <= n_rows; rr_q++) {
		for (qq_q = 1; qq_q <= n_quant; qq_q++) {
			q_n = (rr_q - 1) * n_quant + qq_q
			if (n_rows == 1) {
				row_stripe[q_n, 2] = quant_names[qq_q]
			}
			else {
				row_stripe[q_n, 2] = "r" + strofreal(rr_q) + "_" + quant_names[qq_q]
			}
		}
	}

	st_matrix(stata_matname, summ)

	col_stripe       = J(4, 2, "")
	col_stripe[., 2] = ("estimate" \ "se" \ "ci_lo" \ "ci_hi")
	st_matrixcolstripe(stata_matname, col_stripe)
	st_matrixrowstripe(stata_matname, row_stripe)
}

// Compute bootstrap summary (k^2 x 4) matrix: estimate, SE, BCa CI lo/hi.
// Pushed to Stata as stata_matname with named row/col stripes.
void mvgstudy::push_boot_r_matrix(string scalar eff_name,
                                   string scalar stata_matname,
                                   real scalar   ci_alpha)
{
	real scalar     k_out, v_i, vp_i, pair_idx
	real matrix     emcp_orig, boot_mat, jack_mat, summary
	real colvector  boot_col, jack_col
	real scalar     theta_hat, a, se
	real rowvector  ci
	string matrix   row_stripe, col_stripe

	k_out     = cols(Y_data)
	summary   = J(k_out * k_out, 4, .)
	emcp_orig = asarray(covcomps, eff_name)
	boot_mat  = asarray(boot_store, eff_name)
	jack_mat  = asarray(jack_store, eff_name)

	pair_idx = 0
	for (v_i = 1; v_i <= k_out; v_i++) {
		for (vp_i = 1; vp_i <= k_out; vp_i++) {
			pair_idx++
			theta_hat = emcp_orig[v_i, vp_i]
			boot_col  = select(boot_mat[., pair_idx], !missing(boot_mat[., pair_idx]))
			jack_col  = select(jack_mat[., pair_idx], !missing(jack_mat[., pair_idx]))

			se = (rows(boot_col) > 1 ? sqrt(variance(boot_col)) : .)
			a  = bca_acceleration(jack_col)
			ci = (rows(boot_col) > 0 ? bca_ci(boot_col, theta_hat, a, ci_alpha) : (., .))

			summary[pair_idx, .] = (theta_hat, se, ci)
		}
	}

	st_matrix(stata_matname, summary)

	col_stripe       = J(4, 2, "")
	col_stripe[., 2] = ("estimate" \ "se" \ "ci_lo" \ "ci_hi")
	st_matrixcolstripe(stata_matname, col_stripe)

	row_stripe = J(k_out * k_out, 2, "")
	pair_idx   = 0
	for (v_i = 1; v_i <= k_out; v_i++) {
		for (vp_i = 1; vp_i <= k_out; vp_i++) {
			pair_idx++
			row_stripe[pair_idx, 2] = varlist[v_i] + "_" + varlist[vp_i]
		}
	}
	st_matrixrowstripe(stata_matname, row_stripe)
}

// Phase 11: Fixed-facet D-study augmentation.
// For each effect in the design that contains one or more fixed facets:
//   - add covcomps[effect] / divisor  to covcomps[target_effect], where
//     target_effect is the same effect with fixed-facet tokens removed and
//     divisor is the product of the fixed n values for all fixed facets present.
//   - remove the absorbed effect from the active design.
// Also removes fixed facets from this.facets / this.facetlevels so that
// init_dstudyinputs() validates facetnum() against the correct (reduced) count.
// Originals are saved for restoration by restore_fixed_facets().
void mvgstudy::apply_fixed_facets(string scalar      obj_name,
                                   string rowvector   fixed_names,
                                   real rowvector     fixed_ns)
{
	real scalar      m, nf, nfac, ei, fi, fk, ej, divisor, has_fixed, target_ei, new_idx, found
	real scalar      ntok, nfac_new, tk
	string scalar    eff, target_name, eff_key, target_key, new_key, tok
	string rowvector eff_tokens, remaining_tokens, new_effects_row, new_facets_row
	real rowvector   is_fixed_tok, keep, keep_fac
	real colvector   new_facetlevels_col
	transmorphic     new_covcomps

	m    = length(effects)
	nf   = length(fixed_names)
	nfac = length(facets)

	// --- Validate ---
	for (fi = 1; fi <= nf; fi++) {
		if (fixed_names[fi] == obj_name) {
			errprintf("mvdstudy: fix() facet '%s' is the object of measurement; cannot fix the object\n",
			          fixed_names[fi])
			exit(198)
		}
		found = 0
		for (ei = 1; ei <= m; ei++) {
			if (effects[ei] == fixed_names[fi]) {
				found = 1
				break
			}
		}
		if (!found) {
			errprintf("mvdstudy: fix() facet '%s' not found as a single-facet effect in the design\n",
			          fixed_names[fi])
			exit(198)
		}
	}

	// --- Save originals ---
	orig_effects     = effects
	orig_facets      = facets
	orig_facetlevels = facetlevels
	orig_covcomps    = asarray_create()
	for (ei = 1; ei <= m; ei++) {
		eff_key = "emcp" + strofreal(ei)
		asarray(orig_covcomps, eff_key, asarray(covcomps, eff_key))
	}

	keep = J(1, m, 1)   // 1 = keep, 0 = absorbed into target

	// --- Pass 1: augment targets ---
	for (ei = 1; ei <= m; ei++) {
		eff        = effects[ei]
		tok        = subinstr(subinstr(eff, "#", " ", .), "|", " ", .)
		eff_tokens = tokens(tok)
		ntok       = length(eff_tokens)

		// Determine which fixed facets appear in this effect
		has_fixed    = 0
		divisor      = 1
		is_fixed_tok = J(1, ntok, 0)
		for (fi = 1; fi <= nf; fi++) {
			for (tk = 1; tk <= ntok; tk++) {
				if (eff_tokens[tk] == fixed_names[fi]) {
					has_fixed        = 1
					divisor          = divisor * fixed_ns[fi]
					is_fixed_tok[tk] = 1
				}
			}
		}

		if (!has_fixed) continue

		// Fixed-facet augmentation only supports fully crossed (#) effects.
		// Nested effects (| present) are removed but not augmented.
		if (strpos(eff, "|") > 0) {
			printf("{txt}Note: effect '%s' contains a fixed facet but uses nesting (|); removed from design without augmentation.\n", eff)
			keep[ei] = 0
			continue
		}

		// Build remaining_tokens by loop (avoids select() orientation issues)
		remaining_tokens = J(1, 0, "")
		for (tk = 1; tk <= ntok; tk++) {
			if (!is_fixed_tok[tk]) remaining_tokens = (remaining_tokens, eff_tokens[tk])
		}

		if (length(remaining_tokens) == 0) {
			// Effect consists entirely of fixed facets — discard (fixed offset, not error)
			keep[ei] = 0
			continue
		}

		target_name = invtokens(remaining_tokens, "#")

		target_ei = 0
		for (ej = 1; ej <= m; ej++) {
			if (effects[ej] == target_name) {
				target_ei = ej
				break
			}
		}

		if (target_ei == 0) {
			// Target effect not found (non-full-factorial design) — discard without augmenting
			keep[ei] = 0
			continue
		}

		// Augment: covcomps[target] += covcomps[source] / divisor
		eff_key    = "emcp" + strofreal(ei)
		target_key = "emcp" + strofreal(target_ei)
		asarray(covcomps, target_key,
		        asarray(covcomps, target_key) + asarray(covcomps, eff_key) / divisor)
		keep[ei] = 0
	}

	// --- Pass 2: rebuild effects and covcomps with only kept effects ---
	new_covcomps    = asarray_create()
	new_effects_row = J(1, 0, "")
	new_idx = 0
	for (ei = 1; ei <= m; ei++) {
		if (!keep[ei]) continue
		new_idx++
		eff_key = "emcp" + strofreal(ei)
		new_key = "emcp" + strofreal(new_idx)
		asarray(new_covcomps, new_key, asarray(covcomps, eff_key))
		new_effects_row = (new_effects_row, effects[ei])
	}
	effects  = new_effects_row
	covcomps = new_covcomps

	// --- Remove fixed facets from facets (rowvec) and facetlevels (colvec) ---
	// facets is 1×nfac rowvec; facetlevels is nfac×1 colvec (from st_matrix stacking)
	keep_fac = J(1, nfac, 1)
	for (fi = 1; fi <= nf; fi++) {
		for (fk = 1; fk <= nfac; fk++) {
			if (facets[fk] == fixed_names[fi]) {
				keep_fac[fk] = 0
				break
			}
		}
	}
	nfac_new           = sum(keep_fac)
	new_facets_row     = J(1, nfac_new, "")
	new_facetlevels_col = J(nfac_new, 1, .)
	new_idx = 0
	for (fk = 1; fk <= nfac; fk++) {
		if (!keep_fac[fk]) continue
		new_idx++
		new_facets_row[new_idx]      = facets[fk]
		new_facetlevels_col[new_idx] = facetlevels[fk]
	}
	facets      = new_facets_row
	facetlevels = new_facetlevels_col
}

// Restore effects, covcomps, facets, facetlevels to their pre-augmentation state.
void mvgstudy::restore_fixed_facets()
{
	effects      = orig_effects
	covcomps     = orig_covcomps
	facets       = orig_facets
	facetlevels  = orig_facetlevels
}

// Phase 12: Projected construct-relevant reliability.
// CRR = lambda x Erho2_DCF:PL evaluated at the planned number of lessons nL:
//
//           Jbar^2 s_pi
//   CRR = ----------------------------------------------------------------
//           Jbar^2 s_pi + s_aP + s_jP (mupi^2 + s_pi)
//         + [ s_aL + (Jbar^2 + s_jP) s_piL
//             + s_jL (mupi^2 + s_pi + s_piL) + sigma2_eps ] / nL
//
// (the universe-score variance cancels between lambda and the coefficient).
// Components come from covcomps: emcp1 = person (P) matrix, emcp2 = L:P matrix,
// diagonals at the A, J, prevalence variable positions.  Means come from Y_data.
// sig_mode: 0 = sigma2_eps 0 (default; utterance-sampling noise is already
//               absorbed in the estimated L:P components),
//           1 = user-supplied sigma2_eps (s1),
//           2 = closed form (mupi*nu1^2 + (1-mupi)*nu0^2
//                            + mupi*(1-mupi)*Jbar^2) / n
//               with s1 = nu0^2, s2 = nu1^2, s3 = n (utterances per lesson).
// Negative component estimates are truncated at zero (count reported); if
// truncation drives both person-level DCF components to zero, the ub flag is
// set (lambda = 1 by construction; CRR is an upper bound).
// Results go to the 9 Stata scalars named in scnames_str, in order:
// crr, Abar, lambda, erho2_dcfpl, Jbar, mupi, sigmae, ntrunc, ub;
// the 2x3 post-truncation component matrix (rows P, L_P; cols A, J, prev)
// goes to comps_name.
void mvgstudy::compute_crr(string scalar avn, string scalar jvn,
                           string scalar pvn,
                           real scalar nL, real scalar sig_mode,
                           real scalar s1, real scalar s2, real scalar s3,
                           real scalar pw,
                           string scalar comps_name, string scalar scnames_str)
{
	real rowvector   idx, means, res
	real matrix      comps
	string rowvector scn
	string matrix    rstripe, cstripe

	scn   = tokens(scnames_str)
	idx   = _crr_resolve_idx(avn, jvn, pvn)
	// pw = 1: person-weighted means (mean of within-object means; differs from
	// the observation-weighted grand mean only under unbalance)
	if (pw) means = _pw_means(Y_data, Z_data[., _crr_obj_col()])
	else    means = mean(Y_data)
	res   = _crr_core(asarray(covcomps, "emcp1"), asarray(covcomps, "emcp2"),
	                  means, idx[1], idx[2], idx[3], nL, sig_mode, s1, s2, s3)
	// res: (crr, Abar, lambda, erho2, sigmae, ntrunc, ub,
	//       s_aP, s_jP, s_pi, s_aL, s_jL, s_piL)  -- see _crr_core

	// Push results to Stata
	comps = (res[8], res[9], res[10] \ res[11], res[12], res[13])
	st_matrix(comps_name, comps)
	rstripe        = J(2, 2, "")
	rstripe[., 2]  = ("P" \ "L_P")
	cstripe        = J(3, 2, "")
	cstripe[., 2]  = (avn \ jvn \ pvn)
	st_matrixrowstripe(comps_name, rstripe)
	st_matrixcolstripe(comps_name, cstripe)

	st_numscalar(scn[1], res[1])          // crr
	st_numscalar(scn[2], res[2])          // Abar
	st_numscalar(scn[3], res[3])          // lambda
	st_numscalar(scn[4], res[4])          // erho2_dcfpl
	st_numscalar(scn[5], means[idx[2]])   // Jbar
	st_numscalar(scn[6], means[idx[3]])   // mupi
	st_numscalar(scn[7], res[5])          // sigmae
	st_numscalar(scn[8], res[6])          // ntrunc
	st_numscalar(scn[9], res[7])          // ub
}

// Person-weighted outcome means: mean over groups of within-group means,
// where groups are the levels of zcol (the object-of-measurement column of
// Z_data).  Equals mean(Y) exactly when every group has the same number of
// rows (balanced designs); differs under unbalance, where the grand mean
// overweights objects contributing more rows.
real rowvector mvgstudy::_pw_means(real matrix Y, real colvector zcol)
{
	real colvector uniq, sel
	real matrix    acc
	real scalar    u

	uniq = uniqrows(zcol)
	acc  = J(rows(uniq), cols(Y), .)
	for (u = 1; u <= rows(uniq); u++) {
		sel = selectindex(zcol :== uniq[u])
		acc[u, .] = mean(Y[sel, .])
	}
	return(mean(acc))
}

// Column of Z_data holding the object of measurement, defined as the first
// facet of effects[1] (mvcrr validates object() == effects[1] before calling).
real scalar mvgstudy::_crr_obj_col()
{
	string scalar fac
	real scalar   fc

	fac = tokens(subinstr(subinstr(effects[1], "|", " "), "#", " "))[1]
	for (fc = 1; fc <= length(facets); fc++) {
		if (facets[fc] == fac) return(fc)
	}
	errprintf("mvcrr: object facet '%s' not found among design facets\n", fac)
	exit(198)
}

// Resolve the (A, J, prevalence) variable positions in this.varlist.
// Returns (ai, ji, pvi); errors out if any name is not found.
real rowvector mvgstudy::_crr_resolve_idx(string scalar avn,
                                           string scalar jvn,
                                           string scalar pvn)
{
	real scalar i, ai, ji, pvi

	ai  = 0
	ji  = 0
	pvi = 0
	for (i = 1; i <= length(varlist); i++) {
		if (varlist[i] == avn) ai  = i
		if (varlist[i] == jvn) ji  = i
		if (varlist[i] == pvn) pvi = i
	}
	if (ai == 0 | ji == 0 | pvi == 0) {
		errprintf("mvcrr: could not locate A/J/prevalence variables in the G-study varlist\n")
		exit(198)
	}
	return((ai, ji, pvi))
}

// Shared CRR computation core, used by compute_crr (point estimate) and
// run_crr_bootstrap (per-rep values) so the two can never diverge.
// Inputs: M_P/M_L = person and lesson-within-person component matrices,
// means = row vector of outcome means, ai/ji/pvi = variable positions,
// nL, sig_mode/s1/s2/s3 as documented at compute_crr.
// Returns (crr, Abar, lambda, erho2_dcfpl, sigmae, ntrunc, ub,
//          s_aP, s_jP, s_pi, s_aL, s_jL, s_piL) with components
// post-truncation.
real rowvector mvgstudy::_crr_core(real matrix M_P, real matrix M_L,
                                    real rowvector means,
                                    real scalar ai, real scalar ji,
                                    real scalar pvi,
                                    real scalar nL, real scalar sig_mode,
                                    real scalar s1, real scalar s2,
                                    real scalar s3)
{
	real scalar Abar, Jbar, mupi, s_aP, s_jP, s_pi, s_aL, s_jL, s_piL
	real scalar sigmae, ntrunc, ub, tr_person
	real scalar num_pi, num_tau, err_L, denom, lam, erho2, crr

	Abar  = means[ai]
	Jbar  = means[ji]
	mupi  = means[pvi]

	s_aP  = M_P[ai,  ai]
	s_jP  = M_P[ji,  ji]
	s_pi  = M_P[pvi, pvi]
	s_aL  = M_L[ai,  ai]
	s_jL  = M_L[ji,  ji]
	s_piL = M_L[pvi, pvi]

	// Truncate negative estimates at zero (paper's convention); count them
	ntrunc    = 0
	tr_person = 0
	if (s_aP < 0) {
		s_aP = 0
		ntrunc++
		tr_person++
	}
	if (s_jP < 0) {
		s_jP = 0
		ntrunc++
		tr_person++
	}
	if (s_pi < 0) {
		s_pi = 0
		ntrunc++
	}
	if (s_aL < 0) {
		s_aL = 0
		ntrunc++
	}
	if (s_jL < 0) {
		s_jL = 0
		ntrunc++
	}
	if (s_piL < 0) {
		s_piL = 0
		ntrunc++
	}
	// Upper-bound flag: truncation removed all person-level DCF (lambda = 1)
	ub = (tr_person > 0 & s_aP == 0 & s_jP == 0)

	// Resolve sigma2_eps
	if      (sig_mode == 1) sigmae = s1
	else if (sig_mode == 2) sigmae = (mupi*s2 + (1-mupi)*s1
	                                  + mupi*(1-mupi)*Jbar^2) / s3
	else                    sigmae = 0

	// CRR (chapter eq. line 387, extended to the PL coefficient of line 267)
	num_pi  = Jbar^2 * s_pi
	num_tau = num_pi + s_aP + s_jP*(mupi^2 + s_pi)
	err_L   = s_aL + (Jbar^2 + s_jP)*s_piL + s_jL*(mupi^2 + s_pi + s_piL) + sigmae
	denom   = num_tau + err_L/nL

	lam   = (num_tau > 0 ? num_pi  / num_tau : .)
	erho2 = (denom   > 0 ? num_tau / denom   : .)
	crr   = (denom   > 0 ? num_pi  / denom   : .)

	return((crr, Abar, lam, erho2, sigmae, ntrunc, ub,
	        s_aP, s_jP, s_pi, s_aL, s_jL, s_piL))
}

// Phase 12b: Bootstrap BCa confidence intervals for CRR.
// Propagates the G-study bootstrap reps (boot_store/jack_store, plus the
// per-rep means stored in boot_means/jack_means) through the CRR core and
// summarizes with BCa intervals — the same pattern as run_dstudy_bootstrap.
// Requires a prior mvgstudy run with the bootstrap option.
// Pushes a 4 x 4 summary matrix to table_name:
//   rows crr, Abar, lambda, erho2_dcfpl; cols estimate, se, ci_lo, ci_hi.
void mvgstudy::run_crr_bootstrap(string scalar avn, string scalar jvn,
                                  string scalar pvn,
                                  real scalar nL, real scalar sig_mode,
                                  real scalar s1, real scalar s2,
                                  real scalar s3,
                                  real scalar pw,
                                  real scalar ci_alpha,
                                  string scalar table_name)
{
	real scalar      k_out, b, j_p, q, n_jack, se, a
	real rowvector   idx, th_full, res, ci, means_full, mrow
	real matrix      boot1, boot2, jack1, jack2, bs, jk, summ
	real matrix      bmeans, jmeans
	real colvector   bvals, jvals
	string matrix    rstripe, cstripe

	if (missing(boot_B) | boot_B <= 0) {
		errprintf("mvcrr: bootstrap requires mvgstudy to have been run with the bootstrap option first\n")
		exit(198)
	}
	if (rows(boot_means) != boot_B | (pw & rows(boot_pwmeans) != boot_B)) {
		errprintf("mvcrr: per-rep means not found; rerun mvgstudy with the bootstrap option (this version)\n")
		exit(198)
	}

	k_out = cols(Y_data)
	idx   = _crr_resolve_idx(avn, jvn, pvn)

	// Means matching the pw setting: per-rep stores and full-sample values
	if (pw) {
		bmeans     = boot_pwmeans
		jmeans     = jack_pwmeans
		means_full = _pw_means(Y_data, Z_data[., _crr_obj_col()])
	}
	else {
		bmeans     = boot_means
		jmeans     = jack_means
		means_full = mean(Y_data)
	}

	// Full-sample point estimates (theta-hat for BCa)
	res     = _crr_core(asarray(covcomps, "emcp1"), asarray(covcomps, "emcp2"),
	                    means_full, idx[1], idx[2], idx[3],
	                    nL, sig_mode, s1, s2, s3)
	th_full = res[1..4]   // crr, Abar, lambda, erho2_dcfpl

	boot1  = asarray(boot_store, "emcp1")
	boot2  = asarray(boot_store, "emcp2")
	jack1  = asarray(jack_store, "emcp1")
	jack2  = asarray(jack_store, "emcp2")
	n_jack = rows(jack1)

	// Bootstrap replicates
	bs = J(boot_B, 4, .)
	for (b = 1; b <= boot_B; b++) {
		mrow = bmeans[b, .]
		if (missing(boot1[b, 1]) | missing(mrow[1])) continue
		res = _crr_core(rowshape(boot1[b, .], k_out),
		                rowshape(boot2[b, .], k_out),
		                mrow, idx[1], idx[2], idx[3],
		                nL, sig_mode, s1, s2, s3)
		bs[b, .] = res[1..4]
	}

	// Jackknife replicates (for BCa acceleration)
	jk = J(n_jack, 4, .)
	for (j_p = 1; j_p <= n_jack; j_p++) {
		mrow = jmeans[j_p, .]
		if (missing(jack1[j_p, 1]) | missing(mrow[1])) continue
		res = _crr_core(rowshape(jack1[j_p, .], k_out),
		                rowshape(jack2[j_p, .], k_out),
		                mrow, idx[1], idx[2], idx[3],
		                nL, sig_mode, s1, s2, s3)
		jk[j_p, .] = res[1..4]
	}

	// BCa summary per quantity
	summ = J(4, 4, .)
	for (q = 1; q <= 4; q++) {
		bvals = select(bs[., q], !missing(bs[., q]))
		jvals = select(jk[., q], !missing(jk[., q]))
		se = (rows(bvals) > 1 ? sqrt(variance(bvals)) : .)
		a  = (rows(jvals) >= 2 ? bca_acceleration(jvals) : 0)
		if (rows(bvals) > 0) {
			ci = bca_ci(bvals, th_full[q], a, ci_alpha)
		}
		else {
			ci = (., .)
		}
		summ[q, .] = (th_full[q], se, ci)
	}

	st_matrix(table_name, summ)
	rstripe       = J(4, 2, "")
	rstripe[., 2] = ("crr" \ "Abar" \ "lambda" \ "erho2_dcfpl")
	cstripe       = J(4, 2, "")
	cstripe[., 2] = ("estimate" \ "se" \ "ci_lo" \ "ci_hi")
	st_matrixrowstripe(table_name, rstripe)
	st_matrixcolstripe(table_name, cstripe)
}

// Phase 12d: Single-replication mode.
// For a one-effect G study (A J prevalence = object) — one row per object,
// no lesson replication (e.g., a validation subsample with pooled utterance
// sets per teacher).  Reproduces the paper's person-only CRR
// (lambda x Erho2_DCF:P, chapter eq. line 387):
//   - emcp1 is the across-object sample covariance matrix; its diagonals are
//     DISATTENUATED by subtracting mean per-object sampling variances,
//       SV(A_p)  = nu0^2 * mean(1/n0_p)
//       SV(J_p)  = nu1^2 * mean(1/n1_p) + nu0^2 * mean(1/n0_p)
//       SV(pi_p) = mean(pi_p*(1-pi_p)/(n0_p+n1_p))            (binomial)
//     with n0/n1 = per-object gold-negative/-positive utterance counts read
//     from the dataset in memory (must match the estimation sample rows).
//   - the lesson side enters as D-study parameters:
//     sigma2_pi_LP = rho_lp * max(0, sigma2_pi_disattenuated),
//     sigma2_a_LP = sigma2_j_LP = 0, and sigma2_eps from the nu closed form.
//   - _crr_core does the truncation, flags, and formula; with the lesson-level
//     DCF components zero it reduces exactly to eq. 387.
// Outputs the same 9 scalars and comps matrix as compute_crr, plus the 1x3
// mean sampling variances (sv_name) used for disattenuation.
void mvgstudy::compute_crr_sr(string scalar avn, string scalar jvn,
                               string scalar pvn,
                               real scalar nL,
                               real scalar nu0, real scalar nu1,
                               real scalar nutt, real scalar rho_lp,
                               string scalar n0name, string scalar n1name,
                               string scalar comps_name, string scalar sv_name,
                               string scalar scnames_str)
{
	real rowvector   idx, means, res
	real colvector   n0, n1, piv
	real matrix      M_P, M_L, comps
	real scalar      k, vA, vJ, vP, msvA, msvJ, msvP, d_aP, d_jP, d_pi, s_piL
	string rowvector scn
	string matrix    rstripe, cstripe, svrstripe

	scn   = tokens(scnames_str)
	idx   = _crr_resolve_idx(avn, jvn, pvn)
	means = mean(Y_data)
	k     = cols(Y_data)

	// Per-object counts from the dataset in memory
	if (st_nobs() != rows(Y_data)) {
		errprintf("mvcrr: observations in memory (%g) differ from the mvgstudy estimation sample (%g); re-run mvgstudy on the current data\n",
		          st_nobs(), rows(Y_data))
		exit(198)
	}
	n0 = st_data(., n0name)
	n1 = st_data(., n1name)
	if (missing(n0) > 0 | missing(n1) > 0 | min(n0) <= 0 | min(n1) <= 0) {
		errprintf("mvcrr: n0var()/n1var() must be positive and nonmissing for every object\n")
		exit(198)
	}

	// Raw across-object variances (single-effect emcp1 = sample covariance)
	M_P = asarray(covcomps, "emcp1")
	vA  = M_P[idx[1], idx[1]]
	vJ  = M_P[idx[2], idx[2]]
	vP  = M_P[idx[3], idx[3]]

	// Mean per-object sampling variances (disattenuation corrections)
	piv  = Y_data[., idx[3]]
	msvA = nu0 * mean(1 :/ n0)
	msvJ = nu1 * mean(1 :/ n1) + nu0 * mean(1 :/ n0)
	msvP = mean(piv :* (1 :- piv) :/ (n0 + n1))

	d_aP = vA - msvA
	d_jP = vJ - msvJ
	d_pi = vP - msvP
	// Lesson-side D-study parameter uses the truncated prevalence variance
	s_piL = rho_lp * max((0, d_pi))

	// Build diagonal component matrices; _crr_core truncates negatives (the
	// paper's convention for disattenuated components), sets the lambda = 1
	// upper-bound flag, resolves sigma2_eps by the nu closed form (mode 2),
	// and evaluates the formula.
	M_P = J(k, k, 0)
	M_P[idx[1], idx[1]] = d_aP
	M_P[idx[2], idx[2]] = d_jP
	M_P[idx[3], idx[3]] = d_pi
	M_L = J(k, k, 0)
	M_L[idx[3], idx[3]] = s_piL

	res = _crr_core(M_P, M_L, means, idx[1], idx[2], idx[3],
	                nL, 2, nu0, nu1, nutt)

	// Push results (same contract as compute_crr)
	comps = (res[8], res[9], res[10] \ res[11], res[12], res[13])
	st_matrix(comps_name, comps)
	rstripe        = J(2, 2, "")
	rstripe[., 2]  = ("P" \ "L_P")
	cstripe        = J(3, 2, "")
	cstripe[., 2]  = (avn \ jvn \ pvn)
	st_matrixrowstripe(comps_name, rstripe)
	st_matrixcolstripe(comps_name, cstripe)

	st_matrix(sv_name, (msvA, msvJ, msvP))
	svrstripe        = J(1, 2, "")
	svrstripe[1, 2]  = "sampvar"
	st_matrixrowstripe(sv_name, svrstripe)
	st_matrixcolstripe(sv_name, cstripe)

	st_numscalar(scn[1], res[1])          // crr
	st_numscalar(scn[2], res[2])          // Abar
	st_numscalar(scn[3], res[3])          // lambda
	st_numscalar(scn[4], res[4])          // erho2_dcfp
	st_numscalar(scn[5], means[idx[2]])   // Jbar
	st_numscalar(scn[6], means[idx[3]])   // mupi
	st_numscalar(scn[7], res[5])          // sigmae
	st_numscalar(scn[8], res[6])          // ntrunc
	st_numscalar(scn[9], res[7])          // ub
}

end

//--
