---
name: stata-mata-programmer
description: "Use this skill when asked to write, review, debug, extend, or refactor Stata .ado programs, the Mata code embedded in them (free functions, classes, Stata–Mata hybrids), or the .sthlp help files that accompany them. Triggers include: 'write an ado command for X', 'add an option to this command', 'write a Mata function / class method', 'debug this Mata class', 'write / update the help file', 'add a post-estimation command', 'review this Stata code', or any task touching files ending in .ado, .mata, .sthlp, or .do that contain `mata` blocks."
---

# Stata + Mata Programmer Skill

## Purpose

Produce production-quality Stata packages: an `.ado` wrapper that parses, validates, displays, and returns; a Mata engine with explicit typing, public/private class structure, `asarray` memoization, and clean boundary management; and a SMCL `.sthlp` help file that matches the command's syntax line-for-line. Every convention here survives Stata's batch-mode parser (`stata -b do`).

---

## Part A — Workflow and `.ado` conventions

### Phase 0 — Understand before writing

Before changing anything, answer and report to the user:

1. **New command, new option, or extension of a Mata class?** If extending a class, read the full class declaration first — every field, every method, the public/private split. If adding a Stata option, read the `syntax` line and every validation block that follows it.
2. **What crosses the Stata–Mata boundary?** List every `st_local()`, `st_matrix()`, `st_numscalar()`, `st_data()`, `st_matrixrowstripe()` call needed, and its direction.
3. **Does the command depend on persistent state?** Post-estimation commands typically read a Mata object left behind by the estimation command. Decide whether your command reads it, mutates it (then it needs a restorer), or creates it.
4. **What does the caller get back?** `r()`/`e()` scalars, matrices, macros, displayed tables, side effects on class state. Write these down; they become the help file's "Stored results" section.
5. **Is memoization needed?** If the same subcomputation runs for different inputs in a loop, plan an `asarray` cache before writing the loop.
6. **Which files change?** A new option touches three things: the `.ado` (syntax, validation, display, return), the Mata method, and the `.sthlp` (syntax line, option table, option description, stored results, example). Do all three or say explicitly which you skipped.

### Phase 1 — Package layout: one command per file + a Mata library

**Never ship several user commands in one `.ado`.** Stata autoloads a file only by the name of the command called, and it *auto-drops* autoloaded programs under cache pressure (e.g. after a long bootstrap). A second command living in the first command's file is therefore unrecognized in a fresh session and can vanish mid-script even after it once worked. "It works if you call the main command first" is not a fix.

```
mycmd.ado              program mycmd            (autoloads on `mycmd`)
mypost1.ado            program mypost1          (autoloads on `mypost1`)
_mycmd_parse.ado       program _mycmd_parse     (private helper; package-prefixed name)
lmycmd.mlib            compiled Mata: class + free functions   (shipped)
mycmd.mata             Mata source, `mata … end` block         (repo only, NOT in .pkg)
build_mlib.do          compiles mycmd.mata -> lmycmd.mlib      (developer only)
mycmd.pkg, stata.toc, *.sthlp, example .dta
```

Each `.ado`:

```
//----------------------------------------------------------------------------//
*! mycmd  version 1.4.0  24aug2026
*! v1.4.0: one-line note on this release.  Requires Stata 19.  History: CHANGELOG.md.
//----------------------------------------------------------------------------//
program mycmd, rclass
	version 19
	...
end
```

`build_mlib.do`:

```stata
version 19
clear all
mata: mata clear                        // makes the *() add below exact
do mycmd.mata                           // compiles class + free functions into memory
mata: mata mlib create lmycmd, dir(.) replace
mata: mata mlib add lmycmd *()          // everything in memory == the file's contents
mata: mata mlib index
mata: mata describe using lmycmd        // eyeball: class + N free functions, nothing foreign
```

Rules:
- The `*!` banner is what `which mycmd` reports. Bump it (version **and** date) on every shipped change, in every `.ado`, the `.mata` banner, and every `.sthlp` header (`{* *! version 1.4.0  24aug2026}{...}`). Keep the banner to two lines; release history goes in `CHANGELOG.md`, not repeated in every file.
- An `.mlib` is a binary: rebuild it whenever the `.mata` changes, and run `mata mlib index` (or reinstall) so a stale library does not shadow new code. `f lmycmd.mlib` in the `.pkg` is easy to forget and fatal to omit.
- Persistent state is unaffected by the library: the *class definition* lives in the `.mlib`; the *instance* (`c`) still lives in Mata memory for the session, exactly as with inline Mata.
- Private Stata helpers get a package prefix (`_mycmd_parse`), never a generic name (`parse_equation`) that another package may also define.
- When shipping a rewrite, keep the previous release byte-identical under a `_old` suffix and never edit it.
- Indentation: **match the file**. Do not mix tabs and spaces within a file; use the same scheme in the Stata and Mata portions.
- In `mycmd.mata`, one `mata … end` block. Never split the class definition and its method bodies across separate blocks.

**Version floor.** The `.pkg` `d Requires: Stata version N` line must equal the `version N` statement inside the programs — and once an `.mlib` ships, N is the Stata that compiled it (older Stata cannot read it). State the same N in README and in each help file's Description. A `.pkg` that says 16 while the programs say 19 is a defect.

**No undeclared dependencies.** Before building, `grep -n "mm_" *.mata *.ado` (and any other SSC Mata prefix you know you have installed). A single `mm_cond()` hiding in a helper makes the package fail at runtime for every user without moremata, and your own machine will never show it.

### Phase 2 — Anatomy of a Stata program

Follow this order inside every `program`:

```stata
program mypost1, rclass
	version 19
	syntax, Object(string) ///
	        [NL(integer 0) Avar(string) Weights Bootstrap CI_level(real 95)]

	tempname _nv R_est R_se RT            // every Stata-side scratch name up front

	//---//
	// Validate that prior results exist
	//---//
	capture mata: st_numscalar("`_nv'", length(c.varlist))
	if _rc {
		di as error "mypost1: no mycmd results in memory; run mycmd first"
		exit 301
	}

	//---//
	// Validate options (cheap checks before any heavy compute)
	//---//
	if `nl' < 1 {
		di as error "mypost1: nl() is required and must be a positive integer"
		exit 198
	}
	if ("`weights'" != "" & "`avar'" != "") {
		di as error "mypost1: weights may not be combined with avar()"
		exit 198
	}

	//---//
	// Compute (one or two Mata calls; pass tempnames for the results)
	//---//
	mata c.compute_stat("`avar'", `nl', "`RT'", "`R_est' `R_se'")

	//---//
	// Display
	//---//
	di as text _newline "Projected statistic"
	di as text "{hline 72}"
	di as text "  estimate" _column(39) "= " as result %9.4f scalar(`R_est')
	di as text "  std. err." _column(39) "= " as result %9.4f scalar(`R_se')
	matlist `RT', twidth(8) format(%12.6f) title("Components as used")
	if scalar(`R_tr') > 0 {
		di as text "(note: `=scalar(`R_tr')' negative component(s) truncated at zero)"
	}

	//---//
	// Returned results (last — return matrix clears the local copy)
	//---//
	return scalar est   = scalar(`R_est')
	return scalar nl    = `nl'
	return matrix table = `RT'
end
```

**Syntax line**
- `version 19` (or the minimum you test against) is the first line of every program.
- Capitalize the minimum abbreviation in `syntax` (`Object`, `REPs`, `CI_level`, `BOOTunit`) and use the **same** abbreviation in the help file's `{opt o:bject()}`.
- Typed options (`integer`, `real`, `numlist integer >=1`, `varname`) get type-checking for free; reserve `string` for options you parse yourself (e.g., alternating `name #` pairs).
- Anything that changes defaults or runtime cost is opt-in (`Bootstrap`) — existing output must be byte-identical when the option is omitted.
- `marksample touse` immediately after `syntax` when the command reads data; pass `"`touse'"` to every `st_data()`.

**Validation**
- Every error message starts with the command name and a colon: `"mypost1: nl() must be a positive integer"`.
- Exit codes: `198` invalid syntax/option, `301` required prior results not in memory, `459` object not found, `111` variable not found (let `confirm` raise it).
- Check mutual exclusivity and "must be given together" explicitly, naming both options.
- Use `confirm integer number`, `confirm number`, `confirm numeric variable`, `capture confirm matrix` rather than hand-rolled regexes.
- Validate everything cheap before the first expensive Mata call. Checks that only matter for an opt-in path may live inside that path's block.
- Name validation against a list uses `:list x in y` on whole tokens, never `strpos`.

**Display**
- `di as text` for labels and notes, `di as result` for numbers, `di as error` only for fatal errors.
- Align `=` with `_column(N)`; `%9.4f` for coefficients, `%12.6f` for variance-scale quantities, `%9.0g` for counts.
- Frame a results panel with `{hline 72}`.
- Non-fatal conditions are `(note: ...)` lines in `as text`, lower-case, placed after the table they qualify. Continue long notes with `///`.
- Tables: `matlist M, twidth(20) title("...")`; add `format()` and `names(c)` where defaults are noisy.
- Echo every resolved input that affects the result (source of a default, which mode was chosen) so the output is reproducible from what is printed.

**Returns**
- `return matrix` **moves** the matrix: anything you need from it (sign checks, further display) must happen before that line. Comment the first such site.
- Stata matrix names cannot contain `#` or `|`; sanitize derived names: `local nm_c = subinstr(subinstr("`nm'", "#", "_", .), "|", "_", .)`.
- Return the resolved form of every option that affects the result (`r(nl)`, `r(weights)`), not just the outputs.
- Mode-conditional returns get their own `if` block and a matching line in the help file saying which mode sets them.

### Phase 3 — Stata–Mata boundary

**Reading from Stata**

```stata
mata Y = st_data(., tokens(st_local("varlist")), "`touse'")
mata c = myclass()
mata c.init_inputs(Y, tokens(st_local("groups")), st_matrix("`levels'"))
mata mata drop Y            // drop intermediates immediately after consumption
```

**Handing results back — pass tempnames into Mata; never hard-code global names**

```stata
tempname TAB R_est R_se R_tr
mata c.compute_stat(..., "`TAB'", "`R_est' `R_se' `R_tr'")
return scalar est = scalar(`R_est')
```

```mata
void myclass::compute_stat(..., string scalar tab_name, string scalar scnames_str)
{
	string rowvector scn
	scn = tokens(scnames_str)
	st_matrix(tab_name, tab)
	st_numscalar(scn[1], res[1])          // estimate
	st_numscalar(scn[2], res[2])          // std. err.
	st_numscalar(scn[3], res[3])          // n truncated
}
```

A space-delimited string of tempnames is the cleanest way to pass N result slots through one argument. Comment each `st_numscalar` line with the slot's meaning. Writing to fixed names (`_tmp_1`) leaks scalars into the user's session; don't.

**Label returned matrices in Mata** so `matlist` and users see names without the wrapper re-labeling:

```mata
rstripe        = J(2, 2, "")
rstripe[., 2]  = ("level1" \ "level2")
cstripe        = J(4, 2, "")
cstripe[., 2]  = ("estimate" \ "se" \ "ci_lo" \ "ci_hi")
st_matrixrowstripe(tab_name, rstripe)
st_matrixcolstripe(tab_name, cstripe)
```

Stripes are `n × 2` string matrices; column 1 is the equation name (leave `""`), column 2 the row/column name.

**Iterating an asarray back to r()**

```stata
mata for (loc=asarray_first(c.results); loc!=NULL; loc=asarray_next(c.results, loc)) ///
	st_matrix(asarray_key(c.results, loc), asarray_contents(c.results, loc))
```

**Persistent state**

```stata
// c persists after this program exits so that mypost1 / mypost2 can access it.
// Run "mata drop c" manually when all post-estimation calls are complete.
mata c = myclass()
```

Document persistence at the creation site, in every consumer's header comment, and in every help file (`{stata mata drop c}`). Consumers detect absence with a `capture mata:` probe + `exit 301`, never by letting Mata's `r(3499)` surface.

**Verify the dataset in memory still matches stored state** when a post-estimation command reads new variables from it:

```mata
if (st_nobs() != rows(Y_data)) {
	errprintf("mypost1: observations in memory (%g) differ from the estimation sample (%g)\n",
	          st_nobs(), rows(Y_data))
	exit(198)
}
```

### Phase 4 — State save/restore for mutating operations

When a command mutates class state transiently, save originals first and wrap the Stata body in `capture noisily` so the restorer runs on error:

```stata
local n_mutated 0
if "`fix'" != "" {
	local n_mutated 1
	mata c.apply_transform(tokens("`fix'"))
}

capture noisily {
	mata c.main_routine()
	mata c.export_results()
	... display and return ...
}
local _body_rc = _rc
if `_body_rc' != 0 {
	if `n_mutated' > 0 capture mata c.restore()
	exit `_body_rc'
}
if `n_mutated' > 0 mata c.restore()
```

Each additional opt-in block inside the same program gets its **own** `capture noisily` + restore + `exit` triple. The Mata restorer assigns every saved field back and is idempotent:

```mata
void myclass::restore()
{
	if (orig_state == J(0, 0, .)) return     // nothing to restore
	state      = orig_state
	orig_state = J(0, 0, .)
}
```

Inside Mata the same discipline applies when a method temporarily overwrites class state for a loop (swapping components per replicate): copy originals into a local, loop, restore, then re-run whatever derives downstream state so the object is left exactly as found.

---

## Part B — Mata conventions

### Phase 5 — Class organization

Public section first, private second. **Always.**

```mata
class myclass
{
	//----------//
	public:
	//----------//

	// Estimation routines
	void init_inputs()
	void main_routine()

	// Exporters that push state back to Stata
	void export_results()
	void push_table()

	// Phase 3: projected statistic
	void compute_stat()
	void run_stat_bootstrap()

	// Public state — fields the Stata wrapper reads
	string vector   groups
	real matrix     P
	transmorphic    components
	real scalar     solved
	real scalar     quiet_          // 1 = suppress notes (set on temp instances)

	//----------//
	private:
	//----------//

	// Internal helpers
	void             _cache_cellmeans()
	real rowvector   _stat_core()
	string scalar    _canonical_key()

	// Private state — raw data, caches, intermediate buffers
	real matrix      Y_data
	transmorphic     cmcache
	real scalar      is_bal
}
```

Rules:
- Methods the Stata wrapper calls → `public`. Helpers called only from Mata → `private`, prefixed with `_`.
- State the wrapper reads → `public`. Raw data and caches → `private`.
- Group declarations under `// Phase N: feature` comments so the declaration block doubles as a map of where features live.
- A class may instantiate itself for scratch work (`class myclass scalar tmp` inside a bootstrap loop). Give such instances a `quiet_` flag so their notes don't flood the output; set it **before** each `init`, since `init` may reset state.
- Free functions go above the class when they are used both inside methods and from the wrapper. **No backticks or `${}` in comments above free functions** inside an `.ado` file — Stata's preprocessor expands them before Mata sees the source.

### Phase 6 — Shared computation cores

When a quantity is computed once for the point estimate and again per bootstrap/jackknife replicate (or in two modes), there must be exactly **one** function that evaluates the formula:

```mata
// Shared core, used by compute_stat (point estimate) and run_stat_bootstrap
// (per-replicate values) so the two can never diverge.
// Returns (est, lambda, sigma, ntrunc, ub, c1, c2, c3) with components
// post-truncation.
real rowvector myclass::_stat_core(real matrix M1, real matrix M2,
                                   real rowvector means, real scalar n,
                                   real scalar mode, real scalar s1, real scalar s2)
```

- The core takes plain matrices/scalars, not class state, so it can run on replicate data without mutating `this`.
- It returns a **fixed-layout rowvector**; the layout is documented in the header comment and callers index it with commented positions (`res[1] // est`). Alternative modes build their inputs differently and call the same core — never re-derive the formula.
- Enumerated modes are small integers whose encoding is written at the top of the core's comment (`mode: 0 = ..., 1 = ..., 2 = ...`).

### Phase 7 — Type declarations and style

Declare every local at the top of every function with an explicit type. Mata accepts block-scoped declarations inside `{ }`, but do not use them — they hide variables from the reader.

| Need | Type |
|------|------|
| Loop counter / index / flag / count | `real scalar` |
| Numeric matrix | `real matrix` |
| Numeric vector | `real rowvector` / `real colvector` (match usage; `real vector` only when either orientation is accepted) |
| Single string | `string scalar` |
| Token list | `string rowvector` |
| Stripe matrix | `string matrix` |
| asarray / iteration pointer / genuinely varying | `transmorphic` |
| Scratch class instance | `class myclass scalar tmp` |

Style:

```mata
for (i = 1; i <= length(items); i++) {
	...
}
if (sum2 == 0) return(0)                        // one-line guards are fine
lam = (denom > 0 ? num / denom : .)             // ternary for guarded division
se  = (rows(vals) > 1 ? sqrt(variance(vals)) : .)
```

- Spaces around binary operators; spaces after commas and semicolons.
- Align `=` in runs of ≥3 assignments; align the same in `st_numscalar` blocks.
- `(A, B)` horizontal, `(A \ B)` vertical concatenation.
- Pick `this.` or bare member names per method; don't mix. Init methods use `this.` where a parameter shadows the member.
- Opening brace on the signature line; one-liner methods may use inline `return(...)`.
- Header comment on every non-trivial method: what it computes, the input contract, the output layout, and **why** any non-obvious choice was made.

### Phase 8 — Error handling and messages

```mata
if (ai == 0 | ji == 0) {
	errprintf("mypost1: could not locate the named variables in the stored varlist\n")
	exit(198)
}
```

- `errprintf` + `exit(code)`, message prefixed with the **Stata command name** (not the Mata method) so users can find it.
- Non-fatal notes: `printf("{text}Note: ...{smcl}\n")`; warnings: `printf("{err}Warning: ...{smcl}\n")`. Gate every note on `if (!quiet_)`.
- Optional trailing args: `| real vector weights_in` in the signature, `if (args() > 4)` in the body.

### Phase 9 — Asarray patterns

```mata
cache   = asarray_create()                     // string keys (default)
idx_map = asarray_create("real")               // numeric keys
if (asarray_contains(cache, key)) return       // memoize
asarray(cache, key, value)
for (loc = asarray_first(a); loc != NULL; loc = asarray_next(a, loc)) { ... }
```

- Initialize caches in the init method that owns the data they depend on; **re-create them** whenever that data changes (e.g., after dropping rows with missing values).
- Canonicalize keys through one helper (`_canonical_key`) so `"a#b"` and `"b#a"` never become two entries.
- Matrices stored in an asarray are copied on read; to update, read → modify → write back.
- `""` is a legal key in Mata but not in many porting targets; prefer a sentinel like `"__grand__"` if portability matters.

### Phase 10 — Whole-token string matching

`strpos` is correct for **delimiter detection** (`strpos(s, "|") > 0`) and wrong for **token membership** (`strpos("pr#i", "p")` matches). For membership, tokenize and compare:

```mata
real scalar _name_in_term(string scalar term, string scalar name)
{
	string rowvector toks
	real scalar      t
	toks = tokens(subinstr(subinstr(term, "|", " "), "#", " "))
	for (t = 1; t <= length(toks); t++) {
		if (toks[t] == name) return(1)
	}
	return(0)
}
```

Or on a space-delimited string: `sum(tokens(s) :== name) > 0`. Write one helper and reuse it everywhere. A test suite whose fixtures use only single-letter, non-prefix names cannot detect this bug — always add a fixture with prefix-overlapping names (`p`, `pr`).

### Phase 11 — Matrix algebra and replicate bookkeeping

| Function | When |
|----------|------|
| `luinv(M)` | General square, full rank by construction |
| `invsym(M)` | Symmetric positive-definite |
| `lusolve(K, b)` | Solve `K X = b` — prefer over `luinv(K) * b` |
| `rank(K) < rows(K)` | Check before solving; fall back with a gated warning and a status flag (`solved = 0`) |

- Element-wise ops (`:*`, `:/`, `:==`, `:&`) and `selectindex(mask)` over explicit loops.
- **Flatten/unflatten** replicate matrices: store each k×k replicate as one row of a `B × k²` matrix; rebuild with `rowshape(row, k)`.
- **Skipped replicates stay missing.** Leave the row as `.`, summarize with `select(x, !missing(x))`, count skips, and print one note at the end (`%g of %g reps skipped`), escalating to a warning past 50%.
- Guard every ratio whose denominator can be zero with a ternary returning `.`; clamp BCa-style probabilities to `(0.001, 0.999)` and index with `max((1, floor(...)))` / `min((B, ceil(...)))`.
- Empty-input guards before string ops on possibly-0-row matrices: `if (rows(labels) == 0) { out = J(0, 1, ""); return }`.

### Phase 12 — Batch-mode pitfalls

These work interactively and fail under `stata -b do`:

1. Backticks / `${}` in Mata comments inside `.ado` files → preprocessor expansion. Strip them.
2. Multi-line `mata: { … }` at do-file top level or inside `foreach` → `<istmt> incomplete r(3000)`. Use `mata … end` or a one-line `mata: expr`.
3. `if (cond) { a; b; break }` on one line → spread over lines.
4. `return(...)` directly after a multi-line `if/else if` chain → brace each branch.
5. `mata drop X` when `X` may not exist → `capture mata drop X`.
6. `{stata ...}` links in a help file that reference a relative data path fail when clicked from another directory — ship example data beside the help file or use non-clickable `{cmd:. ...}` lines.
7. `mata drop x` inside an `if { }` / `foreach { }` block in a do-file → `invalid expression r(3000)`. Use `mata: mata drop x`.
8. Any `if cond { a; b }` written on one line inside a loop → `matching close brace not found r(198)`. Always open the brace, put statements on their own lines, close on its own line (generalizes pitfall 3).
9. A test harness that `run`s a package file repeatedly must `capture program drop` **every** program that file (and any split-out `.ado`) defines, not just the main one — a missed `program drop mypost2` surfaces as `program mypost2 already defined r(110)` only once a same-named `.ado` also exists on the adopath.

---

## Part C — Help files (`.sthlp`)

One help file per user-callable command, named exactly `cmdname.sthlp`, next to the `.ado`. Write it in the same change as the code; a command without a help file is not shipped.

### Phase 13 — Skeleton

```smcl
{smcl}
{* *! version 1.3.0  15aug2026}{...}
{viewerdialog mycmd "dialog _mycmd"}{...}
{viewerjumpto "Syntax" "mycmd##syntax"}{...}
{viewerjumpto "Description" "mycmd##description"}{...}
{viewerjumpto "Options" "mycmd##options"}{...}
{viewerjumpto "Remarks" "mycmd##remarks"}{...}
{viewerjumpto "Examples" "mycmd##examples"}{...}
{viewerjumpto "Stored results" "mycmd##results"}{...}
{viewerjumpto "References" "mycmd##reference"}{...}
{vieweralsosee "mypost1" "help mypost1"}{...}

{p2col:{bf:mycmd}} One-line description of what the command does  {p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:mycmd} {cmd:(}{depvarlist} {cmd:=} {it:termlist}{cmd:)} {ifin}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt o:bject(name)}}object of measurement (required){p_end}
{syntab:Bootstrap}
{synopt:{opt boot:strap}}activate bootstrap SEs and BCa confidence intervals{p_end}
{synopt:{opt reps(#)}}number of replications; default is {cmd:reps(1000)}{p_end}
{synoptline}

{p 8 10 2}
where {it:termlist} is ...

{marker description}{...}
{title:Description}

{pstd}{opt mycmd} estimates ... {p_end}

{pstd}Requires Stata 19 or later. No dependencies outside official Stata (the Mata engine ships compiled in {cmd:lmycmd.mlib}).{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt o:bject(name)} specifies ... Required.

{phang}
{opt fix(name # ...)} treats ... e.g., {cmd:fix(i 5)} or {cmd:fix(i 5 h 3)}.

{pmore}
Second paragraph of the same option, indented under it.

{dlgtab:Bootstrap}
...

{marker remarks}{...}
{title:Remarks: <specific topic>}

{pstd}...{p_end}

{marker examples}{...}
{title:Examples}

{pstd}What this example shows{p_end}

{p 8 10 2}{stata use exampledata.dta, clear}{p_end}
{p 8 10 2}{stata mycmd (y1 y2 = p i p#i)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:mycmd} stores the following in {cmd:r()}:

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Scalars}{p_end}
{synopt:{cmd:r(est)}}the projected statistic{p_end}

{p2col 5 26 30 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}k²×4 summary: rows {it:v_vp}, columns estimate, se, ci_lo, ci_hi{p_end}

{marker reference}{...}
{title:References}

{phang}
Author, A. (Year). {it:Title}. Publisher.
{p_end}

{title:Also see}

{psee}
{helpb mypost1}, {helpb mypost2}
```

### Phase 14 — Rules

**Header**
- Line 1 is exactly `{smcl}`. Line 2 is the version stamp, identical to the `.ado` banner.
- One `{viewerjumpto}` per `{marker}`; one `{vieweralsosee}` per related command. Every directive-only line ends with `{...}` so it emits no blank line.

**Syntax section**
- The syntax line must be derivable from the `syntax` command: same option names, same minimum abbreviations (`{opt o:bject()}` ↔ `Object()`), same argument kinds (`#`, `varname`, `numlist`, `matname`, `name # ...`). Required options go outside the `[ ]`; mutually exclusive alternatives are written `{cmd:(} A {cmd:|} B {cmd:)}` with a sentence below saying exactly one is required.
- `{synoptset N tabbed}` where N is the longest option name + 2; `{syntab:}` groups mirror the `{dlgtab:}` groups in Options.
- Under the table, explain grammar terms (`where termlist is ...`) as a bulleted `{p 8 10 2}` list.

**Description**
- First paragraph: what the command does and the citation. Then one paragraph per mode/branch, each opening with a `{bf:label}`.
- Post-estimation commands state preconditions: what must have been run, with what design, and the `{stata mata drop c}` cleanup.

**Options**
- One `{phang}` per option, in syntax-table order. Start with the option in `{opt}` form, then "specifies ...", give the default as `{cmd:opt(value)}`, say "Required." where true, and end with validation behavior ("must be strictly between 1 and 99", "exits with error code 198 if ...").
- Continuation paragraphs of the same option use `{pmore}`. Mode-restricted options say so in italics up front: `({it:single-sample mode only; required there})`.
- Options that interact (mutually exclusive, must be given together, propagate through another option) name the other option with `{opt}`.

**Remarks**
- Separate `{marker}` + `{title:Remarks: topic}` per topic. Use remarks for methodology the user needs to interpret results (why a default is what it is, how a special case is handled). A "Known limitations" section is mandatory if any exist, and must be **re-read and pruned on every release** — a limitation that was fixed two versions ago but still listed is a documentation bug.
- Small tables with `{col N}` and `{hline N}` lines; keep column starts consistent across rows.

**Examples**
- Each example: a `{pstd}` sentence saying what it shows, then commands. Use `{stata ...}` for clickable one-liners when the example data ships beside the help file; use `{phang}{cmd:. ...}{p_end}` otherwise. Separate examples with `{hline}` when there are more than three.
- Include an example for every opt-in option and every mode, and one showing cleanup for commands using persistent state.

**Stored results**
- Every `return scalar/matrix/local` in the `.ado` has a `{synopt:{cmd:r(name)}}` line, grouped as Scalars, Macros, Matrices. Give dimensions and stripes for matrices. Mode-conditional results say which mode sets them and what they replace.

**Cross-checks before you finish a help file**
1. Every option in `syntax` appears in the syntax line, the synopt table, and Options; every `{opt}` in the help file exists in `syntax`.
2. Every `return` in the program appears in Stored results, and vice versa.
3. Every `{marker x}` has a `{viewerjumpto ... "cmd##x"}`, and every `{help cmd##x}` target exists.
4. Version stamp matches the `.ado`; the "Requires Stata N" sentence matches the programs' `version N` and the `.pkg` `Requires:` line.
5. Open it in Stata (`help cmdname`) and read it once; unbalanced braces render as literal text.

---

## Part D — Testing

1. **Smoke test interactively first.** Tiny synthetic dataset, instantiate the class, call each new method, print results before wiring the wrapper.
2. **Validation suite as a do-file** (`test_cmd.do`) that `run`s the `.ado`, counts failures in a local, prints `PASS:`/`FAIL:` per test with the offending value, and errors out at the end if `fails > 0`. Every test has a closed-form or hand-computed target and an explicit tolerance (`1e-10` exact paths; `1e-8` where LU/QR may pivot differently; looser for Monte-Carlo recovery, stated in the test name).
3. **Internal-consistency tests:** recompute the headline statistic from the returned `r()` pieces and compare to the returned headline — catches drift between display and return paths.
4. **Byte-identical regression** for opt-out paths and for behaviour-preserving refactors. Before touching the tree, snapshot the old package **and** its test fixtures (data directories the tests reference by relative path) to a separate directory and run the baselines from there; a baseline that `r(601)`s on a missing fixture still prints its PASS lines for the tests that did run, so check the baseline logs for `^r([0-9]+);` before trusting them. Diff after normalizing paths and Stata tempfile names (`/T//St#####.######`); the only acceptable residue is timestamps and the version banner.
5. **Edge cases:** empty input, single-row, rank-deficient, degenerate option combinations, prefix-overlapping names, missing values, data in memory changed since estimation.
6. **Error paths:** every `exit 198` branch gets a `capture`d call asserting `_rc == 198`.
7. **Cleanup:** after the full suite, `mata describe` shows only the documented persistent object.
8. **Install the way a user does.** Test from `net install pkg, from("<local folder>") replace` — not `run pkg.ado` — so the `.pkg` (including the `f l<pkg>.mlib` line) is exercised. Then `which` every command and `mata: mata mlib query` must list the library.
9. **Cold-start tests, each in a fresh `stata -b` session:** every post-estimation command called with nothing in memory exits 301 with its message (not `r(199)` unrecognized, not a Mata error); the main command followed immediately by a post-estimation command runs; a long bootstrap followed by a post-estimation command runs (this is where auto-drop bites).
10. **Uninstall/reinstall cycle.** `ado uninstall` removes every file including the `.mlib`; reinstall; rerun a cold-start test. Note `ado uninstall <name>` fails with "criterion matches more than one package" once `stata.trk` holds stacked entries from repeated installs — and `capture ado uninstall` hides that silently. Use `ado dir, find(<name>)` and `ado uninstall [#]` by index, highest first.
11. **Post-publish smoke test** from a clean PLUS: `sysdir set PLUS "<tempdir>"`, `net install` from the real URL, run a cold-start test. It is the only test proving the published route serves the binary `.mlib` correctly.

---

## Common pitfalls (quick index)

| Pitfall | Fix |
|---------|-----|
| `strpos` for token membership | Tokenize + whole-token compare (Phase 10); `strpos` only for delimiter presence |
| Formula duplicated in point-estimate and bootstrap code | Single `_core` function returning a fixed-layout rowvector (Phase 6) |
| Results written to fixed global scalar names | Pass `tempname`s into Mata as a space-delimited string (Phase 3) |
| Unlabeled returned matrix | `st_matrixrowstripe` / `st_matrixcolstripe` in Mata (Phase 3) |
| Using a matrix after `return matrix` | Do all checks/display first; comment the reason |
| `#` / `|` in a Stata matrix name | Sanitize with `subinstr` to `_` |
| Post-estimation command run without prior results | `capture mata:` probe + `exit 301` |
| State leak after error | `capture noisily` + paired idempotent restorer per block (Phase 4) |
| Temp class instance floods output | `quiet_` flag set before each `init` |
| Skipped bootstrap reps stored as zeros | Leave rows missing; `select(x, !missing(x))` |
| Cache stale after data filtered | Re-create caches when data changes (Phase 9) |
| Backticks in Mata comments in `.ado` | Strip (Phase 12) |
| Mixed tabs/spaces | Match the file |
| Help-file abbreviation ≠ `syntax` abbreviation | Derive `{opt x:yz()}` from `Xyz()` |
| Help file lists a limitation already fixed | Prune "Known limitations" every release |
| `.ado` and `.sthlp` version stamps differ | Bump both in the same edit |
| Subcommand not found from cold adopath, or "unrecognized" after a long bootstrap | One command per `.ado`; Mata in an `.mlib` (Phase 1) — never a `run` requirement |
| `.pkg` `Requires:` disagrees with `version N` | Set both to the compiling Stata; state it in README and help |
| `mm_*` or other SSC Mata call slips into the package | grep before building; inline the one-liner |
| `mata drop x` inside a do-file `if {}` block | `mata: mata drop x` |
| `ado uninstall name` ambiguous / silently no-op under `capture` | Uninstall by index (Part D §10) |
| Baseline log missing a fixture but still shows PASS lines | Check baselines for `r(#)` before diffing (Part D §4) |
| `luinv(K) * b` | `lusolve(K, b)` |
| `transmorphic` as a shortcut for a known type | Declare the specific type (Phase 7) |
