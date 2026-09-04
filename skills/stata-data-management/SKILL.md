---
name: stata-data-management
description: "Use this skill when asked to write, review, debug, or refactor Stata .do files that import, clean, reshape, merge, append, deduplicate, recode, label, or assemble analytic datasets — everything between raw files and the first estimation command. Triggers include: 'clean this CSV', 'build the analytic file', 'merge these datasets', 'write a cleaning pipeline', 'review my data prep do-file', 'fix this reshape', or any task touching frames, reshape, merge, append, egen, encode, recode, value labels, or dataset saving. Does NOT cover .ado / Mata programming — use stata-mata-programmer for that."
---

# Stata Data Management Skill

## Purpose

Produce cleaning and dataset-assembly do-files that are **reproducible** (run from a clean session on any machine), **self-checking** (every assumption about keys, match rates, and codings is asserted, not eyeballed), and **readable** (a reader can trace every variable in the output back to its source and rule). Every convention here survives `stata -b do`.

Companion skill: `stata-mata-programmer` covers `.ado`/Mata/`.sthlp`. Share its indentation rule (match the file; never mix tabs and spaces) and its voice.

---

## Phase 0 — Understand before writing

Before touching code, answer and report to the user:

1. **Scope: new project, new do-file, or edit to an existing one?** This decides how much structure you impose (see Phase 1 scope rule).
2. **What is the unit of observation of every input and of the output?** Write the key for each (`studentid`, `studentid × year × window`). Every merge and reshape below is checked against these keys.
3. **Which rows are allowed to disappear, and why?** List each filter (`drop if`, `keep if`, unmatched merges, dedup) and its justification. These become `assert`/`count` checks and comments.
4. **What must be true of the output?** Unique key, no missing treatment variable, score ranges, expected N. These become the end-of-file contract (Phase 10).
5. **Which files change?** A new derived variable touches: its `gen` site, its label, the end-of-file contract, and any downstream do-file that lists variables explicitly. Do all or say which you skipped.

---

## Phase 1 — Project structure (scope rule)

**Apply the master structure below only when (a) the user asks for it, or (b) the user is starting a project from scratch.**

When editing or adding to an existing cleaning file, work inside the user's current layout and conventions. Do **not** restructure, rename files, or add `00_config.do` / `00_run_all.do` unprompted. If the existing code would clearly benefit from it (the same path globals repeated in several files, no master run order, value labels redefined per file), say so **once, in one or two sentences, at the end of your response**, and move on.

### Master structure (greenfield or on request)

```
project/
  00_config.do          root path + derived globals, shared value labels, settings
  00_run_all.do         runs every numbered do-file in order; logs to 06_logs/
  01_raw_data/          never written to by code
  02_clean_data/        one .dta per cleaning do-file, plus codebook .txt
  03_cleaning_code/     01_*.do, 02_*.do ... one input family each
  04_analysis_code/
  06_logs/
```

`00_config.do`:

```stata
version 19
clear all
macro drop _all
set varabbrev off
set more off

// Single root; everything else derives from it.  Override by setting
// $PROJECT_ROOT before running, so the code runs on any machine.
if "$PROJECT_ROOT" == "" global PROJECT_ROOT "`c(pwd)'"
global raw    "${PROJECT_ROOT}/01_raw_data"
global clean  "${PROJECT_ROOT}/02_clean_data"
global logs   "${PROJECT_ROOT}/06_logs"

// Shared value labels live in 00_labels.do (see below) -- run it AFTER loading data
do "${PROJECT_ROOT}/00_labels.do"
```

`00_labels.do` — value labels in their **own** file, every one with `, replace`:

```stata
label define yesno     0 "No" 1 "Yes", replace
label define sex       0 "Female" 1 "Male", replace
label define window    1 "Fall" 2 "Winter" 3 "Spring", replace
label define schoolyr  1 "2022-23" 2 "2023-24" 3 "2024-25", replace
```

Why a separate file: **`use`, `import …, clear`, `clear`, and `restore` all discard the value labels in memory.** A label defined at the top of a do-file is gone by the time `encode …, label()` runs after a `use`, and `encode` then silently falls back to alphabetical codes — the exact failure the `label()` option is meant to prevent. So every do-file runs `00_labels.do` *immediately before* its `encode`/`label values` lines, after the last data load, not just at the top. Use `, replace` rather than `capture label drop a b c`: `label drop` aborts at the first name that does not exist and leaves the rest defined.

Every cleaning do-file then begins `do "${PROJECT_ROOT}/00_config.do"` (or `include`), not its own copy of the paths.

---

## Phase 2 — Anatomy of a cleaning do-file

```stata
//----------------------------------------------------------------------------//
// 01_star_panel.do
// Inputs : ${raw}/Historical_Star_{Math,Reading}_*.csv   (6 files)
// Output : ${clean}/star_math_panel.dta, ${clean}/star_read_panel.dta
// Unit   : student x schoolyear x window (one row = last attempt in window)
//----------------------------------------------------------------------------//
version 19
clear all
do "${PROJECT_ROOT}/00_config.do"          // or the file's existing config block
capture log close
log using "${logs}/01_star_panel.log", replace text

//----------------------------------------------------------------------------//
// 1. Import
//----------------------------------------------------------------------------//
...
//----------------------------------------------------------------------------//
// 2. Harmonize and clean
//----------------------------------------------------------------------------//
...
//----------------------------------------------------------------------------//
// 3. Output contract and save
//----------------------------------------------------------------------------//
...
log close
```

Rules:
- Header states inputs, outputs, and the **unit of observation** of the output.
- `version`, `clear all`, logging at the top; `log close` at the bottom. No `set more off` needed in batch mode, harmless otherwise.
- Sections numbered and framed with the same `//---//` bars as the `.ado` skill.
- Comments explain **why** a row is dropped or a value recoded, never restate the command.
- No diagnostic `tab`/`browse`/`list` left in shipped code. Convert each one into an `assert`, a `count`, or delete it.
- No dead code (`gen temp ... drop temp`, repeated `drop` of the same variable).
- Path globals are never hard-coded absolute paths inside a numbered do-file. If the file already has them and you are only editing, leave them — but see the Phase 1 scope rule about suggesting once.

---

## Phase 3 — Import

```stata
import delimited using "${raw}/`f'", varnames(1) case(lower) ///
	bindquote(strict) stringcols(_all) clear
```

- `stringcols(_all)` then `destring` deliberately is safer than letting `import delimited` guess types per file — a column that is numeric in one year's file and has `"K"` in another will otherwise differ in type and break `append`. Destring with explicit handling:

```stata
replace grade = "0" if grade == "K"
destring grade, replace
assert inrange(grade, 0, 12)
```

- Immediately after import, **confirm the columns you depend on exist**, so a changed export fails loudly at the import line, not 80 lines later:

```stata
confirm variable studentidentifier currentgrade unifiedscale activitycompleteddate
```

- Columns that arrive as `v24`, `v30` ...: never rename by position silently. Either import with `varnames(1)` so header text becomes the name, or read the header row and match on its text; at minimum comment the header string next to each positional rename.
- Rename from the **source's structural names**, never from labels or positions. Survey exports already encode structure in the variable name (Qualtrics `Q37#1_5` = block 37, sub-question 1, slot 5 → Stata `q371_5`); parsing question text with regexes or renaming `i1 … i51` by column order is the most fragile thing a cleaning file can do, and it breaks the moment one question is added to one year's form. If names really are uninformative, read the header row once and build the mapping explicitly.
- Hand-typed `input … end` blocks are **data, not code**. Move them to a CSV in a `lookup/` folder, import with `stringcols(_all)`, and assert their format like any raw file. Code that contains data cannot be diffed, validated, or edited by non-programmers.
- Tag the provenance: `gen byte src_file = \`i'` or `gen str src = "`f'"` before appending, and keep it through the pipeline. It is the first thing you need when a merge looks wrong.
- `capture` around `replace`/`destring`/`confirm` hides real errors. Branch explicitly instead:

```stata
capture confirm string variable currentgrade
if !_rc {
	replace currentgrade = "0" if currentgrade == "K"
	destring currentgrade, replace
}
```
(that is, `capture` only on the `confirm`, never on the mutation.)

---

## Phase 4 — Frames and append

Frames are the right tool for holding several inputs at once. Patterns:

```stata
// Import each file into its own frame
forvalues i = 1/6 {
	frame create f`i'
	frame f`i' {
		import delimited using "${raw}/${file`i'}", ...
		gen byte src_file = `i'
	}
}

// Append into a target frame (Stata 18+)
frame create math
frame math: import ... // or frame copy f1 math
forvalues i = 2/3 {
	frame math: frameappend f`i'       // community-contributed; or:
}
```

If `frameappend` (SSC) is not available, or Stata < 18, use tempfiles — but give each its **own name**; reusing one `tempfile` local inside a loop works but is a classic source of confusion when debugging:

```stata
forvalues i = 2/3 {
	tempfile t`i'
	frame f`i': save "`t`i''"
	frame math: append using "`t`i''"
}
```

- After every append: `frame math: assert !missing(src_file)` and `tab src_file` into the log (via `count if src_file == \`i'`) so file-level row counts are recorded.
- Variable types must agree across appended files; `append` promotes silently (string wins), which turns a numeric column into string if one file had `"K"`. Harmonize types per frame **before** appending (Phase 3).
- Do **not** use `preserve`/`restore` inside a `frame X { … }` block or while other frames are in play. `restore` reloads the data, which discards value labels (Phase 1) and can drop rows created in the meantime. Use `frame put` / `frame copy` to carve out a temporary dataset instead.
- Drop scratch frames when done: `frame drop f1` … or `frames reset` at file end — not mid-file, which would also drop the ones you need.

---

## Phase 5 — Keys, duplicates, and "keep the right row"

Every dataset has a declared key. Check it at every stage boundary:

```stata
isid studentidentifier schoolyear window        // errors if not unique
```

`isid` treats `""` in a string key as missing and errors with "should never be missing". When a string component is legitimately blank (a grade suffix that only some rows carry), use `isid …, missok` and say why in a comment.

When the raw data is **not** unique on the intended key, the selection rule must be explicit, ordered, and commented. Use `bysort` with the sort key in parentheses:

```stata
// Rule: within student x year x window keep the LAST attempt by date;
// break remaining ties by highest score.  (Renaissance re-tests within
// a window; the final attempt is the one of record.)
bysort studentidentifier schoolyear window (testdate unifiedscale): ///
	keep if _n == _N
isid studentidentifier schoolyear window
```

Do **not** use:
- `egen dupes = tag(...)` then `drop if dupes == 0` — `tag` marks an *arbitrary* row (first in current sort order), so which row survives depends on sort history. Use `duplicates drop` when rows are truly identical, and `bysort ... (sortkey): keep if _n == _N` when a rule chooses among non-identical rows.
- Two-step `egen max_date` / `egen max_score` / `keep if ==` chains — they are a longer, less transparent way of writing the one `bysort` line above.
- `recode dupes (0=1)(1=0)` style inversions to count duplicates — use `duplicates report` / `duplicates tag, gen()`.

When duplicates are supposed to be byte-identical, prove it rather than assume it:

```stata
duplicates report studentidentifier
duplicates drop                    // drops only fully identical rows
isid studentidentifier             // fails if non-identical dupes remain
```

---

## Phase 6 — Merges

Every merge states what it expects and fails if wrong:

```stata
merge m:1 studentidentifier using "${clean}/state_tests_2425.dta", ///
	keep(master match) gen(_m_state)
// Expected: every state-test student is in the Star panel; a few Star
// students lack a state test (absent / opted out).
assert _m_state != 2                       // already enforced by keep(), kept for the log
count if _m_state == 1
di as text "(note: `r(N)' panel rows without a 2024-25 state test)"
```

Rules:
- Always `gen(_m_<name>)` so multiple merges in one file do not collide and the flag is self-describing. Keep the flag if downstream code needs "has source X"; otherwise drop it **after** the assert.
- Always `keep()` and, where the claim is strong, `assert(match)` / `assert(master match)`. `drop if _merge == 2` after the fact is the same operation with no record of how many rows went where.
- `isid` the key in the **using** file before a `1:1` or `m:1` merge (or `merge` will error late and cryptically). For frames: `frame using_fr: isid key`.
- Merging in "usage"/"participation" files: missing after merge means **not in that file**, not zero. Create the indicator from the merge flag, then fill zeros only where the indicator justifies it:

```stata
gen byte in_tsi_24 = (_m_tsi24 == 3)
foreach v of varlist tsi_attend_24 tsi_engage_24 {
	replace `v' = 0 if in_tsi_24 == 0 & missing(`v')
}
```
  A bare `recode varlist (.=0)` after a merge erases the distinction between "absent from file" and "present with missing value".
- When a merge only adds variables to a subset, `merge` with `update`/`replace` options needs a comment explaining which source wins.
- Prefer `frlink` + `frget` when both sides already live in frames and you need only a few variables; it avoids tempfiles entirely:

```stata
frlink m:1 studentidentifier, frame(demog)
frget gender econdis, from(demog)
```

---

## Phase 7 — Strings, dates, and harmonizing codes

- Trim and normalize once, early: `replace x = strtrim(stritrim(x))`, `ustrlower()` where case varies.
- Chained `subinstr` calls are hard to verify; prefer one regex: `replace schoolyear = ustrregexra(schoolyear, "School Year|\s", "")`.
- Dates: `gen testdate = date(activitycompleteddate, "MDY")`, `format %td`; datetime strings like `2024-09-03T08:15:00Z` → `clock(substr(s,1,19), "YMDhms")` with a comment on the truncation. Always `assert !missing(testdate) if !missing(activitycompleteddate)` so an unexpected format is caught.
- Derive categories from dates with `inrange(month(testdate), 3, 6)` etc.; assert the partition is exhaustive (`assert !missing(window)`), and show the cross-tab against any vendor-supplied equivalent in the log **once** with a comment on which is authoritative and why.
- Map string codes to numeric with an explicit table, not by alphabetic accident:

```stata
gen byte schoolyear = .
replace schoolyear = 1 if schoolyear_str == "2022-2023"
replace schoolyear = 2 if schoolyear_str == "2023-2024"
replace schoolyear = 3 if schoolyear_str == "2024-2025"
assert !missing(schoolyear)
label values schoolyear schoolyr
drop schoolyear_str
```

---

## Phase 8 — Indicators, recodes, and encoding

**Indicators.** Never let missing become 0 or 1 by accident:

```stata
gen byte swd = (swd_str == "Y") if inlist(swd_str, "Y", "N")
assert !missing(swd) | missing(swd_str)
```

- `gen x = (var != "N")` is wrong: blanks and unexpected codes become 1. `gen x = (var == "Y")` is wrong the other way. Test against the **known value set** with `inlist()` and leave everything else missing, then decide explicitly whether blank means "No" (`replace swd = 0 if swd_str == ""` with a comment justifying it).
- Before building indicators, `tab var, missing` once interactively (not in shipped code) to learn the value set, then encode that set in the `inlist()`.

**Encoding.** `encode` assigns codes alphabetically *within the current file*, so `"F"/"M"` gives 1/2 here and something else where a third category appears. Pin the codes:

```stata
do "${PROJECT_ROOT}/00_labels.do"              // AFTER the last use/import/restore (labels are dropped by them)
encode gender_str, gen(gender) label(sex)      // uses existing codes; new values get appended
assert inlist(gender, 0, 1)                    // fail if an unexpected value was appended
drop gender_str
```

- Do not `encode` then `recode (1=0)(2=1)` by inspection — it silently breaks the moment the value set changes.
- Do not `encode` high-cardinality IDs (school, classroom) unless you need a numeric ID; if you do, `encode` from a **sorted, deduplicated** list so the mapping is stable, or use `egen group()` and keep the original string beside it.
- `label variable x "..."` attaches a **variable label**; `label values x lbl` attaches a **value label**. `label variable gender sex` is a bug (it sets the description to the word "sex").
- Multi-level numeric recodes carry their labels in the same statement so code and label cannot drift: `recode frl (2=0 "None") (3=1 "Reduced") (1=2 "Free"), gen(frl_cat)`. For multi-level *string* codes, `label define` the pinned mapping and `encode …, label()` as above; `recode` does not operate on strings.

---

## Phase 9 — Reshape

```stata
isid studentidentifier subject                     // must be unique on i() x j()
reshape wide score classname, i(studentidentifier) j(subject) string
```

- Check `i() × j()` uniqueness with `isid` **before** reshaping; when it fails, apply a documented dedup rule (Phase 5), never a `totdupes` trick.
- Restricting the sample to make a reshape "easier" is a substantive decision — state it as such in the header, not as a trailing comment.
- Name stubs so the result is self-describing: `score_math_2425`, not `statetest_math` renamed later.
- After reshaping, `describe` the new variables into the log and `assert` the expected count of rows (`assert _N == \`n_students'`).
- A two-level grid (e.g. 5 sessions × 26 slots, each with type/recipient/feedback) needs **one** reshape, not two chained ones: encode the index arithmetically, `rename q371_5 otrtype105` (`j = 100*session + slot`), `reshape long otrtype recip fb, i(id) j(j)`, then `gen session = floor(j/100)`, `gen slot = mod(j,100)`. Assert `_N == n_rows * n_cells` before dropping empty cells, and assert the empty cells are empty in *every* stub before dropping them.

---

## Phase 10 — Derived variables, labels, and the output contract

**Derived variables** (treatment, eligibility, composite flags): one block, exhaustive, asserted, labeled at the point of creation.

```stata
// Treatment ladder, 2024-25: 0 none, 1 diagnostic only, 2 + learning path, 3 + TSI.
// TSI students are by construction on a learning path; vendor file sometimes
// codes lp_25 = 0 for them -- correct that first.
replace lp_25 = 1 if tsi_25 == 1
gen byte treat_25 = 0 if diag_25 == 0 & lp_25 == 0 & tsi_25 == 0
replace  treat_25 = 1 if diag_25 == 1 & lp_25 == 0 & tsi_25 == 0
replace  treat_25 = 2 if diag_25 == 1 & lp_25 == 1 & tsi_25 == 0
replace  treat_25 = 3 if diag_25 == 1 & lp_25 == 1 & tsi_25 == 1
assert !missing(treat_25)                     // partition is exhaustive
label define treatlvl 0 "No diagnostic" 1 "Diagnostic" 2 "Learning path" 3 "Learning path + TSI"
label values treat_25 treatlvl
label variable treat_25 "Treatment level, 2024-25"
```

- When the same construction repeats across years or subjects, loop over a local (`foreach y in 24 25`) or write a small `program define` at the top of the file — do not copy-paste the block and edit suffixes by hand. Divergence between the math and reading copies is the most common bug in assembled files.
- Labels go **next to the `gen`**, not in a block at the end; a trailing block of 40 `label variable` lines drifts from the variables it describes. If a project already has such a block, keep it but add new labels at the creation site.
- Variable naming: `stem_subject_yy` or `stem_yy` consistently within a project; `rename (a b c) (x y z)` group syntax; build varlists with `ds`, `unab`, or wildcards (`recode tsi_*_24 (.=0)`) rather than 16 hand-typed names.

**Output contract** — the last block before `save`, every cleaning file:

```stata
//----------------------------------------------------------------------------//
// 3. Output contract
//----------------------------------------------------------------------------//
isid studentidentifier year window
assert inrange(grade, 3, 6)
assert !missing(treat_25, treat_24)
assert inrange(renstar, 0, 1400) if !missing(renstar)
count
assert r(N) > 0
order studentidentifier year window grade school classroom, first
compress
label data "Analytic file, math; built by 05_create_analytic_files.do"
char _dta[built_by]   "05_create_analytic_files.do"
char _dta[built_on]   "`c(current_date)'"
datasignature set, reset saving("${clean}/analytic_math", replace)
save "${clean}/analytic_math.dta", replace
codebook, compact                              // goes to the log
```

- `datasignature set` lets downstream files `datasignature confirm` that the input has not silently changed.
- `tempvar`s live until the **do-file ends**, not until you stop using them, so `save` writes them as `__000000`. Drop them explicitly before the contract and add `capture ds __*` / `assert "\`r(varlist)'" == ""`.
- When asserting that no raw variables remain (`ds q*`), make the pattern tight (`q[0-9]*`): `q*` also matches `quality1`.
- `save , replace` appears **once**, at the end. No intermediate saves over the output path.

---

## Phase 11 — Review checklist (for "review my cleaning code" requests)

Report findings in this order, most consequential first, each with the line and the fix:

1. **Silent data loss** — `drop if _merge`, unasserted `keep if`, `tag`-based dedup, `recode (.=0)` after merges.
2. **Fragile codings** — `encode` + manual `recode`, `!= "N"` indicators, positional `v24` renames, `capture` around mutations.
3. **Unchecked keys** — any `merge`/`reshape`/`bysort … keep` without a preceding or following `isid`.
4. **Repetition** — math/read or year blocks copied instead of looped; path globals duplicated across files.
5. **Reproducibility** — hard-coded absolute paths, no `version`, no log, no output contract.
6. **Hygiene** — dead code, diagnostic `tab`s, labels far from creation, `label variable` used for value labels.

Then, and only then, if the project lacks a master structure and would benefit, one sentence suggesting it (Phase 1 scope rule).

**Rewrites.** When the request is "rewrite/refactor this pipeline", the old output file is the test oracle. Before reporting, `merge 1:1` the new output to the old on the declared key and report: rows only in old, rows only in new, and rows where the values differ — each with a one-line cause. This both proves the refactor is faithful and is the most reliable way to find bugs in the old code (rows the old pipeline silently dropped show up as "only in new").

---

## Common pitfalls (quick index)

| Pitfall | Fix |
|---|---|
| `drop if _merge == 2` with no count | `merge …, keep() assert() gen(_m_x)` + `count`/note |
| `egen tag` → `drop if tag == 0` to dedup | `bysort key (sortvar): keep if _n == _N`, or `duplicates drop` for identical rows |
| `recode varlist (.=0)` after merging a usage file | Indicator from merge flag; zero-fill only where indicator says absent |
| `gen x = (s != "N")` | `gen byte x = (s == "Y") if inlist(s,"Y","N")`; decide blanks explicitly |
| `encode` then `recode (1=0)(2=1)` | `label define` first, `encode …, label()` , `assert` value set |
| `label variable x lbl` for a value label | `label values x lbl` |
| `rename v24 score_ela_2324` uncommented | `varnames(1)` import, or comment the header string |
| `capture replace …` / `capture destring …` | `capture` only on `confirm`; branch explicitly |
| Same block pasted for math and reading | `foreach subj in math read { … }` or a small `program` |
| Absolute Box path globals in every file | One root global in config; derive the rest (suggest only per Phase 1 rule) |
| Reshape on non-unique `i()×j()` fixed by ad-hoc counting | `isid` first; documented dedup rule |
| Labels in a block at file end | Label at the `gen` site |
| No `isid`, `assert`, log, or `datasignature` | Output contract block before `save` |
| `tab`, `browse`, `gen temp … drop temp` left in | Convert to `assert`/`count` or delete |
| Types differ across appended files | `stringcols(_all)` + deliberate `destring` before `append` |
| Imposing a master structure on an existing project unasked | Work in place; suggest once, briefly (Phase 1) |
| `label define` at top of file, then `use` → `encode …, label()` gives alphabetical codes | `use`/`import`/`clear`/`restore` drop labels; run `00_labels.do` right before `encode` |
| `capture label drop a b c` | Aborts at first missing name; use `label define …, replace` |
| `tempvar`s end up in the saved .dta as `__000000` | Drop them before the output contract; `assert` no `__*` remain |
| `preserve`/`restore` inside `frame X { }` | Use `frame put` / `frame copy` |
| `isid` errors "should never be missing" on a blank string key | `isid …, missok` with a comment |
| `ds q*` to assert raw vars are gone | Tight pattern `q[0-9]*`; `q*` matches `quality1` |
| Renaming by label regex or column position | Rename from the source's structural variable names |
| `input … end` data blocks inside do-files | CSV in `lookup/`, imported and asserted like raw data |
| Two chained reshapes for a 2-level grid | One reshape with `j = 100*level1 + level2`, then `floor`/`mod` |
| Refactor reported as "done" without comparing to old output | `merge 1:1` new vs old on the key; report only-old / only-new / differing |
