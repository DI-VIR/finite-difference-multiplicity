# Finite Difference Multiplicity Reproducibility Package v1.5.0

**Authors:** Chao Mu and Qichun Wang.  
**Version:** 1.5.0; prepared on 2026-10-12.  
**Status:** Publicly released on GitHub and archived in Zenodo (version 1.5.0).
**GitHub repository:** https://github.com/DI-VIR/finite-difference-multiplicity
**Zenodo version DOI:** https://doi.org/10.5281/zenodo.23313897


This package is prepared for public distribution under the licenses below.
The authors' licensing authorization has been confirmed. This statement does
not assert journal acceptance, repository publication, DOI registration, or
independent human verification of the manuscript's universal proofs.

## Manuscript and scope

**Extremal Difference Multiplicity in Finite Intervals: Boundary Capacity Bounds, Binary Realizability and Complexity.** Prepared for Discrete Applied Mathematics; no submission or acceptance status is asserted.

This scientific archive contains the finite positive-difference problem, window bounds, binary realizability, the complete finite map, archived medium MILP witnesses, literature comparison, four modular construction families and the explicit-candidate complexity result. Manuscripts are distributed separately. Cryptographic attack data and editorial working files are excluded.

## Contents

- `tables/`: 34 independently numbered frozen CSV tables; no LaTeX exports or separate table catalog are included.
- `data/`: 53 JSON record/protocol files and two CSV inputs, `verified_medium.csv` and `integrated_baselines.csv`. Complete records, witnesses and solver outcomes remain in JSON. `realizability_audit.json` preserves every field of the seven-case diagnostic, including fields omitted from the concise main table.
- `src/`: package-local Python implementations, organized by domain, with compatible flat entry points. Separate organization metadata files are omitted.
- `tests/`: 12 Python files for focused mathematical checks, shared test setup and module-import verification. Test instructions are consolidated in this README. The comprehensive checker is `src/reproduce.py`.
- `manifest/`: software environment, two experiment summaries and data provenance only.
- `MANIFEST.json`: SHA-256 and length of every distributed file, excluding the self-referential manifest itself.
- `LICENSE`: MIT for original software and CC BY 4.0 for original research data and standalone documentation, with explicit third-party exclusions.

## License and reuse

The package is mixed-license; it is not wholly MIT or wholly CC BY 4.0.

| Material | License | Scope |
| --- | --- | --- |
| Python code in `src/` and `tests/`, runnable code examples, `requirements.txt`, and in-code documentation | MIT | Original software contributions by the authors |
| Original records, witnesses, certificates, computed CSV tables, protocols, metadata, and standalone documentation | CC BY 4.0 | Original contributions in `data/`, `tables/`, `manifest/`, `MANIFEST.json`, and the narrative text of this README |
| Published source cells, cited third-party works, and external dependencies | Original rights and license terms | Not relicensed by this package |

The complete MIT text and the binding CC BY 4.0 legal-code link are in
`LICENSE`. Retain the applicable notices when redistributing code, and give
attribution and identify changes when reusing CC BY 4.0 material. Suggested
academic citation is not an additional license restriction.

`data/published_diameters.json` preserves a source URL for every published
numerical cell. `data/baseline_protocol.json` identifies the cited versions
and finite formulas used in comparisons. The authors' compilation,
annotations, and computed comparisons are distinct from the underlying
third-party source content. The package does not claim ownership of those
source works or new exclusive rights over mathematical facts. Check the
original source terms when permission is required for reuse.

External dependencies are installed from `requirements.txt`; their source
distributions are not bundled here. The manuscript, supplementary document,
publication artwork, and private editorial files are outside this release
and outside its license grants.

## Citation and persistent identifiers

This repository contains the reproducibility code, computational data, exact certificates, and verification records associated with the manuscript:

Extremal Difference Multiplicity in Finite Intervals: Boundary Capacity Bounds, Binary Realizability and Complexity

Authors: Chao Mu and Qichun Wang
Release version: 1.5.0
Release year: 2026
GitHub repository: https://github.com/DI-VIR/finite-difference-multiplicity
Zenodo version DOI: https://doi.org/10.5281/zenodo.23313897

The intended record title is:

> Extremal Difference Multiplicity in Finite Intervals: Boundary Capacity
> Bounds, Binary Realizability and Complexity - Reproducibility Code and Data

The citation should contain the authors, this record title, the deposited
version, the actual publication year, Zenodo, and the assigned version DOI.
These are metadata instructions, not a claim that a Zenodo record exists.
The code/data DOI is separate from any future journal article DOI.

## Public release procedure

1. Create a new, independent GitHub repository from this package directory.
   Commit the extracted files, not only a ZIP. Do not publish the entire
   development workspace, its backups, or its previous Git history.
2. Review all 150 files, the applicable licenses, and the source attributions.
   Keep repository descriptions and release text in English. Run the checks
   below in a disposable copy; verification can create new output files.
3. Before creating a GitHub release, connect Zenodo to the intended public
   repository and enable archiving for that repository.
4. Create the release tag `v1.5.0` only when ready to publish. This is the
   proposed tag, not an existing tag. Use the release-note text below.
5. Check the Zenodo archive contents, title, creator order, version, and
   mixed-license metadata. Declare MIT for original software and CC BY 4.0
   for original data/documentation, with the exclusions in `LICENSE`; do not
   describe the complete upload as unconditionally CC BY 4.0.
6. Verify that the assigned DOI resolves and that the files can be downloaded
   without authentication. Only then describe them as publicly available
   in the manuscript and submission system.
7. Record the actual repository URL and version DOI in the maintained README.
   A later repository commit does not change the archived release. Publish
   changed research files as a new archived version and cite the version
   actually used. Updating this README also requires refreshing its SHA-256
   entry in `MANIFEST.json` for any new distributed snapshot.

As an alternative to GitHub integration, deposit the matching ZIP manually
in Zenodo as a Software record and link the GitHub repository in its metadata.
Use one archival route for this version; do not create duplicate DOI records
for the same package. For either route, an unsubmitted draft or reserved DOI
is not evidence that the package is publicly available.

### Release-note text

**Release title:** Reproducibility package v1.5.0

This release provides the reproducibility code and data accompanying
"Extremal Difference Multiplicity in Finite Intervals: Boundary Capacity
Bounds, Binary Realizability and Complexity" by Chao Mu and Qichun Wang.

- Contains 150 files: 55 data files, 41 Python source files, 34 CSV tables,
  12 Python test files, four provenance/environment/experiment JSON files,
  and four root files.
- Provides exact finite records, witnesses, rational Fourier certificates,
  declared baseline comparisons, construction records, and verification code.
- Retains all scientific data, source code, tests, and table values from the
  trimmed v1.5 review snapshot; changes are limited to licensing and release
  documentation/metadata.
- Uses MIT for original software and CC BY 4.0 for original data and
  standalone documentation. Third-party material retains its original terms.
- Separates saved-result verification, fresh enumeration, historical solver
  outcomes, and longer experiment reruns. CSV re-export is not a fresh
  independent recomputation of every numerical result.
- Documents the tested Windows environment and does not claim Linux testing
  or independent human verification of universal mathematical proofs.

This text is prepared for the future release. The local package has not yet
been published, and no public URL or DOI is asserted here.

## Archive cleanup

The initial 2026-10-12 cleanup removed 78 generated reports, historical logs and byte-identical protocol copies. The CSV-only cleanup removed a further 61 files: manuscript reading copies, all manuscript images, LaTeX tables, the table catalog, source organization metadata and metadata outside environment, experiment summaries and data provenance. The final cleanup removes 28 optional data CSV exports and two redundant test entry/documentation files, leaving exactly 150 files. The preceding archives are retained privately.

Every omitted data CSV has a retained JSON counterpart. The two CSV files read by verification commands are retained, as are all 34 CSV tables. `python src/reproduce.py --full` replaces the removed full-map test wrapper. No retained data record, source implementation, test or table value is changed by this final cleanup.

Scientific records, mathematical certificates, witnesses and the retained CSV values are unchanged. File I/O adapters now read and write CSV tables directly, and import checks discover Python modules without an organization index. No numerical bound or certificate-verification formula is changed.

Regenerated outputs in `recomputed/` are not distributed. Running `src/reproduce.py` or `src/regenerate_tables.py` creates this directory again, with CSV table exports only. These exports re-serialize the frozen CSV tables; they do not independently rebuild all numerical values from raw data. Scientific data and certificates are checked separately. Some audit commands likewise recreate report files under `tests/`; these outputs are not required as verification inputs.

The four literature-audit exports previously had byte-identical protocol files. `data/literature_audit_reverse_protocol.json` now supplies their shared protocol; it also applies to `literature_audit_small`, `literature_audit_strata` and `literature_audit_summary`. Their distinct scientific records remain separately archived in JSON.

The cleanup itself did not change licensing. This separately prepared public
version replaces the review-only license with the author-authorized grants
in `LICENSE` and refreshes `MANIFEST.json` for the same 150-file set.

## Independent checks

For this public-release preparation, 14 verification commands passed on
Windows on 2026-10-12 in a disposable copy. These comprised the ten focused
`test_*.py` scripts, module-import checks, `src/reproduce.py --full`, integrated
baseline verification, and CSV re-export. The full enumeration rechecked all
252 nontrivial pairs at lengths 4 through 24 and all 4,046 saved binary-window
minima. This check did not collect a new large solver or stochastic search
population. Generated logs and reports are not part of the 150-file release.

The integrated baseline update adds `python src/integrated_baselines.py` and `python tests/test_integrated_baselines.py`. They recompute the pigeonhole, optimized classical-window, integer-balancing, two explicit Xu finite inequalities, GO lemma/table inversion and boundary bounds on the 702 distinct pairs in the union of the existing domains. `data/baseline_protocol.json` records exact arithmetic, diameter normalization, source versions and coverage limitations; `data/baseline_summary.json` keeps the domains separate. Main Table 3 and Table S2b are generated from these records. Historical solver/search populations and stopping certificates are unchanged. This update evaluates declared finite baselines, not every published LP/SDP or asymptotic result.

The checked Windows environment uses Python 3.12.7, NumPy 2.2.6, SciPy 1.16.3, SymPy 1.14.0 and HiGHS 1.8.0. Use Python 3.12 and install `requirements.txt` into an isolated environment. From the extracted package root:

```sh
python -m pip install -r requirements.txt
python src/reproduce.py
python src/reproduce.py --full
python src/regenerate_tables.py
python tests/test_bounds.py
python tests/test_realizability.py
python tests/test_np_reduction.py
python tests/test_line_endings.py
python tests/test_revision.py
python tests/test_final_details.py
python tests/test_risk_closure.py
python tests/test_extended_revision.py
python tests/test_fourier_packet.py
python tests/test_integrated_baselines.py
python tests/check_imports.py
```

The default command checks all 252 saved witnesses, 21 trivial rows, 714 joint integer vectors, the original five moment-hierarchy certificates and five strengthened vectors, all 126 exceptional walk windows, 1,350 modular construction outputs, the 462-pair literature comparison and reduction fixtures. It independently enumerates intervals through length 12. `--full` re-enumerates every interval length 4 through 24 and checks all 4,046 binary-window minima, extrema, witness choices, optimum counts and translation/reflection classes. These computations do not use the private development workspace. Runtime is machine-dependent.

### Test domains

Run the checks from the extracted package root in an isolated Python environment; no companion or private master folder is required. Canonical implementations are stored in the domain subdirectories of `src/`. Compatible flat entry points load those implementations, and `tests/check_imports.py` discovers and checks the Python modules directly. Revision audits, residual diagnostics and CSV table exporters are standalone modules; modular comparisons are under `src/constructions/`.

`python src/reproduce.py --full` is the complete finite-map verification command. Table re-export is separate from fresh numerical computation and does not recreate historical timings or every large solver run. The reduction tests check 406 finite graph fixtures, not the universal mathematical proof in main Section 7.

`tests/test_final_details.py` checks tent tails, layer remainders, ceiling strictness, the slack identity, worked examples, half-density formulas, algebraic constructors and meeting medium lower certificates. `tests/test_risk_closure.py` checks 14 symbolic identities and repeats the weighted binary-window MILP on 12 small instances against independent pair enumeration. These checks are not evidence of general scalability, a new medium-solver population or independent human proof review.

Runtime categories are distinct and exclude dependency installation:

- Fast verification: `python src/reproduce.py` validates saved certificates and enumerates only lengths through 12. Each run records its measured time in the regenerated `recomputed/verification.json`.
- Complete enumeration: `python src/reproduce.py --full` re-enumerates lengths 4 through 24. Each run records its measured time in `recomputed/verification.json`; this is not a cross-platform guarantee.
- MILP reruns: the four medium feasibility jobs and the original/strengthened residual solves are archived historical experiments. The longest listed medium job took about 287 seconds. Their witnesses are rechecked without representing their historical timings as new solves; larger limits and timeouts are outside the short verification budget.
- Construction experiments: `python src/constructions/modular_comparison.py --seeds 30 --budget 3000` regenerates 1350 streams. The preceding run took approximately 15 seconds. Final closure rechecks the saved outputs and all 36 algebraic constructors; it does not claim a newly collected stochastic population.

Linux execution has not been performed in this environment. The archive does not include a Linux CI recipe. Line-ending tests explicitly compare LF, CRLF, CR and UTF-8 BOM inputs.

`src/milp_difference_optimum.py` retains the binary pair-variable formulation and independent witness checker. Solver status and limits are preserved in `data/milp_complete.json`. A feasible witness meeting an independently valid lower bound proves equality; a timeout does not. `src/np_reduction.py` implements the complete construction and small graph fixtures. Those fixtures are not a substitute for the main-text universal proof.

CSV table re-export compares canonical UTF-8 text after BOM and line-ending normalization. A serialization content change still fails the check. Numerical baseline tests independently compare the relevant table rows with computed values. Table re-export does not rerun every historical medium MILP. Longer solver runs can return different timings or feasible witnesses; unresolved outcomes remain unresolved.

Additional experiment commands are:

```sh
python src/revision_audit.py
python src/residual_diagnostics.py
python src/constructions/modular_comparison.py --seeds 30 --budget 3000
```

The modular comparison uses seeds 2026100200 through 2026100229 and a budget of 3,000 proposals for each method and parameter pair. Exact constructor definitions, finite-field choices, projection losses and all final witnesses are stored in `data/modular_constructors.json` and `data/modular_runs.json`.

`python src/final_details_audit.py` checks 4324 small boundary-window domains, 466 binary sets, four hand examples, 42 half-density instances, all constructor metadata, equal budgets and the meeting lower certificates in main Table 2 (original evidence key 6). It regenerates field trace sets and primitive-order checks. The final Theorem 8 witnesses are `strengthened_witness` in `data/residual_diagnostics.json`; original vectors remain in `data/realizability_audit.json`. The 406 marker-reduction fixtures are finite checks, not independent human proof review.

## Conventions and provenance

The v1.2 scope closure adds Table S15, 14 symbolic identity checks and a weighted binary-window MILP formulation independently compared with enumeration on 12 instances (L=4,6,8,10,11,12). Run `python tests/test_risk_closure.py` to repeat these auxiliary checks. They do not establish scalable performance or replace independent human review. Recorded small-solve outputs are in `data/binary_window_milp_validation.json`; the original medium solver and 1350 construction records are unchanged.

Positions are `0,...,L-1`; ruler diameter is at most `L-1`. The cap constrains positive-difference multiplicity. The historical raw filename `realizability_finite_b2g.json` is retained for traceability, but its values are interpreted as the **difference-bounded** cardinality function, not an automatic identification with the unordered-sum convention of a general B2[g] set. This correction changes no value.

Historical Revision5 data and the preceding submission-revision experiments are retained with their distinct provenance. The final closure adds verified constructor metadata, continuous-k checks in the existing table, revised captions and focused audits without changing the previously observed experimental population. `manifest/revision_provenance.json` distinguishes these sources. Finite computational checks are not a substitute for independent human review of universal mathematical proofs.

The classic-LP discussion does not claim a proved dominance or projection equivalence. The five residual obstructions form a finite result, not an infinite family. Modular methods are the declared prime-field variants, not all possible prime-power constructions. Licensing authorization is distinct from independent mathematical proof review. Repository publication, a release tag, and DOI registration remain pending.

Recent-reference metadata was checked on 2026-10-04. Bailleul--Riblet is cited from the formal European Journal of Combinatorics record, volume 139, article 104443, print year 2027 (DOI registered in 2026); no online publication date is inferred from that year. Bosio--Riblet--Tarr is cited as arXiv:2607.02728v2. Gupta--O'Bryant, Gupta and Croot et al. retain their checked preprint versions. Absence of a matched journal record in the checked sources is not proof that no such record exists; refresh status again at actual submission.

The prepared public snapshot is `DAM_Reproducibility_v1.5_public.zip`, with
package version 1.5.0. The original `DAM_Reproducibility_v1.5_review.zip` and
previous development snapshots are retained privately without modification.
This is a local release preparation, not an assertion of an existing public
repository release. Editorial correspondence, author confirmations, private
proof checklists and document renders are not distributed in this archive.

## New Fourier and experimental comparisons

The new records are separate from the unchanged historical populations. `data/fourier_exact_trees.json` supplies five complete integer branch trees with 758 rationally certified exclusion leaves. `python src/fourier_verify.py` independently reconstructs and checks their inequalities and branch coverage without calling an optimizer. Fourier cuts close the five residual scalar gaps on the stated pairs. This does not prove general exactness, and finite Toeplitz PSD is not treated as the same constraint. No Lasserre hierarchy was executed.

`data/seed_budget_runs.json` contains 1620 streams and checkpoints at 0,3000,30000 proposals, totalling 4860 verified witnesses. Random-SA is an implemented generic annealing control, not a reproduction of a published algorithm. The original 1350 hill-climbing outputs equal the new 3000-proposal prefixes. Pointwise paired bootstrap intervals are exploratory and not used for familywise superiority claims. Initial algebraic sets are deterministic; constructor generation time is not included in the matched proposal budget.

`data/window_scalability.json` records six weighted binary-window MILP runs at lengths 28,32,40,48,64,80, each with a 60-second limit. Two checked objectives meet the exact analytic capacity lower bound. Four are time-limit brackets. Numerical dual bounds are retained but are not represented as exported exact dual certificates. These single-window runs do not compute all-window B_real or Delta.

`data/medium_witness_provenance.json` matches full coordinates to the original historical logs and records source row, source-file hash, seed, method and budget where available. Historical profile objectives and modular direct-multiplicity search are distinct populations. No original seed or budget is invented for solver records. `data/xu_comparison.json` evaluates the printed finite inequalities of Xu (2018) in the same positive-difference model with exact rational ceilings. Reference metadata and author-version theorem numbering were checked on 2026-10-05.

To verify all new records, run `python src/extended_audit.py`. Longer experiment reruns are:

```sh
python src/fourier_comparison.py
python src/fourier_exact_tree.py
python src/seed_budget_ablation.py
python src/window_scalability.py
```

The recorded seed-budget run took about 362 seconds on the tested Windows environment while other solver processes were active; this is not a portable performance estimate. The binary-window limits apply per solve, not to the complete program. Table S16 through Table S21 retain the new comparisons. `python src/review_revision_tables.py` rebuilds these additions from package-local records; the optional `--import-history` is a development import step and is not required by reviewers.

The strong NP-completeness wording concerns explicit sparse candidates with polynomially bounded coordinates. Universal proofs still require independent human review. Public redistribution of the authors' original contributions is authorized under `LICENSE`; third-party rights are unaffected. No existing public repository, DOI, journal acceptance, or independent human proof approval is asserted by this local preparation.

## v1.4 review closure

The Fourier entry point now requires the complete five-case set, the stated caps, rational cosines in [-1,1], complete scalar windows and well-formed integer leaf boxes. `python tests/test_fourier_packet.py` rejects seven packet mutations, including whole-case deletion. Existing leaf-deletion and negative-weight tests remain in `tests/test_extended_revision.py`. These contract checks change neither the 758 certificates nor any recorded scientific result.

Random initial coordinates in `data/seed_budget_runs.json` are stored sorted for set verification, but proposals depend on the internal list order. Regenerate that order with `random.Random(seed).sample(range(L),n)` for Random and Random-SA, and use a separately initialized `random.Random(seed)` for proposals, as in `src/seed_budget_ablation.py`. Do not restart a trajectory from the sorted saved list. Algebraic starts use the constructor's returned order. The full script and saved seed identifiers reproduce the ordering.

The manuscript now attributes the inherited window identity to Carter, Hunter and O'Bryant and the Cauchy--Schwarz window inequality to Caicedo et al., distinguishes Martos et al.'s algebraic constructions from the generic SA control, and credits the known diameter exclusions behind the five Fourier repairs. The LaTeX draft retains the full Section 7 proof. This historical review-closure description is not a current publication-approval record; the present package's licensing and deposit status are stated above.

## v1.5 Main-Text Compression

All raw scientific records are unchanged. The main text keeps the boundary/slack/infinite-family proofs, reconstruction, scalar and Fourier separation, exact classification and the full strong NP-completeness reduction. The original tables 4,6,8 appear as main Tables 1,2,3; the other six main tables move to SI. Original figures 1,4,5 appear as main Figures 1,2,3; figures 2,3 move to SI as S2,S3. New Sections S22--S26 preserve relocated material.

CSV filenames retain their original evidence identifiers for stable reproducibility. `Table_4.csv`, `Table_6.csv` and `Table_8.csv` correspond to current main Tables 1, 2 and 3, respectively. The standalone CSV tables are not the fully typeset article. Manuscript reading copies, images and publication-numbering metadata are not included in this CSV-only archive.
