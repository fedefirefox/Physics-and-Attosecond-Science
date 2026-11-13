# Validation record

Local validation date: **13 November 2026**.
Source inspected: `c85764f07abc484900e5aa2499e3e163d941b3a2`.
Runtime: MATLAB **R2023b Update 5**, `23.2.0.2459199`, Windows x86-64.

## Executed MATLAB tests

Command from the repository root:

```matlab
results = runtests(fullfile(pwd,'SFEXIS','tests'));
assertSuccess(results);
```

**11 passed, 0 failed, 0 incomplete.** Test execution took approximately 9.5 s
in the first local run, excluding MATLAB startup. Runtime depends on the host;
graphics export accounted for most of that time.

| Check | Measured result / acceptance |
| --- | --- |
| Default run | 25×101 finite, nonnegative signal; global normalized maximum 1 |
| Zero-IR analytic Gaussian spectrum | Maximum absolute error / analytic peak: `4.74e-15` (limit `1e-10`) |
| Non-grid-aligned delay translation | Complex phase agrees with `exp(i*(E+Ip)*tau)` within `1e-10` relative error |
| Constant-vector-potential analytic action | Agreement within `1e-11` relative error; checks `p*A` and `A^2/2` |
| Electric field / potential and intensity | Derivative error below `5e-5`; peak intensity agrees within `1e-12` relative error; tails below `1e-9` a.u. |
| Streaking sign | Centroid shift correlation with `-p0*A + A^2/2`: `0.999978`; shift exceeds 1 eV |
| Time-step convergence | Relative Frobenius errors versus 2.5 as: `2.34e-5` at 10 as and `4.68e-6` at 5 as |
| Window extension ±20 to ±25 fs | Raw-signal relative difference below `1e-9` |
| Original rectangular-rule constant-dipole formula | Relative difference from new trapezoidal result: `2.68e-4` at 2.5 as (limit `5e-3`) |
| Invalid inputs | Unknown options, bad axes/intensity, zero step, short window and undersampling rejected |
| Example output | MAT round trip exact, PNG exported; existing outputs rejected |

The legacy-formula test uses the original action/summation formula on the new
test fields. It does **not** execute the original interactive script or claim
identical pulse generation. Constant-vector-potential input is a kernel
mathematical control, not a physically isolated laser pulse.

## Source and dependency inspection

All ten original MATLAB files were inspected, including every difference in
the backup script. All eight original data/output files were inventoried.
MATLAB reports only base MATLAB for the maintained runtime and MATLAB plus
Signal Processing Toolbox and Curve Fitting Toolbox for the legacy main
script. MATLAB Code Analyzer reported no messages for the three new runtime
files. The guide lists unresolved legacy behavior and data provenance.

`verify_legacy.py` compares all 18 current archived files against the original
Git blob hashes. A scoped `.gitattributes` rule disables line-ending conversion
for the archive so the check compares its original bytes on every platform.
`validate_citation.py` checks both CFF files against the official 1.2.0 schema,
checks their equality, and verifies the expected DOI and six-author list.
These are also included in CI. Consult the PR checks for the actual hosted
run and its logs; adding a workflow alone is not evidence that it passed.

Both checks passed locally: **18/18 original blobs match**, and **both CFF files
pass the official CFF 1.2.0 schema and article-metadata checks**. The generated
PNG was also visually inspected: axes, color scale and citation are visible.

## Limits and remaining work

- The original interactive simulation, train and dynamic/tabulated dipole
  branches were not run end to end. They remain archived with known issues.
- No energy-dependent atomic phase, ensemble integration or inverse retrieval
  is implemented by the maintained example. No published figure is reproduced.
- The zero-IR spectrum and other tests validate this implementation, not the
  accuracy of SFA for real argon. Coulomb/multielectron benchmarking remains
  outside this example.
- Convergence is demonstrated for the documented small case, not arbitrary
  configurations. Other parameter regimes require their own grid studies.
- GNU Octave and other MATLAB versions have not been locally tested.
- Atomic-table source/units and third-party redistribution terms remain
  unresolved. No new software license was selected.

The initial restricted MATLAB launch could not access its preferences/licensing
service. The successful tests used the installed licensed MATLAB with a
task-local preferences directory. No repository change or toolbox installation
was needed to resolve that local execution issue.
