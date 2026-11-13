# SFEXIS — Strong Field Electron XUV Interaction Simulator

SFEXIS is Federico Vismarra's MATLAB code for XUV photoemission in an IR field
within the strong-field approximation. This project contains a small,
noninteractive constant-dipole example and the complete original research bundle.

**Credit Federico Vismarra and cite the following paper when using or adapting SFEXIS:**

F. Vismarra, R. Borrego-Varillas, Y. Wu, D. Mocci, M. Nisoli and M. Lucchini,
“Ensemble effects on the reconstruction of attosecond pulses and photoemission
time delays,” *Journal of Physics: Photonics* **4**, 034006 (2022).
[DOI: 10.1088/2515-7647/ac7991](https://doi.org/10.1088/2515-7647/ac7991).
See [CITATION.md](CITATION.md) for BibTeX, attribution and contributor details,
and [CITATION.cff](CITATION.cff) for machine-readable metadata.

## Run the example

The supported reference environment is **MATLAB R2023b, base MATLAB only**.
No downloaded data, additional toolboxes, Python packages or random seed are
needed for the simulation. MATLAB must be installed and licensed. GNU Octave
has not been tested. See [ENVIRONMENT.md](ENVIRONMENT.md) for exact validation
versions and command-line usage.

From the repository root in MATLAB:

```matlab
addpath('SFEXIS', 'SFEXIS/examples');
result = minimal_streaking();
```

This computes 25 delays × 101 electron energies for a 250 as XUV pulse centered
at 50 eV, a 5 fs IR vector-potential envelope at 800 nm, a peak IR intensity of
10^12 W/cm², and the Ar ionization threshold. The dipole is set to one: this is
an Ar-threshold demonstration, not a quantitative argon atomic model.

The example writes `SFEXIS/results/minimal_streaking/streaking.mat` and
`streaking.png`. The MAT file contains parameters, units, fields, the complex
amplitude, raw and normalized signals, MATLAB version and citation. Existing
outputs are protected; choose a new folder for another saved run:

```matlab
result = minimal_streaking(fullfile('SFEXIS','results','my_second_run'));
```

For an in-memory calculation with no figures or file writes:

```matlab
result = run_sfexis();
field_free = run_sfexis(struct('ir_intensity_Wcm2',0));
custom = run_sfexis(struct('delay_fs',-2:0.1:2,'energy_eV',20:0.25:55));
```

`signal` and `normalized_signal` are **delay × energy** arrays. Positive delay
places the XUV pulse later than the IR envelope. The normalized signal uses one
maximum for the entire trace. Raw signals have arbitrary scale, not absolute
ionization probabilities. Finer sampling is needed for convergence studies.

## Run tests

From the repository root:

```matlab
results = runtests(fullfile(pwd,'SFEXIS','tests'));
assertSuccess(results);
```

The 11 tests cover analytic zero-IR and constant-vector-potential spectra,
delay-phase translation, the electric-field/vector-potential relation,
streaking sign, time-step and window convergence, the original constant-dipole
quadrature, invalid inputs, and saved-file round trips. The example test leaves
its outputs in an ignored, uniquely named directory under `SFEXIS/results`.
See [VALIDATION.md](docs/VALIDATION.md) for measured errors and limits.

## Project map

| Path | Purpose |
| --- | --- |
| `run_sfexis.m` | Validated inputs, fields, unit conversions and output metadata |
| `sfexis_defaults.m` | All parameters for the small reference case |
| `sfexis_spectrum.m` | Constant-dipole SFA integral in atomic units |
| `examples/minimal_streaking.m` | Saved example and preview figure |
| `tests/test_sfexis.m` | MATLAB physics and behavior tests |
| `docs/GUIDE.md` | Equations, conventions, implementation, usage and complete source audit |
| `docs/legacy-manifest.json` | Original commit and Git blob hashes for all 18 archived files |
| `tests/verify_legacy.py` | Git/Python check that the archived files are unchanged |
| `legacy/` | Entire original bundle, including historical results |

Add only `SFEXIS` and `SFEXIS/examples` to the MATLAB path for the example.
Do not use `addpath(genpath(...))`: the repository has unrelated functions with
the same name (`EField1D`), and `legacy/` is intentionally isolated.

## Scientific scope

The runnable path evaluates a single spatial point and one detection direction
using a constant dipole and isolated, unchirped Gaussian pulses. The archived
code also contains pulse-train and dipole branches, with unresolved issues
documented in the guide. None of those branches is silently enabled here.

The 2022 paper studies ensemble effects and reconstruction. This repository
bundle does **not** include its full spatial averaging or ePIE retrieval code.
The example is not a reproduction of the paper's figures, a TDSE solver, or a
photoemission-delay retrieval. Numerical checks validate the implemented model
and conventions; they do not establish quantitative atomic accuracy.

## Preservation and reuse

The original `Basic_SFEXIS_SourceCode` tree is moved intact to `legacy/`.
History is retained, and the old path contains a redirect README. No historical
MAT or FIG file is regenerated. See [CITATION.md](CITATION.md) for preserved
third-party credit and the repository's currently unspecified software license.
