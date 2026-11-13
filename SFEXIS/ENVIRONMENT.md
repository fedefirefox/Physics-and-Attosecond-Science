# Environment and reproducibility

## Simulation and tests

Reference runtime: **MATLAB R2023b Update 5**, version `23.2.0.2459199`, on Windows
x86-64. Base MATLAB is sufficient for the new implementation, example, graphics
export and `matlab.unittest` function-based tests. The dependency analyzer reports
only MATLAB for `run_sfexis.m`. No startup script, external dataset, MEX library,
parallel pool or additional MATLAB toolbox is needed.

Run from the repository root in a terminal with MATLAB on PATH:

```text
matlab -batch "addpath('SFEXIS','SFEXIS/examples'); minimal_streaking();"
matlab -batch "results=runtests(fullfile(pwd,'SFEXIS','tests')); assertSuccess(results);"
```

In PowerShell, if MATLAB is not on PATH, prefix the executable explicitly:

```powershell
& 'C:\Program Files\MATLAB\R2023b\bin\matlab.exe' -batch "results=runtests(fullfile(pwd,'SFEXIS','tests')); assertSuccess(results);"
```

The workflow in `../.github/workflows/sfexis.yml` requests MATLAB R2023b on
Ubuntu 22.04 with no extra toolboxes. The hosted installer may provide a later
R2023b maintenance update; its log records the actual version. The workflow
runs the same tests and retains test-generated MAT/PNG outputs. MATLAB's
installation/licensing service and the headless graphics exporter must work
on the host. CI availability is separate from local test success.

No GNU Octave run has been performed. Older MATLAB releases and other operating
systems should not be described as tested until their test logs exist.

## Archive and citation checks

These checks are separate from the MATLAB runtime. Python 3.12 or newer and Git
are sufficient for the archive check (Python 3.14.3 was used locally):

```text
python SFEXIS/tests/verify_legacy.py
```

The citation check needs PyYAML 6.0.3 and jsonschema 4.25.1. Install these in a
virtual environment rather than altering a system Python, then run:

```text
python -m pip install PyYAML==6.0.3 jsonschema==4.25.1
python SFEXIS/tests/validate_citation.py
```

The validator downloads the official CFF 1.2.0 schema over HTTPS. For offline
validation, supply an already downloaded schema as its only argument:

```text
python SFEXIS/tests/validate_citation.py path/to/cff-1.2.0-schema.json
```

Both CFF copies must be identical. DOI, article number and all six author names
are checked in addition to schema validation. The article DOI is placed only
inside `preferred-citation`, never used as a DOI for the software.

## Legacy environment

`legacy/SFEXIS.m` and `backup30Ottobre.m` are archived research scripts with
interactive graphics, relative data loads and Windows-specific paths. Their
static dependency report includes MATLAB, Signal Processing Toolbox
(`envelope`) and Curve Fitting Toolbox (`smooth`). `plotPsi`'s automatic
threshold branch calls a missing `autoThreshold` helper. Installing those
toolboxes does not resolve the other source issues listed in the guide.

The new example does not invoke these scripts or silently depend on their
MAT files. The archived default simulation and its optional branches were not
executed end to end during this reorganization.
