# SFEXIS: physics, implementation and usage

SFEXIS was written by Federico Vismarra. When using or adapting it, cite
F. Vismarra, R. Borrego-Varillas, Y. Wu, D. Mocci, M. Nisoli and M. Lucchini,
“Ensemble effects on the reconstruction of attosecond pulses and photoemission
time delays,” *Journal of Physics: Photonics* **4**, 034006 (2022),
[doi:10.1088/2515-7647/ac7991](https://doi.org/10.1088/2515-7647/ac7991).
See [CITATION.md](../CITATION.md) for provenance and the complete reference.

## 1. What is included

The source inspected for this reorganization is commit
`c85764f07abc484900e5aa2499e3e163d941b3a2`. All 18 files formerly in
`Basic_SFEXIS_SourceCode` are retained in `SFEXIS/legacy`, with their original Git
blob identities. The manifest records each old relative path and blob hash.
No source history was rewritten and no historical simulation was regenerated.

The maintained entry point, `run_sfexis`, extracts the constant-dipole SFA
calculation into a small function with explicit parameters and deterministic
outputs. It uses new, explicitly defined Gaussian fields, a consistent delay
axis and trapezoidal quadrature. The archived files remain available for
comparison and further development, but are not presented as repaired or fully
validated. This separation makes numerical changes visible and preserves the
original research record.

The code is a **forward model**: input pulses produce a photoelectron trace. It
does not solve the inverse problem of reconstructing pulses from a trace.
`ePie_map.mat` is a colormap, not an ePIE implementation.

## 2. Physical model and conventions

### Single-point strong-field integral

All quantities in the integration kernel are in atomic units. Let the final
positive momentum along the polarization/detection axis be

\[
p=\sqrt{2E},\qquad E\geq0.
\]

For the maintained constant-dipole model, the amplitude is

\[
a(p,\tau)=-i\int_{t_0}^{t_f}E_\mathrm{XUV}(t-\tau)
\exp\!\left[i(E+I_p)t-i\int_t^{t_f}
\left(pA_\mathrm{IR}(s)+\frac{A_\mathrm{IR}(s)^2}{2}\right)ds\right]dt.
\]

Here `potential` means the IR **vector potential**, not the atomic binding
potential; the electric field satisfies `E_IR = -dA_IR/dt`. The dipole is one.
The finite endpoints approximate infinity only when both pulse tails and all
delayed XUV pulses are negligible at the edges.

The original main script instead stores the full kinetic action

\[
S_p(t)=\int_t^{t_f}\frac{[p+A(s)-A(t_0)]^2}{2}\,ds
\]

and multiplies `exp(-1i*S_p)` by `exp(1i*Ip*(t-t0))`. When `A(t0)` vanishes,
expanding the square produces the maintained expression up to the constant
phase `exp(-i*E*tf-i*Ip*t0)`. This phase cancels in the signal. The new kernel
uses trapezoidal integrals rather than the original rectangular sums, so exact
floating-point equality to the legacy result is not expected.

The paper writes the shifted wavepacket as `E_WP(t+tau)`. The example uses
`E_XUV(t-tau)`: **positive example delay means the XUV arrives later**, and is
the negative of the paper's delay variable for otherwise identical fields.
The delay-translation test checks this convention using the complex amplitude,
not only an intensity image. No interpolation or integer rounding is used to
shift the analytic XUV field.

### What the plotted signal means

The original code uses `WE = p*abs(a)^2`. The maintained code preserves that
weighting and calls the output `signal`. If `abs(a)^2` is interpreted as a 3D
momentum density, then `d^3p = p^2 dp dOmega = p dE dOmega` gives a directional
energy density proportional to `p*abs(a)^2`. This is **not** the Jacobian for a
strictly 1D momentum density, which would instead involve `1/p` for positive
momentum. Here “1D” refers to a single interaction point and detection direction;
it is not a 1D TDSE wavefunction normalization.

The XUV peak field and dipole are unit amplitudes. The output therefore has an
arbitrary scale. It is not an absolute ionization probability, total yield or
angularly integrated detector count. The raw complex amplitude is returned so
that other explicitly documented measures can be used without changing the
kernel. At zero kinetic energy the retained `p`-weighted signal is zero.

`normalized_signal` divides the complete trace by its single largest value.
It does not normalize each delay independently. When averaging multiple
simulations, combine the raw signals using physically justified weights first;
independent peak normalization would remove relative yields.

### Approximation boundaries

The XUV interaction is treated perturbatively and the IR drives continuum
propagation through the Volkov phase. The example neglects the residual ionic
potential during that propagation, multielectron correlation, depletion,
angular dipole structure, focal geometry and detector response. A constant,
real dipole contains no atomic phase and therefore cannot encode an
energy-dependent atomic photoemission delay. Selecting an Ar ionization
potential only places the threshold; it does not supply an argon dipole.

### Relation to ensemble effects

The cited paper examines how spatially varying fields affect measured traces
and their reconstruction. A spatial ensemble requires an incoherent sum of
local probabilities, schematically

\[
S_\mathrm{ensemble}(E,\tau)=\int n(\mathbf r)\,w(\mathbf r)\,
S(E,\tau;\mathbf r)\,d^3r,
\]

with the local fields, gas density, collection weighting and volume measure
specified consistently. A cylindrically symmetric integration includes the
radial measure `2*pi*r dr dz`. Amplitudes from independent atoms are not summed
coherently for this model. The supplied source has no such spatial integration
and no reconstruction routine. The example and its tests do not reproduce or
validate the paper's ensemble or retrieval results.

## 3. Fields, units and parameters

The conversions retained from the original source are
`1 a.u. of time = 0.02418884326505 fs` and
`1 hartree = 27.21138624 eV`. The new field scaling uses
`I [W/cm²] = 3.5094452e16 * E_peak [a.u.]^2` and
`c = 299.792458 nm/fs`. These constants are returned in every result.

The XUV electric field is

\[
E_X(t-\tau)=\exp[-2\ln2((t-\tau)/T_X)^2]
\cos[\omega_X(t-\tau)+\phi_X],\qquad \omega_X=\hbar^{-1}\mathcal E_X.
\]

`xuv_fwhm_fs` is the FWHM of its **intensity envelope**, not the FWHM of its
electric-field amplitude and not of its oscillating instantaneous intensity.
The same transform-limited pulse is shifted as a whole, including its carrier.

The IR is constructed from a Gaussian vector potential,

\[
A_I(t)=-C\,g(t)\sin(\omega_I t+\phi_I)/\omega_I,
\quad g(t)=\exp[-2\ln2(t/T_I)^2],
\]

and its exact derivative is used for `E_I`. The scale `C` is chosen so that the
sampled maximum of `E_I(t)^2` matches the requested peak intensity.
`ir_fwhm_fs` is the FWHM of the **squared vector-potential envelope**. The
derivative adds an envelope-derivative term to the electric field, so an exact
electric-field-envelope FWHM of 5 fs is not claimed. This distinction becomes
more relevant for very short IR pulses. Both `A` and `E` vanish in the tails
without subtracting an arbitrary electric-field DC component.

| Configuration field | Default | Unit / interpretation |
| --- | --- | --- |
| `dt_fs` | 0.005 | Time step, fs (5 as) |
| `time_half_width_fs` | 20 | Minimum half-window, fs |
| `delay_fs` | `-3:0.25:3` | Increasing XUV-delay samples, fs |
| `energy_eV` | `15:0.5:65` | Increasing final kinetic-energy samples, eV |
| `ionization_potential_eV` | 15.75962 | Ar threshold; constant dipole remains unchanged |
| `xuv_photon_energy_eV` | 50 | Carrier photon energy, eV |
| `xuv_fwhm_fs` | 0.25 | XUV intensity-envelope FWHM, fs |
| `xuv_cep_rad` | 0 | XUV carrier-envelope phase, radians |
| `ir_wavelength_nm` | 800 | IR wavelength, nm |
| `ir_fwhm_fs` | 5 | Squared vector-potential-envelope FWHM, fs |
| `ir_intensity_Wcm2` | `1e12` | Actual sampled peak electric-field intensity |
| `ir_cep_rad` | 0 | IR vector-potential sine phase, radians |

All overrides must be finite real MATLAB doubles. Axes can be scalar or
strictly increasing vectors. Unknown names are rejected to catch misspellings.
Negative energies or intensity, nonpositive widths, short windows and severely
undersampled fields/action are rejected. The sampling guard is a preliminary
check, not proof of convergence for every user-supplied parameter combination.

The time grid is symmetric, contains zero and has spacing exactly `dt_fs`.
The half-window is rounded **up** to an integer number of steps. The returned
`time_fs` records the actual grid. At least four envelope widths are required
around the IR and beyond the most delayed XUV pulse. The default window has
8001 time samples.

The minimal interface does not accept GDD, TOD, pulse trains, custom field
files or tabulated dipoles. These options in the original script are not
silently approximated by the maintained example. To study a new regime, extend
the model deliberately and add a physics control before interpreting results.

## 4. Implementation and reproducible use

`sfexis_defaults.m` returns the reference parameters. `run_sfexis.m` merges
overrides, validates dimensions and units, constructs the fields and calls
`sfexis_spectrum.m`. The kernel loops over final energy; for each energy it
forms the IR phase integral and integrates all delayed XUV fields. It does not
allocate a delay × time × energy tensor. Dominant storage is the time × delay
field matrix and temporary integrands; runtime scales as
`O(N_energy*N_delay*N_time)`.

The kernel accepts nonuniform, increasing time grids and complex XUV fields
with shape `N_time × N_delay`. It requires a real IR potential, increasing
nonnegative energies and positive `Ip`. For this lower-level interface, the
caller is responsible for pulse tails, physical units and sampling; the
higher-level pulse-specific guards in `run_sfexis` are not applied.

The output includes the full effective configuration, atomic-unit constants,
axes, IR fields, complex amplitude, raw signal, globally normalized signal and
citation. `minimal_streaking` adds the MATLAB version and saves one `result`
structure. The simulation function itself has no filesystem, plotting or
working-directory side effects. The example closes only the figure it creates
and refuses to replace an existing MAT or PNG output.

For a field-free control:

```matlab
addpath('SFEXIS');
r = run_sfexis(struct('ir_intensity_Wcm2',0));
plot(r.energy_eV,r.normalized_signal(1,:));
xlabel('Electron kinetic energy (eV)'); ylabel('Normalized signal');
```

The zero-IR spectrum is delay invariant and follows the analytic Fourier
transform of the real Gaussian carrier, multiplied by `p`. Its maximum is
near `photon energy - Ip`; the `p` weighting shifts the maximum slightly, so the
test compares the whole analytic spectrum rather than demanding an exact peak
at that difference.

Before a production calculation, reduce `dt_fs`, enlarge the time window,
refine the energy/delay sampling and widen the energy interval while tracking
the observable of interest. The provided convergence tests establish one
small reference regime. They do not establish arbitrary intensity, bandwidth
or near-threshold accuracy. Independent atomic models and experiment are
needed for physical validation beyond this SFA implementation.

## 5. Complete source and dependency audit

### MATLAB source files (10)

| Archived file | Role and inspected dependencies | Status / material issues |
| --- | --- | --- |
| `SFEXIS.m` | Main script: fields, atomic choice, optional dipoles, action integral, plots and saves | Interactive `pause`, `clear all`, `close all`; Windows paths; very large unused array; see below |
| `backup30Ottobre.m` | Earlier main-script variant | Different defaults, extra factor of 2 in IR field amplitude, mixed dipole data and FFT-based shifting; not an interchangeable reference |
| `EField1D.m` | Gaussian vector-potential spectrum, inverse FFT, finite-difference electric field | Base MATLAB; frequency grid uses `(N-1)*dt`, not the DFT's `N*dt`; duration relates to the constructed potential; positive-peak normalization |
| `EField1D_train.m` | Half-cycle pulse train | Authors F. Vismarra and Yingxuan Wu; calls Signal Processing Toolbox `envelope`; one branch plots and calls `toc`; odd-train shaping is marked approximate in source |
| `DipolestatCalculation1D.m` | Plane-wave dipole integral over a 3D hydrogenic grid | Declared function is named `DipoleCalculation1D`, unlike filename; graphical `waitbar`; returns z-polarized dipole |
| `dipoleShift.m` | Dense interpolation, integer circular shift and smoothing of dipole modulus/phase | Curve Fitting Toolbox `smooth`; positive shift moves samples toward larger p, corresponding to `d(p-A)` rather than its `d(p+A)` comment; edges need validation |
| `AssociatedLaguerre.m` | Polynomial helper for radial hydrogenic state | Includes an additional factorial normalization relative to the usual associated Laguerre polynomial; main script has a compensating factor; do not change one in isolation |
| `AssociatedLegendre.m` | Polynomial helper for angular state | Main script supplies the `(-1)^m` phase; general spherical-harmonic conventions are not newly validated |
| `plotPsi.m` | Inner/outer probability isosurfaces | MathWorks copyright preserved; missing `autoThreshold` for its two-argument form; documentation advertises an output but function declares none; explicit thresholds avoid that missing-helper branch |
| `bluewhitered.m` | Signed blue/white/red colormap | Base MATLAB graphics; not called by the main script; original file preserved |

The complete main script and all helpers were read. The backup was inspected
against the main script to capture all differences. Repository-wide call-site
searches found other, independent `EField1D.m` functions in `Electromagnetism`
and `HHG code/Trajectories Analysis`, with different interfaces; they are not
SFEXIS dependencies and were not moved. No repository `AGENTS.md`, package
manifest, vendored build system or external compiled dependency was present at
the inspected commit.

The MATLAB dependency analyzer confirms **MATLAB only** for the maintained
path and **MATLAB, Signal Processing Toolbox and Curve Fitting Toolbox** for
the original main script's discovered branches. A static dependency list does
not prove that those branches run successfully. Official function references:
[`envelope`](https://www.mathworks.com/help/signal/ref/envelope.html) and
[`smooth`](https://www.mathworks.com/help/curvefit/smooth.html).

### Data files and historical outputs (8)

These files were inventoried with MATLAB. All numeric arrays in the five
top-level MAT files are finite; names, shapes and stored ranges are below.
The energy-axis interpretation in a filename is not treated as validated
provenance.

| File | Stored content | Caution |
| --- | --- | --- |
| `Ar_DipoleT_Khe.mat` | 100×1 `EnergyAr` (1–82.945979), `CrossSectionAr`, `Phase_Ar` | Main custom branch treats the energy as kinetic energy, sets amplitude to one, divides phase by `2*pi`; original phase units/source must be established |
| `Ne_DipoleT_Khe.mat` | 1×100 `Energy` (1–177.25142), `CrossSection_Ne`, `Phase_Ne` | Phase array is entirely zero; not automatically used by the main script's custom branch |
| `Ar_DipoleT_Transition_Mixed.mat` | 1693×1 `Energy` (0.03–80.4), 1×1693 `Amplitude`, 1693×1 `Phase` | Backup multiplies phase by `pi`, a different convention from the main custom branch; provenance and units unresolved |
| `jet_white.mat` | 64×3 `jet_white` RGB colormap in [0,1] | Plotting asset only |
| `ePie_map.mat` | 64×3 `ePie_map` RGB colormap in [0,1] | Plotting asset only; no retrieval code |
| `Simulations/20-Mar-202116_36_29TEST1Ar/Wetot.mat` | `WE`, `fp`: 401×501; `delay`: 1×401; `E_fin`: 1×501 | Historical output lacks a complete saved parameter/environment record; not a newly verified benchmark |
| `Simulations/20-Mar-202116_36_29TEST1Ar/WeFINALE.fig` | MATLAB figure container | Preserved and container-inspected; not regenerated or used as a numerical oracle |
| `Simulations/20-Mar-202116_36_29TEST1Ar/XUVIRField.fig` | MATLAB figure container | Same preservation scope |

No data source or phase convention was inferred from the abbreviation `Khe`.
These tables are not needed by the maintained constant-dipole path.

### Why the legacy default is not the quick start

1. Its initial time axis spans 11,000 fs at 10 as spacing (1.1 million points),
   followed by cropped fields. The unused `EatomicReal` tensor can consume
   roughly 13 GB for a 401×8000×501 nominal grid before useful computation.
2. The hydrogenic dipole branch uses `UR` before defining it. The missing
   `autoThreshold` and helper name mismatch are additional issues, although
   the main call to `plotPsi` supplies explicit thresholds.
3. The delay loop shifts by a value computed from row index and array length,
   not directly from each declared delay; this can introduce an offset and
   relies on integer shifts. Cropping circularly shifted data also needs
   endpoint care.
4. The IR-window check tests a signed sum (`error < 1e-6`) rather than an
   absolute quadrature error or a tail criterion, and omits the time-step
   factor. It is not a convergence certificate.
5. The custom-field branches are placeholders. The custom-dipole branch always
   loads argon data despite the atomic selector, and source comments themselves
   flag phase conversions for checking. A filename does not resolve those units.
6. The backup scales the electric field by an extra factor of two, corresponding
   to four times the intensity for otherwise identical fields. It also has a
   no-dipole branch that references variables defined in another branch.
7. Working-directory changes, relative loads, graphical dialogs, figure saves
   and `clear all` make the scripts unsuitable as library calls or CI tests.

The maintained route removes these obstacles for the **constant-dipole
isolated-pulse case**. It does not assert that the optional legacy modes are
fixed. Restoring tabulated dipoles or quantitative atomic delays requires
source/units verification, a consistent dipole and field convention, and new
reference tests before publication use.

## 6. Validation and citation

[VALIDATION.md](VALIDATION.md) records the actual MATLAB results, convergence
errors, archive check and remaining limitations. The comparison to the legacy
formula tests the same model on shared fields, not an independent physical
theory. The analytic controls are the independent numerical references.

Keep [CITATION.md](../CITATION.md), the code headers and citation metadata with
derivative work. Use the article DOI for the article and the repository commit
for the exact code revision. The article citation should not be used to imply
that this compact example reproduces the paper's ensemble calculations.
