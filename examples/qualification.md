# Native model coverage and numerical settings

The readable `Case` interface has exercised all eight named native surface/water
codes in the small monochromatic matrix below. It uses the same `simulate`
Python/CLI route and durable records as ordinary calculations. These are
specified model connections and finite-resolution calculations; selecting a
model name does not establish accuracy across its parameter or wavelength range.

## Executed model matrix

Each case evaluates **551.25 nm**, with atmospheric aerosol model −1, FMF 0.5,
relative humidity 0.5, optical depth 0.125 at 532 nm, a 15 km Gaussian aerosol
layer with 2 km FWHM, US standard atmospheric profile, and surface pressure
1000 mb. Illumination is explicitly Thuillier. The calculation is elastic vector
RT with an atmosphere, plane-parallel geometry and gases disabled. The ocean
cases enable glint, retain 5 m/s wind, 20 °C water and 37 PSU salinity, and disable
Raman, both fluorescence sources and NPQ. Water-body cases use homogeneous
chlorophyll and 200 m depth. Existing whitecap and bottom closures are retained.

Operating settings are `ncolinput=8`, `nquadainput=12`, `nquadoinput=24`,
`MAXMORDINPUT=3`. They are editable choices, not maxima or accuracy thresholds.
The propagation zenith/relative-azimuth pairs are `(30,63)`, `(30,-63)` and
`(145,120)` degrees. The sun is at zenith 20°, except snow at 60° to exercise its
anisotropic table branch; the native snow model is Lambertian below 51.6° solar
zenith. AirSensor is at 2 km.

| Model | Code | Active model inputs beyond the common state |
| --- | ---: | --- |
| Ocean, no water body | 0 | Interface reflection with glint; no scattering water body |
| Chlorophyll | 1 | Chlorophyll 0.5 mg/m³, Bricaud absorption |
| Chlorophyll/sediment | 2 | Chlorophyll 0.5 mg/m³, Bricaud; particle indices 1.02 and 1.2; size exponents 3.5 and 4; sediment 0.2 g/m³; separate 665 nm mineral normalization |
| Seven-parameter water | 3 | Chlorophyll 0.5 mg/m³, `adg440=0.03 m⁻¹`, `bbp660=0.003 m⁻¹`, `Bp660=0.015`, `Sdg=0.018 nm⁻¹`, `Sbp=1`, `S_Bp=0.1` |
| Lambertian/Fresnel | −200 | Albedo 0.2, Lambertian fraction 0.7, real index 1.5, imaginary index 0 |
| MRPV | −201 | `a=0.08`, `k=0.8`, `b=-0.2`, polarized coefficient `e=0.3`, index `1.5+0i` |
| Snow | −202 | User albedo 0.9, Hudson EOF/BRDF tables at 551.25 nm and 60° sun |
| Ross–Li | −203 | `fiso=0.0399`, `fvol=0.0245`, `fgeo=0.0072`, `Bpol=0.1`, index `1.5+0i` |

The land fixtures are constant spectra over 300–800 nm. The Ross–Li coefficients
are the example coefficients cited to Lorente et al. (2018) in the public native
input examples. The constant spectrum is an illustrative input, not a measured
UV material spectrum. Replace it with the spectrum appropriate to the application.

All eight runs returned finite float64 Stokes I,Q,U,V and irradiances at TOA,
AirSensor and BOA, plus TOO immediately below the interface for ocean. Tests
check native model records, exact evaluated wavelength, dimensions, units,
nonnegative intensity/irradiance, gas optical-depth zeros and the zero downward
diffuse Stokes vector at TOA away from the solar beam. Signed Q/U/V are retained.
The seven non-snow models obey mirror I/Q-even and U/V-odd relations at these angles.
Snow's interpolation exception is described below.

## Reference basis

The retained model reference uses the original public `GLOBAL_DATA`,
`MIE_PHMX_PACE` and `RTSOSINIT` preparation at the pre-extraction fork revision
`9f1bcef2b0bc1d835c1a6a2c67eba6a2dae4d949`, linked to the retained original core.
A separate single-direction driver evaluates the physical wavelength itself,
with exact optical nodes at 551.25, 532 and 665 nm. It performs no instrument
averaging, PACE response processing or double-k approximation. The new Python
adapter and paired-direction/result writer do not generate expected values.
Both native paths consume the same serialized physical inputs, so this comparison
does not independently verify every named parameter's translation. The record
checks and separate surface-equation tests complement that comparison.
The reference values, settings, source identities and selected data hashes are
in `src/rtsos/tests/fixtures/model_reference.json`; input construction is in
`src/rtsos/tests/test_native_models.py`.

All checked Stokes components, detector irradiances, solar irradiance and column
Rayleigh/gas/aerosol optical depths agreed **exactly** at the primary direction
for all eight models in this build. This establishes connection parity with
shared equations and the same underlying solver. It is not independent physical
validation. Regression tolerances `rtol=2e-12, atol=1e-13` reuse the existing
grouped/single native connection allowance for float64 angle/summation roundoff;
they are not application accuracy targets. A single reference direction avoids
the legacy last-zenith attenuation scratch-fit issue, and uses no near-duplicate
cosines. Explicit mode retains its two corrections: exact duplicate removal and
direction-specific zeroth-Fourier attenuation fits. Mapped multi-angle legacy
results cannot be used as an exact reference without accounting for these.

Separate direct-surface tests evaluate the nontrivial MRPV and Ross–Li intensity
equations with the polarized term disabled and the atmospheric optical-depth
floor retained. RossThick/LiSparse-R kernels follow the
[published equations, Appendix B](https://amt.copernicus.org/articles/17/2235/2024/index.html).
MRPV uses the public native exponential phase law and the hotspot factor from
[Rahman et al. (1993)](https://doi.org/10.1029/93JD02072); the original RPV phase
law should not be substituted for the native modified law. Expected radiance is
`E0*cos(solar_zenith)*surface_reflectance/pi`. The existing `E0*2e-9` absolute
surface tolerance accounts for the native disabled-atmosphere Rayleigh floor
(10⁻¹² per layer), not empirical fitting to these outputs. Existing Lambertian/
black-ocean, Rayleigh U-sign, strong-absorption restoration and signed-weight
band-quadrature tests supply their original analytical checks.

Snow needs `SnowBRDF_Data/brdfPaperMaterial/{albedo,eof_A,...,eof_F}.txt` under
the auxiliary directory. The author's
[Hudson et al. supplement](https://www.atmos.washington.edu/~sgwgroup/DC/brdfPaper.html)
provides these tables. The snow test reports a scoped skip if these files are
missing; that skip is not evidence that snow executed. The measured EOF-A table
is exactly mirror symmetric, but the native `interpl2d -> Func_UVIP3P` local
three-node quadratic stencils at azimuths 63 and 297 are not reflected stencils.
The retained snow case has a TOA intensity difference of **0.00444095**
W m⁻² sr⁻¹ µm⁻¹ between those directions, on an intensity of about 254.6
(1.74×10⁻⁵ of I). Both directions agree with original public single-direction
calculations within float64 angular roundoff. Tests retain both reference
directions. This inherited interpolation asymmetry remains for review; no
symmetrization, clipping or replacement model has been introduced.

## Measured numerical sensitivity

Two representative elastic, gas-free cases use the same three directions and
common inputs as the matrix: chlorophyll ocean at **443.125 nm**, sun 20°, and
Ross–Li at **350 nm**, sun 30°. The base is `(ncol,nquad_a,nquad_o,m)=(8,12,24,3)`.
One setting changes at a time:

- Scattering orders: 4, 8, 12, 16.
- Atmospheric quadrature: 8, 12, 20.
- Ocean quadrature: 16, 24, 40, for the ocean case.
- Fourier upper bound: 1, 3, 6, 12.

The table reports the base setting against the highest value **in its own
scan**, with the other settings at base. Entries are maxima across the three
observations at TOA. Radiance differences and the reference maximum |I| are in
W m⁻² sr⁻¹ µm⁻¹. Q/U/V use that I scale, including when a component is zero.
These scan endpoints are comparison values, not a converged truth.

| Case | Change | max abs ΔI | max abs ΔQ | max abs ΔU | max abs ΔV | Reference max abs I |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Ocean 443.125 | Orders 8 → 16 | 0.0991864 | 0.00140669 | 0.00178984 | 0.0000318055 | 75.2025 |
| Ocean 443.125 | Atmosphere 12 → 20 | 0.240657 | 0.00895229 | 0.0365044 | 0.00105243 | 74.8626 |
| Ocean 443.125 | Ocean 24 → 40 | 0.147276 | 0.00745488 | 0.0157474 | 0.0000527334 | 75.2505 |
| Ocean 443.125 | Fourier 3 → 12 | 0.0211485 | 0.00221721 | 0.0280001 | 0.000206557 | 75.0821 |
| Ross–Li 350 | Orders 8 → 16 | 0.519267 | 0.0413155 | 0.00447890 | 0.0000368116 | 69.1564 |
| Ross–Li 350 | Atmosphere 12 → 20 | 0.131319 | 0.0107602 | 0.0244115 | 0.0000412967 | 68.7685 |
| Ross–Li 350 | Fourier 3 → 12 | 0.000755149 | 0.000914534 | 0.00135134 | 0.0000205770 | 68.6364 |

At TOA the 12 → 16 order intensity changes are 0.00717224 (ocean) and 0.0549628
(land); Fourier 6 → 12 changes are 0.0000224382 and 0.00000173816. Sensitivity
varies strongly by detector/direction. For atmospheric quadrature 12 → 20 in the
ocean case, the TOO intensity change is **42.9169**, with reference I=161.433,
at the downward propagation zenith 145°. It is 0.669067 at the two upward
directions, reference I=10.6168. TOA agreement alone does not qualify in-water
angular resolution. All detector/component absolute differences and
differences divided by their detector's reference max |I| were retained locally.

The solver executes its requested scattering-order/Fourier loops; it supplies
no application residual stopping criterion in this path. For `ncol<3` it
disables mapped transport and geometric-series acceleration. All order scans
above stay at or above 4, retaining those algorithms. Their last-order
geometric tail approximation remains active. Changing atmospheric quadrature
also changes the phase expansion degree `max(nquad_a,m)` and the coupled ocean
angular partition: that scan measures the combined native resolution setting.
The Fourier scan holds the phase degree at 12. Ocean quadrature leaves that
degree unchanged. The native requirement `nquad_o >= nquad_a` applies to ocean;
the difference should be even to avoid the core's documented internal adjustment.
All executed settings meet that constraint.

Finite angular/phase/Fourier resolution, the geometric tail, fixed optical-depth
steps, delta fit and absorption-depth reduction are distinct from physical
approximations in the empirical water/land models. Connection-reference parity
does not measure those errors or measurement/model uncertainty. No universal
numerical tolerance or application accuracy target is imposed here. PACE's
existing defaults and internal Fourier bound remain unchanged.

## UV-land example and reproduction

Edit the data paths in `uv_land.toml`, then run from this directory:

```sh
rtsos simulate uv_land.toml uv_land.nc
```

```python
from rtsos import load_case, simulate

case = load_case("uv_land.toml")
result = simulate(case, "uv_land.nc")
print(result.Wavelength_Evaluated.values)  # [350., 380.]
```

The example specifies micrometre wavelengths, a 300–800 nm Ross–Li spectrum,
Thuillier solar illumination, US standard atmosphere and gas columns, aerosol
properties, paired directions and elastic vector/plane-parallel processes.
350 and 380 nm are inside the active aerosol/solar/surface tables. Ozone, NO2
and the public O4 bands contribute according to native coverage; H2O/CO2/CH4
and the O2 line bands are outside their modeled ranges at these wavelengths.
Installed Python and actual CLI products agreed identically, with positive
column gas absorption at both wavelengths. Sensor name is a label and does
not change scene inputs or illumination policy.
Both UV wavelengths also agreed exactly, in all 14 checked fields at the
30°/63° direction, with the separately connected original public preparation.
That is an additional monochromatic connection check with shared physics.

This exercise adds UV land at those two wavelengths. It does not establish
UV water, all aerosol/humidity states, snow's entire parameterization, spherical
shells, scalar RT, inelastic processes or external phase matrices. Existing
single-case and response tests cover their stated domains. Complete supplied
response samples, tails, detector/mirror variants, masks and signed weights
remain intact. Full-response trapezoids integrate the supplied nodes; they do
not establish spectral convergence or turn a band into a monochromatic result.
The [PACE simulator theory](https://www.frontiersin.org/journals/remote-sensing/articles/10.3389/frsen.2022.840188/full)
and [native input manual](../Documentation/Describtion_of_input_files/SOS_io.tex)
describe the underlying physical models and native numerical controls.

With native executables on PATH and the two data locations set, opt in to the
permanent model regressions from the repository root:

```sh
export RTSOS_AUXILIARY=/path/to/Data
export RTSOS_GAS_ABSORPTION=/path/to/Gas_Absorption_Coefficients
PYTHONPATH=src python -m pytest src/rtsos/tests/test_native_models.py --native -q
```

Choose numerical settings for the actual detector, wavelength and application,
then compare higher resolutions against that application's accuracy needs.
