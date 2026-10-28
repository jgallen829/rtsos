# Physical cases, monochromatic evaluation and full-response bands

`pace_ocean.toml` and `pace_land.toml` each specify one physical state. Both select
solar zenith 20 degrees and upper aerosol optical depth 0.125; ocean chlorophyll
is 0.5 mg/m³. These are explicit example selections. The existing PACE initializer
continues to produce its original 20-state sweep and physical defaults.

Set `[data].root` to the directory containing `Data` and
`Gas_Absorption_Coefficients`, then run:

```sh
rtsos simulate pace_ocean.toml ocean.nc
rtsos simulate pace_land.toml land.nc
```

Each command retains the native input text beside the NetCDF result. Add
`--dry-run` to write those records without execution, or `--tmp /path/to/scratch`
to choose the temporary-file location. The native executable must be on `PATH`,
as with the legacy command. Numerics and viewing counts are ordinary adjustable
settings, without a separate workload ceiling. The small example settings do
not establish angular or spectral convergence.

## Case ownership and paths

Every case has eight plain TOML tables (including empty tables when needed):

| Table | Owns |
| --- | --- |
| `sensor` | Sensor identity, channel selection, response source and calibration variants. |
| `scene` | `atmosphere`, `aerosol`, `surface`, and ocean `water` properties. |
| `geometry` | Solar direction, detector altitude, viewing grid or future explicit directions. |
| `physics` | Scattering treatment, spherical approximation and process choices. |
| `numerics` | Numerical resolution and expansion settings. |
| `illumination` | Solar reference policy. |
| `spectral` | Spectral calculation method and requested monochromatic wavelengths. |
| `data` | Auxiliary and gas absorption data locations. |

Relative direct paths use the case file's directory, independent of the process
working directory. This applies to response files, land spectra, aerosol phase
matrices and data locations. Paths may contain spaces. The native text format's
restrictions on `#`, newlines, trailing whitespace and fixed character buffers
still apply. `scene.atmosphere.profile` is a filename under
`gas_absorption/atmmod`, such as `afglus.dat`, rather than a direct path.

Instead of `data.root`, specify both `data.auxiliary` (the `Data` directory) and
`data.gas_absorption` (the `Gas_Absorption_Coefficients` directory). The two forms
are mutually exclusive. Legacy PACE commands still resolve their direct paths
against the process working directory.

Python uses the identical implementation:

```python
from rtsos import Case, load_case, simulate

case = load_case("pace_ocean.toml")
case.scene["water"]["chla"] = 0.7
result = simulate(case, "ocean.nc", tmp="/path/to/scratch")

# For a Python dictionary with the same eight tables:
case = Case.from_mapping(tables, base_dir="/path/to/case_directory")
```

Optional scalar values come from the existing PACE initializer. Solar zenith,
upper aerosol optical depth, and chlorophyll for water-body models must be
explicit; they have no implied unique value from the legacy sweep. Land and
surface-only ocean do not ask for chlorophyll. Unknown or irrelevant settings
are rejected before serialization, so misspelled settings cannot be ignored.

## Named surfaces and water bodies

| `scene.surface.model` | `scene.water.model` | Existing native code |
| --- | --- | --- |
| `ross_li` | Omit water | -203 |
| `snow` | Omit water | -202 |
| `mrpv` | Omit water | -201 |
| `lambertian_fresnel` | Omit water | -200 |
| `ocean` | `none` | 0 |
| `ocean` | `chlorophyll` | 1 |
| `ocean` | `chlorophyll_sediment` | 2 |
| `ocean` | `seven_parameter` | 3 |

Ocean surfaces accept `wind_speed_m_s`. `none` retains the ocean interface,
including surface reflection; the native driver disables water Raman scattering
and fluorescence. It is not an atmosphere-only scene. Codes 1 and 2 use distinct
native particle-optics parameterizations; code 2 also includes sediment.

Water-body models accept `chla` (mg/m³), `water_depth_max` (m) and
`CHLA_HOMOGENEITY` (native 0/1). Codes 1 and 2 also accept `absorption_model`
(`bricaud` or `ciotti`). Case 2 additionally accepts
`phytoplankton_index_refraction`, `phytoplankton_spectral_slope`,
`sediment_index_refraction`, `sediment_spectral_slope`, and
`sediment_concentration` in the native mass concentration convention (g/m³).
Indices and particle-size spectral exponents are dimensionless.

Case 3 requires all seven optical parameters explicitly:

| Parameter | Meaning | Unit |
| --- | --- | --- |
| `chla` | Chlorophyll concentration | mg/m³ |
| `adg440` | Dissolved/detrital absorption at 440 nm | m⁻¹ |
| `bbp660_BackscatterCoeff` | Particle backscatter at 660 nm | m⁻¹ |
| `Bp660_BackscatterFraction` | Backscatter fraction at 660 nm | 1 |
| `Sdg` | Exponential absorption slope | nm⁻¹ |
| `Sbp` | Power-law backscatter exponent | 1 |
| `S_Bp` | Power-law backscatter-fraction exponent | 1 |

The native relations use `(wavelength_nm/660)**(-Sbp)` and
`(wavelength_nm/660)**(-S_Bp)`, so these two exponents are dimensionless. The
adapter retains the seven-parameter model's Bricaud absorption choice.

Land surfaces take `spectrum`, using the established native format: a header, a
wavelength count, then seven numeric columns. The model determines their meaning:

| Model | Columns after wavelength in nm |
| --- | --- |
| `ross_li` | `fiso fvol fgeo Bpol refin_re refin_im` |
| `snow` | `albedo blank blank blank blank blank` |
| `mrpv` | `pBRDFa pBRDFk pBRDFb pBRDFe refin_re refin_im` |
| `lambertian_fresnel` | `albedo fraction_lambertian blank blank refin_re refin_im` |

Ross–Li includes the polarized-reflectance coefficient `Bpol`. Snow uses the
native snow BRDF tables. A Lambertian fraction of one gives pure Lambertian
reflection, as in `land_spectrum.txt`. Fill unused native columns with zero.
No water-body properties belong in a land case.

## Connected PACE settings

`sensor.name = "pace"` connects these complete native channel groups:
`all`, `oci_swir`, `oci_spexone`, and `harp2` (native selectors 0, 2, 3, 4).
Their order and wavelength arrays remain the native output order. The legacy
selector 1 remains available in the legacy command; its output includes retained
SWIR slots even though it skips those integrations, so it is not presented as a
complete-channel selection in this interface. Individual complete OCI channels and detector/mirror variants can instead use
`response_integration` below. The PACE driver retains its fixed OCI mirror side
zero (the first source mirror).

`geometry.mode = "pace_grid"` takes `solar_zenith_deg` in [0, 90),
`sensor_height_km`, `zenith_count` and `azimuth_count`. Counts of one give zero
viewing zenith/azimuth. Even zenith counts generate the existing paired angles
from zero to 180 degrees; odd counts above one are rejected because the native
routine leaves the middle angle unassigned. Azimuths span zero to 180 degrees
in the native solar-plane convention. Native `ThetaV` and `PhiV` coordinates
retain their existing angle definitions. Detector heights that the driver would
silently reset to 2 km are rejected by the new adapter. TOA, interface and the
other native detector products remain available as in the legacy output.

`scene.atmosphere` accepts the column and pressure names in the examples plus
`profile`. `scene.aerosol` accepts the example fields, the equivalent lower-layer
fields and direct `phase_matrix_hi`/`phase_matrix_low` paths when needed by its
native `Aerosol_Model`. `wv_pace_ref` is in nm; optical depths and fractions are
dimensionless; particle centroid heights are in km. Despite the historical
`height_particle_variance_*` names, those values are Gaussian FWHM in km in the
native formula, which squares them in the denominator.

`physics.scattering = "vector"` is the connected PACE treatment.
`spherical_shell` enables both native spherical corrections together.
`atmosphere` and `gas_absorption` are Boolean controls. Ocean surfaces also take
`sun_glint`. Water bodies also take `raman`, `chlorophyll_fluorescence`,
`cdom_fluorescence`, `nonphotochemical_quenching` and
`single_water_phase_matrix`. All these switches are TOML true/false. Fixed native
choices include water temperature 20 °C, salinity 37 PSU, whitecaps enabled and
the existing atmosphere/particle and ocean preparation. Sensor selection does
not change them.

`numerics` accepts `ncolinput`, `nquadainput`, `nquadoinput` and `MAXMORDINPUT`.
These retain their native collocation, atmospheric/ocean quadrature and Fourier
order meanings. Integer inputs must be representable by the native int32 format;
there is no additional fixed limit on the operating settings.

The legacy PACE driver passes `min(30, nquadainput, MAXMORDINPUT)` as the
Fourier-order upper bound to its main per-wavelength RT calls. The requested `MAXMORDINPUT` remains
unchanged in the input coordinates; `pace_rt_fourier_order_limit` records this
native bound in the resolved input attributes and result. It is an upper bound,
not a claim that convergence required every order. Raising the requested value
above 30 does not raise this bound in the current driver. This existing native
restriction remains visible while shared native preparation is developed.

`illumination.source = "pace_native"` explicitly selects the native policy:
use `tsis1_hsrs_binned_fs.csv` when present, otherwise `f0.txt` (Thuillier).
`tsis` and `thuillier` are also accepted when the data directory's native
precedence actually selects that source. A mismatch fails before execution.
Results record the selected file and policy. No campaign-specific solar policy
is imposed on every sensor.

`spectral.method = "pace_native"` retains PACE's response processing and
double-k behavior. Disabling gases or the atmosphere retains the native flag
reset to its compiled-grid monochromatic processing; the effective flag stays
in the output. This does not mean evaluation at arbitrary requested wavelengths.
All existing Stokes, irradiance, wavelength, and detector variables are retained.

## Exact response data and selection

`response.toml` demonstrates unequal-length channels, signed tails, zeros,
separate nominal labels, explicit wavelength units and mirror variants.
`custom_response.toml` shows how to reference and select these data:

```python
from rtsos import load_case
from rtsos.responses import load_oci_rsr, load_responses

selected = load_case("custom_response.toml").load_responses()
blue = load_responses("response.toml").select(["blue"], mirror=1)
oci = load_oci_rsr("/path/to/PACE_OCI_L1B_LUT_RSR_baseline_1.1.1.nc")
channels = oci.select(["oci_red_001", "oci_blue_018"], mirror=0)
```

The OCI loader reads every supplied band and mirror in blue/red/SWIR source
order. IDs are one-based within detectors; mirror indices are zero-based.
Selecting IDs reorders complete channels in requested ID order while preserving
source variant order within each ID. Missing variants raise an error. The loader
keeps calibrated center wavelengths separate from samples and nominal labels;
it does not derive an effective wavelength. The source arrays retain all samples
and their dtype and masks, including signed weights, tails and zeros. A missing
calibration center is `None`, rather than the file's fill number. Neither loader nor
selection trims, clips, normalizes, resamples or interpolates calibration data.
Wavelength units are explicit `nm` or `um`.

## Requested monochromatic wavelengths

`monochromatic_ocean.toml` and `monochromatic_land.toml` evaluate **551.25 nm and
443.125 nm**, in that order. The land example supplies the same wavelengths in
micrometres. Set `data.root` (or both direct data paths) and run:

```sh
rtsos simulate monochromatic_ocean.toml ocean.nc
rtsos simulate monochromatic_land.toml land.nc
```

```python
from rtsos import load_case, simulate

case = load_case("monochromatic_ocean.toml")
case.spectral["wavelengths"] = [551.25, 443.125, 551.25]
result = simulate(case, "ocean.nc")
print(result.Wavelength_Evaluated.values)
```

The explicit `spectral.wavelength_unit` is `nm` or `um`; wavelengths are positive
finite numbers. The `request` coordinate is a zero-based identity, so repeated
wavelengths remain separate requests and the supplied order is preserved.
`wavelength_requested` retains the original numbers/units, `wavelength_nm` is the
conversion passed to native code, and `Wavelength_Evaluated` records the native
value. No PACE group, instrument response file, band center, interval average or
response integration is involved. `[sensor]` can be empty or contain a `name`
used only as a result label. Response files remain independently loadable data.

The native process prepares Mie optics at the actual requested wavelength and at
the separate `wv_pace_ref` aerosol optical-depth reference. Case 2 also prepares
its mineral optics at the model's 665 nm normalization wavelength. These physical
reference wavelengths do not replace the evaluated wavelength. The original
public atmosphere, Rayleigh, built-in aerosol, land, water, phase expansion and
RTSOS solver routines are shared with PACE. Global state is isolated in one
process per scene's wavelength set, through the existing `Program` mechanics.
Atmospheric profiles, gas LUT reads/searches, solar, land/water tables, viewing
grid and aerosol profile are prepared together. Optical storage holds the
current wavelength and physical reference slots; reference cross sections are
reused. Each new wavelength still receives its own Mie, Rayleigh, gas-layer,
phase-mixing and RT calculation. Exact common wavelengths are solved once, then
mapped back to every request in supplied order. Nearby nodes are never merged.
Gas profiles and radiometric output storage scale with unique wavelength count;
Stokes output storage also scales with the requested viewing grid. There is no
hidden count ceiling or automatic thinning.
Native wavelength counts/indices use the compiler's default integer
representation (int32 in the supplied build), as distinct from an operating cap.

All eight named surface/water choices above retain their public equations and
parameters. The ocean example uses the existing chlorophyll model, 20 °C water,
37 PSU salinity, whitecaps, atmosphere, gas absorption, and glint. It explicitly
disables Raman, chlorophyll/CDOM fluorescence and NPQ. The monochromatic route
requires those inelastic excitation/emission connections before they can be
enabled. Requests for them fail before starting the native executable; omitted
water-body process settings retain the existing defaults and consequently also
require an explicit choice. External aerosol/cloud phase-matrix files still use
PACE-segment layouts and need a requested-wavelength matrix loader. Built-in
`Aerosol_Model=-1` or 1–20 are connected. Vector RT and the reviewed `pace_grid`
viewing-count behavior are supported. Explicit direction lists, scalar RT,
and custom illumination retain their explicit capability boundaries.

Select `illumination.source = "tsis"` or `"thuillier"` explicitly. Native code
reads exactly that source, independent of sensor identity or availability of the
other source. TSIS values are converted from W m⁻² nm⁻¹ to W m⁻² µm⁻¹ by 1000;
Thuillier values from mW cm⁻² µm⁻¹ by 10. The existing `FSUNCAL` Akima cubic
interpolator evaluates solar irradiance at the requested wavelength. Elastic
unit-beam results are scaled by `E0 * abs(mu_solar) / native_downwelling_TOA`.
Stokes outputs have units W m⁻² sr⁻¹ µm⁻¹, irradiances W m⁻² µm⁻¹. Stokes I,Q,U,V
and meridian/solar-plane conventions retain native signs. TOA, AirSensor, BOA,
and (ocean only) TOO label the original detector positions. All new native
scientific outputs and explicit real inputs retain double precision; PACE's
existing record precision and result schema are unchanged.

The monochromatic solver receives the requested `MAXMORDINPUT`, without the
legacy PACE `min(30,nquadainput,MAXMORDINPUT)` bound. Phase preparation uses
`max(nquadainput,MAXMORDINPUT)` because spherical harmonics require l ≥ m.
The compiled core stores phase degrees through `MAXLFULL=3000`; exceeding that
representation requires rebuilding the core. There is no wavelength-count,
request-count or workload ceiling. `rt_fourier_order_limit` and
`phase_expansion_order` report the actual bounds; they do not establish
convergence. Existing collocation, quadrature, optical-depth integration steps
(0.05 atmosphere, 0.1 ocean), delta fit and atmospheric/ocean absorption-depth
reduction are retained. As in PACE, the atmospheric absorption removed during
preparation is restored to non-TOA detector products using the exponential of
the preserved column depth. This retains the native approximation.
`physics.atmosphere=false` retains the public 10⁻¹² Rayleigh optical
depth per atmospheric layer and requires `gas_absorption=false`.

### Model and data coverage

Coverage follows each active model's actual tables, rather than the core's
broader UV-to-infrared capability or a sensor's band range. The adapter checks
active supplied table ranges before native execution and never adds extrapolation
or nearest-grid replacement. Data need not include any instrument responses for
monochromatic evaluation.

| Quantity | Established preparation and coverage |
| --- | --- |
| Built-in aerosol refractive indices | Existing SF component tables, 300–2500 nm, including the optical-depth reference; native wavelength/humidity interpolation and mixing, then direct Mie calculations. Public native RH behavior is retained. |
| Land surface spectra | Existing local polynomial interpolation of the supplied seven-column spectrum. Requested wavelengths must lie inside its actual wavelength range. Snow additionally needs the public snow BRDF tables. |
| Pure-water absorption | IOCCG 2018 table through 550 nm and `water_coef.txt` above 550 nm, using existing local quadratic interpolation; actual file coverage is checked. The Lee and Mason files remain inputs of the public reader. |
| Phytoplankton absorption | BIOSOPE UV median shape below 400 nm, normalized at 440 nm; Bricaud (400–700 nm) or Ciotti's supplied pico/micro table in the visible; the public model sets absorption to zero at/above 700 nm. |
| Other water optics | Public chlorophyll/sediment/seven-parameter equations, pure-water power-law scattering, and Quan–Fry refractive index through 1100 nm with the public tabulation above it. Salinity 37 PSU retains the public default; the Quan–Fry comment gives an empirical calibration range below 35 PSU. |
| Solar irradiance | Requested wavelengths must lie within the selected source file's actual nodes. `FSUNCAL` uses its existing cubic interpolation; this changes physical illumination, not response calibration. |
| Gas absorption | Each native species search interpolates the source cross sections directly at the requested wavelength and the atmospheric pressure/temperature/VMR, then integrates layer optical depth. It does not snap to the 37,618-node PACE working grid. |

Native gas spectral interpolation is linear between source LUT nodes (and uses
the existing pressure, temperature and H2O-VMR interpolation). Height integration
is trapezoidal in extinction coefficient; the existing local quadratic
interpolation maps cumulative optical depth onto the RT layers. The public
species cutoffs remain: H2O 396–2300 nm, CO2 781–2300 nm, CH4 1084–2300 nm,
O3 116.65–850 nm, NO2 241.10–759.9 nm. O2 is included inside its seven native
open intervals: 627–637, 686–700, 758–780, 859–872, 1053–1084, 1237–1308 and
1551–1603 nm. The native O4 calculation is added to O2: nine empirical Gaussian bands
centered at 343.4, 360.5, 380.2, 446.7, 477.3, 532.2, 577.2, 630 and 1065.2 nm,
each truncated at twice its tabulated FWHM. These public equations are retained. Outside a species' modeled range its contribution
is zero by the established model, rather than evidence of a measured zero
absorption cross section. Native ozone's negative temperature-fit values are
set to zero by that existing model. These cutoffs, source precision and tabular
interpolation limit gas fidelity; a monochromatic result is not a line-by-line
calculation beyond its supplied data resolution.

The existing water models can set a negative inferred particle-scattering
coefficient to zero, and retain their empirical parameterizations and bottom
closure. Those established model choices are unchanged. Arbitrary wavelengths
do not qualify every model, geometry, parameter range or numerical resolution.
Exact response loading continues to preserve signed samples, tails, masks,
zeros and native calibration precision independently of the selected spectral
calculation method. Full-response integration is described below.


## Full-response band simulation

`response_ocean.toml` and `response_land.toml` run the same physical scene at
**every supplied node** of the selected complete responses. They select red,
then both supplied blue mirror variants from `response.toml`. Grids have unequal
lengths and use both nm and um; the blue signed tail and zero endpoints take part
in the calculation. `custom_response.toml` is a smaller land configuration that
selects mirror zero. Set the data paths and run, from this directory:

```sh
rtsos simulate response_ocean.toml ocean_bands.nc
rtsos simulate response_land.toml land_bands.nc
rtsos simulate custom_response.toml custom_bands.nc
```

```python
from rtsos import load_case, simulate

case = load_case("response_ocean.toml")
case.sensor["channels"] = ["blue", "red"]  # Complete channels, in this order.
result = simulate(case, "ocean_bands.nc")
print(result[["Solar_Irradiance", "Stokes_Vector_TOA"]])
print(result.channel_id.values, result.detector.values, result.mirror.values)
```

The source representation and selection contracts above continue to apply.
`[spectral]` contains only `method = "response_integration"`; `[sensor]` owns
`response_file`, `response_format` (`toml` or `oci`), optional `channels` ID list,
`detector` and `mirror`. Omit `channels` for all supplied complete channels; omit
variant filters to keep all supplied variants. Requested ID order is retained,
with source variant order inside each ID. No absent calibration is interpolated.
Illumination must explicitly select `tsis` or `thuillier`. This route uses the
same supported physical models/process checks and `pace_grid` geometry as the
monochromatic route. Unsupported requested processes fail before execution.

### Integration and data contract

Let x_j be the original sample wavelengths converted to **micrometres**, R_j
the signed source responses, and X_j the evaluated radiance or irradiance per
micrometre. Define composite trapezoids on the sampled product:

```text
T(Y) = sum_j (x[j+1] - x[j]) * (Y[j] + Y[j+1]) / 2
A    = T(R)
band_X = T(R * X) / A
```

This is a normalized spectral-density band product: Stokes radiance remains
W m^-2 sr^-1 um^-1, and solar/detector irradiance remains W m^-2 um^-1.
The response is a relative dimensionless calibration. Normalizing this equation
does not normalize or otherwise change its stored source samples. A negative
signed area is permitted; negative contributions remain signed. A finite zero
band result remains zero. Both products use the same wavelength measure. For
native illumination conversions and unit-beam scaling, see the monochromatic
section above. Solar irradiance is integrated as its own product. No
monochromatic reflectance ratios are averaged, and no derived reflectance is
published. A reflectance subsequently calculated from these products must use
consistently integrated illumination and the applicable geometry definition.

This quadrature uses every supplied response sample, including zero endpoints,
signed tails and overlapping channel nodes. It integrates the *sampled product*
with trapezoids; it does not analytically integrate the product of interpolated
response and interpolated spectra. Constant spectral cases are exact up to
floating-point roundoff; accuracy on variable spectra follows the supplied grid
and this quadrature. Full-sample execution does not establish spectral or angular
convergence and is distinct from PACE double-k processing.

Execution requires at least two finite, positive, strictly increasing wavelengths
and finite unmasked responses in every selected variant. Masked/missing samples,
reversed/duplicate wavelength nodes and zero signed area are rejected without
sorting, filling, clipping or trimming a source calibration. Signed area whose
magnitude is no larger than float64 epsilon times T(abs(R)) is numerically
degenerate and rejected. Data loading still retains such arrays and their masks;
they can be inspected without execution. All nodes must satisfy actual active
model/solar coverage, even where response is zero. There is no out-of-coverage
endpoint exemption. Nonfinite native radiometric values prevent publishing the
affected band product. Independent sweep recovery remains separate work.

### Products and OCI execution

`channel` is a zero-based complete-variant identity, with `channel_id`, `detector`,
`mirror` (text, empty if absent) and `nominal_label`. Calibration centers and their
per-channel `calibration_wavelength_unit` remain separate from the unique
`Wavelength_Evaluated(spectral_node)` array. Missing centers remain missing. No
derived effective wavelength is supplied. `response_area(channel)` is signed,
in um. `Stokes_Vector_*`, `Irradiance_*` and `Solar_Irradiance` are band products;
`Node_` versions keep evaluated spectra. `Tau_*` stay monochromatic diagnostics
on `spectral_node`, with no claim of equivalent physical band optical depth.
Stokes, angle and detector conventions are unchanged. The NetCDF also contains
`calibration/0`, `calibration/1`, ... groups in channel order, each holding exact
source `wavelength` and `response` arrays in their original dtypes and explicit
`wavelength_mask`/`response_mask` arrays. Open those groups with netCDF4 or
`xr.open_dataset(path, group="calibration/0")`.

For a supplied OCI LUT, change the sensor table in either response example:

```toml
[sensor]
name = "oci"
response_format = "oci"
response_file = "/path/to/PACE_OCI_L1B_LUT_RSR_baseline_1.1.1.nc"
channels = ["oci_red_001", "oci_blue_018"]
mirror = 0
```

Use `rtsos simulate your_oci_case.toml oci_bands.nc` or `simulate` as above.
Omit `mirror` to execute every supplied mirror for these IDs. Complete OCI
calibrations may share thousands of nodes and span far beyond a nominal band.
All retained tails are evaluated; all actual coverage requirements still apply.
For example, the baseline 1.1.1 LUT inspected here supplies 20,926 nodes spanning
306–2398.5 nm for every variant. The example land spectrum ends at 2260 nm;
replace it with a supplied spectrum covering the full grid before using that
LUT for land. The adapter rejects this coverage gap rather than extrapolating.
Common wavelengths are evaluated once across channels/mirrors, without dropping
any source sample or variant. This interface does not apply the legacy PACE
response thresholds, band segmentation, mirror default, or double-k method.
Runtime still depends on unique-node count, scene, gas absorption and the
viewing/numerical settings; no dense-OCI performance or convergence guarantee
follows from the small examples.
