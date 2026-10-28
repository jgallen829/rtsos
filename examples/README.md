# Physical cases and exact response data

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
| `spectral` | Spectral calculation method and future requested wavelengths. |
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
complete-channel selection in this interface. Per-channel or detector/mirror
execution needs the next native connection. The PACE driver fixes OCI mirror
side to zero (the first source mirror).

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

## Exact response data and future requests

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

Case loading can also retain requests such as `spectral.method = "monochromatic"`
with `wavelength_unit = "nm"` and an ordered `wavelengths` list. Execution of
those requests, custom response integration, arbitrary viewing lists, scalar RT
and custom illumination require native connections. `simulate` reports the
specific missing operation before running. It never substitutes a band center,
nearest channel or sparse spectral calculation for a requested method.
