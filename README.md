<!--
RTSOS — Radiative Transfer model based on Successive Orders of Scattering
Copyright © 2025 Pengwang Zhai.

Licensed under the Creative Commons Attribution–NonCommercial 4.0
International License (CC BY-NC 4.0).
You may use, modify, and share this code for research and
educational purposes with proper attribution.
Commercial use requires written permission from the author.

Full license: https://creativecommons.org/licenses/by-nc/4.0/
Contact: Pengwang Zhai  |  [pwzhai@gmail.com]
-->

# Radiative Transfer model based on Successive Orders of Scattering (RTSOS)

RTSOS can solve the multiple scattering radiative transfer equation from UV, visible, to
infrared. It can handle atmosphere-land or atmosphere-ocean coupled systems. The
atmosphere can be a mixture of molecules, aerosols, and cloud droplets. The land bottom
can be Lambertian, snow surface, Ross-Li, and a number of other surfaces. The ocean
waters are modeled by a mixture of pure ocean water, phytoplankton, and colored
dissolved organic matter (CDOM), and other hydrosols. The sensors can be placed at
arbitrary levels in the Earth system. The output of the sensor can include the full
polarized Stokes parameters (I, Q, U, V). For more information, see
[References](#references) at the end of this document.

PACE simulator is a wrapper built around the monochromatic RTSOS, which has a list of
built-in aerosol and ocean inherent optical properties. A publication on the PACE
simulator is published on Frontiers in Remote Sensing
(Zhai et al., 2022).

GSFC AC LUT is a wrapper built around the monochromatic RTSOS which builds look up
tables for atmospheric correction for retrievals of surface reflectance from
top-of-atmosphere measurements.

## 1. Installation

Installing the `rtsos` Python package builds both the public Python/CLI interface
and five native commands: `rtsos_Monochromatic`, `rtsos_PACE_Simulator_DoubleK`,
`rtsos_GSFC_AC_LUT`, `rtsos_Aerosol_Phmx_Cal` and `rtsos_rao_dg`.

### 1.1 Source installation

Build prerequisites are [Python] >= 3.12 (CI matrix: 3.12–3.14), C/C++ and
[Fortran] compilers, [HDF5] >= 1.10 with its Fortran and high-level Fortran
libraries/modules, and [LAPACK] >= 3.9 with BLAS. [CMake] >= 3.22 and Make are
required; pip's isolated build installs CMake when needed through
scikit-build-core. Use compatible compiler and HDF5 Fortran toolchains.
For Ubuntu 24.04:

```sh
sudo apt-get update
sudo apt-get install --yes build-essential gfortran cmake libhdf5-dev libblas-dev liblapack-dev python3-venv
git clone https://github.com/jgallen829/rtsos.git
cd rtsos
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
rtsos --help
```

Retain Git history and release tags for setuptools-scm version discovery;
`git fetch --tags` repairs an existing tagless clone. A built source distribution
also carries its version and can be installed with
`python -m pip install /path/to/rtsos-<version>.tar.gz` without Git or
pretend-version overrides. GitHub's automatic source archives do not carry that
generated package metadata; use a Git clone or a built sdist for Python
installation. The [upstream releases][latest release] remain available separately
from this fork's current interface.

With a system package manager or [Conda], provide the same native dependencies
and compilers. If CMake cannot discover HDF5, set `HDF5_ROOT` to the prefix
containing its `lib` and `include` directories. Source installation uses the
existing scikit-build-core/setuptools-scm/CMake backend.

Some macOS HDF5 installations expose the Fortran library but omit the high-level
Fortran library location during CMake discovery. If CMake reports that it cannot
locate `hdf5::hdf5_hl_fortran`, provide the installed library explicitly, using
your HDF5 prefix (and its actual library filename):

```shell
export HDF5_ROOT="$CONDA_PREFIX"
CMAKE_ARGS="-DHDF5_Fortran_HL_LIBRARY=$HDF5_ROOT/lib/libhdf5_hl_fortran.dylib" python -m pip install .
```

The same setting can be passed as `-DHDF5_Fortran_HL_LIBRARY=...` to a direct
CMake configure. Use a compiler and HDF5 Fortran modules from compatible toolchains.

Native executables retain runtime dependencies on the build platform's HDF5,
BLAS/LAPACK and compiler libraries. Keep those compatible libraries installed
and discoverable by the platform loader. A locally built wheel is specific to
its platform and library environment; this build does not bundle or repair
native libraries into a universally portable wheel.

To build distribution artifacts:

```sh
python -m pip install build
python -m build --wheel --sdist
```

The wheel contains Python production modules and native commands. The sdist
also contains examples, regression tests/fixtures and the installed check below;
auxiliary datasets and local audit/build files are excluded.

### 1.2 Auxiliary Data

In addition to the software, simulating the PACE instruments and generating atmospheric
correction look-up tables requires a collection of data files available from a Zenodo
repository:

[RTSOS auxiliary data](https://doi.org/10.5281/zenodo.17410093).

The auxiliary data file is named as `RTSOS_data_v1.0_20251022.tar.gz`, which can be
extracted into a convenient location of your choice. You can extract the files by doing
the following in a Terminal with your chosen location as the working directory.

```shell
$ tar -zxf RTSOS_data_v1.0_20251022.tar.gz
```

You will get two directories. One is `Gas_Absorption_Coefficients`, which contains the
gas absorption cross section lookup table for H2O, O2, CH4, CO2, NO2, ozone. This lookup
table was generated by using the Atmospheric Radiative Transfer
Simulator([ARTS](https://www.radiativetransfer.org)) version 2.3 with the HITRAN2020 ([Gordon et al. 2022](https://doi.org/10.1016/j.jqsrt.2021.107949)) database.

The other is `Data`, which contains the absorption coefficients of pure water, plankton
particles, the bidirectional reflectance distribution functions of snow surface,
scattering matrix of dust particles, the instrument response functions of ocean color
instrument, HARP2, and SPEXone, and a number of other data files used in the PACE
simulator.

Arrange the extracted directories as follows, or use equivalent existing
directories through the two explicit paths:

```text
RTSOS_data/
├── Data/
└── Gas_Absorption_Coefficients/
    └── atmmod/
```

In each example's `[data]` table, set `root = "/absolute/path/to/RTSOS_data"`,
or replace `root` with both `auxiliary = "/path/to/Data"` and
`gas_absorption = "/path/to/Gas_Absorption_Coefficients"`. Relative paths are
resolved against the case file. The shipped `../RTSOS_data` root places these
directories at the checkout root when using cases in `examples/`.
The [coverage guide](examples/qualification.md) describes optional snow tables
and caller-supplied data for scientific checks.

### 1.3 Native executables only

With the same native prerequisites, direct CMake installation remains available:

```sh
cmake -S . -B build --install-prefix "$HOME/.local"
cmake --build build --parallel 2
cmake --install build
```

Add the installation's `bin` directory to `PATH`. The parallel count is an
adjustable build budget. This route installs native commands; use the Python
package installation above for `rtsos` and the public Python API.

### 1.4 Installed-package checks and CI

[GitHub Actions](.github/workflows/installed.yml) builds on Ubuntu 24.04 with
Python 3.12–3.14, then checks fresh non-editable wheel installations outside
the checkout. Python 3.12 also builds and installs the sdist without Git metadata.
The workflow uses system compiler/HDF5/BLAS/LAPACK dependencies and downloads
no auxiliary datasets. Its matrix parallelism and compiler parallel count are
editable resource budgets.

To check either installation yourself, unpack its sdist outside the checkout,
activate a fresh virtual environment, install the wheel or sdist plus pytest,
then run from another directory with source `PYTHONPATH` unset:

```sh
python -m pip install /path/to/rtsos-<version>.tar.gz pytest
unset PYTHONPATH
cd /path/outside/the/checkout
python -I /path/to/unpacked-sdist/validation/check_installed.py /path/to/unpacked-sdist /path/to/new-check-output
```

The check verifies installed module/command ownership, native startup/linking,
CLI help/templates/dry runs and the data-free regression suite, including
synthetic recovery/split-combine and exact response preservation. Native/data
tests retain their existing opt-in `--native` boundary. Build and API checks
are separate from [scientific model qualification](examples/qualification.md).

## 2. Usage

### Quick start: one explicit case

After installation and setting the data paths, run from the checkout or unpacked
sdist root:

```sh
rtsos simulate examples/single_direction.toml single.nc --dry-run
rtsos simulate examples/single_direction.toml single.nc
```

This is a Lambertian land scene at solar zenith 20°, aerosol optical depth
0.125 at 532 nm, with one propagation-zenith/relative-azimuth pair
`(37.25°, 315°)`. It evaluates 551.25 and 443.125 nm, supplied in micrometres.
The product retains the singleton `observation` dimension. A dry run writes
native input records without data access or native execution. Python calls
the same dispatcher:

```python
from rtsos import load_case, simulate

case = load_case("examples/single_direction.toml")
result = simulate(case, "single.nc")
print(result.Wavelength_Evaluated.values)
```

See [the examples guide](examples/README.md#requested-monochromatic-wavelengths)
for the wavelength, model, data, illumination, units and numerical contracts.
The [native coverage guide](examples/qualification.md) records the executed
eight-model matrix, measured numerical sensitivity, reference independence and
the UV-land example.
For full-response band products:

```sh
rtsos simulate examples/response_directions.toml direction_bands.nc
```

These are response-weighted products over every supplied calibration node,
whereas monochromatic products evaluate the requested physical wavelengths.
Full-response integration shares scene
preparation across the wavelength set. See [band simulation](examples/README.md#full-response-band-simulation)
for the trapezoidal method, complete OCI selections and result schema.
For paired observations at requested directions, including a single observation:

```shell
rtsos simulate examples/monochromatic_directions.toml directions.nc
```

Set the example data paths first. The same geometry table serves both spectral
methods, retains order and repetitions, and returns an `observation` axis.
See [ordered viewing directions](examples/README.md#ordered-viewing-directions)
for the propagation-angle, relative-azimuth, detector and Stokes conventions.
The existing low-level `rtsos_rao_dg` command remains available.

An ordered list of existing cases can run independently and resume:

```shell
$ rtsos simulate-cases examples/cases.toml results
$ rtsos simulate-cases examples/cases.toml results --resume
```

Each named case produces `results/<name>.nc`, including its metadata and exact
response calibration groups. Completed native groups survive a later failure;
later independent cases still run. `results/cases.json` reports every requested
case and its failure reason. The Python entry point is `simulate_cases`.
Single-case `simulate` and the existing sweep commands also accept `--resume`.
See [recoverable execution](examples/README.md#recoverable-case-and-sweep-execution)
for the record files, explicit retries, changed requests and combine behavior.

### 2.1A Binary Executables

For the low-level monochromatic text interface, see
[the native input manual](Documentation/Describtion_of_input_files/SOS_io.pdf).
The established PACE binary and NetCDF sweep interfaces below remain available;
the [Python wrapper](#21b-python-wrapper) supports templates and split/combine.

### 2.1A.1 PACE Simulator

#### Input Files 

The current version can conveniently generate synthetic datasets for OCI, HARP, and SPEX
for flexible atmospheric and ocean conditions given as parameters. To run the PACE
simulator from the binary executables, users should become familiar with the parameter
file formats. We provided two Python scripts to generate these input files, one for
ocean and the other for land surfaces, respectively, which are located in the
`test/PACE_Simulator/python_script/` folder of the source code repository:

- [pace_simulator_twolayer_inputfile_land.py](https://github.com/aoog-umbc/rtsos/blob/main/test/PACE_Simulator/python_script/pace_simulator_twolayer_inputfile_land.py)
- [pace_simulator_twolayer_inputfile_ocean.py](https://github.com/aoog-umbc/rtsos/blob/main/test/PACE_Simulator/python_script/pace_simulator_twolayer_inputfile_ocean.py)

For ocean simulations, ocean inherent optical properties are documented in Zhai et al. (2017, 2018, 2022).

Study each file before you execute them. Most importantly, set `aux_dir` to the path
ending in the `Data` folder and `gas_absorption_coeff_dir` to the path ending in the
`Gas_Absorption_Coefficients` folder.

####  Aerosol Models

The `Aerosol_Model[iaerosol]` value in the input file generator scripts is used in the
PACE simulator input file to specify the aerosol models.

In the Python script, `Aerosol_Model` is a array, and `iaerosol` specifies which one of
the aerosol models to use. In the following we use `IAEROSOL` as
`Aerosol_Model[iaerosol]` equivalently.

- The `IAEROSOL` values from 1-10 use the Shettle & Fenn (1979) model.
- The `IAEROSOL` values from 11-20 are the operational aerosol models currently used in
  the atmospheric correction at GSFC (Ahmad et al. 2010). The values span a range of
  fine mode fractions. When using these aerosol models, a value `irh` is also used to
  select relative humidity from the array `RH`.
- The `IAEROSOL` value of `-99` is reserved for the case that a user supplies aerosol
  and cloud phase matrices. See [section 2.2](#21a2-single-scattering-matrix-code) for
  more information.

#### Run the Simulator

Now you may try out the PACE simulator. First generate an input file:

```shell
$ python pace_simulator_twolayer_inputfile_land.py
```

The script will write a parameter file with a long name, which we refer to here as
`input.txt`. Call the executable with the parameter file as the only argument.

```shell
$ rtsos_PACE_Simulator_DoubleK input.txt
```

### 2.1A.2 Single Scattering Matrix Code

If you set `IAEROSOL=-99` in the PACE simulator input file, you will need to prepare the
aerosol phase matrix. A tool was written by Neranga K. Hannadige
(https://scholar.google.com/citations?user=0lAwqNsAAAAJ&hl=en) to generate the input
scattering matrix files for the PACE simulator.

We provided a Python script to generate the input files in the
`test/PACE_Simulator/python_script/` folder of the source code repository:

- [Aerosol_PhaseMatrix_Cal_InputPre.py](https://github.com/aoog-umbc/rtsos/blob/main/test/PACE_Simulator/python_script/Aerosol_PhaseMatrix_Cal_InputPre.py)

In this script, you want to specify `aux_dir` as for the PACE simulator.

```shell
$ python Aerosol_PhaseMatrix_Cal_InputPre.py
```

The script will write a parameter file with a long name, which we refer to here as
`input.txt`. Call the executable with the parameter file as the only argument.

```shell
$ rtsos_Aerosol_Phmx_Cal input.txt
```

The output is an aerosol phase matrix that can be referenced in the parameter file for
the PACE simulator. To do so, set the `aerosol_phasematrix_file_hi` or the
`aerosol_phasematrix_file_low` to the files that are generated with the procedure shown
above.

The single scattering matrix simulation uses the Mie code developed by M. Mishchenko
(Mishchenko et al. 2002). A non spherical dust scattering matrix database is included,
which was provided by Prof. Ping Yang at Texas A&M University (Meng et al. 2010). The
single scattering matrix code assumes either a bimodal (NMODE=2) or trimodal (NMODE=3)
size distribution:

```
----------------------------------------------------------------------
NMODE=2
----------------------------------------------------------------------
Available aerosol modes: 

	*fine mode 		(fmfrac)
		dust like 		(rf1)
		water soluble 	(rf2)
		BrC 			(rf3)
		soot			(1-(rf1+rf2+rf3))
	
	*coarse mode 	(1-fmfrac)
		sea salt		(cmsfrac)
		dust			(1-cmsfrac)
----------------------------------------------------------------------
NMODE=3
----------------------------------------------------------------------
Available aerosol modes:

	*fine mode		(fmfrac)
		water soluble 	(rf2)
		BrC 			(rf3)
		soot			(1-(rf2+rf3))

	*dust mode		(dustfrac)
		dust like		(dsfrac)
		dust			(1-dsfrac)

	*coarse mode 	(1-(fmfrac+dustfrac))
		coarse mode includes sea salt only
----------------------------------------------------------------------

*Note: Different sets of parameters are used for nmode=2 and nmode=3 
		For nmode=2 dustfrac and dsfrac do not exist
		For nmode=3 rf1 and cmsfrac do not exist
```

### 2.1B Python Wrapper

The Python wrapper adds a command line tool `rtsos` with subcommands `GSFC_AC_LUT`,
`PACE_Simulator_DoubleK`, and `Aerosol_Phmx_Cal`. These wrappers presently work by
calling the binary executables `rtsos_GSFC_AC_LUT`, `rtsos_PACE_Simulator_DoubleK`, and
`rtsos_Aerosol_Phmx_Cal` which are included with the install.

The Python wrapper adds the ability to run the programs over arrays of parameters
provided in NetCDF files.  To run the PACE simulator with default values, call
`rtsos PACE_Simulator_DoubleK` followed by two arguments giving, first, the path where
default parameters will be written and second, the path for outputs.

```shell
$ rtsos PACE_Simulator_DoubleK inputs.nc outputs.nc
```

To run the PACE simulator over non-default parameters, you must provide a suitable
NetCDF file.

1. Generate a template for the parameter sets over which you wish to run the
   PACE simulator.
   ```shell
   $ rtsos PACE_Simulator_DoubleK --pre defaults.nc
   ```
1. Using whatever tools you want, generate a new input dataset with the same variables
   as the default dataset. New dimensions are allowed, but every coordinate variable
   must be either scalar or one dimensional. Certain coordinates must be scalar, but
   which ones is not yet presently documented.
1. Choose to run on the outer-product of all parameterizations or of subsets (e.g. for
   parallel execution):
    - To run all parameterizations in one process, provide inputs and outputs file
      paths.
      ```shell
      $ rtsos PACE_Simulator_DoubleK inputs.nc outputs.nc
      ```
    - To run subsets of the parameterizations as separate jobs, add the `--split`
      argument to indicate dimensions and indices over which to slice the inputs.
      ```shell
      $ rtsos PACE_Simulator_DoubleK --split=chla:0-1,Solar_Zenith_Angle inputs.nc outputs.nc
      $ rtsos PACE_Simulator_DoubleK --split=chla:2-4,Solar_Zenith_Angle inputs.nc outputs.nc
      ```
      Note that instead of writing to `outputs.nc`, a tree is written under `outputs/`
      with results from each job. Indices are zero-based and ranges include both
      endpoints. These two commands cover all five `chla` values and both
      `Solar_Zenith_Angle` values in the default template, leaving `tau_ref_hi`
      unsplit. For example, one product is
      `outputs/chla/2/Solar_Zenith_Angle/0/outputs.nc`. Once both jobs are complete,
      use the same dimension order with `--combine` to combine their outputs.
      ```shell
      $ rtsos PACE_Simulator_DoubleK --combine=chla,Solar_Zenith_Angle inputs.nc outputs.nc
      ```
      To aggregate only solar zenith angles while retaining separate `chla`
      directories, use the following instead. Each result is written to
      `outputs/chla/<index>/outputs.nc` and can then be combined across `chla`.
      ```shell
      $ rtsos PACE_Simulator_DoubleK --split=chla --combine=Solar_Zenith_Angle inputs.nc outputs.nc
      $ rtsos PACE_Simulator_DoubleK --combine=chla inputs.nc outputs.nc
      ```

The default PACE template retains the legacy 20-state sweep: five `chla` values,
two solar zenith angles, and two `tau_ref_hi` values. Add `--dry-run` to generate
the native input files without running the simulator.

To use the extracted auxiliary data, pass the directory containing both `Data`
and `Gas_Absorption_Coefficients`:

```shell
$ rtsos PACE_Simulator_DoubleK --data-root /path/to/RTSOS_data --pre inputs.nc
$ rtsos PACE_Simulator_DoubleK --data-root /path/to/RTSOS_data inputs.nc outputs.nc
```

`--data-root` replaces `Aux_Dir` and `Atmos_Dir`, including in an existing input
dataset. Without it, the template retains the legacy `Auxiliary_Files` and
`Gas_Absorption_Coefficients` defaults. The Python initializer uses the same path
handling:

```python
from rtsos.PACE_Simulator_DoubleK import init

pace = init(data_root="/path/to/RTSOS_data")
inputs = pace.defaults
# Select one state from the default sweep when only one simulation is wanted.
inputs = inputs.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
inputs.to_netcdf("inputs.nc")
```

For individual path edits, replace the whole coordinate. In-place assignment such
as `inputs.Aux_Dir[()] = longer_path` can silently truncate fixed-width strings.
Preserve metadata by supplying the coordinate's dimensions and attributes:

```python
inputs = inputs.assign_coords(
    Aux_Dir=((), "/path/to/Data", inputs.Aux_Dir.attrs),
    Atmos_Dir=((), "/path/to/Gas_Absorption_Coefficients", inputs.Atmos_Dir.attrs),
)
```

Direct paths are expanded and resolved relative to the process working directory
before native execution, including aerosol phase-matrix files and the land
spectrum. `CFILE_AP` remains a profile filename such as `afglus.dat`, read under
`Atmos_Dir/atmmod/`. Spaces in paths are supported; `#`, line breaks and trailing
whitespace cannot be represented faithfully by the native input format. The
native readers use fixed character buffers (180 bytes for atmospheric/gas paths,
360 for PACE input and other filenames); full paths, including filenames appended
inside the readers, must fit those buffers.

`OCEAN_CASE_SELECT` is the legacy native selector for both land surfaces and
ocean optical models. Its supported values are:

| Code | Model |
| --- | --- |
| `-203` | Ross–Li land BRDF with a polarized-reflectance term. |
| `-202` | Snow BRDF from the native lookup data, using the supplied spectral albedo. |
| `-201` | Modified Rahman–Pinty–Verstraete (mRPV) land model. |
| `-200` | Lambertian/Fresnel land mixture; a Lambertian fraction of one gives a pure Lambertian surface. |
| `0` | Ocean interface without a scattering water body; surface reflection remains. |
| `1` | Native chlorophyll-based ocean optical model. |
| `2` | Native chlorophyll-and-sediment ocean optical model. |
| `3` | Seven-parameter ocean optical model. |

Case `0` is not an atmosphere-only calculation: it retains the ocean interface
and the applicable surface processes. The native program disables Raman
scattering and both ocean fluorescence flags for this case and for land. Result
coordinates reflect those effective native values; the serialized request keeps
its original flag values. Other documented native resets, such as enabling
monochromatic processing when gas absorption is disabled, are handled similarly.

The seven parameters in case `3` are `chla`, `adg440`,
`bbp660_BackscatterCoeff`, `Bp660_BackscatterFraction`, `Sdg`, `Sbp`, and `S_Bp`.
The integer codes identify the existing implementations; they do not describe
every possible model of a land surface or water body.

For land, choose one of the four negative codes and supply the spectrum file:

```python
import numpy as np

inputs = inputs.assign_coords(
    OCEAN_CASE_SELECT=np.int32(-200),
    land_ref_spectra_filename=(
        (), "/path/to/land_ref_spectra_Lamb", inputs.land_ref_spectra_filename.attrs
    ),
)
inputs.to_netcdf("land.nc")
```

The spectrum has a header, a wavelength count, then seven columns: wavelength in
nm, four model-specific surface coefficients, and the real and imaginary
refractive index. The supplied native land generator documents these columns
and gives examples for each surface type. Negative cases omit the ocean input
records; nonnegative cases, including zero (ocean interface without a water
body), retain them. The land-spectrum record is written only for negative cases.
Ocean coordinates retained in a land dataset do not affect its native input;
select their sweep dimensions explicitly if they would otherwise repeat land work.

### 2.1C Readable physical cases

One-state examples are provided in [examples/pace_ocean.toml](examples/pace_ocean.toml)
and [examples/pace_land.toml](examples/pace_land.toml). Set their `[data]` locations
to your auxiliary data, then use the shared simulation route:

```shell
$ rtsos simulate examples/pace_ocean.toml ocean.nc
$ rtsos simulate examples/pace_land.toml land.nc
$ rtsos simulate --dry-run examples/pace_ocean.toml ocean.nc
```

The same case runs from Python:

```python
from rtsos import load_case, simulate

case = load_case("examples/pace_ocean.toml")
result = simulate(case, "ocean.nc")
```

Case-file paths resolve relative to the case file's directory. Output and `--tmp`
paths resolve relative to the invocation directory. Legacy commands retain their
existing path behavior, 20-state template, native processing and output variables.
New cases choose named surface and water models and supply their active sweep
quantities explicitly. See [the case and response guide](examples/README.md) for
the model choices, units, geometry, exact-response loading and native capabilities.

### Python tests

From the repository root in an environment with the Python dependencies and pytest
installed, run the lightweight checks with
`PYTHONPATH=src python -m pytest src/rtsos/tests`. These include template generation,
dry runs, and synthetic native-output split/combine checks. Native integration
tests are opt-in with `PYTHONPATH=src python -m pytest src/rtsos/tests --native` and require the
installed executables and auxiliary data described above.

## References

Ahmad, Z., Franz, B. A., McClain, C. R., Kwiatkowska, E. J., Werdell, J., Shettle, E.
P., & Holben, B. N. (2010). New aerosol models for the retrieval of aerosol optical
thickness and normalized water-leaving radiances from the SeaWiFS and MODIS sensors over
coastal regions and open oceans. Applied Optics, 49(29), 5545.
https://doi.org/10.1364/ao.49.005545

Gordon, I. E., Rothman, L. S., Hargreaves, R. J., Hashemi, R., Karlovets, E. V.,
Skinner, F. M., Conway, E. K., Hill, C., Kochanov, R. V., Tan, Y., Wcisło, P., Finenko,
A. A., Nelson, K., Bernath, P. F., Birk, M., Boudon, V., Campargue, A., Chance, K. V.,
Coustenis, A., … Yurchenko, S. N. (2022). The HITRAN2020 molecular spectroscopic
database. Journal of Quantitative Spectroscopy and Radiative Transfer, 277, 107949.
https://doi.org/10.1016/j.jqsrt.2021.107949

Shettle, E. P. and R. W. Fenn. (1979). “Models for the aerosols of the lower atmosphere
and the effects of humidity variations on their optical properties,” AFGL-TR 790214, U.
S. Air Force Laboratory, Hanscom Air Force Base, Mass..

Meng, Z., Yang, P., Kattawar, G. W., Bi, L., Liou, K. N., & Laszlo, I. (2010).
Single-scattering properties of tri-axial ellipsoidal mineral dust aerosols: A database
for application to radiative transfer calculations. Journal of Aerosol Science, 41(5),
501–512. https://doi.org/10.1016/j.jaerosci.2010.02.008

Mishchenko, M. I., Travis, L. D., & Lacis, A. A. (2002). Scattering, absorption, and
emission of light by small particles. Cambridge, England: Cambridge University Press.
https://search.worldcat.org/title/48932103

Zhai, P., Hu, Y., Trepte, C. R., & Lucker, P. L. (2009). A vector radiative transfer
model for coupled atmosphere and ocean systems based on successive order of scattering
method. Optics Express, 17(4), 2057. https://doi.org/10.1364/oe.17.002057

Zhai, P., Hu, Y., Chowdhary, J., Trepte, C. R., Lucker, P. L., & Josset, D. B.
(2010). A vector radiative transfer model for coupled atmosphere and ocean systems with
a rough interface. Journal of Quantitative Spectroscopy and Radiative Transfer,
111(7–8), 1025–1040. https://doi.org/10.1016/j.jqsrt.2009.12.005

Zhai, P., Hu, Y., Josset, D. B., Trepte, C. R., Lucker, P. L., & Lin, B. (2013).
Advanced angular interpolation in the vector radiative transfer for coupled atmosphere
and ocean systems. Journal of Quantitative Spectroscopy and Radiative Transfer, 115,
19–27. https://doi.org/10.1016/j.jqsrt.2012.09.018

Zhai, P., Hu, Y., Winker, D. M., Franz, B. A., & Boss, E. (2015). Contribution of
Raman scattering to polarized radiation field in ocean waters. Optics Express, 23(18),
23582. https://doi.org/10.1364/oe.23.023582

Zhai, P., Hu, Y., Winker, D. M., Franz, B. A., Werdell, J., & Boss, E. (2017). Vector
radiative transfer model for coupled atmosphere and ocean systems including inelastic
sources in ocean waters. Optics Express, 25(8), A223.
https://doi.org/10.1364/oe.25.00a223

Zhai, P., Boss, E., Franz, B., Werdell, P. J., & Hu, Y. (2018). Radiative Transfer
Modeling of Phytoplankton Fluorescence Quenching Processes. Remote Sensing, 10(8), 1309.
https://doi.org/10.3390/rs10081309

Zhai, P., & Hu, Y. (2022). An improved pseudo spherical shell algorithm for vector
radiative transfer. Journal of Quantitative Spectroscopy and Radiative Transfer, 282,
108132. https://doi.org/10.1016/j.jqsrt.2022.108132

Zhai, P., Gao, M., Franz, B. A., Werdell, P. J., Ibrahim, A., Hu, Y., & Chowdhary, J.
(2022). A Radiative Transfer Simulator for PACE: Theory and Applications. Frontiers in
Remote Sensing, 3. https://doi.org/10.3389/frsen.2022.840188

[CMake]: https://cmake.org/
[LAPACK]: https://www.netlib.org/lapack/
[Fortran]: https://fortran-lang.org/
[HDF5]: https://www.hdfgroup.org/solutions/hdf5/
[Python]: https://python.org
[latest release]: https://github.com/aoog-umbc/rtsos/releases/latest
[Conda]: https://conda.org/
