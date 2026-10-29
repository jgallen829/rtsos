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

The software is distributed as source code, so the programs must be built locally. There
are two options for building:
- use `cmake` to build and install the binary RTSOS executables
- use `pip` to build and install the binary RTSOS executables along with their Python
   wrappers

### 1.1A First Option

The software is configured for building with `cmake` on a system with the following
prerequistes available:

- [CMake] >= 3.31
- [HDF5] >= 1.10
- [LAPACK] >= 3.9
- a modern [Fortran] compiler (e.g. gfortran)

One option for installing these prerequistes is through [Conda], which is described in
[section 1.1B](#11b-second-option), but installation with a system package manager is
also suitable.

#### Get the Source Code

You can retrieve the [latest release] by downloading the `.zip` or `.tar.gz` assets from
GitHub. Alternatively, to retrieve the development version, use `git`:

```shell
git clone --filter=blob:none https://github.com/aoog-umbc/rtsos
```

#### Build and Install

Open a Terminal with the working directory at the project root (i.e. where this
README.md is located). Choose an unused name for the RTSOS build directory, and
optionally an existing installation prefix, and call `cmake`. For example, to use a
folder named "build" for the build files and the `.local` folder within your home
directory as the prefix for installation:

```shell
$ cmake -B build --install-prefix ~/.local
```

The `--install-prefix` argument is optional. Now change to the build directory and run
cmake twice more, first with `--build` and then with `--install`.

```shell
$ cd build
$ cmake --build .
$ cmake --install .
```

CMake will display the locations of the binary executables installed. If CMake cannot
find the prerequisties, you can provide an environment variable with the prefix, e.g.
`export HDF5_ROOT=...`. The prefix is the path up to the folder containing the "lib"
folder with the HDF5 libraries.

### 1.1B Second Option

The software is configured for building with `pip` on a system with the following
prerequistes available:

- [Python] >= 3.12
- [HDF5] >= 1.10
- [LAPACK] >= 3.9
- a modern [Fortran] compiler (e.g. gfortran)

The installer `pip` calls the Python wrapper's build-system, `scikit-build-core`, which
also uses `cmake` to compile the source code. Those tools are fetched by `pip` as
needed.

One option for installing the prerequistes is through [Conda], which is our recommended
approach for preparing a Python environment for use with the Python wrapper. Choose a
name for the Conda environment and create it with the required prerequistes. For
example, with "my-project" as the environment name:

```shell
$ conda create --yes --name my-project python hdf5 lapack
$ conda activate my-project
```

To install the [latest release] from PyPI:

```shell
(my-project) $ pip install rtsos
```

To install the development version from GitHub:

```shell
(my-project) $ pip install git+https://github.com/aoog-umbc/rtsos
```

If installation fails with the message `ERROR: Could not build wheels for rtsos...`,
then CMake may not have found the prerequisties. You can provide an environment variable
with their prefix, e.g. `export HDF5_ROOT=...`. The prefix is the path up to the folder
containing the "lib" folder with the HDF5 libraries.

For a build from a Git checkout, retain the release tags so the version provider can
derive the package version. The public tag `v1.0-public-2025-10-23` supplies version
`1.0`; subsequent commits receive a development version. No pretend-version
environment variable is needed. To build a wheel locally, run:

```shell
$ python -m pip wheel --no-deps --wheel-dir dist .
```

Some macOS HDF5 installations expose the Fortran library but omit the high-level
Fortran library location during CMake discovery. If CMake reports that it cannot
locate `hdf5::hdf5_hl_fortran`, provide the installed library explicitly, using
your HDF5 prefix (and its actual library filename):

```shell
$ export HDF5_ROOT="$CONDA_PREFIX"
$ CMAKE_ARGS="-DHDF5_Fortran_HL_LIBRARY=$HDF5_ROOT/lib/libhdf5_hl_fortran.dylib" python -m pip wheel --no-deps --wheel-dir dist .
```

The same setting can be passed as `-DHDF5_Fortran_HL_LIBRARY=...` to a direct
CMake configure. Use a compiler and HDF5 Fortran modules from compatible toolchains.

### 1.2 Auxiliary Data

In addition to the software, simulating the PACE instruments and generating atmospheric
correction look-up tables requires a collection of data files available from a Zenodo
repository:

https://doi.org/10.5281/zenodo.17410093

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

## 2. Usage

The first or two ways to use RTSOS is with a single parameterization given to the binary
executables as a plain text file. You can do this regardless of the installation method
above; follow the procedures in [section 2.1A](#21a-binary-executables). The second way
to use RTSOS is with multiple parameterizations given to the Python wrapper in a single
NetCDF file. You can do this if you used the second installation option [section
1.1B](#11b-second-option); follow the procedures
in [section 2.1B](#21b-python-wrapper).

### 2.1A Binary Executables

For monochromatic simulation, users are referred to the documentation located at
"Documentation/Describtion_of_input_files/SOS_io.pdf" in the source code repository.

For readable physical cases at requested wavelengths or through complete response
curves, the installed Python/CLI route is:

```shell
$ rtsos simulate examples/monochromatic_ocean.toml ocean.nc
$ rtsos simulate examples/monochromatic_land.toml land.nc
$ rtsos simulate examples/response_ocean.toml ocean_bands.nc
$ rtsos simulate examples/response_land.toml land_bands.nc
```

Set each example's data locations first. The monochromatic examples evaluate
551.25 and 443.125 nm; the land example supplies micrometres. Python uses the same `simulate` function.
See [the examples guide](examples/README.md#requested-monochromatic-wavelengths)
for the wavelength, model, data, illumination, units and numerical contracts.
Full-response integration uses every supplied calibration node and shares scene
preparation across the wavelength set. See [band simulation](examples/README.md#full-response-band-simulation)
for the trapezoidal method, complete OCI selections and result schema.
For paired observations at requested directions, including a single observation:

```shell
$ rtsos simulate examples/monochromatic_directions.toml directions.nc
$ rtsos simulate examples/response_directions.toml direction_bands.nc
$ rtsos simulate examples/single_direction.toml single.nc
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
