"""Requested-wavelength adapter over the shared public physical preparation."""

from pathlib import Path

import numpy as np
import xarray as xr
from .core import expect_exe
from .parameters import Parameters as P
from .PACE_Simulator_DoubleK import PACE, init
from .physical import _case_keys, physical_case_inputs


def case_inputs(case):
    _case_keys(
        case.spectral,
        {"method", "wavelength_unit", "wavelengths"},
        "spectral",
        required={"method", "wavelength_unit", "wavelengths"},
    )
    _case_keys(case.sensor, {"name"}, "monochromatic sensor")
    unit = case.spectral["wavelength_unit"]
    if unit not in ("nm", "um"):
        raise ValueError("spectral.wavelength_unit must be nm or um")
    wavelengths = np.asarray(case.spectral["wavelengths"])
    if (
        wavelengths.ndim != 1
        or not wavelengths.size
        or wavelengths.dtype.kind not in "iuf"
        or not np.all(np.isfinite(wavelengths))
        or np.any(wavelengths <= 0)
    ):
        raise ValueError("spectral.wavelengths must be a nonempty positive finite list")
    if case.illumination.get("source") not in ("tsis", "thuillier"):
        raise ValueError(
            "Monochromatic illumination.source must explicitly select tsis or thuillier"
        )
    inputs = physical_case_inputs(case, monochromatic=True)
    if inputs.ATMOS_ZERO.item() and inputs.gas_abs_flag.item():
        raise ValueError(
            "Gas absorption requires an atmosphere; set physics.gas_absorption=false when physics.atmosphere=false"
        )
    if inputs.Aerosol_Model.item() < -1:
        raise NotImplementedError(
            "External aerosol/cloud phase matrices need requested-wavelength matrix loading; "
            "the existing files use PACE segments"
        )
    if inputs.Aerosol_Model.item() not in (-1, *range(1, 21)):
        raise ValueError("Aerosol_Model must be -1 or 1 through 20")
    if inputs.OCEAN_CASE_SELECT.item() > 0:
        for flag, operation in (
            ("OCEAN_RAMAN_FLAG", "Raman excitation/emission preparation"),
            (
                "OCEAN_FCHLA_FLAG",
                "chlorophyll fluorescence excitation/emission preparation",
            ),
            ("OCEAN_FCDOM_FLAG", "CDOM fluorescence excitation/emission preparation"),
            ("NonPhotochemicalQuenching_FLAG", "NPQ excitation-field preparation"),
        ):
            if inputs[flag].item():
                raise NotImplementedError(
                    f"Requested-wavelength {operation} is not connected; set the process explicitly false"
                )
    else:
        inputs = inputs.assign_coords(
            OCEAN_RAMAN_FLAG=np.int32(0),
            OCEAN_FCHLA_FLAG=np.int32(0),
            OCEAN_FCDOM_FLAG=np.int32(0),
            NonPhotochemicalQuenching_FLAG=np.int32(0),
            OCEAN_PHMX_ONE=np.int32(1),
        )
    # Spherical harmonics use l >= m. MAXLFULL is the core's compiled storage bound.
    expansion = max(inputs.nquadainput.item(), inputs.MAXMORDINPUT.item())
    if expansion > 3000:
        raise ValueError(
            "Phase expansion exceeds the core's compiled MAXLFULL=3000; rebuild the core to increase it"
        )
    inputs.attrs.update(
        rt_fourier_order_limit=inputs.MAXMORDINPUT.item(),
        phase_expansion_order=expansion,
        solar_policy=case.illumination["source"],
        requested_wavelength_unit=unit,
    )
    wavelengths = wavelengths.astype(np.float64)
    return inputs.assign_coords(
        request=("request", np.arange(wavelengths.size, dtype=np.int64)),
        wavelength_requested=("request", wavelengths, {"units": unit}),
        wavelength_nm=(
            "request",
            wavelengths * (1000 if unit == "um" else 1),
            {"units": "nm"},
        ),
    )


def solar_source(case):
    auxiliary = (
        Path(case.data["root"]) / "Data"
        if "root" in case.data
        else Path(case.data["auxiliary"])
    )
    name = {"tsis": "tsis1_hsrs_binned_fs.csv", "thuillier": "f0.txt"}[
        case.illumination["source"]
    ]
    path = auxiliary / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return path.resolve()


def _coverage(wavelengths, grid, name):
    low, high = np.min(grid), np.max(grid)
    if np.any((wavelengths < low) | (wavelengths > high)):
        raise ValueError(
            f"Requested wavelengths exceed {name} coverage [{low:g}, {high:g}] nm"
        )


def check_data(case, inputs, solar):
    """Check actual active table coverage once before starting native processes."""
    waves = inputs.wavelength_nm.values
    auxiliary = Path(inputs.Aux_Dir.item())
    # All connected built-in aerosol models use the 0.30–2.50 um SF index tables.
    _coverage(
        np.append(waves, inputs.wv_pace_ref.item()),
        [300, 2500],
        "built-in aerosol refractive-index tables",
    )
    if solar.name == "f0.txt":
        solar_waves = np.loadtxt(solar, skiprows=13, usecols=0)
    else:
        solar_waves = np.loadtxt(solar, skiprows=1, delimiter=",", usecols=0)
    _coverage(waves, solar_waves, f"solar file {solar.name}")
    code = inputs.OCEAN_CASE_SELECT.item()
    if code < 0:
        spectrum = np.loadtxt(inputs.land_ref_spectra_filename.item(), skiprows=2)
        _coverage(waves, spectrum[:, 0], "land spectrum")
    elif code > 0:
        short = np.loadtxt(
            auxiliary / "aw_proctol_IOCCG2018.txt", skiprows=13, max_rows=211
        )
        long = np.loadtxt(auxiliary / "water_coef.txt", skiprows=30, max_rows=2250)
        _coverage(waves[waves <= 550], short[:, 0], "IOCCG pure-water absorption")
        _coverage(waves[waves > 550], long[:, 0], "pure-water absorption")
        if np.any(waves < 400):
            uv = np.loadtxt(auxiliary / "ap_UV_biosope.dat", skiprows=5, max_rows=226)
            _coverage(waves[waves < 400], uv[:, 0], "BIOSOPE phytoplankton absorption")
        if inputs.ap_select.item() == 2:
            table = np.loadtxt(auxiliary / "ap_PicoMicro.txt", skiprows=1)
            _coverage(
                np.append(waves[(waves >= 400) & (waves < 700)], [440, 676]),
                table[:, 0],
                "Ciotti phytoplankton absorption",
            )


class Monochromatic(PACE):
    """One isolated native process per scene's wavelength set, through Program."""

    def __init__(self):
        pace = init()
        super().__init__(
            defaults=pace.defaults, params=tuple(getattr(P, n) for n in pace._params)
        )
        self.program = "rtsos_Monochromatic" + (".exe" if expect_exe else "")

    def infile(self, path, dataset):
        infile, outfile = super().infile(path, dataset)
        wavelengths, _ = _evaluation_grid(dataset)
        with (path / infile).open("a") as stream:
            stream.write(f"{dataset.attrs['solar_policy']}\n{len(wavelengths)}\n")
            stream.writelines(f"{wave:.17g}\n" for wave in wavelengths)
        return infile, outfile

    def execution_inputs(self, inputs):
        yield inputs

    def read_output(self, path, inputs):
        _, inverse = _evaluation_grid(inputs)
        result = xr.load_dataset(path)
        return (
            result.isel(evaluation=xr.DataArray(inverse, dims="request"))
            .assign_coords(request=inputs.request)
            .drop_vars("evaluation")
        )

    def result_inputs(self, inputs, outputs):
        if not np.array_equal(outputs.Wavelength_Evaluated, inputs.wavelength_nm):
            raise ValueError(
                "Native evaluated wavelength differs from the requested wavelength"
            )
        # These fields exist only to bridge the established physical record layout.
        return inputs.drop_vars(
            ["MONOCHROMATIC_FLAG", "hyspectral_flag", "WAVEBAND_SEG_FLAG"]
        )


def _evaluation_grid(inputs):
    """Stable exact-value reuse; request identities remain independent of sorting."""
    waves = inputs.wavelength_nm.values.reshape(-1)
    unique, first, inverse = np.unique(waves, return_index=True, return_inverse=True)
    order = np.argsort(first)
    return unique[order], np.argsort(order)[inverse]
