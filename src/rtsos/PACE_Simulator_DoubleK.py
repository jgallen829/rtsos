# ===============================================================
#  RTSOS – Radiative Transfer model based on Successive Orders of Scattering
#  Copyright (C) 2025  Pengwang Zhai

#  Licensed under the Creative Commons Attribution–NonCommercial 4.0
#  International License (CC BY-NC 4.0).
#  You may use, modify, and share this code for research and
#  educational purposes with proper attribution.
#  Commercial use requires written permission from the author.

#  Full license: https://creativecommons.org/licenses/by-nc/4.0/
#  Contact: Pengwang Zhai  |  [pwzhai@gmail.com]
# ===============================================================

from pathlib import Path

import numpy as np
import xarray as xr

from .core import Program
from .parameters import Parameters as P

# Records omitted by the negative-surface branches of the native PACE parser.
OCEAN_PARAMETERS = {
    "albedo_ground",
    "water_depth_max",
    "chla",
    "phytoplankton_index_refraction",
    "phytoplankton_spectral_slope",
    "sediment_index_refraction",
    "sediment_spectral_slope",
    "sediment_concentration",
    "adg440",
    "bbp660_BackscatterCoeff",
    "Bp660_BackscatterFraction",
    "Sdg",
    "Sbp",
    "S_Bp",
    "CHLA_HOMOGENEITY",
    "OCEAN_RAMAN_FLAG",
    "OCEAN_FCHLA_FLAG",
    "OCEAN_FCDOM_FLAG",
    "ap_select",
    "NonPhotochemicalQuenching_FLAG",
    "OCEAN_PHMX_ONE",
    "SUN_GLINT_FLAG",
}


class PACE(Program):
    """Class that drives `rtsos PACE_Simulator_DoubleK`"""

    def __init__(self, *args, data_root=None, **kwargs):
        self.program = __name__.replace(".", "_")
        self.data_root = Path(data_root).expanduser().resolve() if data_root else None
        super().__init__(*args, **kwargs)
        self._selected_params = tuple(
            name for name in self._params if name != "land_ref_spectra_filename"
        )
        if self.data_root is not None:
            self.defaults = self.prepare_inputs(self.defaults)

    def prepare_inputs(self, dataset: xr.Dataset) -> xr.Dataset:
        """Replace path coordinates, resolving direct paths against the working directory.

        CFILE_AP remains a filename relative to Atmos_Dir/atmmod, as read natively.
        Replacing whole coordinates avoids NumPy fixed-width string truncation.
        """
        overrides = {}
        if self.data_root is not None:
            overrides = {
                "Aux_Dir": self.data_root / "Data",
                "Atmos_Dir": self.data_root / "Gas_Absorption_Coefficients",
            }
        for name in (
            "Aux_Dir",
            "Atmos_Dir",
            "Aerosol_Phasematrix_File_Hi",
            "Aerosol_Phasematrix_File_Low",
            "land_ref_spectra_filename",
        ):
            if name not in dataset:
                continue
            coord = dataset[name]
            values = (
                np.broadcast_to(str(overrides[name]), coord.shape)
                if name in overrides
                else coord.values
            )
            resolved = []
            for value in np.asarray(values).flat:
                value = str(value)
                if value:
                    value = str(Path(value).expanduser().resolve())
                    if name == "Aux_Dir":
                        value += "/"  # The native snow reader concatenates without a separator.
                resolved.append(value)
            dataset = dataset.assign_coords(
                {
                    name: (
                        coord.dims,
                        np.array(resolved).reshape(coord.shape),
                        coord.attrs,
                    )
                }
            )
        return dataset

    @property
    def params(self) -> tuple:
        return self._selected_params

    @params.setter
    def params(self, dataset: xr.Dataset) -> None:
        surface = dataset.OCEAN_CASE_SELECT.item()
        if surface not in (-203, -202, -201, -200, 0, 1, 2, 3):
            raise ValueError(
                f"Unsupported OCEAN_CASE_SELECT {surface}; use -203, -202, -201, -200, 0, 1, 2, or 3"
            )
        if surface >= 0:
            self._selected_params = tuple(
                name for name in self._params if name != "land_ref_spectra_filename"
            )
        else:
            name = "land_ref_spectra_filename"
            if name not in dataset or not str(dataset[name].item()).strip():
                raise ValueError(
                    "land_ref_spectra_filename is required for negative OCEAN_CASE_SELECT"
                )
            self._selected_params = tuple(
                name for name in self._params if name not in OCEAN_PARAMETERS
            )
        # Validate the records actually written, after per-state model selection.
        for name in self._selected_params:
            if name not in dataset:
                continue
            value = dataset[name].item()
            if isinstance(value, str) and (
                any(c in value for c in "#\r\n") or value != value.rstrip()
            ):
                raise ValueError(
                    f"{name} cannot contain #, line breaks, or trailing whitespace in native PACE input"
                )

    def result_inputs(self, inputs: xr.Dataset, outputs: xr.Dataset) -> xr.Dataset:
        """Use native values for omitted records and check documented flag resets."""
        replacements = {
            name: outputs[name]
            for name in inputs.coords
            if name in outputs and name not in self.params
        }
        forced = {}
        # main_PACE_Simulator_DoubleK.f90 disables water-body processes here.
        if inputs.OCEAN_CASE_SELECT.item() <= 0:
            forced.update(
                OCEAN_RAMAN_FLAG=0,
                OCEAN_FCHLA_FLAG=0,
                OCEAN_FCDOM_FLAG=0,
                OCEAN_PHMX_ONE=1,
            )
        if inputs.ATMOS_ZERO.item() or not inputs.gas_abs_flag.item():
            forced["MONOCHROMATIC_FLAG"] = 1
        if inputs.ATMOS_ZERO.item():
            forced["SUN_GLINT_FLAG"] = 1
        for name, value in forced.items():
            if name in inputs and name in outputs:
                replacements[name] = xr.full_like(inputs[name], value)
        return inputs.assign_coords(replacements)


def init(*, data_root: str | Path | None = None) -> PACE:
    """Prepare subclass for execution defaults"""

    # values to write, in the order below, to the input file
    values = {
        P.Wind_Speed: 5.0,
        P.Solar_Zenith_Angle: [20.0, 60.0],
        P.wv_pace_ref: 532.0,
        P.tau_ref_hi: [0.1, 0.3],
        P.height_particle_hi: 15.0,
        P.height_particle_variance_hi: 2.0,
        P.tau_ref_low: 0.0,
        P.height_particle_low: -1,
        P.height_particle_variance_low: -1,
        P.Aerosol_Model: -1,
        P.Aerosol_FMF: 0.3,
        P.Relative_Humidity: 0.3,
        P.OCEAN_CASE_SELECT: 1,
        P.albedo_ground: 0.3,
        P.water_depth_max: 200.0,
        P.chla: [0.03, 0.3, 3, 10, 30],
        P.phytoplankton_index_refraction: 1.02,
        P.phytoplankton_spectral_slope: 3.0,
        P.sediment_index_refraction: 1.2,
        P.sediment_spectral_slope: 4.0,
        P.sediment_concentration: 0.0,
        P.adg440: 1.0,
        P.bbp660_BackscatterCoeff: 0.05,
        P.Bp660_BackscatterFraction: 0.01,
        P.Sdg: 0.015,
        P.Sbp: 0.3,
        P.S_Bp: 0.01,
        P.ncolinput: 40,
        P.nquadainput: 40,
        P.nquadoinput: 60,
        P.MAXMORDINPUT: 40,
        P.NTHETAV: 36,
        P.NPHIV: 19,
        P.CHLA_HOMOGENEITY: 1,
        P.OCEAN_RAMAN_FLAG: 1,
        P.OCEAN_FCHLA_FLAG: 1,
        P.OCEAN_FCDOM_FLAG: 1,
        P.ap_select: 1,
        P.MONOCHROMATIC_FLAG: 0,
        P.hyspectral_flag: 0,
        P.NonPhotochemicalQuenching_FLAG: 0,
        P.OCEAN_PHMX_ONE: 1,
        P.ATMOS_ZERO: 0,
        P.SUN_GLINT_FLAG: 0,
        P.gas_abs_flag: 1,
        P.OZONE_COLUMN_DobsonUnit: 345.66,
        P.H2O_COLUMN: 1.4387,
        P.NO2_COLUMN_DobsonUnit: 0.2075,
        P.Pressure_Surface_mb: 1013.0,
        P.WAVEBAND_SEG_FLAG: 0,
        P.AirSensor_Height: 2.2,
        P.pss_flag: 0,
        P.Aux_Dir: "Auxiliary_Files",
        P.Atmos_Dir: "Gas_Absorption_Coefficients",
        P.CFILE_AP: "afglus.dat",
        P.Aerosol_Phasematrix_File_Hi: "",
        P.Aerosol_Phasematrix_File_Low: "",
        P.land_ref_spectra_filename: "",
    }

    # create a dataset to hold the inputs as coordinates
    params = P()
    dataset = params.make_dataset(values)

    # the callable object that runs the given program
    return PACE(defaults=dataset, params=tuple(values), data_root=data_root)
