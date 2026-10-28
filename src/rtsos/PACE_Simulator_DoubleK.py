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


# Common-case translation belongs to this scientific adapter, not Program.
LAND_MODELS = {"ross_li": -203, "snow": -202, "mrpv": -201, "lambertian_fresnel": -200}
WATER_MODELS = {
    "none": 0,
    "chlorophyll": 1,
    "chlorophyll_sediment": 2,
    "seven_parameter": 3,
}
PACE_CHANNELS = {
    "all": 0,
    "oci_swir": 2,
    "oci_spexone": 3,
    "harp2": 4,
}
WATER_INPUTS = {
    "chla",
    "water_depth_max",
    "CHLA_HOMOGENEITY",
}
CASE2_INPUTS = {
    "phytoplankton_index_refraction",
    "phytoplankton_spectral_slope",
    "sediment_index_refraction",
    "sediment_spectral_slope",
    "sediment_concentration",
}
CASE3_INPUTS = {
    "adg440",
    "bbp660_BackscatterCoeff",
    "Bp660_BackscatterFraction",
    "Sdg",
    "Sbp",
    "S_Bp",
}


def _case_keys(table, allowed, section, required=()):
    extra = table.keys() - set(allowed)
    if extra:
        raise ValueError(f"Unsupported {section} settings: {', '.join(sorted(extra))}")
    missing = set(required) - table.keys()
    if missing:
        raise ValueError(f"{section} requires {', '.join(sorted(missing))}")


def case_inputs(case) -> xr.Dataset:
    """Resolve one physical case to the existing native PACE input layout.

    Scalar defaults have one owner, init(). Every sweep quantity that is active
    in the case must be supplied explicitly. Inactive native record placeholders
    stay internal; a user never needs ocean coefficients to describe land.
    """
    if case.spectral.get("method") != "pace_native":
        method = case.spectral.get("method")
        raise NotImplementedError(
            f"Spectral method {method!r} needs native preparation/evaluation at requested "
            "wavelengths and, for responses, custom response integration. The PACE "
            "driver only connects its existing band/double-k processing; its "
            "MONOCHROMATIC_FLAG does not accept arbitrary wavelengths."
        )
    _case_keys(case.spectral, {"method"}, "spectral")
    if case.sensor.get("name") != "pace" or "response_file" in case.sensor:
        raise NotImplementedError(
            "Custom sensors/responses need native requested-wavelength evaluation and "
            "response integration. Exact response data can be loaded/selected, but "
            "cannot be passed to the existing PACE driver."
        )
    channels = case.sensor.get("channels")
    if not isinstance(channels, str) or channels not in PACE_CHANNELS:
        raise NotImplementedError(
            "PACE connects whole native channel groups: "
            + ", ".join(PACE_CHANNELS)
            + ". Individual channel/variant execution needs native channel selection."
        )
    _case_keys(case.sensor, {"name", "channels"}, "sensor")
    if case.geometry.get("mode") != "pace_grid":
        raise NotImplementedError(
            "Explicit viewing directions need a native geometry input connection; PACE generates its viewing grid from counts"
        )
    _case_keys(
        case.geometry,
        {
            "mode",
            "solar_zenith_deg",
            "sensor_height_km",
            "zenith_count",
            "azimuth_count",
        },
        "geometry",
        required={"solar_zenith_deg"},
    )
    _case_keys(
        case.scene,
        {"atmosphere", "aerosol", "surface", "water"},
        "scene",
        required={"aerosol", "surface"},
    )
    _case_keys(
        case.numerics,
        {"ncolinput", "nquadainput", "nquadoinput", "MAXMORDINPUT"},
        "numerics",
    )
    _case_keys(case.illumination, {"source"}, "illumination")
    if case.illumination.get("source") not in ("pace_native", "tsis", "thuillier"):
        raise NotImplementedError(
            "Custom illumination needs a native solar-spectrum input connection; use pace_native, tsis, or thuillier"
        )
    _case_keys(case.data, {"root", "auxiliary", "gas_absorption"}, "data")
    if "root" in case.data:
        if len(case.data) != 1:
            raise ValueError("Use data.root or data.auxiliary and data.gas_absorption")
        auxiliary = str(Path(case.data["root"]) / "Data")
        gases = str(Path(case.data["root"]) / "Gas_Absorption_Coefficients")
    else:
        _case_keys(
            case.data,
            {"auxiliary", "gas_absorption"},
            "data",
            required={"auxiliary", "gas_absorption"},
        )
        auxiliary = case.data["auxiliary"]
        gases = case.data["gas_absorption"]

    defaults = init().defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    values = {name: coord.item() for name, coord in defaults.coords.items()}
    values.update(
        Aux_Dir=auxiliary, Atmos_Dir=gases, WAVEBAND_SEG_FLAG=PACE_CHANNELS[channels]
    )
    values.update(
        Solar_Zenith_Angle=case.geometry["solar_zenith_deg"],
        AirSensor_Height=case.geometry.get(
            "sensor_height_km", values["AirSensor_Height"]
        ),
        NTHETAV=case.geometry.get("zenith_count", values["NTHETAV"]),
        NPHIV=case.geometry.get("azimuth_count", values["NPHIV"]),
    )
    atmosphere = case.scene.get("atmosphere", {})
    _case_keys(
        atmosphere,
        {
            "OZONE_COLUMN_DobsonUnit",
            "H2O_COLUMN",
            "NO2_COLUMN_DobsonUnit",
            "Pressure_Surface_mb",
            "profile",
        },
        "scene.atmosphere",
    )
    values.update({k: v for k, v in atmosphere.items() if k != "profile"})
    values["CFILE_AP"] = atmosphere.get("profile", values["CFILE_AP"])
    if Path(values["CFILE_AP"]).name != values["CFILE_AP"]:
        raise ValueError("atmosphere.profile is a filename under gas_absorption/atmmod")
    aerosol = case.scene["aerosol"]
    _case_keys(
        aerosol,
        {
            "wv_pace_ref",
            "tau_ref_hi",
            "height_particle_hi",
            "height_particle_variance_hi",
            "tau_ref_low",
            "height_particle_low",
            "height_particle_variance_low",
            "Aerosol_Model",
            "Aerosol_FMF",
            "Relative_Humidity",
            "phase_matrix_hi",
            "phase_matrix_low",
        },
        "scene.aerosol",
        required={"tau_ref_hi"},
    )
    values["tau_ref_hi"] = aerosol["tau_ref_hi"]
    for key, value in aerosol.items():
        name = {
            "phase_matrix_hi": "Aerosol_Phasematrix_File_Hi",
            "phase_matrix_low": "Aerosol_Phasematrix_File_Low",
        }.get(key, key)
        values[name] = value
    surface = case.scene["surface"]
    model = surface.get("model")
    if model == "ocean":
        _case_keys(surface, {"model", "wind_speed_m_s"}, "scene.surface")
        values["Wind_Speed"] = surface.get("wind_speed_m_s", values["Wind_Speed"])
        if "water" not in case.scene:
            raise ValueError(
                "Ocean surfaces require scene.water.model, including none for surface-only ocean"
            )
        water = case.scene["water"]
        water_model = water.get("model")
        if water_model not in WATER_MODELS:
            raise ValueError(
                f"Unknown water model {water_model!r}; use {', '.join(WATER_MODELS)}"
            )
        code = WATER_MODELS[water_model]
        allowed = (
            (WATER_INPUTS if code else set())
            | (CASE2_INPUTS if code == 2 else set())
            | (CASE3_INPUTS if code == 3 else set())
        )
        if code in (1, 2):
            allowed = allowed | {"absorption_model"}
        required = ({"chla"} if code else set()) | (
            CASE3_INPUTS if code == 3 else set()
        )
        _case_keys(water, allowed | {"model"}, "scene.water", required=required)
        if code:
            values["chla"] = water["chla"]
        values.update(
            {k: v for k, v in water.items() if k not in ("model", "absorption_model")}
        )
        if "absorption_model" in water:
            absorption = {"bricaud": 1, "ciotti": 2}
            if water["absorption_model"] not in absorption:
                raise ValueError("absorption_model must be bricaud or ciotti")
            values["ap_select"] = absorption[water["absorption_model"]]
    elif model in LAND_MODELS:
        _case_keys(
            surface, {"model", "spectrum"}, "scene.surface", required={"spectrum"}
        )
        if "water" in case.scene:
            raise ValueError("Land surfaces do not take scene.water inputs")
        code = LAND_MODELS[model]
        values["land_ref_spectra_filename"] = surface["spectrum"]
    else:
        raise ValueError(
            f"Unknown surface model {model!r}; use ocean or {', '.join(LAND_MODELS)}"
        )
    values["OCEAN_CASE_SELECT"] = code

    physics = case.physics
    common_physics = {"scattering", "spherical_shell", "atmosphere", "gas_absorption"}
    ocean_physics = {"sun_glint"} if code >= 0 else set()
    water_physics = (
        {
            "raman",
            "chlorophyll_fluorescence",
            "cdom_fluorescence",
            "nonphotochemical_quenching",
            "single_water_phase_matrix",
        }
        if code > 0
        else set()
    )
    _case_keys(physics, common_physics | ocean_physics | water_physics, "physics")
    if physics.get("scattering", "vector") != "vector":
        raise NotImplementedError(
            "Scalar scattering needs a native PACE control connection; this driver uses vector RT"
        )
    switches = {
        "spherical_shell": ("pss_flag", False),
        "atmosphere": ("ATMOS_ZERO", True),
        "gas_absorption": ("gas_abs_flag", False),
        "sun_glint": ("SUN_GLINT_FLAG", True),
        "raman": ("OCEAN_RAMAN_FLAG", False),
        "chlorophyll_fluorescence": ("OCEAN_FCHLA_FLAG", False),
        "cdom_fluorescence": ("OCEAN_FCDOM_FLAG", False),
        "nonphotochemical_quenching": ("NonPhotochemicalQuenching_FLAG", False),
        "single_water_phase_matrix": ("OCEAN_PHMX_ONE", False),
    }
    for key, (name, invert) in switches.items():
        if key in physics:
            if not isinstance(physics[key], bool):
                raise ValueError(f"physics.{key} must be true or false")
            values[name] = int(not physics[key] if invert else physics[key])
    values.update(case.numerics)
    # Protect the native integer/real scalar boundary, without coercing bad
    # shapes, fractional counts, nonfinite values or integer overflow.
    for name, value in values.items():
        dtype = defaults[name].dtype
        if dtype.kind not in "if":
            continue
        array = np.asarray(value)
        if array.ndim or array.dtype.kind not in "ifb" or not np.isfinite(array):
            raise ValueError(f"{name} must be a finite numeric scalar")
        if dtype.kind == "i" and (
            value != int(value)
            or not np.iinfo(dtype).min <= value <= np.iinfo(dtype).max
        ):
            raise ValueError(f"{name} must be an integer representable by {dtype}")
        if dtype.kind == "f" and abs(value) > np.finfo(dtype).max:
            raise ValueError(f"{name} is outside native {dtype} precision")
    for name in (
        "ncolinput",
        "nquadainput",
        "nquadoinput",
        "MAXMORDINPUT",
        "NTHETAV",
        "NPHIV",
    ):
        if values[name] < 1:
            raise ValueError(f"{name} must be positive")
    if not 0 <= values["Solar_Zenith_Angle"] < 90:
        raise ValueError(
            "solar_zenith_deg must be in [0, 90) for downward solar illumination"
        )
    if not 1e-6 <= abs(values["AirSensor_Height"]) < 100:
        raise ValueError(
            "sensor_height_km would be reset by PACE; supply a nonzero height with magnitude below 100 km"
        )
    ntheta = values["NTHETAV"]
    if ntheta != 1 and ntheta % 2:
        raise ValueError(
            "PACE zenith_count must be 1 or even; its odd-grid middle angle is uninitialized"
        )
    if values["height_particle_hi"] < values["height_particle_low"]:
        raise ValueError(
            "Upper aerosol height must be at least the lower aerosol height"
        )
    inputs = P.make_dataset({getattr(P, name): value for name, value in values.items()})
    # The legacy driver passes this upper bound to RTSOS at each wavelength.
    # Keep the requested MAXMORDINPUT intact and expose its internal restriction.
    inputs.attrs["pace_rt_fourier_order_limit"] = min(
        30, values["nquadainput"], values["MAXMORDINPUT"]
    )
    return inputs


def solar_source(case) -> Path:
    """Check explicit illumination against the native reader's actual precedence."""
    data = case.data
    auxiliary = (
        Path(data["root"]) / "Data" if "root" in data else Path(data["auxiliary"])
    )
    tsis = auxiliary / "tsis1_hsrs_binned_fs.csv"
    chosen = tsis if tsis.is_file() else auxiliary / "f0.txt"
    source = case.illumination["source"]
    if source != "pace_native" and source != (
        "tsis" if chosen == tsis else "thuillier"
    ):
        raise ValueError(
            f"Requested illumination {source!r}, but the native solar reader selects {chosen}; explicit solar-file selection needs a native connection"
        )
    if not chosen.is_file():
        raise FileNotFoundError(chosen)
    return chosen.resolve()
