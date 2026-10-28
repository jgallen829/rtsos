"""Shared physical-case translation; scalar defaults retain their existing owner."""

from pathlib import Path

import numpy as np
import xarray as xr

from .parameters import Parameters as P

# Common-case translation belongs to this scientific adapter, not Program.
LAND_MODELS = {"ross_li": -203, "snow": -202, "mrpv": -201, "lambertian_fresnel": -200}
WATER_MODELS = {
    "none": 0,
    "chlorophyll": 1,
    "chlorophyll_sediment": 2,
    "seven_parameter": 3,
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


def physical_case_inputs(case, *, channel_group=0, monochromatic=False):
    """Shared public scene/geometry/process translation with PACE scalar defaults."""
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

    from .PACE_Simulator_DoubleK import init

    defaults = init().defaults.isel(Solar_Zenith_Angle=0, tau_ref_hi=0, chla=0)
    values = {name: coord.item() for name, coord in defaults.coords.items()}
    values.update(Aux_Dir=auxiliary, Atmos_Dir=gases, WAVEBAND_SEG_FLAG=channel_group)
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
            "Scalar scattering needs a native control connection; this path uses vector RT"
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
        if monochromatic and dtype.kind == "f":
            dtype = np.dtype("float64")
        if dtype.kind not in "if":
            continue
        array = np.asarray(value)
        if array.ndim or array.dtype.kind not in "iufb" or not np.isfinite(array):
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
    if monochromatic:
        # Native physical records are double precision; preserve explicit case values.
        return xr.Dataset(
            coords={
                name: (
                    (),
                    np.asarray(
                        value,
                        dtype=(
                            "float64"
                            if defaults[name].dtype.kind == "f"
                            else (
                                defaults[name].dtype
                                if defaults[name].dtype.kind == "i"
                                else None
                            )
                        ),
                    ),
                    defaults[name].attrs,
                )
                for name, value in values.items()
            }
        )
    inputs = P.make_dataset({getattr(P, name): value for name, value in values.items()})
    # The legacy driver passes this upper bound to RTSOS at each wavelength.
    # Keep the requested MAXMORDINPUT intact and expose its internal restriction.
    inputs.attrs["pace_rt_fourier_order_limit"] = min(
        30, values["nquadainput"], values["MAXMORDINPUT"]
    )
    return inputs
