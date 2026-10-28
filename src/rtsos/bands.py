"""Full-sample response integration of physical spectral radiance/irradiance.

Composite trapezoids act on the sampled product R(lambda)*X(lambda), with
wavelength in micrometres. Band X = trapz(R*X, lambda)/trapz(R, lambda).
The source calibration is never modified by this mathematical normalization.
"""

from dataclasses import replace

import numpy as np

from . import monochromatic
from .physical import _case_keys


class ResponseIntegration(monochromatic.Monochromatic):
    """Consume the shared native node results before Program publishes a band file."""

    def __init__(self, responses):
        super().__init__()
        self.responses = responses

    def result_dataset(self, dataset):
        return integrate(dataset, self.responses)


def integration_grid(channel):
    """Validate an executable complete calibration, without repairing its data."""
    for name, values in (
        ("wavelength", channel.wavelength),
        ("response", channel.response),
    ):
        if (
            np.any(np.ma.getmaskarray(values))
            or values.dtype.kind not in "iuf"
            or not np.all(np.isfinite(values))
        ):
            raise ValueError(
                f"{channel.id}: {name} samples must be numeric, finite and unmasked"
            )
    wave_um = np.asarray(channel.wavelength, dtype=np.float64) * (
        0.001 if channel.wavelength_unit == "nm" else 1.0
    )
    if wave_um.size < 2 or np.any(wave_um <= 0) or np.any(np.diff(wave_um) <= 0):
        raise ValueError(
            f"{channel.id}: at least two positive strictly increasing wavelengths are required"
        )
    response = np.asarray(channel.response, dtype=np.float64)
    area = np.trapezoid(response, wave_um)
    absolute_area = np.trapezoid(np.abs(response), wave_um)
    if not np.isfinite(area) or abs(area) <= np.finfo(np.float64).eps * absolute_area:
        raise ValueError(
            f"{channel.id}: zero or numerically degenerate signed response area"
        )
    return wave_um, area


def case_inputs(case, responses):
    _case_keys(case.spectral, {"method"}, "response spectral", required={"method"})
    _case_keys(
        case.sensor,
        {"name", "response_file", "response_format", "channels", "detector", "mirror"},
        "response sensor",
        required={"response_file"},
    )
    for channel in responses.channels:
        integration_grid(channel)
    # Convert directly to native nm, independently of the quadrature's um measure.
    # Reuse only equal float64 physical wavelengths, never nearby nodes.
    wavelengths = np.unique(
        np.concatenate(
            [
                np.asarray(c.wavelength, dtype=np.float64)
                * (1000 if c.wavelength_unit == "um" else 1)
                for c in responses.channels
            ]
        )
    )
    physical = replace(
        case,
        sensor={"name": case.sensor.get("name", "custom")},
        spectral={
            "method": "monochromatic",
            "wavelength_unit": "nm",
            "wavelengths": wavelengths,
        },
    )
    return monochromatic.case_inputs(physical)


def integrate(nodes, responses):
    """Keep unique-node diagnostics; publish band products on ordered variants."""
    result = nodes.rename(request="spectral_node")
    channels = responses.channels
    result = result.assign_coords(
        channel=np.arange(len(channels), dtype=np.int64),
        channel_id=("channel", [c.id for c in channels]),
        detector=("channel", [c.detector for c in channels]),
        mirror=(
            "channel",
            ["" if c.mirror is None else str(c.mirror) for c in channels],
        ),
        nominal_label=("channel", [c.nominal_label for c in channels]),
        calibration_center_wavelength=(
            "channel",
            [
                (
                    np.nan
                    if c.calibration_center_wavelength is None
                    else c.calibration_center_wavelength
                )
                for c in channels
            ],
        ),
        calibration_wavelength_unit=("channel", [c.wavelength_unit for c in channels]),
    )
    fields = [
        name
        for name in nodes.data_vars
        if name.startswith(("Stokes_Vector_", "Irradiance_"))
        or name == "Solar_Irradiance"
    ]
    band_fields = {name: [] for name in fields}
    areas = []
    for channel in channels:
        wave_um, area = integration_grid(channel)
        wave_nm = np.asarray(channel.wavelength, dtype=np.float64) * (
            1000 if channel.wavelength_unit == "um" else 1
        )
        indices = nodes.wavelength_nm.to_index().get_indexer(wave_nm)
        if np.any(indices < 0):
            raise ValueError(
                f"{channel.id}: a response node has no evaluated physical result"
            )
        areas.append(area)
        for name in fields:
            values = nodes[name].isel(request=indices).transpose("request", ...)
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{channel.id}: nonfinite native {name}")
            weights = np.asarray(channel.response, dtype=np.float64).reshape(
                (-1,) + (1,) * (values.ndim - 1)
            )
            band = np.trapezoid(values.values * weights, wave_um, axis=0) / area
            if not np.all(np.isfinite(band)):
                raise ValueError(f"{channel.id}: nonfinite integrated {name}")
            band_fields[name].append(band)
    for name in fields:
        # Node optical depths and radiometric products remain explicit diagnostics.
        result = result.rename({name: "Node_" + name})
        dims = ("channel",) + tuple(d for d in nodes[name].dims if d != "request")
        result[name] = (dims, np.stack(band_fields[name]), nodes[name].attrs)
    result["response_area"] = (
        "channel",
        areas,
        {"units": "um", "long_name": "signed trapezoidal response area"},
    )
    result.attrs.update(
        response_source=responses.source,
        spectral_quadrature="composite trapezoid of sampled R*X in um, divided by signed trapezoid of R in um; every source sample",
        native_unique_wavelength_count=nodes.sizes["request"],
        response_sample_count=sum(c.wavelength.size for c in channels),
        optical_depth_convention="Tau_* are monochromatic spectral_node diagnostics; no equivalent band optical depth is defined",
    )
    return result


def write_calibrations(output, responses):
    """Retain exact source arrays/dtypes and explicit masks in NetCDF groups."""
    group = output.createGroup("calibration")
    group.source = responses.source
    for index, channel in enumerate(responses.channels):
        variant = group.createGroup(str(index))
        variant.setncatts(
            dict(
                channel_id=channel.id,
                detector=channel.detector,
                nominal_label=channel.nominal_label,
                wavelength_unit=channel.wavelength_unit,
            )
        )
        if channel.mirror is not None:
            variant.mirror = channel.mirror
        if channel.calibration_center_wavelength is not None:
            variant.calibration_center_wavelength = (
                channel.calibration_center_wavelength
            )
        variant.createDimension("sample", channel.wavelength.size)
        for name, values in (
            ("wavelength", channel.wavelength),
            ("response", channel.response),
        ):
            variant.createVariable(name, values.dtype, ("sample",), fill_value=False)[
                :
            ] = np.ma.getdata(values)
            variant.createVariable(name + "_mask", "u1", ("sample",))[:] = (
                np.ma.getmaskarray(values)
            )
