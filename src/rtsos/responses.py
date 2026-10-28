"""Exact source response samples, separate from spectral calculation methods."""

from dataclasses import dataclass
from pathlib import Path
import tomllib

from netCDF4 import Dataset
import numpy as np


@dataclass
class ChannelResponse:
    """One complete channel calibration variant, without sample modification."""

    id: str
    wavelength: np.ndarray
    response: np.ndarray
    wavelength_unit: str
    nominal_label: str = ""
    detector: str = ""
    mirror: int | str | None = None
    calibration_center_wavelength: float | None = None

    def __post_init__(self):
        self.wavelength = np.asanyarray(self.wavelength).copy()
        self.response = np.asanyarray(self.response).copy()
        if self.wavelength_unit not in ("nm", "um"):
            raise ValueError("wavelength_unit must be explicit: nm or um")
        if (
            not self.id
            or self.wavelength.ndim != 1
            or self.wavelength.size == 0
            or self.response.shape != self.wavelength.shape
        ):
            raise ValueError("A channel needs an ID and paired nonempty 1-D arrays")


@dataclass
class Responses:
    channels: tuple[ChannelResponse, ...]
    source: str

    def __post_init__(self):
        identities = [(c.id, c.detector, c.mirror) for c in self.channels]
        if not identities or len(set(identities)) != len(identities):
            raise ValueError(
                "Response channel/detector/mirror identities must be unique"
            )

    def select(self, ids=None, *, detector=None, mirror=None) -> "Responses":
        """Select whole channels in requested ID order and source variant order.

        With no IDs, preserve source order. A missing calibration is an error;
        another detector or mirror is never substituted or interpolated.
        """
        available = [
            c
            for c in self.channels
            if (detector is None or c.detector == detector)
            and (mirror is None or c.mirror == mirror)
        ]
        if ids is None:
            selected = available
        else:
            if isinstance(ids, str) or not ids or len(set(ids)) != len(ids):
                raise ValueError("channels must be a nonempty list of unique IDs")
            selected = []
            for identifier in ids:
                variants = [c for c in available if c.id == identifier]
                if not variants:
                    raise ValueError(
                        f"No supplied calibration for channel {identifier!r} and selected variants"
                    )
                selected.extend(variants)
        if not selected:
            raise ValueError("No supplied calibration for the selected detector/mirror")
        return Responses(tuple(selected), self.source)


def load_responses(path: str | Path) -> Responses:
    """Read [[channels]] TOML arrays, preserving sample values and ordering."""
    path = Path(path).expanduser().resolve()
    with path.open("rb") as stream:
        tables = tomllib.load(stream)
    return Responses(
        tuple(ChannelResponse(**channel) for channel in tables["channels"]), str(path)
    )


def load_oci_rsr(path: str | Path) -> Responses:
    """Read every OCI LUT detector, band and mirror without native trimming.

    IDs are one-based within detectors, e.g. oci_blue_018, oci_swir_001.
    Mirror indices are zero-based source indices. Calibration centers are kept
    separately from sample wavelengths; they are not effective wavelengths.
    Source masks are retained; a missing calibration center is represented by None.
    """
    path = Path(path).expanduser().resolve()
    channels = []
    with Dataset(path) as source:
        for detector in ("blue", "red", "SWIR"):
            group = source.groups[detector]
            unit = group["sampled_wavelength"].getncattr("units")
            wavelength = group["sampled_wavelength"][:]
            centers = group["center_wavelength"][:]
            for index in range(centers.shape[0]):
                for mirror in range(centers.shape[1]):
                    center = centers[index, mirror]
                    channels.append(
                        ChannelResponse(
                            id=f"oci_{detector.lower()}_{index + 1:03d}",
                            wavelength=wavelength,
                            response=group["RSR"][index, :, mirror],
                            wavelength_unit=unit,
                            detector=detector.lower(),
                            mirror=mirror,
                            calibration_center_wavelength=(
                                None if np.ma.is_masked(center) else center.item()
                            ),
                        )
                    )
    return Responses(tuple(channels), str(path))
