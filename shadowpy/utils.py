"""
Module containing utility functions for the ShadowPy package.
"""

import numpy as np
import Shadow
from attr import dataclass

# from .optical_elements import OpticalElement
# ↓ Moving away from mixing data acquistion and processing
# from caxscripts.image_statistics import Histogram2DAnalyzer


def rotation_matrix(axis, angle):
    """
    Generate a rotation matrix for a given axis and angle.
    
    Parameters:
        axis (str): The axis of rotation ('x', 'y', or 'z').
        angle (float): The angle of rotation in degrees.
    
    Returns:
        np.ndarray: The 3x3 rotation matrix.
    """
    angle = np.radians(angle)
    c = np.cos(angle)
    s = np.sin(angle)

    if axis == 'x':
        return np.array([[1, 0, 0],
                         [0, c, -s],
                         [0, s, c]])
    elif axis == 'y':
        return np.array([[c, 0, s],
                         [0, 1, 0],
                         [-s, 0, c]])
    elif axis == 'z':
        return np.array([[c, -s, 0],
                         [s, c, 0],
                         [0, 0, 1]])
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'.")
    

class ReferenceFrame:
    """
    Right-handed orthonormal reference frame.

    Convention:
        v_lab = R @ v_local

    Columns of R are the local basis vectors expressed
    in lab coordinates.
    """

    def __init__(self, origin=None, orientation=None, name=None):

        self.name = name or "frame"

        self.origin = np.array(
            origin if origin is not None else [0., 0., 0.],
            dtype=float
        )

        self.R = np.array(
            orientation if orientation is not None else np.eye(3),
            dtype=float
        )

        self._validate_rotation_matrix()

    # -------------------------------------------------
    # Validation
    # -------------------------------------------------

    def _validate_rotation_matrix(self, atol=1e-10):

        should_be_identity = self.R.T @ self.R

        if not np.allclose(should_be_identity, np.eye(3), atol=atol):
            raise ValueError("Orientation matrix is not orthogonal")

        det = np.linalg.det(self.R)

        if not np.isclose(det, 1.0, atol=atol):
            raise ValueError("Orientation matrix must have determinant +1")

    # -------------------------------------------------
    # Vector transforms
    # -------------------------------------------------

    def vector_to_lab(self, v_local):
        v_local = np.asarray(v_local, dtype=float)
        return self.R @ v_local

    def vector_from_lab(self, v_lab):
        v_lab = np.asarray(v_lab, dtype=float)
        return self.R.T @ v_lab

    # -------------------------------------------------
    # Point transforms
    # -------------------------------------------------

    def point_to_lab(self, p_local):
        p_local = np.asarray(p_local, dtype=float)
        return self.origin + self.R @ p_local

    def point_from_lab(self, p_lab):
        p_lab = np.asarray(p_lab, dtype=float)
        return self.R.T @ (p_lab - self.origin)

    # -------------------------------------------------
    # Relative frame construction
    # -------------------------------------------------

    def child_frame(self, relative_origin, relative_rotation, name=None):
        """
        Create a new frame defined relative to THIS frame.
        """

        relative_origin = np.asarray(relative_origin, dtype=float)
        relative_rotation = np.asarray(relative_rotation, dtype=float)

        new_origin = self.point_to_lab(relative_origin)
        new_rotation = self.R @ relative_rotation

        return ReferenceFrame(
            origin=new_origin,
            orientation=new_rotation,
            name=name
        )

    # -------------------------------------------------
    # Relations between frames
    # -------------------------------------------------

    def rotation_to(self, other):
        """
        Rotation matrix mapping vectors from this frame
        into another frame.
        """

        return other.R.T @ self.R

    def __repr__(self):
        return (
            f"ReferenceFrame(name={self.name}, "
            f"origin={self.origin})"
        )
    
# ==== BEAM IMAGE SAVING FUNCTION ====

def _calc_sigma(bin_edges, histogram, axis=0):
    """
    Calculate the weighted standard deviation of a histogram projection.
    
    Parameters:
        bin_edges (np.ndarray): The edges of the bins.
        histogram (np.ndarray): The histogram counts.
        axis (int): The axis along which to calculate the standard deviation. 
                    Default is 0.
    
    Returns:
        float: The standard deviation of the histogram.
    """
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    proj = np.sum(histogram, axis=axis)
    total = np.sum(proj)

    if total == 0:
        return 0.0

    mask = proj > proj.max()/2
    lengts = bin_centers[mask]
    fwhm = np.ptp(lengts)  # Full width at half maximum (FWHM)
    return fwhm / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to standard deviation

@dataclass
class Image:
    histogram: np.ndarray
    bin_h_edges: np.ndarray
    bin_v_edges: np.ndarray
    bin_h_centers: np.ndarray
    bin_v_centers: np.ndarray
    sigma_h: float
    sigma_v: float
    centroid_h: float
    centroid_v: float

    @classmethod
    def from_ticket(cls, ticket: dict):
        histogram = ticket.get("histogram", np.empty((0, 0)))
        bin_h_edges = ticket.get("bin_h_edges", np.array([]))
        bin_v_edges = ticket.get("bin_v_edges", np.array([]))
        bin_h_centers = ticket.get("bin_h_center", np.array([]))
        bin_v_centers = ticket.get("bin_v_center", np.array([]))

        fwhm_h = ticket.get("fwhm_h")
        fwhm_v = ticket.get("fwhm_v")
        coordinates_h = ticket.get("fwhm_coordinates_h")
        coordinates_v = ticket.get("fwhm_coordinates_v")

        sigma_h = np.nan if fwhm_h is None else float(fwhm_h) / (2 * np.sqrt(2 * np.log(2)))
        sigma_v = np.nan if fwhm_v is None else float(fwhm_v) / (2 * np.sqrt(2 * np.log(2)))

        centroid_h = np.nan if coordinates_h is None or len(coordinates_h) == 0 else float(np.mean(coordinates_h))
        centroid_v = np.nan if coordinates_v is None or len(coordinates_v) == 0 else float(np.mean(coordinates_v))

        return cls(histogram=np.asarray(histogram),
                   bin_h_edges=np.asarray(bin_h_edges, dtype=float),
                   bin_v_edges=np.asarray(bin_v_edges, dtype=float),
                   bin_h_centers=np.asarray(bin_h_centers, dtype=float),
                   bin_v_centers=np.asarray(bin_v_centers, dtype=float),
                   sigma_h=sigma_h,
                   sigma_v=sigma_v,
                   centroid_h=centroid_h,
                   centroid_v=centroid_v)

    @property
    def bin_area(self) -> float:
        """Physical area of one histogram bin."""
        dh = np.diff(self.bin_h_edges)
        dv = np.diff(self.bin_v_edges)
        if dh.size == 0 or dv.size == 0:
            return 0.0
        return float(np.median(dh) * np.median(dv))

    @property
    def beam_area(self) -> float:
        """Gaussian-equivalent beam area, 2*pi*sigma_h*sigma_v."""
        return float(2.0 * np.pi * self.sigma_h * self.sigma_v)

    @property
    def peak(self) -> float:
        """Maximum histogram occupancy."""
        if self.histogram.size == 0:
            return 0.0
        return float(np.max(self.histogram))

    @property
    def peak_density(self) -> float:
        """Maximum histogram occupancy per unit area."""
        area = self.bin_area
        if area <= 0:
            return 0.0
        return self.peak / area

    # @property
    # def sigma_h(self):
    #     return _calc_sigma(self.bin_h_edges, self.histogram, axis=1)

    # @property
    # def sigma_v(self):
    #     return _calc_sigma(self.bin_v_edges, self.histogram, axis=0)


def _robust_beam_range(beam: Shadow.Beam, col: int,
                       quantiles: tuple = (0.001, 0.999),
                       pad_fraction: float = 0.05,
                       nolost: int = 1):
    """
    Return a robust, translation-invariant coordinate range for a beam column.

    The range follows the requested quantiles of the selected rays and is padded
    by a fraction of the robust span. It does not depend on the absolute beam
    position, so translating the beam does not change the window width.
    """
    if not 0.0 <= quantiles[0] < quantiles[1] <= 1.0:
        raise ValueError("quantiles must satisfy 0 <= low < high <= 1")
    if pad_fraction < 0.0:
        raise ValueError("pad_fraction must be non-negative")

    values = np.asarray(beam.getshonecol(col, nolost=nolost), dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return [-1.0, 1.0]

    lo, hi = np.quantile(values, quantiles)
    span = hi - lo

    if span <= 0.0:
        scale = float(np.std(values))
        span = 10.0 * scale if scale > 0.0 else 2.0e-9
        center = float(np.median(values))
        lo = center - span / 2.0
        hi = center + span / 2.0
        span = hi - lo

    margin = pad_fraction * span
    return [float(lo - margin), float(hi + margin)]


def image_from_beam(beam: Shadow.Beam, frame: ReferenceFrame,
                   image_size: tuple = (210, 280),
                   range_mode: str = "robust",
                   quantiles: tuple = (0.0001, 0.9999),
                   pad_fraction: float = 0.10):
        """Build an :class:`Image` from a beam and its reference frame.

        Parameters
        ----------
        beam : Shadow.Beam
            Beam to histogram.
        frame : ReferenceFrame
            Reference frame defining the beam image axes.
        image_size : tuple
            Number of bins as ``(nbins_v, nbins_h)``.
        range_mode : str
            ``"robust"`` uses quantile-based, beam-relative ranges. ``"shadow"``
            preserves the legacy ``Beam.get_good_range`` behavior.
        quantiles : tuple
            Low and high quantiles used by ``"robust"`` mode.
        pad_fraction : float
            Fraction of the robust span added to each side in ``"robust"`` mode.
        """
        if frame is None:
            raise ValueError("frame must be provided to derive the beam image axes.")

        lab_x = np.array([1, 0, 0])
        lab_y = np.array([0, 1, 0])

        local_x = frame.vector_from_lab(lab_x)
        local_y = frame.vector_from_lab(lab_y)

        x_index = np.argmax(np.abs(local_x)) + 1
        y_index = np.argmax(np.abs(local_y)) + 1

        if range_mode == "robust":
            x_range = _robust_beam_range(beam, x_index, quantiles=quantiles,
                                         pad_fraction=pad_fraction, nolost=1)
            y_range = _robust_beam_range(beam, y_index, quantiles=quantiles,
                                         pad_fraction=pad_fraction, nolost=1)
            histo_range_kwargs = {"xrange": x_range, "yrange": y_range}
        elif range_mode == "shadow":
            x_range = beam.get_good_range(x_index)
            y_range = beam.get_good_range(y_index)
            histo_range_kwargs = {}
        else:
            raise ValueError("range_mode must be 'robust' or 'shadow'")

        nbins_v, nbins_h = image_size
        x_span = x_range[1] - x_range[0]
        y_span = y_range[1] - y_range[0]
        if y_span > x_span:
            nbins_h, nbins_v = nbins_v, nbins_h

        image_ticket = beam.histo2(col_h=x_index, col_v=y_index,
                                  nbins_h=nbins_h, nbins_v=nbins_v, nolost=1,
                                  **histo_range_kwargs)

        return Image.from_ticket(image_ticket)


def save_image(element, beam: Shadow.Beam, image_size: tuple = (210, 280),
               range_mode: str = "robust",
               quantiles: tuple = (0.0001, 0.9999),
               pad_fraction: float = 0.10):
        """Backward-compatible wrapper around :func:`image_from_beam`."""
        if element is None:
            raise ValueError("element must be provided to save an image.")
        if getattr(element, "frame", None) is None:
            raise ValueError("element.frame must be defined to save a beam image.")
        return image_from_beam(beam, element.frame, image_size=image_size,
                               range_mode=range_mode, quantiles=quantiles,
                               pad_fraction=pad_fraction)