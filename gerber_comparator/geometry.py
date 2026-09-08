"""Physical-geometry resolution and controlled GEOS validity recovery."""
from __future__ import annotations
from math import cos, pi, sin
from typing import Iterable
from .model import Aperture, ParsedLayer

class GeometryError(RuntimeError):
    """A vector geometry could not be safely made valid for comparison."""

def _shapely():
    try:
        from shapely.geometry import LineString, Point, Polygon
        from shapely.ops import unary_union
        return LineString, Point, Polygon, unary_union
    except ImportError as exc:
        raise RuntimeError("Geometry support requires Shapely. Install project dependencies before comparing.") from exc

def _make_valid(geometry):
    """Use the strongest available Shapely repair API, with a documented fallback."""
    try:
        from shapely import make_valid
        return make_valid(geometry), "make_valid"
    except ImportError:
        try:
            from shapely.validation import make_valid
            return make_valid(geometry), "make_valid"
        except ImportError:
            # GEOS buffer(0) is only a compatibility fallback for old Shapely.
            return geometry.buffer(0), "buffer(0) compatibility fallback"

def _polygonal_components(geometry):
    """Retain every polygonal component; lines/points have no copper area."""
    if geometry.is_empty: return []
    if geometry.geom_type == "Polygon": return [geometry]
    if geometry.geom_type in {"MultiPolygon", "GeometryCollection"}:
        parts = []
        for member in geometry.geoms: parts.extend(_polygonal_components(member))
        return parts
    return []


def _polygonal_geometry(geometry, *, label: str):
    """Return all polygonal content from a validity-repair result.

    ``make_valid`` is allowed to return a GeometryCollection.  Gerber copper is
    area geometry, so retaining its polygonal members is lossless for the
    comparison model while avoiding lines left behind by a self-intersection.
    """
    _, _, Polygon, unary_union = _shapely()
    polygons = _polygonal_components(geometry)
    if not polygons:
        if geometry.is_empty:
            return Polygon()
        raise GeometryError(f"{label} contains no polygonal copper geometry.")
    try:
        return unary_union(polygons)
    except Exception as exc:
        raise GeometryError(f"{label} polygonal components could not be unioned safely.") from exc


def _precision_normalize(geometry, *, label: str):
    """Apply the last-resort 1 nm precision grid used for GEOS recovery only."""
    try:
        from shapely import set_precision
    except ImportError as exc:
        raise GeometryError(
            f"{label} needs precision recovery, but this Shapely version has no set_precision support."
        ) from exc
    try:
        # Coordinates are normalized to millimetres.  One nanometre is far
        # below the supported manufacturing tolerances and is never a normal
        # comparison path.
        return set_precision(geometry, 0.000001)
    except Exception as exc:
        raise GeometryError(f"{label} precision normalization failed.") from exc

def normalize_geometry(geometry, *, label: str = "Geometry", warnings: list[str] | None = None):
    """Return valid polygonal geometry or raise a descriptive ``GeometryError``.

    No simplification, scaling, or area threshold is applied. GeometryCollection
    repairs are flattened only to their polygonal copper components.
    """
    _, _, Polygon, _ = _shapely()
    if geometry is None or geometry.is_empty: return Polygon()
    valid_before = geometry.is_valid
    if valid_before:
        return geometry
    repaired, method = _make_valid(geometry)
    normalized = _polygonal_geometry(repaired, label=f"{label} repair")
    if normalized.is_empty or not normalized.is_valid:
        raise GeometryError(f"{label} remains invalid after {method} repair.")
    if warnings is not None:
        warnings.append(f"{label} was invalid and required GEOS {method} repair.")
    return normalized

def safe_union(geometries: Iterable, *, label: str = "Geometry", warnings: list[str] | None = None):
    """Union vector copper safely, retrying normalized inputs after a GEOS error."""
    _, _, Polygon, unary_union = _shapely()
    inputs = [geometry for geometry in geometries if geometry is not None and not geometry.is_empty]
    if not inputs: return Polygon()

    # Validate every primitive before asking GEOS to overlay them.  In
    # particular, this keeps a single bow-tie region from poisoning the union.
    normalized_inputs = [
        normalize_geometry(item, label=f"{label} primitive {index}", warnings=warnings)
        for index, item in enumerate(inputs, 1)
    ]
    try:
        return normalize_geometry(unary_union(normalized_inputs), label=label, warnings=warnings)
    except Exception:
        # A second pass is intentional: GEOS can still reject two individually
        # valid polygons with nearly coincident boundaries.
        repaired = [normalize_geometry(item, label=f"{label} primitive {index}", warnings=warnings) for index, item in enumerate(normalized_inputs, 1)]
        try:
            result = normalize_geometry(unary_union(repaired), label=label, warnings=warnings)
        except Exception as exc:
            raise GeometryError(f"{label} union failed after individual geometry repair.") from exc
        if warnings is not None:
            warnings.append(f"{label} union required retry after a GEOS topology failure.")
        return result

def safe_symmetric_difference(original, working, *, warnings: list[str] | None = None):
    """Calculate complete vector XOR after validation, with precision retry only on GEOS failure."""
    original = normalize_geometry(original, label="Original geometry", warnings=warnings)
    working = normalize_geometry(working, label="Working geometry", warnings=warnings)
    try:
        return normalize_geometry(original.symmetric_difference(working), label="Raw XOR geometry", warnings=warnings)
    except Exception:
        precision_original = _precision_normalize(original, label="Original geometry")
        precision_working = _precision_normalize(working, label="Working geometry")
        try:
            xor = precision_original.symmetric_difference(precision_working)
            result = normalize_geometry(xor, label="Raw XOR geometry", warnings=warnings)
        except Exception as exc:
            raise GeometryError("Raw XOR failed after repair and controlled precision normalization.") from exc
        if warnings is not None:
            warnings.append("XOR required 0.000001 mm precision normalization after a GEOS topology exception.")
        return result

def aperture_shape(aperture: Aperture):
    LineString, Point, Polygon, _ = _shapely()
    if aperture.geometry_type == "circle": return Point(0, 0).buffer(aperture.parameters[0] / 2)
    if aperture.geometry_type == "rectangle":
        w, h = aperture.parameters[:2]; return Polygon([(-w/2,-h/2),(w/2,-h/2),(w/2,h/2),(-w/2,h/2)])
    if aperture.geometry_type == "oblong":
        w, h = aperture.parameters[:2]
        if w >= h: return LineString([(-((w-h)/2),0),((w-h)/2,0)]).buffer(h/2)
        return LineString([(0,-((h-w)/2)),(0,((h-w)/2))]).buffer(w/2)
    if aperture.geometry_type == "polygon":
        diameter, vertices = aperture.parameters[:2]; r = diameter / 2
        return Polygon([(r*cos(2*pi*i/int(vertices)), r*sin(2*pi*i/int(vertices))) for i in range(int(vertices))])
    raise ValueError(f"No geometry provider available for aperture {aperture.number} ({aperture.name})")

def resolve_geometry(layer: ParsedLayer):
    """Resolve primitives then validate/repair their final normalized-mm union."""
    LineString, _, Polygon, _ = _shapely(); dark = []
    for primitive in layer.primitives:
        kind = primitive[0]
        if kind == "region": dark.append(Polygon(primitive[1])); continue
        aperture = layer.apertures.get(primitive[1])
        if aperture is None:
            layer.warnings.append(f"Primitive uses undefined aperture D{primitive[1]}"); continue
        if aperture.geometry_type is None:
            if aperture not in layer.unresolved_apertures: layer.unresolved_apertures.append(aperture)
            continue
        shape = aperture_shape(aperture)
        if kind == "flash":
            from shapely import affinity
            dark.append(affinity.translate(shape, primitive[2], primitive[3]))
        else: dark.append(LineString([(primitive[2], primitive[3]), (primitive[4], primitive[5])]).buffer(aperture.parameters[0]/2))
    layer.geometry_valid_before = all(geometry.is_valid for geometry in dark)
    result = safe_union(dark, label=f"{layer.filename} geometry", warnings=layer.warnings)
    layer.geometry_valid_after = result.is_valid
    layer.geometry_type = result.geom_type
    return result
