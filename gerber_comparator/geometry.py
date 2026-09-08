"""Physical-geometry resolution using Shapely; no raster comparison is used."""
from __future__ import annotations
from math import cos, pi, sin
from .model import Aperture, ParsedLayer

def _shapely():
    try:
        from shapely.geometry import LineString, Point, Polygon
        from shapely.ops import unary_union
        return LineString, Point, Polygon, unary_union
    except ImportError as exc:
        raise RuntimeError("Geometry support requires Shapely. Install project dependencies before comparing.") from exc

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
    """Return normalized mm geometry or explicit unresolved-aperture diagnostics."""
    LineString, Point, Polygon, unary_union = _shapely(); dark = []
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
        else:
            width = aperture.parameters[0]
            dark.append(LineString([(primitive[2], primitive[3]), (primitive[4], primitive[5])]).buffer(width/2))
    return unary_union(dark) if dark else Polygon()
