import pytest

shapely = pytest.importorskip("shapely")
from shapely.geometry import LineString, Polygon
from gerber_comparator.geometry import normalize_geometry, safe_symmetric_difference, safe_union

def test_invalid_polygon_is_repaired_to_valid_polygonal_geometry():
    bow_tie = Polygon([(0, 0), (2, 2), (0, 2), (2, 0), (0, 0)])
    repaired = normalize_geometry(bow_tie)
    assert repaired.is_valid
    assert repaired.area > 0

def test_touching_polygons_and_overlapping_strokes_union_safely():
    touching = safe_union([Polygon([(0,0),(1,0),(1,1),(0,1)]), Polygon([(1,0),(2,0),(2,1),(1,1)])])
    strokes = safe_union([LineString([(0,0),(2,0)]).buffer(.25), LineString([(1,0),(3,0)]).buffer(.25)])
    assert touching.is_valid and touching.area == pytest.approx(2)
    assert strokes.is_valid and strokes.area > 0

def test_invalid_geometries_produce_valid_vector_xor():
    original = Polygon([(0,0),(2,2),(0,2),(2,0),(0,0)])
    working = Polygon([(0,0),(3,3),(0,3),(3,0),(0,0)])
    xor = safe_symmetric_difference(original, working)
    assert xor.is_valid
    assert xor.area > 0
