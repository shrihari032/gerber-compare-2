import pytest

shapely = pytest.importorskip("shapely")
from shapely.geometry import GeometryCollection, LineString, MultiLineString, Point, Polygon
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


def test_polygonal_members_of_a_geometry_collection_are_retained():
    warnings = []
    geometry = GeometryCollection([
        Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]),
        LineString([(0, 0), (2, 2)]),
    ])

    normalized = normalize_geometry(geometry, label="mixed repair", warnings=warnings)

    assert normalized.is_valid
    assert normalized.area == pytest.approx(4)
    assert any("polygonal copper components were retained" in warning for warning in warnings)


@pytest.mark.parametrize("geometry", [
    GeometryCollection([LineString([(0, 0), (1, 1)])]),
    GeometryCollection([Point(1, 1)]),
    MultiLineString([[(0, 0), (1, 1)], [(1, 1), (2, 1)]]),
])
def test_non_area_geometry_is_explicitly_normalized_to_empty_polygon(geometry):
    warnings = []

    normalized = normalize_geometry(geometry, label="primitive 509 repair", warnings=warnings)

    assert normalized.is_empty
    assert normalized.geom_type == "Polygon"
    assert any("primitive 509 repair" in warning for warning in warnings)
    assert any("zero-area geometry" in warning for warning in warnings)
