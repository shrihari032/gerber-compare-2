import pytest

pytest.importorskip("shapely")
from gerber_comparator.core import _local_deviation, compare_gerbers
from gerber_comparator.model import ComparisonConfig

def gerber(aperture="R,2X2", flashes=((1, 1),)):
    commands = ["%MOMM*%", "%FSLAX24Y24*%", f"%ADD10{aperture}*%", "D10*"]
    commands.extend(f"X{round(x * 10000):06d}Y{round(y * 10000):06d}D03*" for x, y in flashes)
    return "".join(commands) + "M02*"

def compare(original, working, tolerance=.05):
    return compare_gerbers(original, working, ComparisonConfig(geometric_tolerance_mm=tolerance, auto_alignment=False, snapshot_generation=False))

def test_identical_gerbers_pass_with_empty_raw_xor():
    result = compare(gerber(), gerber())
    assert result.overall_result == "PASS"
    assert result.raw_xor.area == pytest.approx(0)

def test_small_change_is_ignored_and_large_change_is_flagged():
    ignored = compare(gerber(), gerber("R,2.02X2"))
    flagged = compare(gerber(), gerber("R,2.2X2"))
    assert {region.classification for region in ignored.regions} == {"IGNORE"}
    assert {region.classification for region in flagged.regions} == {"FLAG"}

def test_component_count_topology_change_is_flagged():
    result = compare(gerber(), gerber(flashes=((1, 1), (5, 1))))
    assert result.regions
    assert {region.classification_reason for region in result.regions} == {"TOPOLOGY_CHANGE"}


def test_complete_object_removal_is_directional():
    result = compare(gerber(), gerber(flashes=()))

    assert not result.missing_geometry.is_empty
    assert result.added_geometry.is_empty
    assert {region.classification_reason for region in result.regions} == {"MISSING_FROM_WORKING"}
    assert result.statistics()["missing_from_working"] == 1


def test_complete_object_addition_is_directional():
    result = compare(gerber(flashes=()), gerber())

    assert result.missing_geometry.is_empty
    assert not result.added_geometry.is_empty
    assert {region.classification_reason for region in result.regions} == {"ADDED_IN_WORKING"}
    assert result.statistics()["added_in_working"] == 1


def test_partial_removal_reports_only_original_minus_working():
    result = compare(gerber("R,2X2"), gerber("R,1X2"))

    assert result.missing_geometry.area == pytest.approx(2)
    assert result.added_geometry.is_empty
    assert {region.classification_reason for region in result.regions} == {"MISSING_FROM_WORKING"}


def test_local_deviation_handles_identical_and_empty_neighbourhoods():
    from shapely.geometry import Polygon

    geometry = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    assert _local_deviation(geometry, geometry, geometry, 0.1) == pytest.approx(0)
    assert _local_deviation(geometry, geometry, Polygon(), 0.1) == float("inf")
