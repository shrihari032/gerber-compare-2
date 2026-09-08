"""Comparison orchestration: parse -> normalized vectors -> XOR -> classification."""
from __future__ import annotations
from dataclasses import asdict, dataclass, field
from pathlib import Path
import time
from .geometry import resolve_geometry
from .model import ComparisonConfig, DifferenceRegion, ParsedLayer
from .parser import parse_gerber

@dataclass
class ComparisonResult:
    original: ParsedLayer
    working: ParsedLayer
    config: ComparisonConfig
    original_geometry: object
    working_geometry: object
    raw_xor: object
    flagged_xor: object
    regions: list[DifferenceRegion]
    alignment: dict
    warnings: list[str] = field(default_factory=list)
    timings: dict[str, float] = field(default_factory=dict)
    @property
    def overall_result(self) -> str:
        if any(r.classification == "UNRESOLVED" for r in self.regions) or self.original.unresolved_apertures or self.working.unresolved_apertures: return "UNRESOLVED"
        if any(r.classification == "FLAG" for r in self.regions): return "FAIL"
        return "PASS WITH WARNINGS" if self.warnings else "PASS"
    def statistics(self):
        return {"total_difference_regions": len(self.regions), "flagged_regions": sum(r.classification == "FLAG" for r in self.regions), "ignored_regions": sum(r.classification == "IGNORE" for r in self.regions), "unresolved_regions": sum(r.classification == "UNRESOLVED" for r in self.regions), "total_xor_area_mm2": self.raw_xor.area, "flagged_xor_area_mm2": self.flagged_xor.area}

def _parts(geometry): return list(geometry.geoms) if geometry.geom_type.startswith("Multi") or geometry.geom_type == "GeometryCollection" else ([geometry] if not geometry.is_empty else [])
def compare_gerbers(original: str | Path, working: str | Path, config: ComparisonConfig | None = None) -> ComparisonResult:
    from shapely import affinity
    config = config or ComparisonConfig(); started = time.perf_counter()
    original_layer, working_layer = parse_gerber(original), parse_gerber(working); parsed = time.perf_counter()
    original_geometry, working_geometry = resolve_geometry(original_layer), resolve_geometry(working_layer); geometry_time = time.perf_counter()
    dx, dy, rotation = config.dx_mm, config.dy_mm, config.rotation_deg
    if config.auto_alignment and not original_geometry.is_empty and not working_geometry.is_empty:
        oc, wc = original_geometry.centroid, working_geometry.centroid; dx += oc.x - wc.x; dy += oc.y - wc.y
    aligned_working = affinity.rotate(affinity.translate(working_geometry, dx, dy), rotation, origin="centroid")
    raw_xor = original_geometry.symmetric_difference(aligned_working); xor_time = time.perf_counter()
    regions = []; flagged = []
    topology_changed = len(_parts(original_geometry)) != len(_parts(aligned_working))
    for index, part in enumerate(_parts(raw_xor), 1):
        # Symmetric Hausdorff distance measures physical boundary displacement;
        # XOR area is retained only as a reporting metric.
        deviation = original_geometry.boundary.hausdorff_distance(aligned_working.boundary)
        translation = (dx*dx + dy*dy) ** .5
        if topology_changed and config.topology_check: classification, reason = "FLAG", "TOPOLOGY_CHANGE"
        elif deviation > config.geometric_tolerance_mm: classification, reason = "FLAG", "GEOMETRY_CHANGE"
        elif translation > config.translation_tolerance_mm: classification, reason = "FLAG", "TRANSLATION"
        else: classification, reason = "IGNORE", "BELOW_TOLERANCE"
        if classification == "FLAG": flagged.append(part)
        b = part.bounds; c = part.centroid
        regions.append(DifferenceRegion(f"F-{index:03d}", part, part.area, part.length, b, (c.x,c.y), deviation, translation, topology_changed, classification, reason))
    from shapely.ops import unary_union
    result = ComparisonResult(original_layer, working_layer, config, original_geometry, aligned_working, raw_xor, unary_union(flagged), regions, {"translation_x_mm": dx, "translation_y_mm": dy, "rotation_deg": rotation, "method": "centroid" if config.auto_alignment else "manual"}, original_layer.warnings + working_layer.warnings, {"parse_seconds": parsed-started, "geometry_seconds": geometry_time-parsed, "xor_and_analysis_seconds": time.perf_counter()-xor_time, "total_seconds": time.perf_counter()-started})
    return result
