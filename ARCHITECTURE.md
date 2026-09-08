# Architecture

The pipeline is `RS-274X input -> ParsedLayer -> Shapely normalized geometry (mm) -> alignment -> raw symmetric difference -> region analysis -> export adapters`. Raw XOR is always generated before tolerance classification. Region significance uses symmetric boundary Hausdorff distance, translation magnitude, and component-count topology checks; XOR area is only reported.

`Aperture` stores identity, source name, numeric parameters, source definition, resolved type/provider, and macro linkage. A resolver is intentionally conservative: unsupported vendors remain unresolved and affect the final result rather than receiving invented geometry.

The core has no Colab imports. V2 can directly call `compare_gerbers(original, working, config)` and render the structured `ComparisonResult` in a desktop workspace.
