# Development

Run `pytest -q` and `python -m compileall gerber_comparator` before committing. Tests currently cover parser normalization and safe handling of vendor apertures. Extend the suite with physical geometry fixtures for every newly supported command.

The project depends on Shapely for robust vector operations and Matplotlib only for snapshots. Do not replace vector comparison with pixel differences. Do not add unknown-aperture fallbacks: implement a named resolver/provider with tests instead.
