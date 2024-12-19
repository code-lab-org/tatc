# TAT-C Change Log

## 3.3.0

Added:
- Schema for streets-of-coverage constellation `SOCConstellation`.
- Schemas for eccentric high-altitude orbits: `MolniyaOrbit` and `TundraOrbit`.
- Optional dependency `cartopy` for examples.

Changed:
- Signature for functions `collect_orbit_track`, `collect_ground_track`, `compute_ground_track`, and `collect_observations` passes in `instrument_index` rather than the `instrument` object, with a default value of 0 (first instrument).
- Signature for function `collect_ro_observations` replaces `max_azimuth` with `max_yaw` constraint.
- Signature for private function `_collect_ro_series` passes in receiver vertical, normal, binormal (VNB) unit vectors.
- Examples use `cartopy` for coastlines and stock imagery.