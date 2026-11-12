# scenario-shaking-exposure

Ensemble regression of seismic shaking: from rupture geometry to exposure, by tree and by forest — a weighted ground-motion logic tree, lithology-conditioned site response and spatially blocked Random Forest residual learning over an intracontinental collision domain.

Four ruptures — two crustal strike-slip, one blind oblique thrust, one intermediate-depth intraslab — are propagated through a weighted ground-motion logic tree, conditioned on a slope-proxy $V_{S30}$ field tested against mapped Quaternary lithology, confronted with recorded ShakeMaps, interrogated by spatially blocked Random Forest residual learning, and projected onto gridded population as macroseismic exposure with aleatory bounds.

## Architecture

| Stage | Scripts |
|---|---|
| Rupture kinematics | `fig_ruptures.sh` → `ruptures.csv` |
| Ground motion (logic tree) | `scenario_gm.py`, `city_gm.py`, `gmm_spread.py` |
| Site conditioning | `site_effect.py`, `fig_vs30_violin.py` |
| Validation and residual learning | `validation.py`, `residual_analysis.py` |
| Exposure | `city_exposure.py` |
| Hazard reconciliation | `psha_compare.py` |
| Cartography and graphics | `fig_*.sh` (GMT), `fig_*.py` (Matplotlib) |

Execution is stage-ordered; each script reads the artefacts of its predecessors from `data/` or `out/` and emits its own. Graphics accept `IN`, `OUT` and `WORK` from the environment.

## Stack

Python ≥ 3.12 · openquake.hazardlib 3.26.2 · NumPy · SciPy · pandas · scikit-learn · rasterio · GeoPandas · Shapely · pyproj · Matplotlib  
GMT 6.6 · GDAL · Ghostscript · curl

```bash
pip install -r requirements.txt
```

Ground-motion models: BSSA14, CZ15, ASB14 (active shallow crust); PEA20, KBCG20 (intraslab). Basin terms via CY14 and CB14 $z_{1.0}$/$z_{2.5}$.

## Data

Open sources only, retrieved on demand: USGS global $V_{S30}$, ComCat, ShakeMap v4 and W-phase solutions; GHSL GHS-POP / GHS-BUILT-S R2023A; GEM GSHM v2023.1; SRTM15+ via GMT remote datasets; GSHHG; PB2002. Third-party licences apply to each.

## Citation

See `CITATION.cff`.

## Licence

MIT.
