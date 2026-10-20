# Building the TOPMODEL Watershed Explorer

`index.html` at the repository root is a single self-contained file: the model, the data for
11 CAMELS watersheds and the page itself. This folder has everything needed to rebuild it.

| Path | What it is |
|---|---|
| `src/model.js` | Enhanced TOPMODEL core, a line-for-line port of `topmodel_core.R` (plus metrics and SCE-UA) |
| `src/setup.js` | Parameter table, watershed descriptions, starting parameters, forcing and topography helpers |
| `src/app.html` | Page template (layout, charts, calibration worker, land-use scenarios) |
| `data/basins.json` | Daymet forcing, USGS discharge (mm/d) and CAMELS attributes, WY2005–2014 |
| `data/twi.json` | TWI classes, channel-distance histogram and map layers from the 30 m DEMs |
| `data/map.json` | US state outlines and gauge locations for the watershed picker |
| `data/presets/` | SCE-UA results; the best per gauge becomes the "Pre-calibrated" set |
| `scripts/` | The pipeline below |
| `matlab/topmodel_explorer_twi_batch.m` | Runs `twi_watershed_usgs.m` for each gauge at 30 m and writes slim `.mat` files |
| `tests/` | Checks that the JavaScript model reproduces the R model |

## Watersheds

Virginia: 02017500 Johns Creek, 02028500 Rockfish River, 02038850 Holiday Creek, 01669000 Piscataway Creek.
Other climates: 01013500 Fish River (ME), 12040500 Queets River (WA), 11143000 Big Sur River (CA),
11266500 Merced River (CA), 09352900 Vallecito Creek (CO), 09512280 Cave Creek (AZ), 08202700 Seco Creek (TX).
The list lives in `scripts/gauges.json`.

## Pipeline

Run from this `build/` folder. Python needs `numpy`, `pandas` and `scipy`; JavaScript needs Node 18+.

1. **Climate and discharge** (only if the period or gauges change):
   `python scripts/prep_camels.py` clones the needed CAMELS files from three public GitHub
   repositories into `cache/` and writes `data/basins.json`. Discharge in cfs is converted to mm/d
   with the GAGES-II drainage area.
2. **Topography** (in MATLAB, with `twi_watershed_usgs.m` on the path):
   `topmodel_explorer_twi_batch` writes `twi_slim_30m/twi_slim_<gauge>.mat` and
   `topmodel_twi_30m.zip`. Then `python scripts/twi_prep.py <folder with the .mat files>`
   writes `data/twi.json`: 30 equal-width TWI classes (as `ti_classes()` in R), the 200-bin
   channel-distance histogram for routing, and the downsampled hillshade and saturation layers.
3. **Map** (only if gauges change): `cd scripts && npm install && node us_map.mjs`.
4. **Pre-calibration**: `sh scripts/run_calibration.sh` runs SCE-UA on KGE for every gauge
   (calibration = first 60% after a 365-day warm-up, 3,000 runs; the two dryland basins get
   6,000 runs, two seeds and infiltration excess). Results go to `data/presets/`.
5. **Assemble**: `python scripts/assemble.py` writes `../index.html` (for GitHub Pages) and
   `out/topmodel_explorer.html` (for a claude.ai artifact). The output is pure ASCII, so special
   characters display correctly whatever encoding a server sends.

## Check against the R model

```sh
node tests/compare_with_R.js
Rscript tests/compare_with_R.R path/to/topmodel_r/R/topmodel_core.R
```

This runs five configurations (base model; stormflow + infiltration excess + impervious +
routing; canopy on pervious area only; two land-use change scenarios). The largest
differences from R should be about 1e-13 mm/d.

## Data sources

- CAMELS (Newman et al. 2015; Addor et al. 2017), from the GitHub copies in `scripts/gauges.json`.
- USGS 3DEP 1 arc-second elevation, resampled to 30 m by `twi_watershed_usgs.m`.
- US state outlines: `us-atlas` (US Census Bureau cartographic boundaries).
