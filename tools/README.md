# Adding basins to the page

`add_sites.py` downloads everything a basin needs, calibrates a preset and writes it into `index.html`.

```
pip install -r tools/requirements.txt
python tools/add_sites.py index.html            # the six basins listed in NEW_SITES
python tools/add_sites.py index.html --sites 01654000 --iters 400   # one basin, quicker calibration
```

It needs internet access to USGS (NLDI, 3DEP, NWIS, NHD), ORNL Daymet and MRLC NLCD. Downloads are cached in
`tools/cache/<gauge>/`, so a re-run only redoes the processing and calibration. Each basin takes a few minutes,
mostly the two calibrations (2010–2013, and the 2010–11 fit scored on 2012–13). The page file is saved after each
basin. To add other gauges, add an entry to `NEW_SITES` (id, name, loc, state, climate, setting).

`--refresh-nhd` recomputes the NHD perennial stream length for every basin in the page by clipping the
high-resolution flowlines to the basin outline (the six original basins currently carry tile-based estimates).

Files:

- `topmodel.py` – Python copy of the page's model, metrics and random-search calibration
- `terrain.py` – depression filling, FD8 flow accumulation, topographic index classes and the map layer
- `sitedata.py` – reads and writes the page's embedded `DATA` and its packed daily series
- `tests/check_port.py index.html` – checks `topmodel.py` against the page's JavaScript (needs Node)
