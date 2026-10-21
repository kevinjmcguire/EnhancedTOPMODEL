# Enhanced TOPMODEL: Watershed Explorer

A browser-based version of Enhanced TOPMODEL built for the classroom. Students run a daily
watershed model on eleven real gauged watersheds, switch hydrologic processes on and off,
adjust parameters, calibrate the model against observed streamflow, and test how land-use
change alters the water balance. Everything runs in the browser: there is nothing to install,
no data to download, and no programming required.

**Open the app:** [`index.html`](index.html). If GitHub Pages is enabled for this repository,
the app is at <https://kevinjmcguire.github.io/EnhancedTOPMODEL/>; otherwise download
`index.html` and open it in any modern browser (an internet connection is needed once, to load
the charting library and fonts).

## Teaching goals

The app is designed to help students:

- **Compare process representations.** See how a hydrograph changes when the model includes,
  or leaves out, saturation-excess overland flow, shallow subsurface stormflow,
  infiltration-excess runoff, impervious-surface runoff, snow, canopy interception and channel
  routing.
- **Connect topography to runoff generation.** Watch saturated areas expand and contract across
  the topographic wetness index (TWI) distribution and on a map of each watershed during storms
  and dry spells.
- **Understand parameters and calibration.** Explore parameter sensitivity, calibrate with SCE-UA,
  compare calibration and validation performance, and discuss equifinality and model structural
  error, especially in the dryland basins where TOPMODEL struggles.
- **Evaluate land-use change.** Run forest harvest, urban development, afforestation, conifer
  conversion, insect defoliation, tree mortality, wildfire and agricultural conversion scenarios
  against a baseline, and examine changes in annual water yield, evapotranspiration, peak flows
  and low flows.
- **Contrast climates.** The same model structure is applied from Coastal Plain to alpine
  snowmelt to desert watersheds, highlighting where its assumptions hold and where they break down.

## The watersheds

Ten water years (October 2004 – September 2014) of Daymet climate and USGS daily discharge from the
CAMELS dataset, with topography from 30 m USGS 3DEP elevation data.

| Watershed | USGS gauge | Setting | Climate |
|---|---|---|---|
| Johns Creek, VA | 02017500 | Valley & Ridge | Humid temperate |
| Rockfish River, VA | 02028500 | Blue Ridge | Humid subtropical |
| Holiday Creek, VA | 02038850 | Piedmont | Humid subtropical |
| Piscataway Creek, VA | 01669000 | Coastal Plain | Humid subtropical |
| Fish River, ME | 01013500 | Northern Appalachians | Cold continental, snowy |
| Queets River, WA | 12040500 | Olympic Mountains | Temperate rainforest |
| Big Sur River, CA | 11143000 | Santa Lucia Range | Mediterranean |
| Merced River, CA | 11266500 | Sierra Nevada (Yosemite) | Snowmelt-dominated Mediterranean |
| Vallecito Creek, CO | 09352900 | San Juan Mountains | Subalpine continental |
| Cave Creek, AZ | 09512280 | Sonoran Desert uplands | Hot semi-arid / desert |
| Seco Creek, TX | 08202700 | Edwards Plateau | Semi-arid subtropical, karst |

## Using the app

The pages follow a modelling workflow:

1. **Watershed** – pick a watershed from the map or the cards.
2. **Data** – climate forcing, observed discharge, monthly climate, watershed attributes and the
   TWI distribution. A precipitation multiplier lets students correct gauge undercatch
   (try it on the Queets River).
3. **Model** – turn components on and off, adjust parameters with sliders, and view the
   hydrograph (with flow components and a flow duration curve), water balance by water year,
   ET and canopy fluxes, storage and snow, the animated saturated-area map, and fit statistics
   (KGE, NSE, log NSE, bias). Each watershed starts from attribute-based values; a
   **Pre-calibrated** set is also available.
4. **Calibrate** – choose an objective function, the calibration period and the parameters, then
   run SCE-UA. The search runs in the background and can be stopped at any time.
5. **Scenarios** – choose a land-use change preset, where it happens (whole watershed, patches,
   driest uplands, or all but a riparian buffer), when it happens and how it recovers.
6. **Export & about** – copy the daily results, parameters or scenario table, and read the
   model formulation.

### Suggested exercises

- Calibrate Rockfish River with stormflow off, then on. Which flows improve, and why?
- Halve and double *m* in Piscataway Creek and Johns Creek. How do recessions and the
  saturated area respond in flat versus steep terrain?
- Compare snow accumulation and melt timing in Fish River, Merced River and Vallecito Creek.
  What does the degree-day factor control?
- Clearcut Johns Creek and follow the change in water yield and ET through the recovery. How
  does that compare with paired-watershed studies?
- Calibrate Seco Creek or Cave Creek. Why does a saturation-excess model struggle here, and
  does adding infiltration-excess runoff help?
- Develop 30% of Holiday Creek, with and without stormwater detention. What happens to peak
  flows versus annual yield?

## The model

Enhanced TOPMODEL is a daily version of TOPMODEL (Beven & Kirkby, 1979) with Penman–Monteith
evapotranspiration, canopy rain and snow interception, a degree-day snowpack, optional shallow
subsurface stormflow (Scanlon et al., 2000), infiltration-excess and impervious-area runoff,
a two-tile land-use change scheme and channel routing. The JavaScript model in this app is a
line-for-line port of the R/Shiny version and reproduces its results to about 10⁻¹³ mm/d. The
formulation is listed on the app's **Export & about** page.

## Limitations

- Watersheds are limited to the eleven examples; delineating new watersheds and downloading
  data require the full R/Shiny or MATLAB tools.
- The TWI and channel distances come from 30 m DEMs; finer DEMs give different distributions.
- The saturated-area map is downsampled for display (one map pixel covers up to 8 × 8 DEM cells).
- Daymet precipitation can be biased; the Queets River is the clearest example.
- DREAMzs uncertainty analysis from the R version is not included.

## Repository contents

| Path | What it is |
|---|---|
| `index.html` | The complete app in one file (model, data and interface) |
| `build/` | Sources, data and scripts to rebuild `index.html`; see [`build/README.md`](build/README.md) |
| `tools/` | Python scripts for the earlier version of the page; they do not work with the current `index.html` |

## Data sources and references

- CAMELS: Newman, A. J., et al. (2015), *Hydrology and Earth System Sciences*, 19, 209–223;
  Addor, N., et al. (2017), *Hydrology and Earth System Sciences*, 21, 5293–5313.
- Daymet daily surface weather (Thornton et al.), via CAMELS basin means.
- USGS daily streamflow, via CAMELS.
- USGS 3D Elevation Program (3DEP) 1 arc-second DEM, processed to MD-infinity TWI
  (Seibert & McGlynn, 2007).
- Beven, K. J., & Kirkby, M. J. (1979). A physically based, variable contributing area model of
  basin hydrology. *Hydrological Sciences Bulletin*, 24(1), 43–69.
- Scanlon, T. M., Raffensperger, J. P., Hornberger, G. M., & Clapp, R. B. (2000). Shallow
  subsurface storm flow in a forested headwater catchment. *Water Resources Research*, 36(9), 2575–2586.
- Duan, Q., Sorooshian, S., & Gupta, V. (1992). Effective and efficient global optimization for
  conceptual rainfall-runoff models. *Water Resources Research*, 28(4), 1015–1031.
