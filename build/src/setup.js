// Basin setup shared by the app and the offline pre-calibration script.
(function (root) {
"use strict";
// key, label, group, min, max, default, isLog, digits, tip   (from param_table() in model_utils.R)
const PARAMS = [
  ["m", "m · transmissivity decay", "top", 0.001, 0.2, 0.02, true, 4, "m", "Exponential decay of transmissivity with deficit; small m = flashy recessions"],
  ["lnTe", "ln(Te) · saturated transmissivity", "top", -2, 12, 5, false, 2, "ln m²/d", "Areal-average log transmissivity when the soil is saturated to the surface"],
  ["Srmax", "Srmax · root-zone capacity", "top", 0.005, 0.3, 0.08, true, 4, "m", "Plant-available water in the root zone"],
  ["td", "td · unsaturated-zone delay", "top", 0.01, 50, 1, true, 3, "d/m", "Drainage from the unsaturated zone = SUZ / (S_local · td)"],
  ["SR0", "SR0 · initial root-zone deficit", "init", 0, 0.2, 0.02, false, 3, "m", "Root-zone deficit on day 1"],
  ["Q0", "Q0 · initial discharge", "init", 0.01, 20, 1, true, 3, "mm/d", "Sets the initial mean saturated-zone deficit"],
  ["Mssf", "Mssf · stormflow-zone capacity", "ssf", 0.002, 0.3, 0.05, true, 4, "m", "Drainable water the shallow stormflow zone can hold (porosity × thickness)"],
  ["kssf", "kssf · stormflow recession rate", "ssf", 0.01, 5, 0.5, true, 3, "1/d", "Qssf = kssf × W; linear transmissivity profile (Scanlon et al. 2000)"],
  ["kv", "kv · vertical drainage to groundwater", "ssf", 0.001, 2, 0.1, true, 3, "1/d", "Leakage kv × W from the stormflow zone recharges groundwater"],
  ["Ks", "Ks · mean infiltration capacity", "ie", 1, 500, 100, true, 3, "mm/h", "Mean of an exponential distribution; forest soils ~50–300 mm/h, compacted or fine soils ~1–20"],
  ["ieMoist", "Dry-soil factor", "ie", 0, 5, 1, false, 2, "–", "Capacity × (1 + factor × deficit/Srmax): a dry root zone absorbs more"],
  ["stormSum", "Summer storm duration", "ie", 0.5, 24, 4, false, 1, "h", "Hours over which a summer rain day falls (convective storms are short)"],
  ["stormWin", "Winter storm duration", "ie", 0.5, 24, 10, false, 1, "h", "Hours over which a winter rain day falls (frontal storms are long)"],
  ["ieConn", "Share of infiltration excess reaching stream", "ie", 0, 1, 1, false, 2, "–", "The rest re-infiltrates downslope (run-on)"],
  ["LAImax", "LAI max · leaf-on", "can", 0.5, 10, 5, false, 2, "m²/m²", "Peak (summer) leaf area index"],
  ["LAImin", "LAI min · leaf-off", "can", 0, 6, 1, false, 2, "m²/m²", "Winter plant area (stems, conifers); = LAI max for evergreen"],
  ["leafOn", "Leaf-out start", "can", 60, 200, 130, false, 0, "day", "Start of leaf-out (LAI ramps up over the transition)"],
  ["leafOff", "Leaf-fall start", "can", 200, 345, 275, false, 0, "day", "Start of leaf fall (LAI ramps down)"],
  ["ramp", "Leaf transition length", "can", 1, 60, 20, false, 0, "days", "Length of the seasonal LAI ramps"],
  ["rl", "rl · mean leaf resistance", "can", 30, 800, 250, false, 0, "s/m", "Canopy resistance rc = rl / (0.5 LAI)"],
  ["hc", "Canopy height", "can", 0.1, 40, 20, false, 1, "m", "Sets aerodynamic resistance: d = 2/3 hc, z0m = 0.123 hc"],
  ["u", "Wind speed above canopy", "can", 0.3, 8, 2, false, 2, "m/s", "Used for aerodynamic resistance (no wind in the forcing)"],
  ["cint", "Interception capacity", "can", 0, 1, 0.3, false, 2, "mm/LAI", "Canopy storage = cint × LAI; evaporates at the wet-canopy rate"],
  ["albedo", "Snow-free albedo (constant)", "alb", 0.05, 0.5, 0.15, false, 2, "–", "Shortwave reflectance without snow (forest ~0.12–0.18)"],
  ["albLeafOn", "Snow-free albedo, leaf-on", "alb", 0.05, 0.5, 0.17, false, 3, "–", "Albedo at the seasonal LAI high"],
  ["albLeafOff", "Snow-free albedo, leaf-off", "alb", 0.05, 0.5, 0.13, false, 3, "–", "Albedo at the seasonal LAI low"],
  ["Tsnow", "Rain/snow threshold", "snow", -4, 4, 0, false, 1, "°C", "Precipitation is snow when daily mean T ≤ this"],
  ["Tmelt", "Melt base temperature", "snow", -4, 4, 0, false, 1, "°C", "Melt = DDF × (Tmean − Tmelt)"],
  ["DDF", "Degree-day melt factor", "snow", 0.5, 10, 3, false, 2, "mm/°C/d", "Forest ~2–3, open ~4–6"],
  ["albSnow", "Open-snow albedo", "snow", 0.4, 0.9, 0.75, false, 2, "–", "Albedo of snow with no canopy"],
  ["kAlb", "Canopy masking of snow albedo", "snow", 0, 2, 0.75, false, 2, "1/LAI", "Snow albedo under canopy decays with exp(−k LAI)"],
  ["sweFull", "SWE for full snow cover", "snow", 1, 100, 20, false, 0, "mm", "Snow-cover fraction = min(SWE / this, 1)"],
  ["fEv", "Evergreen fraction of leaf-off LAI", "snow", 0, 1, 0, false, 2, "–", "0 = deciduous (bare branches), 1 = evergreen"],
  ["cSnow", "Foliage snow capacity", "snow", 0, 8, 3.5, false, 2, "mm/LAI", "Snow load per unit leaf area (Hedstrom & Pomeroy 1998)"],
  ["cSnowBr", "Branch snow capacity", "snow", 0, 4, 1.0, false, 2, "mm", "Snow load per unit leafless branch area"],
  ["unload", "Canopy snow unloading rate", "snow", 0, 3, 0.5, false, 2, "1/d", "Fraction 1 − exp(−rate) falls to the snowpack each day"],
  ["chv", "Channel flow velocity", "rout", 0.05, 3, 0.7, false, 2, "m/s", "Travel time = channel distance / velocity"],
  ["impPct", "Impervious area", "imp", 0, 100, 0, false, 2, "%", "Share of the watershed that is impervious"],
  ["fConn", "Directly connected impervious fraction", "imp", 0, 1, 1, false, 2, "–", "Share draining straight to the stream; the rest runs on to pervious soil"],
  ["Sdep", "Impervious depression storage", "imp", 0, 10, 1.5, false, 2, "mm", "Water held on impervious surfaces before runoff starts"],
].map(r => ({ key: r[0], label: r[1], group: r[2], min: r[3], max: r[4], def: r[5], isLog: r[6], dig: r[7], unit: r[8], tip: r[9] }));
const PBYKEY = Object.fromEntries(PARAMS.map(p => [p.key, p]));

// Basin descriptions (climate / setting) and canopy starting values per basin
const META = {
  "01669000": { short: "Piscataway Creek", state: "VA", region: "Coastal Plain", climate: "Humid subtropical (Cfa)",
    blurb: "Low-relief Coastal Plain stream on unconsolidated sands and clays, with broad wet bottomlands. Deep soils, slow groundwater-fed recessions, and large saturated areas.",
    fEv: 0.3, hc: 20, on: 110, off: 290 },
  "02038850": { short: "Holiday Creek", state: "VA", region: "Piedmont", climate: "Humid subtropical (Cfa)",
    blurb: "Small forested Piedmont catchment in the Appomattox–Buckingham State Forest on deeply weathered saprolite. A long-term USGS small-watershed research site.",
    fEv: 0.3, hc: 22, on: 105, off: 290 },
  "02028500": { short: "Rockfish River", state: "VA", region: "Blue Ridge", climate: "Humid subtropical (Cfa)",
    blurb: "Steep, deciduous-forested Blue Ridge basin on metamorphic rock (Nelson County). Thin soils give a flashy response to tropical and frontal storms.",
    fEv: 0, hc: 24, on: 110, off: 285 },
  "02017500": { short: "Johns Creek", state: "VA", region: "Valley & Ridge", climate: "Humid temperate (Cfb)",
    blurb: "Forested Valley & Ridge basin in the Jefferson National Forest near New Castle. Sandstone ridges and shale valleys, with more winter snow than the eastern Virginia basins.",
    fEv: 0, hc: 24, on: 120, off: 280 },
  "01013500": { short: "Fish River", state: "ME", region: "Northern Appalachians", climate: "Cold humid continental (Dfb)",
    blurb: "Large boreal–northern hardwood basin at the Canadian border. About a third of precipitation falls as snow; a big spring freshet follows months of frozen winter.",
    fEv: 0.4, hc: 18, on: 140, off: 270 },
  "12040500": { short: "Queets River", state: "WA", region: "Olympic Mountains", climate: "Temperate maritime rainforest (Cfb)",
    blurb: "Old-growth conifer rainforest draining the Olympic Mountains: almost 3 m of precipitation a year, mostly from winter atmospheric rivers. Daymet under-catches here, so observed Q exceeds P.",
    fEv: 1, hc: 35, on: 130, off: 275 },
  "11143000": { short: "Big Sur River", state: "CA", region: "Santa Lucia Range", climate: "Mediterranean (Csb)",
    blurb: "Steep coastal mountains with redwood canyons and chaparral. Wet winters and rainless summers, so the root zone dries out completely every year.",
    fEv: 0.6, hc: 15, on: 60, off: 300 },
  "11266500": { short: "Merced River", state: "CA", region: "Sierra Nevada (Yosemite)", climate: "Mediterranean montane, snow-dominated (Dsb)",
    blurb: "High Sierra granite above Yosemite Valley. Nearly 90% of precipitation falls as snow; the hydrograph peaks with May–June snowmelt.",
    fEv: 0.8, hc: 15, on: 150, off: 270 },
  "08202700": { short: "Seco Creek", state: "TX", region: "Edwards Plateau", climate: "Semi-arid subtropical (BSh)",
    blurb: "Karst limestone of the Balcones escarpment. High PET and thin soils: only about 3% of rain becomes streamflow, almost all of it in a few intense storms. A hard test for TOPMODEL.",
    fEv: 0.5, hc: 5, on: 80, off: 310 },
  "09512280": { short: "Cave Creek", state: "AZ", region: "Sonoran Desert uplands", climate: "Hot semi-arid / desert (BSh)",
    blurb: "Ephemeral desert stream north of Phoenix with winter frontal storms and summer monsoon thunderstorms. Aridity index near 3; the channel is dry most days.",
    fEv: 0.5, hc: 1.5, on: 75, off: 320 },
  "09352900": { short: "Vallecito Creek", state: "CO", region: "San Juan Mountains", climate: "Subalpine continental (Dfc)",
    blurb: "High-elevation wilderness basin (mean 3,280 m) of spruce–fir and alpine tundra. Snowmelt dominates, with a smaller late-summer monsoon signal.",
    fEv: 0.8, hc: 10, on: 155, off: 255 },
};
const ORDER = ["02017500", "02028500", "02038850", "01669000", "01013500", "12040500", "11143000", "11266500", "09352900", "09512280", "08202700"];

// Topographic index distribution: shifted gamma with mean lambda estimated from mean slope
function tiLambda(b) { return 4.6 - Math.log(b.slope / 1000); }
function tiDist(lam, k = 4, th = 0.75, nc = 30) {
  const t0 = lam - k * th, hi = k * th + 6 * Math.sqrt(k) * th, N = 3000, pdf = [];
  for (let i = 0; i < N; i++) { const x = (i + 0.5) * hi / N; pdf.push(Math.pow(x, k - 1) * Math.exp(-x / th)); }
  const val = [], frac = [], per = N / nc;
  for (let c = 0; c < nc; c++) {
    let s = 0, sx = 0;
    for (let i = c * per; i < (c + 1) * per; i++) { const x = (i + 0.5) * hi / N; s += pdf[i]; sx += pdf[i] * x; }
    if (s > 0) { val.push(t0 + sx / s); frac.push(s); }
  }
  const S = frac.reduce((a, b) => a + b, 0);
  return { val, frac: frac.map(f => f / S), lambda: lam };
}
function dates(b) {
  const d0 = Date.parse(b.start + "T00:00:00Z"), n = b.P.length, out = new Array(n);
  for (let t = 0; t < n; t++) out[t] = new Date(d0 + t * 864e5);
  return out;
}
function doyOf(d) { return Math.round((d - Date.UTC(d.getUTCFullYear(), 0, 1)) / 864e5) + 1; }
// forcing between indices i0..i1 (inclusive), with precipitation multiplier
function forcing(b, i0, i1, pMult = 1) {
  const D = dates(b), n = i1 - i0 + 1, F = { P: new Float64Array(n), Tmax: new Float64Array(n), Tmin: new Float64Array(n),
    Rs: new Float64Array(n), ea: new Float64Array(n), doy: new Float64Array(n), lat: b.latF, elev: b.elev };
  for (let j = 0; j < n; j++) { const t = i0 + j; F.P[j] = b.P[t] * pMult; F.Tmax[j] = b.Tx[t]; F.Tmin[j] = b.Tn[t]; F.Rs[j] = b.Rs[t]; F.ea[j] = b.ea[t]; F.doy[j] = doyOf(D[t]); }
  return F;
}
// Channel routing without a DEM: Hack's law main channel length, uniform distances (as routing_distances() in R)
function hackRouting(area_km2) {
  const ha = area_km2 * 100, Lmax = 1.4 * Math.pow(ha / 259, 0.6) * 1609.34, dist = [], frac = [];
  for (let i = 0; i < 100; i++) { dist.push(i * Lmax / 99); frac.push(0.01); }
  return { on: true, dist, frac, Lmax };
}
function defaults() {
  const p = {}; for (const r of PARAMS) p[r.key] = r.def;
  return Object.assign(p, { ssfOn: false, ieOn: false, impOn: false, albMode: 0, canopyPervOnly: false });
}
function startParams(b) {
  const m = META[b.id], p = defaults();
  p.LAImax = +Math.max(b.lai_max, 0.5).toFixed(2);
  p.LAImin = +Math.max(b.lai_max - b.lai_diff, 0.2).toFixed(2);
  if (m.fEv >= 1) p.LAImin = +(0.85 * p.LAImax).toFixed(2);
  p.fEv = m.fEv; p.hc = m.hc; p.leafOn = m.on; p.leafOff = m.off;
  const q = b.Q.find(v => v !== null && v > 0);
  p.Q0 = +Math.max(q || 1, 0.01).toFixed(3);
  return p;
}
// Topography from the 30 m 3DEP DEM (MATLAB twi_watershed_usgs.m): 30 equal-width TWI classes
// as ti_classes() in model_utils.R, and the 200-bin channel-distance histogram used for routing.
function realTI(tw) { return { val: tw.ti.val, frac: tw.ti.frac, lnLambdaLin: tw.ti.lnLambdaLin, lambda: tw.ti.lambda_ }; }
function demRouting(tw) { return { on: true, dist: tw.route.dist, frac: tw.route.frac, Lmax: tw.route.max }; }
const CAL_KEYS = ["m", "lnTe", "Srmax", "td", "DDF", "rl", "chv"];
const api = { realTI, demRouting, PARAMS, PBYKEY, META, ORDER, tiLambda, tiDist, dates, doyOf, forcing, hackRouting, defaults, startParams, CAL_KEYS };
if (typeof module !== "undefined" && module.exports) module.exports = api;
root.SETUP = api;
})(typeof self !== "undefined" ? self : this);
