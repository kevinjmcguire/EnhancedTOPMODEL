"""Add USGS-gauged basins to index.html (the Enhanced TOPMODEL page).

For each gauge this downloads the basin outline (NLDI), a 30 m 3DEP DEM, Daymet v4 daily forcing and USGS daily
discharge for 2009-2013, NLCD 2011 land cover and NHD high-resolution perennial streams; builds the topographic
index classes and the map layer; calibrates a parameter preset with the same random search as the page's
Auto-calibrate (plus a 2010-11 / 2012-13 split-sample test); and writes the basin into the page's DATA object.

Usage (from the repository folder):
    pip install -r tools/requirements.txt
    python tools/add_sites.py index.html                 # adds the six basins in NEW_SITES
    python tools/add_sites.py index.html --sites 01654000 --iters 400
    python tools/add_sites.py index.html --refresh-nhd   # also recompute NHD perennial length for every basin

Downloads are cached in tools/cache/<gauge>/ so a re-run only recalibrates. Runs about 5-10 minutes per basin,
mostly calibration (lower --iters for a quicker, rougher preset).
"""
import argparse
import json
import math
import os
import pickle
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import terrain  # noqa: E402
import topmodel as tm  # noqa: E402
from sitedata import N, load_html, save_html, encode_pack, deflate_b64  # noqa: E402

DATES = ('2009-01-01', '2013-12-31')
SPIN = 365

# Three Virginia settings that differ from the forested Blue Ridge (Piney River, already in the page),
# plus three regions the page does not cover yet.
NEW_SITES = [
    dict(id='01654000', name='Accotink Creek', loc='near Annandale, Virginia', state='VA',
         climate='Humid subtropical', setting='Urban Piedmont'),
    dict(id='01632900', name='Smith Creek', loc='near New Market, Virginia', state='VA',
         climate='Humid subtropical', setting='Agricultural karst valley'),
    dict(id='01669000', name='Piscataway Creek', loc='near Tappahannock, Virginia', state='VA',
         climate='Humid subtropical', setting='Low-relief Coastal Plain'),
    dict(id='06879650', name='Kings Creek', loc='near Manhattan, Kansas', state='KS',
         climate='Humid continental', setting='Tallgrass prairie'),
    dict(id='05454220', name='Clear Creek', loc='near Oxford, Iowa', state='IA',
         climate='Humid continental', setting='Row-crop agriculture'),
    dict(id='08155200', name='Barton Creek', loc='at SH 71 near Oak Hill, Texas', state='TX',
         climate='Subhumid subtropical', setting='Karst Hill Country'),
]


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------ downloads (HyRiver stack)
def _need(mod):
    try:
        return __import__(mod)
    except ImportError:
        raise SystemExit(f'Missing package {mod!r}. Run: pip install -r tools/requirements.txt')


def cached(site_id, name, fn):
    d = os.path.join(HERE, 'cache', site_id)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, name + '.pkl')
    if os.path.exists(path):
        with open(path, 'rb') as f:
            return pickle.load(f)
    t = time.time()
    val = fn()
    with open(path, 'wb') as f:
        pickle.dump(val, f)
    log(f'    {name}: downloaded in {time.time() - t:.0f} s')
    return val


def fetch_basin(site_id):
    pynhd = _need('pynhd')
    nldi = pynhd.NLDI()
    basin = nldi.get_basins(site_id)
    geom = basin.to_crs(4326).geometry.iloc[0]
    st = nldi.getfeature_byid('nwissite', f'USGS-{site_id}').to_crs(4326).geometry.iloc[0]
    return dict(geom_wkt=geom.wkt, gauge=(st.x, st.y))


def fetch_dem(geom_wkt, buffer_m=1000, res=30):
    py3dep = _need('py3dep')
    shapely = _need('shapely')
    from shapely import wkt
    from shapely.ops import transform as stransform
    from pyproj import Transformer
    g = wkt.loads(geom_wkt)
    to5070 = Transformer.from_crs(4326, 5070, always_xy=True).transform
    to4326 = Transformer.from_crs(5070, 4326, always_xy=True).transform
    box = shapely.geometry.box(*stransform(to5070, g).buffer(buffer_m).bounds)
    box4326 = stransform(to4326, box)
    try:
        dem = py3dep.get_map('DEM', box4326, resolution=res, geo_crs=4326, crs=5070)
    except Exception:  # older/newer py3dep
        dem = py3dep.get_dem(box4326, res, crs=5070)
    dem = dem.squeeze()
    x = dem.x.values
    y = dem.y.values
    z = dem.values.astype(float)
    if y[0] < y[-1]:
        z = z[::-1]
        y = y[::-1]
    cell = float(abs(x[1] - x[0]))
    return dict(z=z, x0=float(x[0] - cell / 2), y0=float(y[0] + cell / 2), cell=cell)


def sample_points(geom_wkt, n=16):
    """About n points spread evenly inside the basin (lon, lat), like the 16 Daymet pixels of the original basins."""
    from shapely import wkt
    from shapely.geometry import Point
    g = wkt.loads(geom_wkt)
    x0, y0, x1, y1 = g.bounds
    for k in range(4, 40):
        xs = np.linspace(x0, x1, k + 2)[1:-1]
        ys = np.linspace(y0, y1, k + 2)[1:-1]
        pts = [(x, y) for y in ys for x in xs if g.contains(Point(x, y))]
        if len(pts) >= n:
            idx = np.linspace(0, len(pts) - 1, n).round().astype(int)
            return [pts[i] for i in idx]
    c = g.representative_point()
    return [(c.x, c.y)]


DAYMET_VARS = ['prcp', 'tmax', 'tmin', 'srad', 'vp', 'dayl']


def _daymet_point(lon, lat):
    """Daymet single-pixel API (no login). Daymet's calendar has 365 days; in leap years Dec 31 is dropped."""
    import urllib.request
    url = ('https://daymet.ornl.gov/single-pixel/api/data?lat=%.5f&lon=%.5f&vars=%s&start=%s&end=%s'
           % (lat, lon, ','.join(DAYMET_VARS), DATES[0], DATES[1]))
    req = urllib.request.Request(url, headers={'User-Agent': 'EnhancedTOPMODEL-add_sites'})
    txt = urllib.request.urlopen(req, timeout=120).read().decode()
    lines = txt.splitlines()
    h = next(i for i, ln in enumerate(lines) if ln.startswith('year,yday'))
    cols = [c.split(' (')[0] for c in lines[h].split(',')]
    rows = [ln.split(',') for ln in lines[h + 1:] if ln.strip()]
    tab = {c: np.array([float(r[i]) for r in rows]) for i, c in enumerate(cols)}
    import datetime
    dates = [(datetime.date(int(y), 1, 1) + datetime.timedelta(days=int(d) - 1)).isoformat() for y, d in zip(tab['year'], tab['yday'])]
    return dates, {v: tab[v] for v in DAYMET_VARS}


def fetch_daymet_points(geom_wkt):
    pts = sample_points(geom_wkt)
    acc, dates = None, None
    for lon, lat in pts:
        dates, t = _daymet_point(lon, lat)
        acc = t if acc is None else {v: acc[v] + t[v] for v in DAYMET_VARS}
    out = {v: acc[v] / len(pts) for v in DAYMET_VARS}
    out.update(dates=dates, npix=len(pts), source=f'Daymet v4, mean of {len(pts)} pixels')
    return out


def fetch_daymet_grid(geom_wkt):
    """Gridded Daymet through pydaymet (ORNL THREDDS). Needs a NASA Earthdata login in ~/.netrc."""
    pydaymet = _need('pydaymet')
    from shapely import wkt
    g = wkt.loads(geom_wkt)
    ds = pydaymet.get_bygeom(g, DATES, variables=DAYMET_VARS, crs=4326)
    out = {}
    for v in DAYMET_VARS:
        a = ds[v]
        out[v] = a.mean(dim=[d for d in a.dims if d != 'time'], skipna=True).values.astype(float)
    out['dates'] = [str(t)[:10] for t in ds['time'].values]
    out['npix'] = int(np.isfinite(ds['prcp'].isel(time=0).values).sum())
    out['source'] = f"Daymet v4, mean of {out['npix']} pixels"
    return out


def fetch_gridmet(geom_wkt, zmean):
    """gridMET 4 km (no login), put on Daymet's calendar (Dec 31 dropped in leap years)."""
    pygridmet = _need('pygridmet')
    from shapely import wkt
    g = wkt.loads(geom_wkt)
    ds = pygridmet.get_bygeom(g, DATES, variables=['pr', 'tmmn', 'tmmx', 'srad', 'sph'], crs=4326)

    def pick(*names):
        for k in ds.data_vars:
            if any(n in k.lower() for n in names):
                a = ds[k]
                return a.mean(dim=[d for d in a.dims if d != 'time'], skipna=True).values.astype(float)
        raise KeyError(f'gridMET variable {names} not found in {list(ds.data_vars)}')
    pr = pick('pr', 'precip'); tmn = pick('tmmn', 'minimum_temp'); tmx = pick('tmmx', 'maximum_temp')
    srad = pick('srad', 'shortwave', 'radiation'); sph = pick('sph', 'specific_hum')
    if np.nanmean(tmx) > 100:  # Kelvin
        tmn, tmx = tmn - 273.15, tmx - 273.15
    pres = 101325 * (1 - 2.25577e-5 * zmean) ** 5.25588
    vp = sph * pres / (0.622 + 0.378 * sph)
    dates = [str(t)[:10] for t in ds['time'].values]
    keep = [i for i, d in enumerate(dates) if not (d.endswith('-12-31') and int(d[:4]) % 4 == 0)]
    sel = lambda a: np.asarray(a)[keep]  # noqa: E731
    return dict(prcp=sel(pr), tmax=sel(tmx), tmin=sel(tmn), srad=sel(srad), dayl=np.full(len(keep), 86400.0),
                vp=sel(vp), dates=[dates[i] for i in keep], npix=0, source='gridMET (4 km), basin mean')


def fetch_forcing(geom_wkt, zmean):
    errors = []
    for name, fn in [('Daymet point service', lambda: fetch_daymet_points(geom_wkt)),
                     ('Daymet gridded (Earthdata login)', lambda: fetch_daymet_grid(geom_wkt)),
                     ('gridMET', lambda: fetch_gridmet(geom_wkt, zmean))]:
        try:
            out = fn()
            if len(out['dates']) != N:
                raise ValueError(f"{len(out['dates'])} days returned, expected {N}")
            log(f'    forcing: {out["source"]}')
            return out
        except BaseException as e:  # noqa: BLE001 - try the next source
            if isinstance(e, KeyboardInterrupt):
                raise
            msg = str(e).splitlines()[0][:160] if str(e) else type(e).__name__
            log(f'    {name} failed: {msg}')
            errors.append(f'{name}: {msg}')
    raise SystemExit('No forcing source worked:\n  ' + '\n  '.join(errors))


def fetch_q_cfs(site_id):
    pygeohydro = _need('pygeohydro')
    nwis = pygeohydro.NWIS()
    df = nwis.get_streamflow(site_id, DATES, mmd=False)
    s = df.iloc[:, 0]
    idx = s.index
    if getattr(idx, 'tz', None) is not None:
        idx = idx.tz_convert(None) if hasattr(idx, 'tz_convert') else idx
    return {str(d)[:10]: float(v) for d, v in zip(idx, s.values)}


def fetch_nlcd(geom_wkt):
    pygeohydro = _need('pygeohydro')
    gpd = _need('geopandas')
    from shapely import wkt
    gdf = gpd.GeoDataFrame(geometry=[wkt.loads(geom_wkt)], crs=4326)
    res = pygeohydro.nlcd_bygeom(gdf, resolution=30, years={'cover': [2011], 'impervious': [2011]})
    ds = res[list(res)[0]]
    cover = ds['cover_2011'].values.astype(float)
    imp = ds['impervious_2011'].values.astype(float)
    return dict(cover=cover[np.isfinite(cover) & (cover > 0) & (cover < 255)].astype(int),
                imp=imp[np.isfinite(imp) & (imp <= 100)])


def fetch_nhd_km(geom_wkt):
    pynhd = _need('pynhd')
    gpd = _need('geopandas')
    from shapely import wkt
    g = wkt.loads(geom_wkt)
    fl = pynhd.NHD('flowline_hr').bygeom(g, geo_crs=4326)
    col = 'fcode' if 'fcode' in fl.columns else 'FCODE'
    per = fl[fl[col].astype(int) == 46006]
    if per.empty:
        return 0.0
    basin = gpd.GeoDataFrame(geometry=[g], crs=4326).to_crs(5070)
    return float(per.to_crs(5070).clip(basin).length.sum() / 1000)


# ------------------------------------------------------------------ geometry helpers
def to5070(lon, lat):
    from pyproj import Transformer
    return Transformer.from_crs(4326, 5070, always_xy=True).transform(lon, lat)


def rings_5070(geom_wkt):
    from shapely import wkt
    from shapely.ops import transform as stransform
    from pyproj import Transformer
    g = stransform(Transformer.from_crs(4326, 5070, always_xy=True).transform, wkt.loads(geom_wkt))
    polys = list(g.geoms) if g.geom_type == 'MultiPolygon' else [g]
    rings = [np.asarray(p.exterior.coords) for p in polys]
    holes = [np.asarray(h.coords) for p in polys for h in p.interiors]
    return rings, holes, g.area / 1e6, g.centroid


def rasterize(rings, holes, shape, x0, y0, cell):
    """Even-odd scanline fill of polygon rings onto the DEM grid (cell centres)."""
    H, W = shape
    mask = np.zeros(shape, bool)
    xc = x0 + (np.arange(W) + 0.5) * cell
    for ring_set, val in ((rings, True), (holes, False)):
        for r in ring_set:
            m = np.zeros(shape, bool)
            xs, ys = r[:, 0], r[:, 1]
            for row in range(H):
                yc = y0 - (row + 0.5) * cell
                y1, y2 = ys[:-1], ys[1:]
                cross = (y1 <= yc) != (y2 <= yc)
                if not cross.any():
                    continue
                xi = xs[:-1][cross] + (yc - y1[cross]) * (xs[1:][cross] - xs[:-1][cross]) / (y2[cross] - y1[cross])
                xi.sort()
                for a, b in zip(xi[0::2], xi[1::2]):
                    m[row, (xc >= a) & (xc < b)] = True
            mask = (mask | m) if val else (mask & ~m)
    return mask


def simplify(pts, tol):
    pts = [tuple(p) for p in pts]

    def dp(seg):
        if len(seg) < 3:
            return seg
        a, b = np.array(seg[0]), np.array(seg[-1])
        ab = b - a
        L = np.hypot(*ab) or 1e-12
        d = [abs(ab[0] * (a[1] - p[1]) - (a[0] - p[0]) * ab[1]) / L for p in seg[1:-1]]
        i = int(np.argmax(d)) + 1
        if d[i - 1] > tol:
            return dp(seg[:i + 1])[:-1] + dp(seg[i:])
        return [seg[0], seg[-1]]
    far = max(range(len(pts)), key=lambda i: (pts[i][0] - pts[0][0]) ** 2 + (pts[i][1] - pts[0][1]) ** 2)
    return dp(pts[:far + 1])[:-1] + dp(pts[far:])[:-1]


def albers_fit(sites):
    """Affine map from Albers equal-area (29.5/45.5) to the page's locator-map coordinates, fitted on existing dots."""
    def albers(lon, lat, lat1=29.5, lat2=45.5, lat0=37.5, lon0=-96):
        r = math.radians
        n = (math.sin(r(lat1)) + math.sin(r(lat2))) / 2
        C = math.cos(r(lat1)) ** 2 + 2 * n * math.sin(r(lat1))
        rho = lambda p: math.sqrt(C - 2 * n * math.sin(r(p))) / n  # noqa: E731
        th = n * r(lon - lon0)
        return rho(lat) * math.sin(th), rho(lat0) - rho(lat) * math.cos(th)
    A = np.array([[*albers(s['lon'], s['lat']), 1] for s in sites])
    U = np.array([s['usPt'] for s in sites])
    cx = np.linalg.lstsq(A, U[:, 0], rcond=None)[0]
    cy = np.linalg.lstsq(A, U[:, 1], rcond=None)[0]
    return lambda lon, lat: [round(float(np.dot([*albers(lon, lat), 1], cx)), 1), round(float(np.dot([*albers(lon, lat), 1], cy)), 1)]


# ------------------------------------------------------------------ land cover -> canopy parameters
LC_GROUPS = dict(water=[11, 12], developed=[21, 22, 23, 24], barren=[31], forest=[41, 42, 43], shrub=[51, 52],
                 grass=[71, 72, 73, 74], pasture=[81], crops=[82], wetland=[90, 95])
HEIGHT = dict(water=0.1, developed=8, barren=0.1, forest=20, shrub=2, grass=0.5, pasture=0.5, crops=1.0, wetland=5)


def land_cover(nlcd):
    cov = nlcd['cover']
    n = max(len(cov), 1)
    fr = {g: float(np.isin(cov, codes).sum() / n) for g, codes in LC_GROUPS.items()}
    fr['evergreen'] = float(((cov == 42).sum() + 0.5 * (cov == 43).sum()) / n)
    fr['imp'] = float(nlcd['imp'].mean()) if len(nlcd['imp']) else 0.0
    return fr


def canopy_preset(fr):
    """Starting values for the canopy and impervious parameters from NLCD fractions (not calibrated)."""
    tot = sum(fr[g] for g in HEIGHT) or 1
    hc = sum(fr[g] * HEIGHT[g] for g in HEIGHT) / tot
    p = dict(hc=round(min(max(hc, 0.3), 30), 1))
    if fr['forest'] > 0.05:
        p['fEv'] = round(min(fr['evergreen'] / fr['forest'], 1), 2)
    herb = fr['grass'] + fr['pasture'] + fr['crops']
    if herb > 0.5:  # grassland / cropland phenology: low winter LAI, later green-up for crops
        p['LAImin'] = 0.3
        p['leafOn'] = 150 if fr['crops'] > fr['grass'] + fr['pasture'] else 115
        p['leafOff'] = 265 if fr['crops'] > fr['grass'] + fr['pasture'] else 285
    if fr['imp'] >= 5:
        p['impOn'] = True
        p['impPct'] = round(fr['imp'], 1)
        p['fConn'] = round(min(1.0, 0.1 * math.sqrt(fr['imp'])), 2)  # Sutherland: EIA = 0.1 TIA^1.5
    return p


# ------------------------------------------------------------------ build one site record
def build_site(meta, basin, dem, met, q_cfs, nlcd, nhd_km, us_proj, iters=1200, seed=7):
    rings, holes, area_km2, cen = rings_5070(basin['geom_wkt'])
    z, x0, y0, cell = dem['z'], dem['x0'], dem['y0'], dem['cell']
    valid = np.isfinite(z) & (z > -1000)
    inside = rasterize(rings, holes, z.shape, x0, y0, cell) & valid
    log(f'    terrain: {inside.sum()} cells of {cell:.0f} m in the basin')
    ti_all, acc = terrain.topographic_index(z, valid, cell)
    ti = np.where(inside, ti_all, np.nan)
    val, frac, lam, lnL = terrain.ti_classes(ti)
    zin = z[inside]
    # display map, cropped to the basin
    ys, xs = np.nonzero(inside)
    r0, r1, c0, c1 = max(ys.min() - 2, 0), ys.max() + 3, max(xs.min() - 2, 0), xs.max() + 3
    q, cellM, k = terrain.display_map(ti[r0:r1, c0:c1], inside[r0:r1, c0:c1], cell)
    mx0, my0 = x0 + c0 * cell, y0 - r0 * cell
    gx, gy = to5070(*basin['gauge'])
    outline = []
    for r in rings:
        g = [((X - mx0) / cellM, (my0 - Y) / cellM) for X, Y in r]
        outline.append([[round(a, 1), round(b, 1)] for a, b in simplify(g, 0.35)])
    lon_c, lat_c = to4326(cen.x, cen.y)
    mp = dict(w=int(q.shape[1]), h=int(q.shape[0]), cellM=round(float(cellM), 1), q=deflate_b64(q.tobytes()),
              outline=outline, gauge=[round((gx - mx0) / cellM, 1), round((my0 - gy) / cellM, 1), basin['gauge'][0], basin['gauge'][1]])

    # forcing and discharge on Daymet's 365-day calendar
    dates = met['dates']
    if len(dates) != N:
        raise SystemExit(f"{meta['id']}: Daymet returned {len(dates)} days, expected {N}")
    P = met['prcp']; Tx = met['tmax']; Tn = met['tmin']
    Rs = met['srad'] * met['dayl'] / 1e6
    ea = met['vp'] / 1000
    Q = np.array([q_cfs.get(d, np.nan) for d in dates], float)
    Q = np.where(Q >= 0, Q * 0.0283168 * 86400 / (area_km2 * 1e6) * 1000, np.nan)
    log(f'    discharge: {np.isfinite(Q).sum()} of {N} days, mean {np.nanmean(Q[SPIN:]) * 365:.0f} mm/yr; '
        f'precip {P[SPIN:].mean() * 365:.0f} mm/yr ({met.get("source", "Daymet")})')
    pack = encode_pack(P, Tx, Tn, Rs, ea, Q)

    # decode exactly as the page does (rounded values) so the preset is fitted to what the page will run
    from sitedata import decode_pack
    site = dict(meta)
    site.update(A=round(area_km2, 1), zmean=int(round(zin.mean())), zmin=int(round(zin.min())), zmax=int(round(zin.max())),
                lat=round(lat_c, 4), lon=round(lon_c, 4), lam=lam, lnL=lnL, val=val, frac=frac, pack=pack, map=mp,
                usPt=us_proj(lon_c, lat_c))
    F, TI, Qd = decode_pack(site)

    fr = land_cover(nlcd)
    base = dict(tm.defaults(), **canopy_preset(fr))
    pct = lambda x: round(100 * x)  # noqa: E731
    ag = fr['pasture'] + fr['crops']
    site['lc'] = dict(forest=pct(fr['forest']), ag=pct(ag), grass=pct(fr['grass'] + fr['shrub']), developed=pct(fr['developed']),
                      water=pct(fr['water'] + fr['wetland']), imp=round(fr['imp'], 1))
    site['note'] = (f"{meta['setting'].lower()}; {pct(fr['forest'])}% forest, {pct(ag)}% agriculture, "
                    f"{pct(fr['developed'])}% developed ({fr['imp']:.0f}% impervious), NLCD 2011"
                    + ('' if met.get('source', 'Daymet').startswith('Daymet') else f"; forcing from {met['source']}"))
    if nhd_km is not None:
        site['nhdKm'] = round(nhd_km, 1)

    log(f'    calibrating 2010-2013 ({iters} runs)...')
    full, _ = tm.calibrate(F, TI, Qd, base, (SPIN, N), iters, seed, log)
    log(f'    split-sample: fitting 2010-11, scoring 2012-13 ({iters} runs)...')
    split, _ = tm.calibrate(F, TI, Qd, base, (SPIN, 3 * 365), iters, seed + 1, log)

    def scores(p, a, b):
        M = tm.perf_metrics(Qd[a:b], tm.run(F, TI, p)['Qsim'][a:b])
        return dict(NSE=round(M['NSE'], 2), KGE=round(M['KGE'], 2), PBIAS=round(M['PBIAS'], 1))
    site['fullCal'] = scores(full, SPIN, N)
    site['split'] = dict(cal=scores(split, SPIN, 3 * 365), val=scores(split, 3 * 365, N))
    calib = {k: float(f'{full[k]:.4g}') for k in tm.CAL_KEYS}
    calib.update(canopy_preset(fr))
    site['calib'] = calib
    log(f"    preset: NSE {site['fullCal']['NSE']}, KGE {site['fullCal']['KGE']}, PBIAS {site['fullCal']['PBIAS']}%; "
        f"split-sample 2012-13 NSE {site['split']['val']['NSE']}")
    return site


def to4326(x, y):
    from pyproj import Transformer
    return Transformer.from_crs(5070, 4326, always_xy=True).transform(x, y)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('html', help='path to index.html (updated in place)')
    ap.add_argument('--sites', help='comma-separated gauge IDs from NEW_SITES (default: all)')
    ap.add_argument('--iters', type=int, default=1200, help='calibration runs per fit (default 1200)')
    ap.add_argument('--refresh-nhd', action='store_true', help='recompute NHD perennial length for all basins in the page')
    a = ap.parse_args()
    html, data = load_html(a.html)
    want = set(a.sites.split(',')) if a.sites else None
    us_proj = albers_fit([s for s in data['sites'] if 'usPt' in s][:6])
    todo = [m for m in NEW_SITES if want is None or m['id'] in want]
    for meta in todo:
        sid = meta['id']
        log(f"{sid} {meta['name']} {meta['loc']}")
        basin = cached(sid, 'basin', lambda: fetch_basin(sid))
        dem = cached(sid, 'dem', lambda: fetch_dem(basin['geom_wkt']))
        zmean = float(np.nanmean(np.where(dem['z'] > -1000, dem['z'], np.nan)))
        met = cached(sid, 'forcing', lambda: fetch_forcing(basin['geom_wkt'], zmean))
        q = cached(sid, 'nwis', lambda: fetch_q_cfs(sid))
        nlcd = cached(sid, 'nlcd', lambda: fetch_nlcd(basin['geom_wkt']))
        try:
            nhd = cached(sid, 'nhd', lambda: fetch_nhd_km(basin['geom_wkt']))
        except Exception as e:  # the page falls back to a fixed channel threshold
            log(f'    NHD perennial length not available ({e}); the map will use a fixed channel threshold')
            nhd = None
        site = build_site(meta, basin, dem, met, q, nlcd, nhd, us_proj, a.iters)
        data['sites'] = [s for s in data['sites'] if s['id'] != sid] + [site]
        save_html(a.html, html, data)  # save after each basin so a long run can be interrupted
        html, data = load_html(a.html)
        log(f'    added to {a.html}')
    if a.refresh_nhd:
        for s in data['sites']:
            try:
                b = cached(s['id'], 'basin', lambda: fetch_basin(s['id']))
                s['nhdKm'] = round(cached(s['id'], 'nhd', lambda: fetch_nhd_km(b['geom_wkt'])), 1)
                log(f"{s['id']} NHD perennial {s['nhdKm']} km")
            except Exception as e:
                log(f"{s['id']} NHD failed: {e}")
        save_html(a.html, html, data)
    log('done')


if __name__ == '__main__':
    main()
