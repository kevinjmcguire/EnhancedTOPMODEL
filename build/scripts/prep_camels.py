"""Extract Daymet forcing, USGS discharge and CAMELS attributes for the 11 gauges.

    python scripts/prep_camels.py          (run from build/)

Clones three public GitHub repos (sparse, only the needed files) into
build/cache/ and writes data/basins.json:
  P (mm/d), Tx, Tn (C), Rs = srad x dayl / 1e6 (MJ m-2 d-1), ea = vp / 1000 (kPa),
  Q (mm/d, cfs converted with the GAGES-II area), plus attributes.
"""
import json, os, subprocess
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
CFG = json.load(open(os.path.join(HERE, "gauges.json")))
CACHE = os.path.join(ROOT, "cache"); os.makedirs(CACHE, exist_ok=True)

def sparse_clone(url, paths):
    d = os.path.join(CACHE, url.rstrip("/").split("/")[-1].replace(".git", ""))
    if not os.path.isdir(os.path.join(d, ".git")):
        subprocess.check_call(["git", "clone", "-q", "--depth", "1", "--filter=blob:none", "--no-checkout", url, d])
    subprocess.check_call(["git", "-C", d, "checkout", "-q", "HEAD", "--", *paths])
    return d

ids = CFG["virginia"] + CFG["other"]
src = {}
for grp in ("virginia", "other"):
    s = CFG["sources"][grp]
    paths = [s[k].format(id=i) for i in CFG[grp] for k in ("forcing", "flow")]
    d = sparse_clone(s["repo"], paths)
    for i in CFG[grp]:
        src[i] = (os.path.join(d, s["forcing"].format(id=i)), os.path.join(d, s["flow"].format(id=i)))
sa = CFG["sources"]["attributes"]
da = sparse_clone(sa["repo"], [sa["dir"]])
att = pd.concat([pd.read_csv(os.path.join(da, sa["dir"], f"camels_{n}.txt"), sep=";", dtype={"gauge_id": str}).set_index("gauge_id")
                 for n in ["name", "topo", "clim", "vege", "soil", "hydro", "geol"]], axis=1)

t0, t1 = (pd.Timestamp(x) for x in CFG["period"])
r = lambda s, d: [None if np.isnan(v) else round(float(v), d) for v in s]
out = {}
for b in ids:
    ff, fq = src[b]
    L = open(ff).read().splitlines(); latF = float(L[0])
    F = pd.read_csv(ff, sep=r"\s+", skiprows=4, header=None,
                    names=["Y", "M", "D", "Hr", "dayl", "prcp", "srad", "swe", "tmax", "tmin", "vp"])
    F.index = pd.to_datetime(dict(year=F.Y, month=F.M, day=F.D))
    Q = pd.read_csv(fq, sep=r"\s+", header=None, names=["id", "Y", "M", "D", "Q", "flag"])
    Q.index = pd.to_datetime(dict(year=Q.Y, month=Q.M, day=Q.D))
    F = F.loc[t0:t1]; Q = Q["Q"].reindex(F.index); Q[Q < 0] = np.nan
    a = att.loc[b]
    qmm = Q * 0.0283168 * 86400 / (a.area_gages2 * 1e6) * 1000
    out[b] = dict(id=b, name=a.gauge_name, lat=float(a.gauge_lat), lon=float(a.gauge_lon), latF=latF, elev=float(a.elev_mean),
                  area_km2=float(a.area_gages2), slope=float(a.slope_mean), p_mean=float(a.p_mean), pet_mean=float(a.pet_mean),
                  aridity=float(a.aridity), frac_snow=float(a.frac_snow), frac_forest=float(a.frac_forest), lai_max=float(a.lai_max),
                  lai_diff=float(a.lai_diff), land=a.dom_land_cover.strip(), geol=a.geol_1st_class, soil_depth=float(a.soil_depth_pelletier),
                  runoff_ratio=float(a.runoff_ratio), bfi=float(a.baseflow_index), start=str(F.index[0].date()),
                  P=r(F.prcp, 1), Tx=r(F.tmax, 1), Tn=r(F.tmin, 1), Rs=r(F.srad * F.dayl / 1e6, 2), ea=r(F.vp / 1000, 3), Q=r(qmm, 3))
    print(b, len(F), "days, missing Q:", int(qmm.isna().sum()))
json.dump(out, open(os.path.join(ROOT, "data", "basins.json"), "w"), separators=(",", ":"))
print("wrote data/basins.json")
