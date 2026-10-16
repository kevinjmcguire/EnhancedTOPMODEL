"""Python port of the page's Enhanced TOPMODEL (index.html, `TM.run` and `TM.perfMetrics`).

Kept line-for-line with the JavaScript so parameter presets calibrated here behave the same in the page.
tests/check_port.py compares the two on the bundled basins.
"""
import math
import numpy as np

LV_LS = 2.501 / 2.834


def _clamp(x, a, b):
    return min(max(x, a), b)


def _esf(T):
    return 0.6108 * np.exp(17.27 * T / (T + 237.3))


def pm_pet(Ta, Tx, Tn, Rs, ea, u, LAI, doy, lat, elev, p, alb):
    uu = max(u, 0.1)
    Pa = 101.3 * ((293 - 0.0065 * elev) / 293) ** 5.26
    gam = 0.000665 * Pa
    cp = 1.013e-3
    phi = lat * math.pi / 180
    hc = max(p['hc'], 0.1)
    d = 2 / 3 * hc
    z0m = 0.123 * hc
    z0h = 0.1 * z0m
    zr = hc + 2
    ra = math.log((zr - d) / z0m) * math.log((zr - d) / z0h) / (0.41 * 0.41 * uu)
    es = (_esf(Tx) + _esf(Tn)) / 2
    vpd = np.maximum(es - ea, 0)
    Dl = 4098 * _esf(Ta) / (Ta + 237.3) ** 2
    rho = 3.486 * Pa / (1.01 * (Ta + 273.15))
    lv = 2.501 - 0.002361 * Ta
    dr = 1 + 0.033 * np.cos(2 * math.pi * doy / 365)
    dec = 0.409 * np.sin(2 * math.pi * doy / 365 - 1.39)
    ws = np.arccos(np.clip(-math.tan(phi) * np.tan(dec), -1, 1))
    Ra = 24 * 60 / math.pi * 0.0820 * dr * (ws * math.sin(phi) * np.sin(dec) + math.cos(phi) * np.cos(dec) * np.sin(ws))
    Rso = (0.75 + 2e-5 * elev) * Ra
    rr = np.minimum(Rs / np.maximum(Rso, 1e-6), 1)
    Rnl = 4.903e-9 * ((Tx + 273.16) ** 4 + (Tn + 273.16) ** 4) / 2 * (0.34 - 0.14 * np.sqrt(np.maximum(ea, 0))) * (1.35 * rr - 0.35)
    rn = (1 - alb) * Rs - Rnl
    rc = p['rl'] / (0.5 * LAI)
    aer = 86400 * rho * cp * vpd / ra
    PET = np.maximum((Dl * rn + aer) / (Dl + gam * (1 + rc / ra)) / lv, 0)
    Ewet = np.maximum((Dl * rn + aer) / (Dl + gam) / lv, 0)
    return PET, Ewet, rn


def _lai_phenology(LAI):
    v = np.sort(LAI)
    n = len(v)
    lo = v[max(1, round(0.05 * n)) - 1]
    hi = v[max(1, round(0.95 * n)) - 1]
    if hi - lo < 1e-3:
        return np.ones_like(LAI)
    return np.clip((LAI - lo) / (hi - lo), 0, 1)


def _snowfree_albedo(p, LAI):
    if p.get('albMode', 0) >= 1:
        ph = _lai_phenology(LAI)
        a = p['albLeafOff'] + (p['albLeafOn'] - p['albLeafOff']) * ph
    else:
        a = np.full(len(LAI), p['albedo'])
    return np.clip(a, 0.02, 0.95)


def _ie_runoff(P, D, K):
    i = P / D
    return P * i / (i + max(K, 1e-6))


def run(F, TI, p):
    """F: dict P,Tmax,Tmin,Rs,ea,doy (arrays), lat, elev, u(optional); TI: dict val,frac,lnLambdaLin."""
    P = np.asarray(F['P'], float); Tx = np.asarray(F['Tmax'], float); Tn = np.asarray(F['Tmin'], float)
    n = len(P)
    Ta = (Tx + Tn) / 2
    doy = np.asarray(F['doy'], float)
    fOn = np.clip((doy - p['leafOn']) / max(p['ramp'], 1), 0, 1)
    fOff = np.clip((doy - p['leafOff']) / max(p['ramp'], 1), 0, 1)
    LAI = np.maximum(p['LAImin'] + (p['LAImax'] - p['LAImin']) * fOn * (1 - fOff), 0.01)
    alb0 = _snowfree_albedo(p, LAI)
    albS = alb0 + (p['albSnow'] - alb0) * np.exp(-p['kAlb'] * LAI)
    uu = F['u'] if F.get('u') is not None else p['u']
    PET0, Ewet0, Rn0 = pm_pet(Ta, Tx, Tn, F['Rs'], F['ea'], uu, LAI, doy, F['lat'], F['elev'], p, alb0)
    PETS, EwetS, RnS = pm_pet(Ta, Tx, Tn, F['Rs'], F['ea'], uu, LAI, doy, F['lat'], F['elev'], p, albS)

    useImp = bool(p.get('impOn'))
    fI = _clamp(p['impPct'] / 100, 0, 0.99) if useImp else 0.0
    fP = 1 - fI
    fC = _clamp(p['fConn'], 0, 1) if useImp else 1.0
    canPerv = useImp and fI > 0 and bool(p.get('canopyPervOnly'))
    gC = fP if canPerv else 1.0
    Dimp = 0.0

    frac = np.asarray(TI['frac'], float)
    ti = np.asarray(TI['val'], float)
    af = frac / frac.sum()
    nc = len(ti)
    lam = float((af * ti).sum())
    szq = math.exp(p['lnTe'] - lam)
    q0 = max(p['Q0'], 1e-6) / 1000
    Sbar = max(-p['m'] * math.log(q0 / szq), 0)
    SRZ = np.full(nc, float(p['SR0']))
    SUZ = np.zeros(nc)
    SWE = C = Cs = SWEi = 0.0

    useSSF = bool(p.get('ssfOn'))
    lnLam = Mssf = fOut = shareQ = 0.0
    if useSSF:
        lnLam = TI['lnLambdaLin']
        Mssf = max(p['Mssf'], 1e-6)
        fOut = 1 - math.exp(-(p['kssf'] + p['kv']))
        shareQ = p['kssf'] / max(p['kssf'] + p['kv'], 2.220446e-16)
    W = 0.0
    useIE = bool(p.get('ieOn'))
    Dh = np.maximum(p['stormWin'] + (p['stormSum'] - p['stormWin']) * 0.5 * (1 - np.cos(2 * math.pi * (doy - 15) / 365)), 0.25)

    keys = ['rain', 'snow', 'melt', 'SWE', 'Ei', 'Ea', 'Es', 'Eimp', 'Qof', 'Qb', 'Quz', 'Qssf', 'Qimp', 'Qie', 'satFrac', 'storage']
    o = {k: np.zeros(n) for k in keys}
    Srmax, mm, td, DDF, Tmelt = p['Srmax'], p['m'], p['td'], p['DDF'], p['Tmelt']

    for t in range(n):
        PETt, Ewt = PET0[t], Ewet0[t]
        if Ta[t] <= p['Tsnow']:
            sn, rn = P[t], 0.0
        else:
            sn, rn = 0.0, P[t]
        L = LAI[t]
        woody = (1 - p['fEv']) * min(p['LAImin'], L)
        CsMax = p['cSnow'] * (L - woody) + p['cSnowBr'] * woody
        Is = 0.0
        if sn > 0 and CsMax > Cs:
            Cc = 1 - math.exp(-0.5 * L)
            Is = min(sn, (CsMax - Cs) * (1 - math.exp(-Cc * sn / CsMax)))
        Cs += Is
        unl = Cs * (1 - math.exp(-p['unload'])); Cs -= unl
        dmelt = max(DDF * (Ta[t] - Tmelt), 0)
        drip = min(Cs, dmelt); Cs -= drip
        if Cs > CsMax:
            unl += Cs - CsMax; Cs = CsMax
        SWE += (sn - Is) + unl
        fsc = min(SWE / max(p['sweFull'], 1e-6), 1)
        if fsc > 0:
            PETt = PET0[t] + fsc * (PETS[t] - PET0[t])
            Ewt = Ewet0[t] + fsc * (EwetS[t] - Ewet0[t])
        mlt = min(SWE, dmelt); SWE -= mlt
        mltI = 0.0
        if canPerv:
            SWEi += sn; mltI = min(SWEi, dmelt); SWEi -= mltI
        Cmax = p['cint'] * L
        I = min(rn, max(Cmax - C, 0))
        C += I; thr = rn - I
        Ei = min(C, Ewt); C -= Ei
        Es = min(Cs, (Ewt - Ei) * LV_LS); Cs -= Es
        fw = min((Ei + Es / LV_LS) / Ewt, 1) if Ewt > 0 else 0
        Edry = PETt * (1 - fw)

        win = thr + mlt + drip
        qimp = eimp = ron = xs = 0.0
        if useImp and fI > 0:
            if canPerv:
                Dimp += rn + mltI; eDem = Ewt
            else:
                Dimp += win; eDem = max(Ewt - Ei - Es / LV_LS, 0)
            xs = max(Dimp - p['Sdep'], 0); Dimp -= xs
            eimp = min(Dimp, eDem); Dimp -= eimp
            qimp = fC * xs
            ron = (1 - fC) * xs * fI / fP

        SDraw = Sbar + mm * (lam - ti)
        ieK = np.zeros(nc)
        if useIE and thr > 0:
            K = p['Ks'] * (1 + p['ieMoist'] * np.minimum(SRZ / Srmax, 1))
            i_ = thr / Dh[t]
            ieK = np.where(SDraw > 0, p['ieConn'] * thr * i_ / (i_ + np.maximum(K, 1e-6)), 0.0)
        ep = Edry / 1000
        SD = np.maximum(SDraw, 0)
        if useSSF and W > 0:
            satSSF = ti >= lnLam + math.log(Mssf / W)
        else:
            satSSF = np.zeros(nc, bool)
        pin = (win + ron - ieK) / 1000
        SRZ -= pin
        neg = SRZ < 0
        SUZ[neg] -= SRZ[neg]; SRZ[neg] = 0
        SDe = np.where(satSSF, 0, SD)
        ex = np.maximum(SUZ - SDe, 0); SUZ -= ex
        with np.errstate(divide='ignore', invalid='ignore'):
            uz = np.where(SDe > 0, np.minimum(SUZ / (SD * td), SUZ), 0.0)
        SUZ -= uz
        ea_ = ep * (1 - SRZ / Srmax)
        ea_ = np.minimum(np.maximum(ea_, 0), Srmax - SRZ)
        SRZ += ea_
        sumEx = float((af * ex).sum()); quz = float((af * uz).sum()); sumEa = float((af * ea_).sum()); sumIe = float((af * ieK).sum())
        qb = szq * math.exp(-Sbar / mm)
        qssf = ovf = 0.0
        if useSSF:
            W += quz
            if W > Mssf:
                ovf = W - Mssf; W = Mssf
            dW = W * fOut
            qssf = dW * shareQ; qv = dW - qssf; W -= dW
            Sbar = Sbar - qv + qb
        else:
            Sbar = Sbar - quz + qb
        sat = float(af[(SD <= 0) | satSSF].sum())
        o['rain'][t] = rn; o['snow'][t] = sn
        o['melt'][t] = gC * mlt + (1 - gC) * mltI; o['SWE'][t] = gC * SWE + (1 - gC) * SWEi
        o['Ei'][t] = gC * Ei; o['Es'][t] = gC * Es; o['Ea'][t] = fP * 1000 * sumEa; o['Eimp'][t] = fI * eimp
        o['Qof'][t] = fP * 1000 * (sumEx + ovf); o['Qb'][t] = fP * 1000 * qb; o['Quz'][t] = fP * 1000 * quz
        o['Qssf'][t] = fP * 1000 * qssf; o['Qimp'][t] = fI * qimp; o['Qie'][t] = fP * sumIe
        o['satFrac'][t] = fP * sat
    o['Qsim'] = o['Qof'] + o['Qb'] + o['Qssf'] + o['Qimp'] + o['Qie']
    o['ET'] = o['Ei'] + o['Es'] + o['Ea'] + o['Eimp']
    o['LAI'] = LAI
    return o


def perf_metrics(obs, sim):
    obs = np.asarray(obs, float); sim = np.asarray(sim, float)
    ok = np.isfinite(obs) & np.isfinite(sim) & (obs >= 0)
    oo, ss = obs[ok], sim[ok]
    n = len(oo)
    mo, ms = oo.mean(), ss.mean()

    def nse(a, b):
        return 1 - ((b - a) ** 2).sum() / ((a - a.mean()) ** 2).sum()
    e = 0.01 * mo
    M = {'n': n, 'mean_obs': mo, 'mean_sim': ms}
    M['NSE'] = nse(oo, ss)
    M['NSE_log'] = nse(np.log(np.maximum(oo, 0) + e), np.log(np.maximum(ss, 0) + e))
    M['NSE_sqrt'] = nse(np.sqrt(np.maximum(oo, 0)), np.sqrt(np.maximum(ss, 0)))
    so, sS = oo.std(ddof=1), ss.std(ddof=1)
    cov = ((oo - mo) * (ss - ms)).sum() / (n - 1)
    r = cov / (so * sS) if so > 0 and sS > 0 else float('nan')
    M['KGE_r'] = r; M['KGE_alpha'] = sS / so; M['KGE_beta'] = ms / mo
    M['KGE'] = 1 - math.sqrt((r - 1) ** 2 + (sS / so - 1) ** 2 + (ms / mo - 1) ** 2)
    M['PBIAS'] = 100 * (ms - mo) / mo
    return M


# ---- parameter table (same as the page's PT) ----
PT = [
    ('m', 0.001, 0.2, 0.02, 1), ('lnTe', -2, 12, 5, 0), ('Srmax', 0.005, 0.3, 0.08, 1), ('td', 0.01, 50, 1, 1),
    ('SR0', 0, 0.2, 0.02, 0), ('Q0', 0.01, 20, 1, 1), ('Mssf', 0.002, 0.3, 0.05, 1), ('kssf', 0.01, 5, 0.5, 1),
    ('kv', 0.001, 2, 0.1, 1), ('Ks', 1, 500, 100, 1), ('ieMoist', 0, 5, 1, 0), ('stormSum', 0.5, 24, 4, 0),
    ('stormWin', 0.5, 24, 10, 0), ('ieConn', 0, 1, 1, 0), ('LAImax', 0.5, 10, 5, 0), ('LAImin', 0, 6, 1, 0),
    ('leafOn', 60, 200, 130, 0), ('leafOff', 200, 345, 275, 0), ('ramp', 1, 60, 20, 0), ('rl', 30, 800, 250, 0),
    ('hc', 0.1, 40, 20, 0), ('u', 0.3, 8, 2, 0), ('cint', 0, 1, 0.3, 0), ('albedo', 0.05, 0.5, 0.15, 0),
    ('Tsnow', -4, 4, 0, 0), ('Tmelt', -4, 4, 0, 0), ('DDF', 0.5, 10, 3, 0), ('albSnow', 0.4, 0.9, 0.75, 0),
    ('kAlb', 0, 2, 0.75, 0), ('sweFull', 1, 100, 20, 0), ('fEv', 0, 1, 0, 0), ('cSnow', 0, 8, 3.5, 0),
    ('cSnowBr', 0, 4, 1, 0), ('unload', 0, 3, 0.5, 0), ('impPct', 0, 100, 0, 0), ('fConn', 0, 1, 1, 0),
    ('Sdep', 0, 10, 1.5, 0)]
PM = {k: dict(min=a, max=b, default=c, log=bool(lg)) for k, a, b, c, lg in PT}
CAL_KEYS = ['m', 'lnTe', 'Srmax', 'td', 'Mssf', 'kssf', 'kv', 'DDF', 'rl', 'LAImax', 'Ks']


def defaults():
    p = {k: v['default'] for k, v in PM.items()}
    p.update(ssfOn=False, ieOn=False, impOn=False, albMode=0, canopyPervOnly=False)
    return p


def to_unit(k, v):
    m = PM[k]
    if m['log']:
        return (math.log(v) - math.log(m['min'])) / (math.log(m['max']) - math.log(m['min']))
    return (v - m['min']) / (m['max'] - m['min'])


def from_unit(k, x):
    m = PM[k]
    if m['log']:
        return math.exp(math.log(m['min']) + x * (math.log(m['max']) - math.log(m['min'])))
    return m['min'] + x * (m['max'] - m['min'])


def objective(M):
    v = 0.5 * M['KGE'] + 0.5 * M['NSE_sqrt'] - 0.25 * abs(M['PBIAS']) / 100
    return v if math.isfinite(v) else -9


def calibrate(F, TI, Q, base, window, n_iter=1200, seed=7, log=None):
    """DDS search over CAL_KEYS in unit space (same scheme as the page's Auto-calibrate)."""
    rng = np.random.default_rng(seed)
    a, b = window
    base = dict(base, ssfOn=True, ieOn=True)

    def params(u):
        p = dict(base)
        for k, x in zip(CAL_KEYS, u):
            p[k] = from_unit(k, min(max(x, 0), 1))
        return p

    def f(u):
        o = run(F, TI, params(u))
        return objective(perf_metrics(Q[a:b], o['Qsim'][a:b]))
    best = np.clip([to_unit(k, base[k]) for k in CAL_KEYS], 0, 1)
    bf = f(best)
    for it in range(1, n_iter + 1):
        pr = 1 - math.log(it) / math.log(n_iter)
        cand = best.copy()
        sel = rng.random(len(CAL_KEYS)) < pr
        if not sel.any():
            sel[rng.integers(len(CAL_KEYS))] = True
        for j in np.nonzero(sel)[0]:
            x = cand[j] + 0.2 * rng.standard_normal()
            if x < 0: x = -x
            if x > 1: x = 2 - x
            cand[j] = min(max(x, 0), 1)
        v = f(cand)
        if v > bf:
            bf, best = v, cand
        if log and it % 200 == 0:
            log(f'    {it}/{n_iter} objective {bf:.3f}')
    return params(best), bf
