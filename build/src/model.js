// =============================================================================
// Enhanced TOPMODEL core - JavaScript port of R/topmodel_core.R (non-compiled
// loop + land-use-change two-tile loop), plus routing, metrics and SCE-UA.
// Works in the browser, in a Web Worker and in Node (module.exports).
// =============================================================================
(function (root) {
"use strict";
const LV_LS = 2.501 / 2.834;
const clamp = (x, a, b) => Math.min(Math.max(x, a), b);
const dflt = (v, d) => (v === undefined || v === null || (typeof v === "number" && Number.isNaN(v))) ? d : v;

// ---------------------------------------------------------------- Penman-Monteith
function esf(T) { return 0.6108 * Math.exp(17.27 * T / (T + 237.3)); }
function pmPet(Ta, Tx, Tn, Rs, ea, uIn, LAI, doy, lat, elev, p, alb, hcArr) {
  const n = Ta.length;
  const PET = new Float64Array(n), Ewet = new Float64Array(n), ra = new Float64Array(n),
        rc = new Float64Array(n), Rn = new Float64Array(n);
  const Pa = 101.3 * Math.pow((293 - 0.0065 * elev) / 293, 5.26);
  const gam = 0.000665 * Pa, cp = 1.013e-3, phi = lat * Math.PI / 180;
  for (let t = 0; t < n; t++) {
    const u = Math.max(uIn ? uIn[t] : p.u, 0.1);
    const es = (esf(Tx[t]) + esf(Tn[t])) / 2;
    const vpd = Math.max(es - ea[t], 0);
    const Dl = 4098 * esf(Ta[t]) / Math.pow(Ta[t] + 237.3, 2);
    const rho = 3.486 * Pa / (1.01 * (Ta[t] + 273.15));
    const lv = 2.501 - 0.002361 * Ta[t];
    const dr = 1 + 0.033 * Math.cos(2 * Math.PI * doy[t] / 365);
    const dec = 0.409 * Math.sin(2 * Math.PI * doy[t] / 365 - 1.39);
    const ws = Math.acos(clamp(-Math.tan(phi) * Math.tan(dec), -1, 1));
    const Ra = 24 * 60 / Math.PI * 0.0820 * dr * (ws * Math.sin(phi) * Math.sin(dec) + Math.cos(phi) * Math.cos(dec) * Math.sin(ws));
    const Rso = (0.75 + 2e-5 * elev) * Ra;
    const Rns = (1 - alb[t]) * Rs[t];
    const rr = Math.min(Rs[t] / Math.max(Rso, 1e-6), 1);
    const Rnl = 4.903e-9 * (Math.pow(Tx[t] + 273.16, 4) + Math.pow(Tn[t] + 273.16, 4)) / 2 *
                (0.34 - 0.14 * Math.sqrt(Math.max(ea[t], 0))) * (1.35 * rr - 0.35);
    Rn[t] = Rns - Rnl;
    const hc = Math.max(hcArr ? hcArr[t] : p.hc, 0.1);
    const d = 2 / 3 * hc, z0m = 0.123 * hc, z0h = 0.1 * z0m, zr = hc + 2;
    ra[t] = Math.log((zr - d) / z0m) * Math.log((zr - d) / z0h) / (0.41 * 0.41 * u);
    rc[t] = p.rl / (0.5 * LAI[t]);
    const aer = 86400 * rho * cp * vpd / ra[t];
    PET[t] = Math.max((Dl * Rn[t] + aer) / (Dl + gam * (1 + rc[t] / ra[t])) / lv, 0);
    Ewet[t] = Math.max((Dl * Rn[t] + aer) / (Dl + gam) / lv, 0);
  }
  return { PET, Ewet, ra, rc, Rn };
}

function laiPhenology(LAI) {
  const v = Array.from(LAI).filter(x => !Number.isNaN(x)).sort((a, b) => a - b);
  const n = LAI.length, out = new Float64Array(n).fill(1);
  if (!v.length) return out;
  const lo = v[Math.max(1, Math.round(0.05 * v.length)) - 1], hi = v[Math.max(1, Math.round(0.95 * v.length)) - 1];
  if (hi - lo < 1e-3) return out;
  for (let t = 0; t < n; t++) out[t] = clamp((LAI[t] - lo) / (hi - lo), 0, 1);
  return out;
}
function snowfreeAlbedo(p, LAI, n) {
  const a = new Float64Array(n).fill(p.albedo);
  if (p.albMode >= 1) {
    const on = dflt(p.albLeafOn, p.albedo), off = dflt(p.albLeafOff, p.albedo), ph = laiPhenology(LAI);
    for (let t = 0; t < n; t++) a[t] = off + (on - off) * ph[t];
  }
  for (let t = 0; t < n; t++) a[t] = clamp(a[t], 0.02, 0.95);
  return a;
}

// infiltration excess: P i/(i+K), i = P/D; K per class; inactive (saturated) classes 0
function ieRunoff(P, D, K, active, out) {
  const i = P / D;
  for (let c = 0; c < out.length; c++) out[c] = active[c] ? P * i / (i + Math.max(K[c], 1e-6)) : 0;
  return out;
}
function ieRunoffScalar(P, D, K) { const i = P / D; return P * i / (i + Math.max(K, 1e-6)); }

// [laiMode laiFactor LAImax LAImin leafOn leafOff ramp hc albedo DDF srFrac fEv compFrac compK]
function lucState(st, LAIb, alb0, p, doy) {
  st = st.slice(); while (st.length < 14) st.push(NaN);
  const n = LAIb.length, isna = v => v === null || v === undefined || Number.isNaN(v);
  const pick = (a, d) => isna(a) ? d : a;
  const X = {};
  X.LAI = new Float64Array(n);
  if (isna(st[0]) || st[0] === 0) {
    const f = pick(st[1], 1);
    for (let t = 0; t < n; t++) X.LAI[t] = f * LAIb[t];
  } else {
    const dq = [p.LAImax, p.LAImin, p.leafOn, p.leafOff, p.ramp];
    const q = [2, 3, 4, 5, 6].map((j, i) => isna(st[j]) ? dq[i] : st[j]);
    for (let t = 0; t < n; t++) {
      const fOn = clamp((doy[t] - q[2]) / Math.max(q[4], 1), 0, 1);
      const fOff = clamp((doy[t] - q[3]) / Math.max(q[4], 1), 0, 1);
      X.LAI[t] = q[1] + (q[0] - q[1]) * fOn * (1 - fOff);
    }
  }
  const rep = v => new Float64Array(n).fill(v);
  X.hc = rep(pick(st[7], p.hc));
  X.alb = isna(st[8]) ? Float64Array.from(alb0) : rep(st[8]);
  X.ddf = rep(pick(st[9], p.DDF));
  X.sr = rep(pick(st[10], 1));
  X.fEv = rep(pick(st[11], p.fEv));
  X.comp = rep(pick(st[12], 0));
  X.compK = rep(pick(st[13], 25));
  return X;
}

function canopyStep(C, Cs, SWE, sn, rn, LAIt, woody, p, DDF, Ta, PET0, PETs, Ewet0, EwetS, Rn0, RnS, alb0, albS) {
  const CsMax = p.cSnow * (LAIt - woody) + p.cSnowBr * woody;
  let Is = 0;
  if (sn > 0 && CsMax > Cs) {
    const Cc = 1 - Math.exp(-0.5 * LAIt);
    Is = Math.min(sn, (CsMax - Cs) * (1 - Math.exp(-Cc * sn / CsMax)));
  }
  Cs += Is;
  let unl = Cs * (1 - Math.exp(-p.unload)); Cs -= unl;
  const dm = Math.max(DDF * (Ta - p.Tmelt), 0);
  const drip = Math.min(Cs, dm); Cs -= drip;
  if (Cs > CsMax) { unl += Cs - CsMax; Cs = CsMax; }
  SWE += (sn - Is) + unl;
  const fsc = Math.min(SWE / Math.max(p.sweFull, 1e-6), 1);
  const PET = PET0 + fsc * (PETs - PET0), Ewet = Ewet0 + fsc * (EwetS - Ewet0);
  const Rn = Rn0 + fsc * (RnS - Rn0), alb = alb0 + fsc * (albS - alb0);
  const mlt = Math.min(SWE, dm); SWE -= mlt;
  const Cmax = p.cint * LAIt;
  const I = Math.min(rn, Math.max(Cmax - C, 0));
  C += I; const thr = rn - I;
  const Ei = Math.min(C, Ewet); C -= Ei;
  const Es = Math.min(Cs, (Ewet - Ei) * LV_LS); Cs -= Es;
  const fw = Ewet > 0 ? Math.min((Ei + Es / LV_LS) / Ewet, 1) : 0;
  return { C, Cs, SWE, Is, unl, drip, mlt, thr, Ei, Es, Edry: PET * (1 - fw), PET, Ewet, Rn, alb, fsc };
}

const OUT_KEYS = ["rain", "snow", "melt", "SWE", "Ei", "Ea", "Qssf", "Wssf", "qv", "satFracSSF", "Isnow", "Es",
  "Csnow", "unloaded", "drip", "Ccan", "Qof", "Qb", "Quz", "Qsim", "Sbar", "satFrac", "SRZ", "SUZ", "storage",
  "Qimp", "Eimp", "Dimp", "runon", "Qinf", "Qcomp", "LAIharv", "Qdev", "Sdet", "impChanged", "Qie"];

// F: {P,Tmax,Tmin,Rs,ea,doy,u?,lat,elev}; TI: {val,frac,lnLambdaLin?}; p: params; R: {on,dist,frac}
function core(F, TI, p, R) {
  p = Object.assign({}, p);
  const n = F.P.length, P = F.P, Tx = F.Tmax, Tn = F.Tmin, doy = F.doy;
  const Ta = new Float64Array(n);
  for (let t = 0; t < n; t++) Ta[t] = (Tx[t] + Tn[t]) / 2;

  // seasonal LAI
  let LAI = new Float64Array(n);
  for (let t = 0; t < n; t++) {
    const fOn = clamp((doy[t] - p.leafOn) / Math.max(p.ramp, 1), 0, 1);
    const fOff = clamp((doy[t] - p.leafOff) / Math.max(p.ramp, 1), 0, 1);
    LAI[t] = Math.max(p.LAImin + (p.LAImax - p.LAImin) * fOn * (1 - fOff), 0.01);
  }

  p.albSnow = dflt(p.albSnow, p.albedo); p.kAlb = dflt(p.kAlb, 0.75); p.sweFull = dflt(p.sweFull, 20);
  p.cSnow = dflt(p.cSnow, 3.5); p.cSnowBr = dflt(p.cSnowBr, 1.0); p.fEv = dflt(p.fEv, 0);
  p.unload = dflt(p.unload, 0.5); p.albMode = dflt(p.albMode, 0);
  let alb0 = snowfreeAlbedo(p, LAI, n);
  const albS = new Float64Array(n);
  for (let t = 0; t < n; t++) albS[t] = alb0[t] + (p.albSnow - alb0[t]) * Math.exp(-p.kAlb * LAI[t]);
  const pm0 = pmPet(Ta, Tx, Tn, F.Rs, F.ea, F.u, LAI, doy, F.lat, F.elev, p, alb0);
  const pmS = pmPet(Ta, Tx, Tn, F.Rs, F.ea, F.u, LAI, doy, F.lat, F.elev, p, albS);
  const PET0 = pm0.PET, Ewet0 = pm0.Ewet, Rn0 = pm0.Rn, PETs = pmS.PET, EwetS = pmS.Ewet, RnS = pmS.Rn;
  const PET = Float64Array.from(PET0), Ewet = Float64Array.from(Ewet0), Rn = Float64Array.from(Rn0);
  const alb = Float64Array.from(alb0), fscv = new Float64Array(n);

  // impervious area
  const useImp = !!p.impOn;
  let fI = 0;
  if (useImp) { p.fConn = dflt(p.fConn, 1); p.Sdep = dflt(p.Sdep, 1.5); fI = clamp(p.impPct / 100, 0, 0.99); }
  const fP = 1 - fI;
  const fC = useImp ? clamp(p.fConn, 0, 1) : 1;
  const canPerv = useImp && fI > 0 && !!p.canopyPervOnly;
  const gC = canPerv ? fP : 1;
  let Dimp = 0;

  // TOPMODEL setup
  const nc = TI.val.length, ti = Float64Array.from(TI.val), af = new Float64Array(nc);
  let sf = 0; for (let c = 0; c < nc; c++) sf += TI.frac[c];
  let lam = 0; for (let c = 0; c < nc; c++) { af[c] = TI.frac[c] / sf; lam += af[c] * ti[c]; }
  const szq = Math.exp(p.lnTe - lam);
  const q0 = Math.max(p.Q0, 1e-6) / 1000;
  let Sbar = Math.max(-p.m * Math.log(q0 / szq), 0);
  let SRZ = new Float64Array(nc).fill(p.SR0), SUZ = new Float64Array(nc);
  let SWE = 0, C = 0, Cs = 0, SWEi = 0;

  const useSSF = !!p.ssfOn;
  let lnLam = NaN, Mssf = NaN, fOut = NaN, shareQ = NaN;
  if (useSSF) {
    if (TI.lnLambdaLin !== undefined && TI.lnLambdaLin !== null) lnLam = TI.lnLambdaLin;
    else { let mx = -Infinity; for (const v of ti) mx = Math.max(mx, v); let s = 0; for (let c = 0; c < nc; c++) s += af[c] * Math.exp(ti[c] - mx); lnLam = mx + Math.log(s); }
    Mssf = Math.max(p.Mssf, 1e-6);
    fOut = 1 - Math.exp(-(p.kssf + p.kv));
    shareQ = p.kssf / Math.max(p.kssf + p.kv, Number.EPSILON);
  }
  let W = 0;

  const out = {};
  for (const k of OUT_KEYS) out[k] = new Float64Array(n);

  // infiltration excess
  const useIE = !!p.ieOn;
  p.Ks = dflt(p.Ks, 20); p.ieMoist = dflt(p.ieMoist, 1); p.stormSum = dflt(p.stormSum, 4);
  p.stormWin = dflt(p.stormWin, 10); p.ieConn = dflt(p.ieConn, 1);
  const Dh = new Float64Array(n);
  for (let t = 0; t < n; t++) Dh[t] = Math.max(p.stormWin + (p.stormSum - p.stormWin) * 0.5 * (1 - Math.cos(2 * Math.PI * (doy[t] - 15) / 365)), 0.25);

  const useHarv = p.luc && p.luc.on;
  let aH = 0;
  const DDF = p.DDF, Tmelt = p.Tmelt, td = p.td, mm = p.m, Srmax = p.Srmax;
  const SD = new Float64Array(nc), SDe = new Float64Array(nc), ex = new Float64Array(nc), uz = new Float64Array(nc),
        ea = new Float64Array(nc), satSSF = new Uint8Array(nc), unsat = new Uint8Array(nc), ieV = new Float64Array(nc),
        Kv = new Float64Array(nc);

  // shared lower part of the step (unsaturated zone, root zone ET, baseflow, SSF); returns fluxes
  function lowerStep(srzEaFn) {
    for (let c = 0; c < nc; c++) SD[c] = Math.max(Sbar + mm * (lam - ti[c]), 0);
    satSSF.fill(0);
    if (useSSF && W > 0) { const thr = lnLam + Math.log(Mssf / W); for (let c = 0; c < nc; c++) satSSF[c] = ti[c] >= thr ? 1 : 0; }
    let sEx = 0, quz = 0, sat = 0, satS = 0;
    for (let c = 0; c < nc; c++) {
      SDe[c] = satSSF[c] ? 0 : SD[c];
      ex[c] = Math.max(SUZ[c] - SDe[c], 0); SUZ[c] -= ex[c];
      uz[c] = SDe[c] > 0 ? Math.min(SUZ[c] / (SD[c] * td), SUZ[c]) : 0;
      SUZ[c] -= uz[c];
      sEx += af[c] * ex[c]; quz += af[c] * uz[c];
      if (SD[c] <= 0 || satSSF[c]) sat += af[c];
      if (satSSF[c] && SD[c] > 0) satS += af[c];
    }
    const sEa = srzEaFn();
    const qb = szq * Math.exp(-Sbar / mm);
    let qssf = 0, qv = 0, ovf = 0;
    if (useSSF) {
      W += quz;
      if (W > Mssf) { ovf = W - Mssf; W = Mssf; }
      const dW = W * fOut; qssf = dW * shareQ; qv = dW - qssf; W -= dW;
      Sbar = Sbar - qv + qb;
    } else Sbar = Sbar - quz + qb;
    return { sEx, quz, sat, satS, sEa, qb, qssf, qv, ovf };
  }

  if (!useHarv) {
    for (let t = 0; t < n; t++) {
      let sn, rn;
      if (Ta[t] <= p.Tsnow) { sn = P[t]; rn = 0; } else { sn = 0; rn = P[t]; }
      const woody = (1 - p.fEv) * Math.min(p.LAImin, LAI[t]);
      const CsMax = p.cSnow * (LAI[t] - woody) + p.cSnowBr * woody;
      let Is = 0;
      if (sn > 0 && CsMax > Cs) { const Cc = 1 - Math.exp(-0.5 * LAI[t]); Is = Math.min(sn, (CsMax - Cs) * (1 - Math.exp(-Cc * sn / CsMax))); }
      Cs += Is;
      let unl = Cs * (1 - Math.exp(-p.unload)); Cs -= unl;
      const dmelt = Math.max(DDF * (Ta[t] - Tmelt), 0);
      const drip = Math.min(Cs, dmelt); Cs -= drip;
      if (Cs > CsMax) { unl += Cs - CsMax; Cs = CsMax; }
      SWE += (sn - Is) + unl;
      const fsc = Math.min(SWE / Math.max(p.sweFull, 1e-6), 1);
      if (fsc > 0) {
        PET[t] = PET0[t] + fsc * (PETs[t] - PET0[t]); Ewet[t] = Ewet0[t] + fsc * (EwetS[t] - Ewet0[t]);
        Rn[t] = Rn0[t] + fsc * (RnS[t] - Rn0[t]); alb[t] = alb0[t] + fsc * (albS[t] - alb0[t]); fscv[t] = fsc;
      }
      const mlt = Math.min(SWE, dmelt); SWE -= mlt;
      let mltI = 0;
      if (canPerv) { SWEi += sn; mltI = Math.min(SWEi, dmelt); SWEi -= mltI; }
      const Cmax = p.cint * LAI[t];
      const I = Math.min(rn, Math.max(Cmax - C, 0));
      C += I; const thr = rn - I;
      const Ei = Math.min(C, Ewet[t]); C -= Ei;
      const Es = Math.min(Cs, (Ewet[t] - Ei) * LV_LS); Cs -= Es;
      const fw = Ewet[t] > 0 ? Math.min((Ei + Es / LV_LS) / Ewet[t], 1) : 0;
      const Edry = PET[t] * (1 - fw);

      const win = thr + mlt + drip;
      let qimp = 0, eimp = 0, ron = 0, xs = 0;
      if (useImp && fI > 0) {
        let eDem;
        if (canPerv) { Dimp += rn + mltI; eDem = Ewet[t]; }
        else { Dimp += win; eDem = Math.max(Ewet[t] - Ei - Es / LV_LS, 0); }
        xs = Math.max(Dimp - p.Sdep, 0); Dimp -= xs;
        eimp = Math.min(Dimp, eDem); Dimp -= eimp;
        qimp = fC * xs; ron = (1 - fC) * xs * fI / fP;
      }
      let sIe = 0;
      ieV.fill(0);
      if (useIE && thr > 0) {
        for (let c = 0; c < nc; c++) { Kv[c] = p.Ks * (1 + p.ieMoist * Math.min(SRZ[c] / Srmax, 1)); unsat[c] = (Sbar + mm * (lam - ti[c]) > 0) ? 1 : 0; }
        ieRunoff(thr, Dh[t], Kv, unsat, ieV);
        for (let c = 0; c < nc; c++) { ieV[c] *= p.ieConn; sIe += af[c] * ieV[c]; }
      }
      const ep = Edry / 1000;
      for (let c = 0; c < nc; c++) {
        SRZ[c] -= (win + ron - ieV[c]) / 1000;
        if (SRZ[c] < 0) { SUZ[c] -= SRZ[c]; SRZ[c] = 0; }
      }
      const L = lowerStep(() => {
        let s = 0;
        for (let c = 0; c < nc; c++) {
          let e = ep * (1 - SRZ[c] / Srmax); e = Math.min(Math.max(e, 0), Srmax - SRZ[c]);
          SRZ[c] += e; s += af[c] * e;
        }
        return s;
      });
      out.rain[t] = rn; out.snow[t] = sn;
      out.melt[t] = gC * mlt + (1 - gC) * mltI; out.SWE[t] = gC * SWE + (1 - gC) * SWEi;
      out.Isnow[t] = gC * Is; out.Es[t] = gC * Es; out.Csnow[t] = gC * Cs;
      out.unloaded[t] = gC * unl; out.drip[t] = gC * drip;
      out.Ei[t] = gC * Ei; out.Ea[t] = fP * 1000 * L.sEa; out.Ccan[t] = gC * C;
      out.Qof[t] = fP * 1000 * (L.sEx + L.ovf); out.Qb[t] = fP * 1000 * L.qb; out.Quz[t] = fP * 1000 * L.quz;
      out.Qssf[t] = fP * 1000 * L.qssf; out.Wssf[t] = 1000 * W; out.qv[t] = fP * 1000 * L.qv;
      out.Qimp[t] = fI * qimp; out.Eimp[t] = fI * eimp; out.Dimp[t] = fI * Dimp;
      out.runon[t] = fP * ron; out.Qinf[t] = fI * xs; out.Qie[t] = fP * sIe;
      out.Sbar[t] = 1000 * Sbar; out.satFrac[t] = fP * L.sat; out.satFracSSF[t] = fP * L.satS;
      let sR = 0, sU = 0; for (let c = 0; c < nc; c++) { sR += af[c] * SRZ[c]; sU += af[c] * SUZ[c]; }
      out.SRZ[t] = 1000 * sR; out.SUZ[t] = 1000 * sU;
      out.storage[t] = gC * SWE + (1 - gC) * SWEi + gC * C + gC * Cs + fI * Dimp + fP * 1000 * (sU - sR - Sbar + W);
    }
  } else {
    // ---------------------------------------------- land-use change: two canopy tiles
    const lu = p.luc;
    let h = Float64Array.from(lu.h.length === nc ? lu.h : new Array(nc).fill(lu.h.reduce((a, b) => a + b, 0) / lu.h.length));
    for (let c = 0; c < nc; c++) h[c] = clamp(h[c], 0, 1);
    aH = 0; for (let c = 0; c < nc; c++) aH += af[c] * h[c];
    const aF = 1 - aH;
    const tdays = new Float64Array(n), chg = new Uint8Array(n), rr = new Float64Array(n);
    const curve = (lu.curve || "none").toLowerCase();
    for (let t = 0; t < n; t++) {
      tdays[t] = (t + 1) - lu.tIdx;
      chg[t] = (lu.change && tdays[t] >= 0) ? 1 : 0;
      rr[t] = curve === "linear" ? clamp(tdays[t] / Math.max(lu.Tdays, 1), 0, 1)
            : curve === "exponential" ? 1 - Math.exp(-3 * Math.max(tdays[t], 0) / Math.max(lu.Tdays, 1))
            : curve === "step" ? 1 : 0;
    }
    const X0 = lucState(lu.pre, LAI, alb0, p, doy), Xs = lucState(lu.start, LAI, alb0, p, doy), Xe = lucState(lu.final, LAI, alb0, p, doy);
    const XH = {};
    for (const nm of Object.keys(X0)) {
      const v = Float64Array.from(X0[nm]);
      for (let t = 0; t < n; t++) if (chg[t]) v[t] = Xs[nm][t] + (Xe[nm][t] - Xs[nm][t]) * rr[t];
      XH[nm] = v;
    }
    const finiteRec = Number.isFinite(lu.compRecDays);
    for (let t = 0; t < n; t++) if (chg[t]) {
      const rc = finiteRec ? clamp(tdays[t] / Math.max(lu.compRecDays, 1), 0, 1) : 0;
      XH.comp[t] = Xs.comp[t] + (Xe.comp[t] - Xs.comp[t]) * rc;
      XH.compK[t] = Xs.compK[t];
    }
    const LAIH = Float64Array.from(XH.LAI);
    for (let t = 0; t < n; t++) { if (lu.laiMult && lu.laiMult.length === n) LAIH[t] *= lu.laiMult[t]; LAIH[t] = Math.max(LAIH[t], 0.01); }
    const fImpH = new Float64Array(n);
    if (lu.impFrac > 0) for (let t = 0; t < n; t++) if (chg[t]) fImpH[t] = lu.impFrac * clamp((tdays[t] + 1) / Math.max(lu.impRampDays, 1), 0, 1);
    const alb0H = XH.alb, albSH = new Float64Array(n);
    for (let t = 0; t < n; t++) albSH[t] = alb0H[t] + (p.albSnow - alb0H[t]) * Math.exp(-p.kAlb * LAIH[t]);
    const pmH0 = pmPet(Ta, Tx, Tn, F.Rs, F.ea, F.u, LAIH, doy, F.lat, F.elev, p, alb0H, XH.hc);
    const pmHS = pmPet(Ta, Tx, Tn, F.Rs, F.ea, F.u, LAIH, doy, F.lat, F.elev, p, albSH, XH.hc);
    const LAIminH = new Float64Array(n), SrH = new Float64Array(n);
    for (let t = 0; t < n; t++) { LAIminH[t] = Math.min(p.LAImin * LAIH[t] / LAI[t], LAIH[t]); SrH[t] = Srmax * XH.sr[t]; }
    const SRZf = Float64Array.from(SRZ), SRZh = Float64Array.from(SRZ);
    let SWEf = 0, SWEh = 0, Cf = 0, Ch = 0, Csf = 0, Csh = 0, Sdet = 0; SWEi = 0;
    const useDet = !!lu.detention;
    const kdet = useDet ? 1 - Math.exp(-Math.max(lu.kdet, 1e-6)) : 0;
    const ieF = new Float64Array(nc), oF = new Float64Array(nc), oH = new Float64Array(nc), KF = new Float64Array(nc), KH = new Float64Array(nc), ieHv = new Float64Array(nc);
    const LAIout = new Float64Array(n), alb0out = new Float64Array(n);
    for (let t = 0; t < n; t++) {
      let sn, rn;
      if (Ta[t] <= p.Tsnow) { sn = P[t]; rn = 0; } else { sn = 0; rn = P[t]; }
      const wF = (1 - p.fEv) * Math.min(p.LAImin, LAI[t]);
      const sF = canopyStep(Cf, Csf, SWEf, sn, rn, LAI[t], wF, p, DDF, Ta[t], PET0[t], PETs[t], Ewet0[t], EwetS[t], Rn0[t], RnS[t], alb0[t], albS[t]);
      Cf = sF.C; Csf = sF.Cs; SWEf = sF.SWE;
      const wH = (1 - XH.fEv[t]) * LAIminH[t];
      const sH = canopyStep(Ch, Csh, SWEh, sn, rn, LAIH[t], wH, p, XH.ddf[t], Ta[t], pmH0.PET[t], pmHS.PET[t], pmH0.Ewet[t], pmHS.Ewet[t], pmH0.Rn[t], pmHS.Rn[t], alb0H[t], albSH[t]);
      Ch = sH.C; Csh = sH.Cs; SWEh = sH.SWE;
      const winF = sF.thr + sF.mlt + sF.drip, winH = sH.thr + sH.mlt + sH.drip;
      const qI = fImpH[t] * Math.max(winH - lu.impSdep, 0);
      const qIc = lu.impConn * qI;
      const fCmp = (1 - fImpH[t]) * XH.comp[t];
      let xsC = 0;
      if (fCmp > 0 && sH.thr > 0) xsC = fCmp * ieRunoffScalar(sH.thr, Dh[t], XH.compK[t]);
      const winHs = winH - qIc - xsC;
      let qDev;
      if (useDet) { Sdet += qIc; qDev = Sdet * kdet; Sdet -= qDev; } else qDev = qIc;
      let qimp = 0, eimp = 0, ron = 0, xs = 0, mltI = 0;
      if (canPerv && fI > 0) { SWEi += sn; mltI = Math.min(SWEi, Math.max(DDF * (Ta[t] - Tmelt), 0)); SWEi -= mltI; }
      if (useImp && fI > 0) {
        let eDem;
        if (canPerv) { Dimp += rn + mltI; eDem = sF.Ewet; }
        else {
          Dimp += aF * winF + aH * winH;
          eDem = aF * Math.max(sF.Ewet - sF.Ei - sF.Es / LV_LS, 0) + aH * Math.max(sH.Ewet - sH.Ei - sH.Es / LV_LS, 0);
        }
        xs = Math.max(Dimp - p.Sdep, 0); Dimp -= xs;
        eimp = Math.min(Dimp, eDem); Dimp -= eimp;
        qimp = fC * xs; ron = (1 - fC) * xs * fI / fP;
      }
      ieF.fill(0); ieHv.fill(0);
      if (useIE) {
        for (let c = 0; c < nc; c++) unsat[c] = (Sbar + mm * (lam - ti[c]) > 0) ? 1 : 0;
        if (sF.thr > 0) {
          for (let c = 0; c < nc; c++) KF[c] = p.Ks * (1 + p.ieMoist * Math.min(SRZf[c] / Srmax, 1));
          ieRunoff(sF.thr, Dh[t], KF, unsat, ieF); for (let c = 0; c < nc; c++) ieF[c] *= p.ieConn;
        }
        const fNat = 1 - fImpH[t] - fCmp;
        if (sH.thr > 0 && fNat > 0) {
          for (let c = 0; c < nc; c++) KH[c] = p.Ks * (1 + p.ieMoist * Math.min(SRZh[c] / SrH[t], 1));
          ieRunoff(sH.thr, Dh[t], KH, unsat, ieHv); for (let c = 0; c < nc; c++) ieHv[c] *= fNat * p.ieConn;
        }
      }
      for (let c = 0; c < nc; c++) {
        SRZf[c] -= (winF + ron - ieF[c]) / 1000; SRZh[c] -= (winHs + ron - ieHv[c]) / 1000;
        oF[c] = Math.max(-SRZf[c], 0); SRZf[c] = Math.max(SRZf[c], 0);
        oH[c] = Math.max(-SRZh[c], 0); SRZh[c] = Math.max(SRZh[c], 0);
        SUZ[c] += (1 - h[c]) * oF[c] + h[c] * oH[c];
      }
      const L = lowerStep(() => {
        let s = 0;
        for (let c = 0; c < nc; c++) {
          let eF = (sF.Edry / 1000) * (1 - SRZf[c] / Srmax); eF = Math.min(Math.max(eF, 0), Math.max(Srmax - SRZf[c], 0));
          let eH = (sH.Edry * (1 - fImpH[t]) / 1000) * (1 - SRZh[c] / SrH[t]); eH = Math.min(Math.max(eH, 0), Math.max(SrH[t] - SRZh[c], 0));
          SRZf[c] += eF; SRZh[c] += eH;
          s += af[c] * ((1 - h[c]) * eF + h[c] * eH);
        }
        return s;
      });
      const wa = (xf, xh) => gC * (aF * xf + aH * xh);
      const SWEt = wa(SWEf, SWEh) + (1 - gC) * SWEi;
      out.rain[t] = rn; out.snow[t] = sn; out.SWE[t] = SWEt;
      out.melt[t] = wa(sF.mlt, sH.mlt) + (1 - gC) * mltI;
      out.Isnow[t] = wa(sF.Is, sH.Is); out.Es[t] = wa(sF.Es, sH.Es); out.Csnow[t] = wa(Csf, Csh);
      out.unloaded[t] = wa(sF.unl, sH.unl); out.drip[t] = wa(sF.drip, sH.drip);
      out.Ei[t] = wa(sF.Ei, sH.Ei); out.Ccan[t] = wa(Cf, Ch);
      out.Ea[t] = fP * 1000 * L.sEa;
      out.Qof[t] = fP * 1000 * (L.sEx + L.ovf); out.Qb[t] = fP * 1000 * L.qb; out.Quz[t] = fP * 1000 * L.quz;
      out.Qssf[t] = fP * 1000 * L.qssf; out.Wssf[t] = 1000 * W; out.qv[t] = fP * 1000 * L.qv;
      out.Qimp[t] = fI * qimp; out.Eimp[t] = fI * eimp; out.Dimp[t] = fI * Dimp;
      out.runon[t] = fP * ron; out.Qinf[t] = fI * xs;
      out.Qcomp[t] = fP * aH * xsC; out.Qdev[t] = fP * aH * qDev; out.Sdet[t] = fP * aH * Sdet;
      let sIe = 0; for (let c = 0; c < nc; c++) sIe += af[c] * ((1 - h[c]) * ieF[c] + h[c] * ieHv[c]);
      out.Qie[t] = fP * sIe;
      out.Sbar[t] = 1000 * Sbar; out.satFrac[t] = fP * L.sat; out.satFracSSF[t] = fP * L.satS;
      let sR = 0, sU = 0; for (let c = 0; c < nc; c++) { sR += af[c] * ((1 - h[c]) * SRZf[c] + h[c] * SRZh[c]); sU += af[c] * SUZ[c]; }
      out.SRZ[t] = 1000 * sR; out.SUZ[t] = 1000 * sU;
      out.storage[t] = SWEt + wa(Cf + Csf, Ch + Csh) + fI * Dimp + fP * aH * Sdet + fP * 1000 * (sU - sR - Sbar + W);
      PET[t] = aF * sF.PET + aH * sH.PET; Ewet[t] = aF * sF.Ewet + aH * sH.Ewet;
      Rn[t] = aF * sF.Rn + aH * sH.Rn; alb[t] = aF * sF.alb + aH * sH.alb; fscv[t] = aF * sF.fsc + aH * sH.fsc;
      LAIout[t] = aF * LAI[t] + aH * LAIH[t]; alb0out[t] = aF * alb0[t] + aH * alb0H[t];
    }
    out.LAIharv = LAIH; out.impChanged = fImpH;
    alb0 = alb0out; LAI = LAIout;
  }
  out.lucOn = !!useHarv; out.lucFrac = aH;
  out.Qtot = new Float64Array(n);
  for (let t = 0; t < n; t++) out.Qtot[t] = out.Qof[t] + out.Qb[t] + out.Qssf[t] + out.Qimp[t] + out.Qcomp[t] + out.Qdev[t] + out.Qie[t];
  out.ieOn = useIE; out.stormHours = Dh; out.ssfOn = useSSF;
  out.impOn = useImp && fI > 0; out.impFrac = fI; out.canopyPervOnly = canPerv;

  // channel routing (time-delay histogram)
  out.routed = false; out.uh = [1];
  if (R && R.on && R.dist && R.dist.length) {
    const tt = R.dist.map(d => d / Math.max(p.chv, 1e-3) / 86400);
    let Lmax = 0; for (const v of tt) Lmax = Math.max(Lmax, Math.floor(v));
    const uh = new Float64Array(Lmax + 2);
    for (let i = 0; i < tt.length; i++) { const L = Math.floor(tt[i]), fr = tt[i] - L; uh[L] += R.frac[i] * (1 - fr); uh[L + 1] += R.frac[i] * fr; }
    let s = 0; for (const v of uh) s += v; for (let i = 0; i < uh.length; i++) uh[i] /= s;
    const q = new Float64Array(n);
    for (let k = 0; k < uh.length; k++) { if (uh[k] === 0 || k >= n) continue; for (let t = k; t < n; t++) q[t] += uh[k] * out.Qtot[t - k]; }
    out.Qsim = q; out.uh = Array.from(uh); out.routed = true;
    out.chanStore = new Float64Array(n); let a = 0, b = 0;
    for (let t = 0; t < n; t++) { a += out.Qtot[t]; b += q[t]; out.chanStore[t] = a - b; }
  } else { out.Qsim = Float64Array.from(out.Qtot); out.chanStore = new Float64Array(n); }
  out.PET = PET; out.Ewet = Ewet; out.LAI = LAI; out.ra = pm0.ra; out.rc = pm0.rc; out.Rn = Rn;
  out.ET = new Float64Array(n);
  for (let t = 0; t < n; t++) out.ET[t] = out.Ei[t] + out.Es[t] + out.Ea[t] + out.Eimp[t];
  out.lambda = lam; out.albedo = alb; out.snowCover = fscv; out.albedo0 = alb0; out.albMode = p.albMode;
  out.S0 = fP * (p.SR0 * (-1000) - 1000 * Math.max(-p.m * Math.log(q0 / szq), 0));
  return out;
}

// ---------------------------------------------------------------- metrics
function mean(a) { let s = 0; for (const v of a) s += v; return s / a.length; }
function sd(a) { const m = mean(a); let s = 0; for (const v of a) s += (v - m) * (v - m); return Math.sqrt(s / (a.length - 1)); }
function cor(a, b) {
  const ma = mean(a), mb = mean(b); let sab = 0, saa = 0, sbb = 0;
  for (let i = 0; i < a.length; i++) { const x = a[i] - ma, y = b[i] - mb; sab += x * y; saa += x * x; sbb += y * y; }
  return sab / Math.sqrt(saa * sbb);
}
function nse(o, s) { const m = mean(o); let a = 0, b = 0; for (let i = 0; i < o.length; i++) { a += (s[i] - o[i]) ** 2; b += (o[i] - m) ** 2; } return 1 - a / b; }
function perfMetrics(oIn, sIn) {
  const o = [], s = [];
  for (let i = 0; i < oIn.length; i++) if (oIn[i] !== null && !Number.isNaN(oIn[i]) && !Number.isNaN(sIn[i])) { o.push(oIn[i]); s.push(sIn[i]); }
  const M = { n: o.length, mean_obs: mean(o), mean_sim: mean(s) };
  M.NSE = nse(o, s);
  const e = 0.01 * M.mean_obs;
  M.NSE_log = nse(o.map(v => Math.log(Math.max(v, 0) + e)), s.map(v => Math.log(Math.max(v, 0) + e)));
  M.NSE_sqrt = nse(o.map(v => Math.sqrt(Math.max(v, 0))), s.map(v => Math.sqrt(Math.max(v, 0))));
  const so = sd(o), ss = sd(s);
  M.KGE_r = (M.n > 2 && so > 0 && ss > 0) ? cor(o, s) : NaN;
  M.KGE_alpha = ss / so; M.KGE_beta = M.mean_sim / M.mean_obs;
  M.KGE = 1 - Math.sqrt((M.KGE_r - 1) ** 2 + (M.KGE_alpha - 1) ** 2 + (M.KGE_beta - 1) ** 2);
  let so2 = 0, ss2 = 0, se = 0, sa = 0;
  for (let i = 0; i < o.length; i++) { so2 += o[i]; ss2 += s[i]; se += (s[i] - o[i]) ** 2; sa += Math.abs(s[i] - o[i]); }
  M.PBIAS = 100 * (ss2 - so2) / so2; M.RMSE = Math.sqrt(se / o.length); M.MAE = sa / o.length; M.R2 = M.KGE_r * M.KGE_r;
  return M;
}
function objective(name, o, s) {
  const M = perfMetrics(o, s);
  switch (name) {
    case "NSE": return M.NSE;
    case "NSE of log Q": return M.NSE_log;
    case "Mean of NSE and log NSE": return (M.NSE + M.NSE_log) / 2;
    case "KGE of sqrt Q": return perfMetrics(o.map(v => Math.sqrt(v)), Array.from(s, v => Math.sqrt(Math.max(v, 0)))).KGE;
    default: return M.KGE;
  }
}

// ---------------------------------------------------------------- SCE-UA
function mulberry32(a) { return function () { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
function sceua(fun, x0, bl, bu, opts) {
  opts = Object.assign({ ngs: 2, maxn: 3000, kstop: 5, pcento: 0.1, peps: 1e-3, seed: 1, includeX0: true, outFcn: null }, opts || {});
  const rnd = mulberry32(opts.seed >>> 0 || 1);
  const n = bl.length, npg = 2 * n + 1, nps = n + 1, nspl = npg, ngs = Math.max(1, Math.round(opts.ngs)), npt = npg * ngs;
  const bound = bu.map((u, i) => u - bl[i]);
  const feval = x => { let f; try { f = fun(x); } catch (e) { f = Infinity; } if (typeof f !== "number" || Number.isNaN(f)) f = Infinity; return f; };
  const geoRange = X => { let s = 0; for (let j = 0; j < n; j++) { let mx = -Infinity, mn = Infinity; for (const r of X) { mx = Math.max(mx, r[j]); mn = Math.min(mn, r[j]); } s += Math.log(Math.max((mx - mn) / bound[j], Number.MIN_VALUE)); } return Math.exp(s / n); };
  const randPt = () => bl.map((b, j) => b + rnd() * bound[j]);
  function cce(s, sf) {
    const sw = s[s.length - 1], fw = sf[sf.length - 1];
    const ce = new Array(n).fill(0);
    for (let i = 0; i < s.length - 1; i++) for (let j = 0; j < n; j++) ce[j] += s[i][j] / (s.length - 1);
    let snew = ce.map((c, j) => c + (c - sw[j]));
    if (snew.some((v, j) => v < bl[j] || v > bu[j])) snew = randPt();
    let fnew = feval(snew), ne = 1;
    if (fnew > fw) {
      snew = ce.map((c, j) => c - 0.5 * (c - sw[j])); fnew = feval(snew); ne++;
      if (fnew > fw) { snew = randPt(); fnew = feval(snew); ne++; }
    }
    return { x: snew, f: fnew, ne };
  }
  let x = [], xf = [];
  for (let i = 0; i < npt; i++) x.push(randPt());
  if (opts.includeX0 && x0 && x0.length) x[0] = x0.map((v, j) => clamp(v, bl[j], bu[j]));
  for (let i = 0; i < npt; i++) xf.push(feval(x[i]));
  let icall = npt;
  const sortPop = () => { const o = xf.map((f, i) => i).sort((a, b) => xf[a] - xf[b]); x = o.map(i => x[i]); xf = o.map(i => xf[i]); };
  sortPop();
  let loop = 0, gnrng = geoRange(x), change = Infinity, exitflag = "";
  const hist = [[loop, icall, xf[0], xf[npt - 1], gnrng]];
  const callOut = () => opts.outFcn ? !!opts.outFcn({ nevals: icall, loop, bestf: xf[0], bestx: x[0], gnrng, change }) : false;
  if (callOut()) exitflag = "stopped";
  while (exitflag === "") {
    loop++;
    for (let igs = 0; igs < ngs; igs++) {
      const k2 = []; for (let k = 0; k < npg; k++) k2.push(k * ngs + igs);
      let cx = k2.map(i => x[i].slice()), cf = k2.map(i => xf[i]);
      for (let it = 0; it < nspl; it++) {
        const lcs = [0];
        for (let k3 = 1; k3 < nps; k3++) {
          let lpos = 0;
          for (let tries = 0; tries < 1000; tries++) {
            lpos = Math.floor(npg + 0.5 - Math.sqrt((npg + 0.5) ** 2 - npg * (npg + 1) * rnd()));
            if (!lcs.includes(lpos)) break;
          }
          lcs.push(lpos);
        }
        lcs.sort((a, b) => a - b);
        const s = lcs.map(i => cx[i]), sfv = lcs.map(i => cf[i]);
        const r = cce(s, sfv); icall += r.ne;
        cx[lcs[lcs.length - 1]] = r.x; cf[lcs[lcs.length - 1]] = r.f;
        const o = cf.map((f, i) => i).sort((a, b) => cf[a] - cf[b]); cx = o.map(i => cx[i]); cf = o.map(i => cf[i]);
        if (icall >= opts.maxn) break;
      }
      k2.forEach((i, j) => { x[i] = cx[j]; xf[i] = cf[j]; });
      if (icall >= opts.maxn) break;
    }
    sortPop();
    gnrng = geoRange(x);
    hist.push([loop, icall, xf[0], xf[npt - 1], gnrng]);
    if (loop >= opts.kstop) {
      const a = hist.slice(hist.length - opts.kstop).map(h => h[2]);
      let den = mean(a.map(Math.abs)); if (den === 0) den = 1;
      change = Math.abs(a[0] - a[a.length - 1]) * 100 / den;
    }
    if (callOut()) exitflag = "stopped";
    else if (icall >= opts.maxn) exitflag = "maxn";
    else if (gnrng < opts.peps) exitflag = "population";
    else if (change < opts.pcento) exitflag = "converged";
  }
  const msgs = { maxn: "maximum number of model runs reached", converged: `best value changed < ${opts.pcento}% over ${opts.kstop} shuffles`,
                 population: "population converged", stopped: "stopped by the user" };
  return { bestx: x[0], bestf: xf[0], nevals: icall, loops: loop, exitflag, message: msgs[exitflag], history: hist };
}

const api = { core, pmPet, perfMetrics, objective, sceua, mulberry32, lucState };
if (typeof module !== "undefined" && module.exports) module.exports = api;
root.TM = api;
})(typeof self !== "undefined" ? self : this);
