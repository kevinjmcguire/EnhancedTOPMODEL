"""Topographic index and display map from a projected DEM (no GIS dependencies beyond numpy).

Steps: priority-flood depression filling with a small gradient across flats (Barnes et al. 2014),
multiple-flow-direction accumulation with Quinn et al. (1991) slope weighting (exponent 1.1, Freeman 1991),
TI = ln(a / tan b) with a = upslope area per unit contour width and tan b the flow-weighted downslope gradient.
"""
import heapq
import math
import numpy as np

NB = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
# contour-length weights for cardinal / diagonal neighbours (Quinn et al. 1991)
CL = [0.354, 0.5, 0.354, 0.5, 0.5, 0.354, 0.5, 0.354]


def fill_depressions(z, valid, eps=1e-4):
    """Priority-flood fill from the edge of `valid`; every cell gets a strictly downhill path to the edge."""
    H, W = z.shape
    out = np.where(valid, z, np.nan).astype(float)
    done = ~valid.copy()
    pq = []
    cnt = 0
    for y in range(H):
        for x in range(W):
            if not valid[y, x]:
                continue
            edge = y == 0 or x == 0 or y == H - 1 or x == W - 1
            if not edge:
                for dy, dx in NB:
                    if not valid[y + dy, x + dx]:
                        edge = True
                        break
            if edge:
                heapq.heappush(pq, (out[y, x], cnt, y, x)); cnt += 1
                done[y, x] = True
    while pq:
        e, _, y, x = heapq.heappop(pq)
        for dy, dx in NB:
            yy, xx = y + dy, x + dx
            if 0 <= yy < H and 0 <= xx < W and not done[yy, xx]:
                done[yy, xx] = True
                if out[yy, xx] <= e:
                    out[yy, xx] = e + eps
                heapq.heappush(pq, (out[yy, xx], cnt, yy, xx)); cnt += 1
    return out


def topographic_index(z, valid, cell, p=1.1, min_tan=1e-3):
    """Returns (ti, acc_cells) arrays (NaN / 0 outside `valid`)."""
    H, W = z.shape
    zf = fill_depressions(z, valid)
    acc = np.where(valid, 1.0, 0.0)
    tanb = np.full((H, W), np.nan)
    order = np.argsort(np.where(valid, zf, -np.inf), axis=None)[::-1]
    dist = [math.hypot(dy, dx) * cell for dy, dx in NB]
    for idx in order:
        y, x = divmod(int(idx), W)
        if not valid[y, x]:
            break
        z0 = zf[y, x]
        ws, js, ts = [], [], []
        for k, (dy, dx) in enumerate(NB):
            yy, xx = y + dy, x + dx
            if 0 <= yy < H and 0 <= xx < W and valid[yy, xx]:
                g = (z0 - zf[yy, xx]) / dist[k]
                if g > 0:
                    ws.append((g ** p) * CL[k]); js.append((yy, xx)); ts.append(g)
        if ws:
            s = sum(ws)
            tanb[y, x] = sum(w * t for w, t in zip(ws, ts)) / s
            a = acc[y, x]
            for w, (yy, xx) in zip(ws, js):
                acc[yy, xx] += a * w / s
        else:
            # basin edge outlet: use the steepest gradient towards any neighbour (or the floor value)
            tanb[y, x] = min_tan
    tanb = np.where(valid, np.maximum(tanb, min_tan), np.nan)
    a_spec = acc * cell  # upslope area per unit contour width, m
    ti = np.where(valid, np.log(a_spec / tanb), np.nan)
    return ti, acc


def ti_classes(ti, n_cls=40):
    v = ti[np.isfinite(ti)]
    edges = np.linspace(v.min(), v.max() + 1e-9, n_cls + 1)
    k = np.clip(np.digitize(v, edges) - 1, 0, n_cls - 1)
    val, frac = [], []
    for c in range(n_cls):
        sel = v[k == c]
        if len(sel):
            val.append(round(float(sel.mean()), 3)); frac.append(round(len(sel) / len(v), 5))
    lam = float(np.dot(val, frac) / sum(frac))
    lnL = float(np.log(np.mean(np.exp(v))))
    return val, frac, round(lam, 3), round(lnL, 3)


def display_map(ti, valid, cell, max_side=150):
    """Coarsen to <= max_side cells per side: mean TI of in-basin cells per block, quantized as the page expects
    (v in 0..30 with TI = 2 + v/30*20; 255 = outside)."""
    H, W = ti.shape
    k = max(1, math.ceil(max(H, W) / max_side))
    h, w = math.ceil(H / k), math.ceil(W / k)
    q = np.full((h, w), 255, np.uint8)
    for by in range(h):
        for bx in range(w):
            blk = ti[by * k:(by + 1) * k, bx * k:(bx + 1) * k]
            m = valid[by * k:(by + 1) * k, bx * k:(bx + 1) * k]
            if m.sum() * 2 >= m.size and np.isfinite(blk[m]).any():
                t = float(np.nanmean(blk[m]))
                q[by, bx] = int(min(max(round((t - 2) / 20 * 30), 0), 30))
    return q, k * cell, k
