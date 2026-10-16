"""Read and write the page's embedded DATA object (the `const DATA=...;` line in index.html)."""
import base64
import json
import re
import zlib
import numpy as np

N = 1825
SCALES = [10, 10, 10, 10, 100, 100]  # P, Tmax, Tmin, Rs, ea, Q
DATA_RE = re.compile(r'^const DATA=(\{.*\});\s*$', re.M)


def load_html(path):
    html = open(path, encoding='utf-8').read()
    m = DATA_RE.search(html)
    if not m:
        raise SystemExit('Could not find the `const DATA=...;` line in ' + path)
    return html, json.loads(m.group(1))


def save_html(path, html, data):
    line = 'const DATA=' + json.dumps(data, separators=(',', ':'), ensure_ascii=False) + ';'
    html = DATA_RE.sub(lambda _: line, html, count=1)
    open(path, 'w', encoding='utf-8').write(html)


def inflate(b64):
    return zlib.decompress(base64.b64decode(b64), -15)


def deflate_b64(raw):
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    return base64.b64encode(c.compress(raw) + c.flush()).decode()


def decode_pack(site):
    """Same as the page's decode(): returns forcing dict and Q (mm/d, NaN = missing)."""
    b = np.frombuffer(inflate(site['pack']), np.uint8)
    T = 6 * N
    v = (b[:T].astype(np.int32) | (b[T:2 * T].astype(np.int32) << 8))
    v = np.where(v > 32767, v - 65536, v).astype(np.int64)
    S = []
    for k in range(6):
        a = v[k * N:(k + 1) * N].astype(float)
        if k > 0:
            a = np.cumsum(a)
        S.append(a)
    doy = np.arange(N) % 365 + 1
    F = dict(P=S[0] / 10, Tmax=S[1] / 10, Tmin=S[2] / 10, Rs=S[3] / 10, ea=S[4] / 100, doy=doy,
             lat=site['lat'], elev=site['zmean'])
    Q = np.where(S[5] < 0, np.nan, S[5] / 100)
    TI = dict(val=site['val'], frac=site['frac'], lnLambdaLin=site['lnL'])
    return F, TI, Q


def encode_pack(P, Tmax, Tmin, Rs, ea, Q):
    """Inverse of decode_pack. Q may contain NaN (stored as -1 after rounding)."""
    cols = []
    for k, (a, sc) in enumerate(zip([P, Tmax, Tmin, Rs, ea, Q], SCALES)):
        a = np.asarray(a, float)
        if len(a) != N:
            raise ValueError(f'series {k} has {len(a)} values, expected {N}')
        if k == 5:
            q = np.where(np.isfinite(a), np.round(a * sc), -1).astype(np.int64)
        else:
            if not np.isfinite(a).all():
                raise ValueError(f'forcing series {k} has gaps')
            q = np.round(a * sc).astype(np.int64)
        if q.max() > 32767 or q.min() < -32768:
            raise ValueError(f'series {k} out of int16 range after scaling')
        if k > 0:
            q = np.diff(q, prepend=0)
            if q.max() > 32767 or q.min() < -32768:
                raise ValueError(f'series {k} steps out of int16 range')
        cols.append(q)
    v = np.concatenate(cols) & 0xFFFF
    raw = np.concatenate([(v & 0xFF).astype(np.uint8), (v >> 8).astype(np.uint8)]).tobytes()
    return deflate_b64(raw)
