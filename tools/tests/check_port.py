"""Compare tools/topmodel.py with the page's JavaScript model (needs node).

    python tools/tests/check_port.py index.html
"""
import json, os, re, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import topmodel as tm
from sitedata import load_html, decode_pack

html, data = load_html(sys.argv[1] if len(sys.argv) > 1 else 'index.html')
js = re.search(r'<script>\n(// Enhanced TOPMODEL.*?)</script>', html, re.S).group(1)
tmp = tempfile.mkdtemp()
open(os.path.join(tmp, 'tm.js'), 'w').write(js)
worst = 0
for s in data['sites']:
    F, TI, Q = decode_pack(s)
    for label, p in [('defaults', tm.defaults()), ('preset', dict(tm.defaults(), **s['calib'], ssfOn=True, ieOn=True)),
                     ('preset+imp', dict(tm.defaults(), **s['calib'], ssfOn=True, ieOn=True, impOn=True, impPct=25, fConn=0.6))]:
        Fj = {k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in F.items()}
        open(os.path.join(tmp, 'in.json'), 'w').write(json.dumps(dict(F=Fj, TI=TI, p=p)))
        out = subprocess.run(['node', '-e', f"""
const TM=require({json.dumps(os.path.join(tmp,'tm.js'))});const I=require({json.dumps(os.path.join(tmp,'in.json'))});
const F=I.F;for(const k of ['P','Tmax','Tmin','Rs','ea','doy'])F[k]=Float64Array.from(F[k]);
const o=TM.run(F,I.TI,I.p);process.stdout.write(JSON.stringify({{Q:Array.from(o.Qsim),ET:Array.from(o.ET),sat:Array.from(o.satFrac)}}));"""],
                             capture_output=True, text=True, check=True)
        J = json.loads(out.stdout)
        o = tm.run(F, TI, p)
        d = max(np.abs(np.array(J['Q']) - o['Qsim']).max(), np.abs(np.array(J['ET']) - o['ET']).max(), np.abs(np.array(J['sat']) - o['satFrac']).max())
        mj = tm.perf_metrics(Q[365:], np.array(J['Q'])[365:]); mp = tm.perf_metrics(Q[365:], o['Qsim'][365:])
        worst = max(worst, d)
        print(f"{s['id']:9s} {label:10s} max|diff| {d:.2e}  NSE js {mj['NSE']:.4f} py {mp['NSE']:.4f}  calib-in-page NSE {s['fullCal']['NSE']}")
print('worst difference', worst)
sys.exit(0 if worst < 1e-6 else 1)
