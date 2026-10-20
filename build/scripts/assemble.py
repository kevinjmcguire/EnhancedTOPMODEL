"""Build the single-file app.

    python scripts/assemble.py        (run from build/)

Inlines src/model.js, src/setup.js and the data/*.json files into src/app.html, keeps the best
pre-calibration per gauge (highest calibration KGE over seeds), makes the file pure ASCII, and writes
  out/topmodel_explorer.html   page body (for the claude.ai artifact, which adds its own <head>)
  ../index.html                full HTML document for GitHub Pages
"""
import json,glob,re,os
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P=lambda *a: os.path.join(ROOT,*a)
pres={}
for f in glob.glob(P('data','presets','cal_*.json')):
    r=json.load(open(f)); i=r['id']
    if i not in pres or r['cal']>pres[i]['cal']: pres[i]={'cal':r['cal'],'val':r['val'],'best':r['best']}
s=open(P('src','app.html'),encoding='utf-8').read()
s=s.replace('/*MODEL*/',open(P('src','model.js'),encoding='utf-8').read()).replace('/*SETUP*/',open(P('src','setup.js'),encoding='utf-8').read())
s=s.replace('/*DATA*/',open(P('data','basins.json')).read()).replace('/*MAP*/',open(P('data','map.json')).read()).replace('/*PRESETS*/',json.dumps(pres)).replace('/*TWI*/',open(P('data','twi.json')).read())
print('presets', {k:(round(v['cal'],3),round(v['val'],3)) for k,v in sorted(pres.items())})

# --- make the file encoding-proof: ASCII only, plus a charset declaration
import re
def js_esc(m):
    body=m.group(2)
    body=''.join(c if ord(c)<128 else ('\\u%04x'%ord(c) if ord(c)<0x10000 else ''.join('\\u%04x'%u for u in c.encode('utf-16-be') and [int.from_bytes(c.encode('utf-16-be')[i:i+2],'big') for i in (0,2)])) for c in body)
    return m.group(1)+body+m.group(3)
s=re.sub(r'(<script[^>]*>)(.*?)(</script>)',js_esc,s,flags=re.S)
s=re.sub(r'(<style>)(.*?)(</style>)',lambda m:m.group(1)+''.join(c if ord(c)<128 else '\\%06x'%ord(c) for c in m.group(2))+m.group(3),s,flags=re.S)
s=''.join(c if ord(c)<128 else '&#x%x;'%ord(c) for c in s)
os.makedirs(P('out'),exist_ok=True)
open(P('out','topmodel_explorer.html'),'w',encoding='ascii').write('<meta charset="utf-8">\n'+s)
i=s.index('</style>')+len('</style>')
doc=('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
     '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
     + s[:i] + '\n<style>html,body{margin:0} :root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}</style>\n</head>\n<body>'
     + s[i:] + '\n</body>\n</html>\n')
open(os.path.join(os.path.dirname(ROOT),'index.html'),'w',encoding='ascii').write(doc)
print('wrote out/topmodel_explorer.html and ../index.html (%.1f MB)' % (len(doc)/1e6))
print('non-ascii left:',sum(1 for c in s if ord(c)>127))
