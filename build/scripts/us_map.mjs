// US state outlines (us-atlas, Albers USA, simplified) + gauge points -> data/map.json.  Run from build/scripts: npm install && node us_map.mjs
import fs from 'fs'; import * as topo from 'topojson-client'; import * as simp from 'topojson-simplify'; import {geoPath, geoAlbersUsa, geoIdentity} from 'd3-geo';
let us=JSON.parse(fs.readFileSync('node_modules/us-atlas/states-albers-10m.json'));
us=simp.presimplify(us); us=simp.simplify(us, 1.5);
const path=geoPath(); // already projected
const st=topo.feature(us,us.objects.states).features;
const r=x=>x.replace(/(\d+\.\d)\d+/g,'$1');
const out={states:st.map(f=>({n:f.properties.name,d:r(path(f)||'')})), border:r(path(topo.mesh(us,us.objects.states,(a,b)=>a!==b)))};
const proj=geoAlbersUsa().scale(1300).translate([487.5,305]);
const B=JSON.parse(fs.readFileSync('../data/basins.json')); out.pts={};
for(const [k,b] of Object.entries(B)){const p=proj([b.lon,b.lat]); out.pts[k]=p.map(v=>+v.toFixed(1));}
// Virginia inset: project VA with its own identity fit
const va=st.find(f=>f.properties.name==='Virginia'); out.vaD=r(path(va)); out.vaBounds=path.bounds(va);
fs.writeFileSync('../data/map.json',JSON.stringify(out)); console.log(fs.statSync('../data/map.json').size, out.pts, out.vaBounds);
