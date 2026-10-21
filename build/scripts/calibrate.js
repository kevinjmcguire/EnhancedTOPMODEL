// SCE-UA pre-calibration of one gauge on KGE (calibration = first 60% after a 365-day warm-up).
//   node scripts/calibrate.js <gauge> [seed]      (run from build/)  -> data/presets/cal_<gauge>_<seed>.json
const fs=require('fs'), path=require('path'); const ROOT=path.join(__dirname,'..');
const TM=require(path.join(ROOT,'src','model.js')); const S=require(path.join(ROOT,'src','setup.js'));
const id=process.argv[2]; const LONG=['08202700','09512280'].includes(id); const SEED=+(process.argv[3]||1); const B=JSON.parse(fs.readFileSync(path.join(ROOT,'data','basins.json')))[id];
const TW=JSON.parse(fs.readFileSync(path.join(ROOT,'data','twi.json')))[id]; const n=B.P.length, warm=365, F=S.forcing(B,0,n-1), TI=S.realTI(TW), R=S.demRouting(TW);
const kC=[]; const nCal=Math.round(0.6*(n-1-warm)); for(let t=warm;t<=warm+nCal;t++) if(B.Q[t]!==null) kC.push(t);
const kV=[]; for(let t=warm+nCal+1;t<n;t++) if(B.Q[t]!==null) kV.push(t);
const p0=S.startParams(B); const arid=['08202700','09512280'].includes(id);
if(arid){p0.ieOn=true;}
const keys=S.CAL_KEYS.concat(arid?['Ks']:[]);
const P=keys.map(k=>S.PBYKEY[k]); const bl=P.map(r=>r.isLog?Math.log10(r.min):r.min), bu=P.map(r=>r.isLog?Math.log10(r.max):r.max);
const toP=x=>{const p={...p0}; keys.forEach((k,i)=>p[k]=P[i].isLog?10**x[i]:x[i]); return S.runParams(p);};
const evalK=(q,k)=>TM.objective('KGE',k.map(t=>B.Q[t]),k.map(t=>q[t]));
const x0=keys.map((k,i)=>P[i].isLog?Math.log10(p0[k]):p0[k]);
const q0=TM.core(F,TI,S.runParams(p0),R).Qsim; const start={cal:evalK(q0,kC),val:evalK(q0,kV)};
const t0=Date.now();
const r=TM.sceua(x=>1-evalK(TM.core(F,TI,toP(x),R).Qsim,kC),x0,bl,bu,LONG?{ngs:6,maxn:6000,kstop:10,pcento:0.01,seed:SEED}:{ngs:3,maxn:3000,kstop:5,pcento:0.1,seed:SEED});
const pb=toP(r.bestx); const qb=TM.core(F,TI,pb,R).Qsim; const M=TM.perfMetrics(kV.map(t=>B.Q[t]),kV.map(t=>qb[t]));
const best={}; keys.forEach(k=>best[k]=+pb[k].toPrecision(4)); if(arid) best.ieOn=true;
const res={id,start,cal:1-r.bestf,val:M.KGE,valNSE:M.NSE,nevals:r.nevals,sec:(Date.now()-t0)/1000,msg:r.message,best};
const od=path.join(ROOT,'data','presets'); fs.mkdirSync(od,{recursive:true}); fs.writeFileSync(path.join(od,`cal_${id}_${SEED}.json`),JSON.stringify(res)); console.log(JSON.stringify(res));
