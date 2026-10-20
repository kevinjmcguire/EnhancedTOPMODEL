// Writes tests/out/in.json (inputs) and tests/out/js.json (JS results) for compare_with_R.R.
//   node tests/compare_with_R.js      (run from build/)
const fs=require('fs'), path=require('path'); const ROOT=path.join(__dirname,'..'); const TM=require(path.join(ROOT,'src','model.js'));
const OUT=path.join(__dirname,'out'); fs.mkdirSync(OUT,{recursive:true});
const B=JSON.parse(fs.readFileSync(path.join(ROOT,'data','basins.json')));
function tiDist(lam, k=4, th=0.75, nc=30){ // shifted gamma, mean lam
  const t0=lam-k*th, hi=k*th+6*Math.sqrt(k)*th, N=3000; const pdf=[];
  for(let i=0;i<N;i++){const x=(i+0.5)*hi/N; pdf.push(Math.pow(x,k-1)*Math.exp(-x/th));}
  const val=[],frac=[]; const per=N/nc;
  for(let c=0;c<nc;c++){let s=0,sx=0;for(let i=c*per;i<(c+1)*per;i++){const x=(i+0.5)*hi/N;s+=pdf[i];sx+=pdf[i]*x;} if(s>0){val.push(t0+sx/s);frac.push(s);} }
  const S=frac.reduce((a,b)=>a+b,0); return {val, frac:frac.map(f=>f/S)};
}
function forcing(b){ const n=b.P.length, d0=new Date(b.start+'T00:00:00Z'); const doy=[];
  for(let t=0;t<n;t++){const d=new Date(d0.getTime()+t*864e5); const y0=Date.UTC(d.getUTCFullYear(),0,1); doy.push(Math.round((d-y0)/864e5)+1);}
  return {P:b.P,Tmax:b.Tx,Tmin:b.Tn,Rs:b.Rs,ea:b.ea,doy,lat:b.latF,elev:b.elev}; }
const base={m:0.02,lnTe:5,Srmax:0.08,td:1,SR0:0.02,Q0:1,Mssf:0.05,kssf:0.5,kv:0.1,Ks:100,ieMoist:1,stormSum:4,stormWin:10,ieConn:1,LAImax:5,LAImin:1,leafOn:130,leafOff:275,ramp:20,rl:250,hc:20,u:2,albedo:0.15,albLeafOn:0.17,albLeafOff:0.13,cint:0.3,Tsnow:0,Tmelt:0,DDF:3,albSnow:0.75,kAlb:0.75,sweFull:20,fEv:0,cSnow:3.5,cSnowBr:1,unload:0.5,chv:0.7,impPct:0,fConn:1,Sdep:1.5,ssfOn:false,ieOn:false,impOn:false,albMode:0,canopyPervOnly:false};
const b=B['01013500']; const F=forcing(b); const TI=tiDist(7.5);
const area_ha=b.area_km2*100, Lmax=1.4*Math.pow(area_ha/259,0.6)*1609.34; const dist=[],frac=[];for(let i=0;i<100;i++){dist.push(i*Lmax/99);frac.push(0.01);}
const R={on:true,dist,frac};
const cases={
 base:{p:base,R:{on:false}},
 all:{p:{...base,ssfOn:true,ieOn:true,impOn:true,impPct:8,fConn:0.6,albMode:1,Ks:15},R},
 canperv:{p:{...base,impOn:true,impPct:12,canopyPervOnly:true,ieOn:true,Ks:10},R:{on:false}},
 harv:{p:{...base,ssfOn:true,ieOn:true,Ks:30,luc:{on:true,change:true,h:TI.val.map((v,i)=>i<10?1:0.3),tIdx:900,
   pre:[0,1,null,null,null,null,null,null,null,null,1,null,0,null],start:[0,0.2,null,null,null,null,null,1,0.19,5,0.5,0,0.1,5],
   final:[0,1,null,null,null,null,null,null,null,null,1,null,0,null],curve:'exponential',Tdays:365.25*8,compRecDays:365.25*5,
   impFrac:0.2,impRampDays:365,impConn:0.6,impSdep:1.5,detention:true,kdet:0.3}},R},
 harvImp:{p:{...base,impOn:true,impPct:5,ieOn:true,luc:{on:true,change:true,h:TI.val.map(()=>0.4),tIdx:400,
   pre:[1,null,2.5,0.6,100,290,30,0.4,0.2,4.5,0.5,0,0.1,15],start:[1,null,2.5,0.6,110,290,25,4,0.18,4,0.4,0,0.3,10],
   final:[1,null,2.5,0.6,110,290,25,4,0.18,4,0.4,0,0.3,10],curve:'step',Tdays:365,compRecDays:Infinity,
   impFrac:0.35,impRampDays:1826,impConn:0.6,impSdep:1.5,detention:false,kdet:0.3,laiMult:F.P.map((v,i)=>(i>500&&i<600)?0.3:1)}},R:{on:false}},
};
const res={};
for(const [k,c] of Object.entries(cases)){ const t0=Date.now(); const o=TM.core(F,TI,c.p,c.R); res[k]={ms:Date.now()-t0}; for(const f of ['Qsim','ET','SWE','Sbar','satFrac','Qie','Qimp','Qssf','storage','Qcomp','Qdev','PET','Ea'])res[k][f]=Array.from(o[f]); }
fs.writeFileSync(path.join(OUT,'in.json'),JSON.stringify({F,TI,cases:Object.fromEntries(Object.entries(cases).map(([k,c])=>[k,{p:{...c.p,luc:c.p.luc?{...c.p.luc,compRecDays:Number.isFinite(c.p.luc.compRecDays)?c.p.luc.compRecDays:'Inf'}:undefined},R:c.R}]))}));
fs.writeFileSync(path.join(OUT,'js.json'),JSON.stringify(res));
for(const k in res) console.log(k,res[k].ms+' ms');
