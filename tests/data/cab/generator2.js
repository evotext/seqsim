let seed=777;const rnd=()=>{seed=(seed*1103515245+12345)%2147483648;return seed/2147483648};
const W=['A','B','C','D','E','F','G','H'];
// transliteration mode: aligned vectors of codes, 0 = missing
const L=600, vectors={};
for(const w of W){vectors[w]=new Uint32Array(L)}
for(let i=0;i<L;i++){const nr=1+Math.floor(rnd()*4);for(const w of W){vectors[w][i]=rnd()<0.1?0:1+Math.floor(Math.pow(rnd(),2)*nr)}}
const prof=A.buildMsProfiles({vectors},W);const ms={};
for(let i=0;i<W.length;i++)for(let j=i+1;j<W.length;j++){const r=A.pairMsStats(W[i],W[j],{vectors},prof);ms[W[i]+'|'+W[j]]=r}
// apparatus mode with a subset of witnesses selected
const loci=[];
for(let i=0;i<300;i++){const nr=1+Math.floor(rnd()*4);const rows=[];const assign={};
 for(const w of W){if(rnd()<0.1)continue;assign[w]=Math.floor(Math.pow(rnd(),2)*nr)}
 for(let k=0;k<nr;k++){const wits=W.filter(w=>assign[w]===k);if(!wits.length)continue;
  const t=k===0?'lem':(rnd()<0.7?'var':'reading');rows.push(['r'+k,[t],wits])}
 loci.push({r:rows})}
(async()=>{const sub=['A','C','D','F','H'];const m=await A.appPairMatrix(sub,{loci},s=>s,'varreading',null);
 require('fs').writeFileSync('cab_out2.json',JSON.stringify({vectors:Object.fromEntries(W.map(w=>[w,[...vectors[w]]])),ms,loci,sub,app:Object.fromEntries([...m])}))})();
