// Synthetic apparatus: random loci, types lem/var/reading, some missing witnesses
let seed=12345;const rnd=()=>{seed=(seed*1103515245+12345)%2147483648;return seed/2147483648};
const W=['A','B','C','D','E','F','G','H'];const loci=[];
for(let i=0;i<400;i++){const nr=1+Math.floor(rnd()*4);const rows=[];const assign={};
 for(const w of W){if(rnd()<0.1)continue;assign[w]=Math.floor(Math.pow(rnd(),2)*nr)}
 for(let k=0;k<nr;k++){const wits=W.filter(w=>assign[w]===k);if(!wits.length)continue;
  const t=k===0?'lem':(rnd()<0.7?'var':'reading');rows.push(['r'+k,[t],wits])}
 loci.push({r:rows})}
(async()=>{const out={loci,results:{}};
 for(const preset of ['var','varreading','all']){const m=await A.appPairMatrix(W,{loci},s=>s,preset,null);
  out.results[preset]=Object.fromEntries([...m].map(([k,v])=>[k,v]))}
 const tests=[];for(let t=0;t<50;t++){const qs=[];for(let j=0;j<1+Math.floor(rnd()*6);j++)qs.push([rnd(),rnd()*20]);const mm=rnd()<0.3?rnd()*3000:rnd()*30;tests.push({qs,mm,post:A.bayesPosterior(qs,mm)})}
 out.tests=tests;require('fs').writeFileSync('cab_out.json',JSON.stringify(out))})();
