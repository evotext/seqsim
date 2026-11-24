let seed=4242;const rnd=()=>{seed=(seed*1103515245+12345)%2147483648;return seed/2147483648};
const W=['A','B','C','D','E','F','G','H'];const T=['var','reading','orth','om'];const loci=[];
for(let i=0;i<300;i++){const nr=1+Math.floor(rnd()*4);const rows=[];const assign={};
 for(const w of W){if(rnd()<0.1)continue;assign[w]=Math.floor(Math.pow(rnd(),2)*nr)}
 for(let k=0;k<nr;k++){const wits=W.filter(w=>assign[w]===k);if(!wits.length)continue;
  const t=k===0?'lem':T[Math.floor(rnd()*T.length)];
  // sometimes split a reading into two rows of different types (merged by CAB)
  if(wits.length>1&&rnd()<0.3){const t2=T[Math.floor(rnd()*T.length)];rows.push(['r'+k,[t],wits.slice(0,1)]);rows.push(['r'+k,[t2],wits.slice(1)])}
  else rows.push(['r'+k,[t],wits])}
 loci.push({r:rows})}
(async()=>{const m=await A.appPairMatrix(W,{loci},s=>s,'all',null);
 require('fs').writeFileSync('cab_out3.json',JSON.stringify({loci,app:Object.fromEntries([...m])}))})();
