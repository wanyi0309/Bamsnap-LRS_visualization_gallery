(() => {
const DATA=window.GALLERY_DATA;
const order=['INS','DEL','INV','DUP','RNA','PHASE'];
const toolDefs=[['bamsnap','Bamsnap-LRS'],['wally','Wally'],['svhawkeye','SVhawkeye']];
const qs=new URLSearchParams(location.search);
let type=DATA[qs.get('type')]?qs.get('type'):'INS';
let idx=0; let stacked=false;
const $=s=>document.querySelector(s);
const tabs=$('#tabs'), select=$('#locusSelect'), comp=$('#comparison');
function summary(){const total=order.reduce((n,k)=>n+DATA[k].records.length,0);$('#summary').innerHTML=`<span>${total} matched loci</span><span>${order.length} visualization classes</span><span>${toolDefs.length} tools represented</span><span>GitHub Pages ready</span>`}
function buildTabs(){tabs.innerHTML='';order.forEach(k=>{let b=document.createElement('button');b.className='tab'+(k===type?' active':'');b.textContent=`${DATA[k].title} · ${DATA[k].records.length}`;b.onclick=()=>{type=k;idx=0;renderAll(true)};tabs.appendChild(b)})}
function fillSelect(){select.innerHTML='';DATA[type].records.forEach((r,i)=>{const o=document.createElement('option');o.value=i;o.textContent=`${r.id}  ·  ${r.region}`;select.appendChild(o)});select.value=idx}
function esc(s){return String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}
function renderCards(){const r=DATA[type].records[idx];comp.innerHTML='';const visibleTools=toolDefs.filter(([key])=>DATA[type].records.some(x=>x.images&&x.images[key]));comp.classList.toggle('two-tools',visibleTools.length===2);visibleTools.forEach(([key,label])=>{const path=r.images[key], source=r.source_files[key];const card=document.createElement('section');card.className='card';let body=path?`<div class="image-wrap"><img loading="eager" decoding="async" src="${esc(path)}" alt="${esc(label)} visualization for ${esc(r.id)}" data-tool="${esc(label)}"></div>`:`<div class="image-wrap"><div class="missing-box"><div><strong>Not available</strong><br>No corresponding ${esc(label)} image was present for ${esc(r.id)}.</div></div></div>`;card.innerHTML=`<div class="card-head"><span class="tool-name">${esc(label)}</span><span class="status ${path?'':'missing'}">${path?'Available':'Not available'}</span></div>${body}<div class="card-foot">${path?`<a href="${esc(path)}" target="_blank" rel="noopener">Open original</a> · ${esc(source)}`:'No source image'}</div>`;comp.appendChild(card)});comp.querySelectorAll('img').forEach(im=>im.onclick=()=>openModal(im.src,`${im.dataset.tool} · ${r.id} · ${r.region}`))}
function updateMeta(){const r=DATA[type].records[idx];$('#locusTitle').textContent=`${r.id} · ${DATA[type].title}`;$('#locusMeta').textContent=`${r.region}  ·  displayed interval ${r.span.toLocaleString()} bp`;$('#prevBtn').disabled=idx===0;$('#nextBtn').disabled=idx===DATA[type].records.length-1;select.value=idx;const u=new URL(location.href);u.search='';u.searchParams.set('type',type);u.searchParams.set('id',r.id);history.replaceState(null,'',u)}
function renderAll(resetTabs=false){if(resetTabs)buildTabs();fillSelect();renderCards();updateMeta();$('#searchInput').value='';$('#searchNote').hidden=true}
select.onchange=()=>{idx=Number(select.value);renderCards();updateMeta()};
$('#prevBtn').onclick=()=>{if(idx>0){idx--;renderCards();updateMeta()}};$('#nextBtn').onclick=()=>{if(idx<DATA[type].records.length-1){idx++;renderCards();updateMeta()}};
$('#layoutBtn').onclick=()=>{stacked=!stacked;comp.classList.toggle('stacked',stacked);$('#layoutBtn').textContent=stacked?'Use columns':'Stack vertically'};
$('#copyBtn').onclick=async()=>{try{await navigator.clipboard.writeText(location.href);$('#copyBtn').textContent='Copied';setTimeout(()=>$('#copyBtn').textContent='Copy direct link',1200)}catch{prompt('Copy this link:',location.href)}};
$('#searchInput').addEventListener('input',e=>{const q=e.target.value.trim().toLowerCase();const note=$('#searchNote');if(!q){note.hidden=true;return}const records=DATA[type].records;const hit=records.findIndex(r=>r.id.toLowerCase()===q||r.region.toLowerCase().replaceAll(',','').includes(q.replaceAll(',',''))||`${r.chrom}:${r.start}-${r.end}`.toLowerCase().includes(q.replaceAll(',','')));if(hit>=0){idx=hit;renderCards();updateMeta();note.hidden=true}else{note.textContent='No matching locus in this category.';note.hidden=false}});
const modal=$('#modal'), modalBody=$('#modalBody'), modalImg=$('#modalImg'), zoomLabel=$('#zoomLabel');
let zoom=1, fitZoom=1, dragging=false, dragX=0, dragY=0, dragLeft=0, dragTop=0;
const clamp=(v,min,max)=>Math.min(max,Math.max(min,v));
function calcFitZoom(){
  if(!modalImg.naturalWidth)return 1;
  const usableWidth=Math.max(160,modalBody.clientWidth-36);
  return Math.min(1,usableWidth/modalImg.naturalWidth);
}
function applyZoom(next,anchorX=null,anchorY=null){
  if(!modalImg.naturalWidth)return;
  next=clamp(next,0.1,5);
  const rect=modalBody.getBoundingClientRect();
  const px=anchorX==null?rect.left+modalBody.clientWidth/2:anchorX;
  const py=anchorY==null?rect.top+modalBody.clientHeight/2:anchorY;
  const ox=modalBody.scrollLeft+(px-rect.left);
  const oy=modalBody.scrollTop+(py-rect.top);
  const ratio=next/zoom;
  zoom=next;
  modalImg.style.width=`${Math.round(modalImg.naturalWidth*zoom)}px`;
  modalImg.style.height='auto';
  zoomLabel.textContent=`${Math.round(zoom*100)}%`;
  requestAnimationFrame(()=>{
    modalBody.scrollLeft=ox*ratio-(px-rect.left);
    modalBody.scrollTop=oy*ratio-(py-rect.top);
  });
}
function fitImage(){fitZoom=calcFitZoom();zoom=fitZoom;modalImg.style.width=`${Math.round(modalImg.naturalWidth*zoom)}px`;modalImg.style.height='auto';zoomLabel.textContent=`${Math.round(zoom*100)}%`;modalBody.scrollLeft=0;modalBody.scrollTop=0}
function openModal(src,title){
  modalImg.onload=()=>{fitImage()};
  modalImg.src=src;$('#modalTitle').textContent=title;modal.classList.add('open');document.body.style.overflow='hidden';
}
function closeModal(){modal.classList.remove('open');modalImg.onload=null;modalImg.src='';modalImg.style.width='';document.body.style.overflow='';dragging=false;modalBody.classList.remove('dragging')}
$('#zoomIn').onclick=()=>applyZoom(zoom*1.25);
$('#zoomOut').onclick=()=>applyZoom(zoom/1.25);
$('#zoomFit').onclick=fitImage;
$('#zoomActual').onclick=()=>applyZoom(1);
$('#closeModal').onclick=closeModal;
modalBody.addEventListener('wheel',e=>{if(!modal.classList.contains('open'))return;e.preventDefault();applyZoom(zoom*(e.deltaY<0?1.15:1/1.15),e.clientX,e.clientY)},{passive:false});
modalImg.addEventListener('dblclick',e=>{e.preventDefault();const nearActual=Math.abs(zoom-1)<0.05;applyZoom(nearActual?fitZoom:1,e.clientX,e.clientY)});
modalBody.addEventListener('mousedown',e=>{if(e.button!==0)return;dragging=true;dragX=e.clientX;dragY=e.clientY;dragLeft=modalBody.scrollLeft;dragTop=modalBody.scrollTop;modalBody.classList.add('dragging');e.preventDefault()});
window.addEventListener('mousemove',e=>{if(!dragging)return;modalBody.scrollLeft=dragLeft-(e.clientX-dragX);modalBody.scrollTop=dragTop-(e.clientY-dragY)});
window.addEventListener('mouseup',()=>{dragging=false;modalBody.classList.remove('dragging')});
window.addEventListener('resize',()=>{if(modal.classList.contains('open'))fitZoom=calcFitZoom()});
document.addEventListener('keydown',e=>{if($('#modal').classList.contains('open')){if(e.key==='Escape')closeModal();return}if(e.target.matches('input,select'))return;if(e.key==='ArrowLeft')$('#prevBtn').click();if(e.key==='ArrowRight')$('#nextBtn').click()});
const requestedId=qs.get('id');if(requestedId){const n=DATA[type].records.findIndex(r=>r.id===requestedId);if(n>=0)idx=n}
summary();buildTabs();renderAll(false);
})();