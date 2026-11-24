// The home page: every receptor, grouped by class, ligand type and family, with a
// search box that matches gene, receptor name, UniProt name, accession or family.

const escape = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

export function renderIndex(index, { list, search, summary }) {
  const haystack = r => [r.gene, r.name, r.entry, r.accession, r.family, r.ligand_type].join(' ').toLowerCase();
  const prepared = index.map(r => ({ r, text: haystack(r) }));
  const CLASS_NAMES = { A: 'Class A (Rhodopsin)', B1: 'Class B1 (Secretin)',
    B2: 'Class B2 (Adhesion)', C: 'Class C (Glutamate)', F: 'Class F (Frizzled)', T2: 'Class T2 (Taste 2)' };
  const families = new Set(index.map(r => r.family)).size;
  const perClass = d3.rollups(index, v => v.length, r => r.class)
    .map(([c, n]) => `${n} Class ${c}`).join(', ');
  summary.textContent = `${index.length} human receptors (${perClass}) in ${families} families. ` +
    'Each page shows a template-based snake plot next to the AlphaFold model.';

  function draw() {
    const words = search.value.toLowerCase().split(/\s+/).filter(Boolean);
    const hits = prepared.filter(p => words.every(w => p.text.includes(w))).map(p => p.r);
    if (!hits.length) { list.innerHTML = `<p class="empty">No receptor matches “${escape(search.value)}”.</p>`; return; }
    const byClass = d3.groups(hits, r => r.class, r => r.ligand_type || 'Other', r => r.family);
    list.innerHTML = byClass.map(([cls, types]) => `
      <h2 class="class-head">${escape(CLASS_NAMES[cls] || 'Class ' + cls)}</h2>
      ${types.map(([type, fams]) => `
      <h3 class="ltype">${escape(type)}</h3>
      ${fams.map(([fam, rs]) => `
        <h4 class="family">${escape(fam)} <span class="count">${rs.length}</span></h4>
        <div class="rec-grid">
          ${rs.map(r => `<a class="rec" href="#${r.entry}">
              <span class="g">${escape(r.gene || r.entry)}</span>
              <span class="p" title="Mean AlphaFold pLDDT">pLDDT ${r.mean_plddt.toFixed(0)}</span>
              <span class="n">${escape(r.name)} · ${escape(r.entry)}</span>
            </a>`).join('')}
        </div>`).join('')}`).join('')}`).join('');
  }
  search.addEventListener('input', draw);
  draw();
}
