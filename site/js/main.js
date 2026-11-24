// Page entry point.
//
// Two views, chosen by the address hash:
//   (no hash)       the index: every receptor, searchable
//   #adrb2_human    one receptor page: snake plot(s) + 3D viewer
//
// View modes on a receptor page
//   single   one receptor (A)
//   side     A and a partner B as two plots
//   overlay  A and B drawn into the same template slots in one plot
// In both compare modes the 3D viewer shows the two models superposed,
// A in blue and B in orange.
import { createSnakePlot, bounds } from './snakeplot.js';
import { createViewer3D } from './viewer3d.js';
import { layoutReceptor } from './layout.js';
import { PLDDT_BANDS, COMPARE_CLASSES, RECEPTOR_COLORS, CHARGE_CLASSES, TYPE_CLASSES, hydroScale } from './colors.js';
import { parseSelection, resolveSelection } from './selection.js';
import { buildSvg, svgToPng, toCsv, saveFile } from './export.js';
import { motifsForClass } from './motifs.js';
import { renderIndex } from './index-page.js';

const hex = c => parseInt(c.slice(1), 16);
const getJSON = async url => {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  return r.json();
};

// ------------------------------------------------------------------- data
const index = await getJSON('data/index.json');
const byEntry = new Map(index.map(r => [r.entry, r]));
const templates = {};
const receptorCache = new Map();

// Receptor data, one JSON file per receptor. A deployment can instead ship
// them in one file per class by setting window.RECEPTOR_BUNDLE to a URL
// with {class} in it, e.g. "data/bundle-{class}.json".
const bundles = {};
async function loadReceptor(entry) {
  if (!receptorCache.has(entry)) {
    receptorCache.set(entry, (async () => {
      if (window.RECEPTOR_BUNDLE) {
        const cls = byEntry.get(entry)?.class;
        bundles[cls] ??= getJSON(window.RECEPTOR_BUNDLE.replace('{class}', cls));
        return (await bundles[cls])[entry];
      }
      return getJSON(`data/receptors/${entry}.json`);
    })());
  }
  return receptorCache.get(entry);
}
async function loadTemplate(cls) {
  templates[cls] ??= await getJSON(`data/templates/class${cls}.json`);
  return templates[cls];
}
const short = rec => rec.gene || rec.entry_name;

// --------------------------------------------------------------- state
const view = { mode: 'single', color: 'segment', label: 'aa', motifs: new Set(), custom: [] };   // toolbar choices
// view.custom: the user's highlight groups, [{id, text, color, tokens}]. They
// are kept as typed and resolved again for each receptor shown.
let A = null, B = null, template = null, frame = null;
let plotA = null, plotB = null, viewer = null;

// Residue of `rec` with the same generic number as residue r (or null).
const counterpart = (r, rec) => (r && r.generic && rec ? rec.residues.find(x => x.generic === r.generic) || null : null);
const compareMode = () => view.mode !== 'single';

// The plot frame: the reference receptor's box, grown to fit A and B. Using
// the same frame for every receptor keeps residues with the same generic
// number on the same spot on screen when you move between receptors.
async function computeFrame() {
  const ref = await loadReceptor(template.reference);
  const recs = [ref, A, ...(compareMode() ? [B] : [])];
  return recs.map(rec => bounds(layoutReceptor(rec, template, new Set()), template.membrane_half_px))
    .reduce((a, b) => ({ x0: Math.min(a.x0, b.x0), x1: Math.max(a.x1, b.x1),
                         y0: Math.min(a.y0, b.y0), y1: Math.max(a.y1, b.y1) }));
}

// ---------------------------------------------------------------- tooltip
const tip = document.getElementById('tooltip');
function tipRow(r, rec, color) {
  const dot = color ? `<span class="dot" style="background:${color}"></span>` : '';
  if (!r) return `<div class="row">${dot}<span class="muted">${short(rec)}: no residue here</span></div>`;
  return `<div class="row">${dot}<b>${r.aa}${r.seq}</b> <span class="muted">${r.display || 'no generic number'}</span></div>
    <div class="muted">${compareMode() ? short(rec) + ' · ' : ''}${r.segment} · pLDDT ${r.plddt ?? '–'}</div>`;
}
function showTip(ra, rb, event) {
  if (!ra && !rb) { tip.hidden = true; return; }
  tip.innerHTML = compareMode() ? tipRow(ra, A, RECEPTOR_COLORS.a) + tipRow(rb, B, RECEPTOR_COLORS.b) : tipRow(ra, A);
  tip.hidden = false;
  const pad = 14;
  let x = event.clientX + pad, y = event.clientY + pad;
  const box = tip.getBoundingClientRect();
  if (x + box.width > innerWidth - 8) x = event.clientX - box.width - pad;
  if (y + box.height > innerHeight - 8) y = event.clientY - box.height - pad;
  tip.style.left = `${x}px`; tip.style.top = `${y}px`;
}

// -------------------------------------------------------------- selection
const selectionLabel = document.getElementById('selection-label');

// Select residue ra of A and/or rb of B everywhere: plots, 3D, label.
function selectPair(ra, rb, { focus3d = true } = {}) {
  if (compareMode()) { ra = ra || counterpart(rb, A); rb = rb || counterpart(ra, B); }
  const fmt = r => (r ? `${r.aa}${r.seq}` : '–');
  selectionLabel.textContent = compareMode()
    ? `${short(A)} ${fmt(ra)} · ${short(B)} ${fmt(rb)} · ${(ra || rb).display || (ra || rb).segment}`
    : ra ? `${fmt(ra)} · ${ra.display || ra.segment}` : '';
  if (view.mode === 'side') { plotA.select(ra?.seq ?? null); plotB.select(rb?.seq ?? null); }
  else plotA.select(ra?.seq ?? null, rb?.seq ?? null);
  if (focus3d) viewer?.focus([{ index: 0, seq: ra?.seq }, { index: 1, seq: rb?.seq }]);
}

function hoverPair(ra, rb, event) {
  if (compareMode() && (ra || rb)) { ra = ra || counterpart(rb, A); rb = rb || counterpart(ra, B); }
  showTip(ra, rb, event);
  viewer?.highlight(ra || rb ? [{ index: 0, seq: ra?.seq }, { index: 1, seq: rb?.seq }] : null);
}

// ------------------------------------------------------------------ render
function renderHeader() {
  document.title = `${A.name} · GPCR Snake Plot`;
  document.getElementById('receptor-name').textContent = A.name;
  document.getElementById('receptor-gene').textContent = A.gene;
  const meta = byEntry.get(A.entry_name);
  document.getElementById('crumb-path').textContent =
    meta ? `· ${meta.ligand_type} › ${meta.family}` : '';
  document.getElementById('receptor-meta').innerHTML = [
    A.class,
    `<a href="https://gpcrdb.org/protein/${A.entry_name}/" target="_blank" rel="noopener">GPCRdb</a>`,
    `<a href="https://www.uniprot.org/uniprotkb/${A.accession}" target="_blank" rel="noopener">UniProt ${A.accession}</a>`,
    `<a href="https://alphafold.ebi.ac.uk/entry/${A.accession}" target="_blank" rel="noopener">AlphaFold ${A.model.id}</a>`,
  ].join(' · ');
  const m = rec => `${rec.model.id} (mean pLDDT ${rec.model.mean_plddt.toFixed(1)})`;
  document.getElementById('model-label').textContent = compareMode()
    ? `AlphaFold ${m(A)} and ${m(B)}, both superposed on ${short(byEntry.get(template.reference) || { entry_name: template.reference })} TM C-alphas`
    : `AlphaFold model ${A.model.id} v${A.model.version} · mean pLDDT ${A.model.mean_plddt.toFixed(1)}`;
  const cap = (rec, c) => `<span class="dot" style="background:${c}"></span>${rec.name} (${short(rec)})`;
  document.getElementById('caption-a').innerHTML = view.mode === 'overlay'
    ? `${cap(A, RECEPTOR_COLORS.a)} &nbsp; ${cap(B, RECEPTOR_COLORS.b)}`
    : cap(A, RECEPTOR_COLORS.a);
  document.getElementById('caption-b').innerHTML = cap(B, RECEPTOR_COLORS.b);
}

async function renderPlots() {
  document.body.className = `mode-${view.mode}`;
  frame = await computeFrame();
  const opts = (rec, partner, overlay) => ({
    frame, color: view.color, label: view.label, partner, overlay, motifs: activeMotifs(), custom: resolvedCustom(),
    onHover(r, event, rOther) {
      if (!r && !rOther) return showTip(null);
      rec === A ? hoverPair(r, rOther, event) : hoverPair(rOther, r, event);
    },
    onClick(r, rOther) { rec === A ? selectPair(r, rOther) : selectPair(rOther, r); },
  });
  plotA = createSnakePlot(document.getElementById('snakeplot'), A, template,
    opts(A, compareMode() ? B : null, view.mode === 'overlay'));
  plotB = view.mode === 'side'
    ? createSnakePlot(document.getElementById('snakeplot-b'), B, template, opts(B, A, false))
    : null;
  if (view.mode !== 'side') d3.select('#snakeplot-b').selectAll('*').remove();
  selectionLabel.textContent = '';
  renderTemplateNote();
  renderHeader();
  renderLegend();
}

async function render3D() {
  if (!viewer) {
    viewer = await createViewer3D(document.getElementById('viewer3d'), [A.model], {
      onPick(index, seq) {
        const rec = index === 1 && compareMode() ? B : A;
        const r = rec.residues.find(x => x.seq === seq);
        rec === A ? selectPair(r, null, { focus3d: false }) : selectPair(null, r, { focus3d: false });
      },
    });
    if (!compareMode()) { viewer.mark(motifTargets()); return; }
  }
  if (compareMode()) await viewer.load([A.model, B.model], [hex(RECEPTOR_COLORS.a), hex(RECEPTOR_COLORS.b)]);
  else await viewer.load([A.model]);
  viewer.mark(motifTargets());
}

// Where the class template's helix positions come from.
function renderTemplateNote() {
  const el = document.getElementById('template-note');
  const n = template.structures?.length;
  if (template.method === 'experimental') {
    el.textContent = `${A.class} template: helix heights and phases averaged over ${n} inactive-state ` +
      `structures placed in the membrane by OPM (one per receptor).`;
    el.title = template.structures.map(s => `${s.pdb} ${s.receptor}`).join(', ');
  } else {
    el.textContent = `${A.class} template: approximate. There are too few inactive-state structures ` +
      `for an average, so it is measured on one AlphaFold model (${template.reference}).`;
    el.title = '';
  }
}

function renderLegend() {
  const el = document.getElementById('legend');
  const sw = (c, b) => `<span><span class="swatch" style="background:${c};${b ? `border-color:${b}` : ''}"></span>`;
  let entries;
  if (view.color === 'plddt') entries = PLDDT_BANDS.map(b => `${sw(b.color)}${b.label}</span>`);
  else if (view.color === 'hydro') entries = [`<span>Hydrophilic −4.5<span class="gradient" style="background:linear-gradient(90deg,${
    [-4.5, -2.25, 0, 2.25, 4.5].map(hydroScale).join(',')})"></span>+4.5 Hydrophobic (Kyte–Doolittle)</span>`];
  else if (view.color === 'charge') entries = CHARGE_CLASSES.map(c => `${sw(c.color)}${c.label}</span>`);
  else if (view.color === 'type') entries = TYPE_CLASSES.map(c => `${sw(c.color)}${c.label}</span>`);
  else if (view.color === 'compare') entries = COMPARE_CLASSES.map(c => `${sw(c.color, c.edge)}${c.label}</span>`);
  else entries = [`${sw('var(--tm)', 'var(--tm-edge)')}TM helix</span>`,
                  `${sw('var(--h8)', 'var(--h8-edge)')}Helix 8</span>`,
                  `${sw('var(--loop)', 'var(--loop-edge)')}Loop / terminus</span>`];
  if (view.mode === 'overlay') entries.push(
    `<span>Split bead = residues differ (<b style="color:${RECEPTOR_COLORS.a}">${short(A)}</b> left, ` +
    `<b style="color:${RECEPTOR_COLORS.b}">${short(B)}</b> right); dashed = ${short(B)} only</span>`);
  activeMotifs().forEach(m => entries.push(
    `<span><span class="ring" style="border-color:${m.color}"></span>${m.name} (${m.generic.join(', ')})</span>`));
  view.custom.forEach(g => entries.push(`<span>${sw(g.color)}${escapeHtml(g.text)}</span>`));
  entries.push('<span>Faded = faces the lipids, solid = faces the helix bundle</span>');
  el.innerHTML = entries.join('');
}

// ------------------------------------------------------------------ motifs
// Toggle buttons, one per motif of the receptor's class. Several can be on.
const activeMotifs = () => motifsForClass(A?.class_template).filter(m => view.motifs.has(m.id));
const motifGroup = document.getElementById('motif-group');
function renderMotifButtons() {
  motifGroup.querySelectorAll('button').forEach(b => b.remove());
  const motifs = motifsForClass(A.class_template);
  motifGroup.hidden = motifs.length === 0;   // no presets for this class yet
  motifs.forEach(m => {
    const b = document.createElement('button');
    b.dataset.motif = m.id;
    b.style.setProperty('--m', m.color);
    b.title = m.generic.join(', ');
    b.innerHTML = `<span class="ring"></span>${m.name}`;
    b.classList.toggle('on', view.motifs.has(m.id));
    b.setAttribute('aria-pressed', String(view.motifs.has(m.id)));
    motifGroup.appendChild(b);
  });
}
// The motif residues of every receptor on screen, as 3D targets.
// (Also the residues of the user's highlight groups.)
function motifTargets() {
  const generics = new Set(activeMotifs().flatMap(m => m.generic));
  const recs = compareMode() ? [A, B] : [A];
  const custom = resolvedCustom();
  return recs.flatMap((rec, index) => rec.residues
    .filter(r => generics.has(r.generic) || custom.some(g => g.seqs.get(rec.entry_name)?.has(r.seq)))
    .map(r => ({ index, seq: r.seq })));
}

// ------------------------------------------------------- custom highlights
const escapeHtml = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const hlForm = document.getElementById('highlight-form');
const hlInput = document.getElementById('highlight-input');
const hlColor = document.getElementById('highlight-color');
const hlStatus = document.getElementById('highlight-status');
const hlChips = document.getElementById('highlight-chips');
const PALETTE = ['#e6375a', '#2a9d8f', '#f4a20c', '#7b4fd6', '#1f8bd6', '#7a8a1f'];
let nextId = 1;

// The groups as residue sets for each receptor on screen.
function resolvedCustom() {
  if (!A) return [];
  return view.custom.map(g => {
    const seqs = new Map([[A.entry_name, resolveSelection(g.tokens, A, true).seqs]]);
    if (compareMode() && B) seqs.set(B.entry_name, resolveSelection(g.tokens, B, false).seqs);
    return { color: g.color, seqs };
  });
}

function renderChips() {
  hlChips.innerHTML = '';
  view.custom.forEach(g => {
    const n = resolveSelection(g.tokens, A, true).seqs.size;
    const chip = document.createElement('span');
    chip.className = 'chip';
    chip.innerHTML = `<span class="swatch" style="background:${g.color}"></span><span></span>` +
      `<span class="muted"></span><button type="button" aria-label="Remove">×</button>`;
    chip.children[1].textContent = g.text;
    chip.children[2].textContent = `${n} residue${n === 1 ? '' : 's'}`;
    chip.querySelector('button').addEventListener('click', () => {
      view.custom = view.custom.filter(x => x.id !== g.id);
      customChanged();
    });
    hlChips.appendChild(chip);
  });
}

function customChanged() {
  const groups = resolvedCustom();
  plotA?.set('custom', groups); plotB?.set('custom', groups);
  renderChips(); renderLegend();
  viewer?.mark(motifTargets());
}

hlForm.addEventListener('submit', e => {
  e.preventDefault();
  const text = hlInput.value.trim();
  if (!text) return;
  const tokens = parseSelection(text);
  const bad = tokens.filter(t => t.kind === 'bad').map(t => t.raw);
  const good = tokens.filter(t => t.kind !== 'bad');
  if (!good.length) { hlStatus.textContent = `Not understood: ${bad.join(', ')}`; return; }
  const { seqs, missing } = resolveSelection(good, A, true);
  view.custom.push({ id: nextId++, text: good.map(t => t.raw).join(', '), color: hlColor.value, tokens: good });
  const notes = [];
  if (bad.length) notes.push(`not understood: ${bad.join(', ')}`);
  if (missing.length) notes.push(`not in ${short(A)}: ${missing.join(', ')}`);
  hlStatus.textContent = `${seqs.size} residue${seqs.size === 1 ? '' : 's'} highlighted` + (notes.length ? ` (${notes.join('; ')})` : '');
  hlInput.value = '';
  hlColor.value = PALETTE[view.custom.length % PALETTE.length];   // next group gets a new colour
  customChanged();
});
document.getElementById('highlight-clear').addEventListener('click', () => {
  view.custom = []; hlStatus.textContent = '';
  customChanged();
});

// ------------------------------------------------------------------ export
const COLOR_NAMES = { segment: 'segment', plddt: 'pLDDT', compare: 'compare', hydro: 'hydrophobicity (Kyte–Doolittle)', charge: 'charge', type: 'residue type' };
function exportName(ext) {
  const parts = [short(A)];
  if (compareMode()) parts.push(view.mode === 'overlay' ? 'overlay' : 'vs', short(B));
  return `${parts.join('_')}_snakeplot.${ext}`.replace(/[^\w.+-]+/g, '_');
}
function currentSvg() {
  const els = [document.getElementById('snakeplot')];
  if (view.mode === 'side') els.push(document.getElementById('snakeplot-b'));
  const title = compareMode()
    ? `${A.name} (${short(A)}) ${view.mode === 'overlay' ? 'overlaid with' : 'and'} ${B.name} (${short(B)})`
    : `${A.name} (${short(A)})`;
  const subtitle = [A.class, `colour: ${COLOR_NAMES[view.color]}`,
    ...view.custom.map(g => `highlight ${g.text}`),
    ...activeMotifs().map(m => `motif ${m.name}`)].join(' · ');
  return buildSvg(els, { title, subtitle });
}
async function runExport(btn, job) {
  const label = btn.textContent;
  btn.disabled = true; btn.textContent = '…';
  try {
    const result = await job();
    btn.textContent = result === 'saved' ? '✓' : result === 'declined' ? label : '!';
    if (result !== 'saved' && result !== 'declined') btn.title = String(result);
  } catch (e) { btn.textContent = '!'; btn.title = e.message; }
  setTimeout(() => { btn.textContent = label; btn.disabled = false; }, 1200);
}
document.getElementById('export-svg').addEventListener('click', e =>
  runExport(e.currentTarget, () => saveFile(exportName('svg'), currentSvg().text)));
document.getElementById('export-png').addEventListener('click', e =>
  runExport(e.currentTarget, async () => saveFile(exportName('png'), await svgToPng(currentSvg(), 3))));
document.getElementById('export-list').addEventListener('click', e => runExport(e.currentTarget, () => {
  // One row per residue: each highlight group and each active motif, for every receptor on screen.
  const recs = compareMode() ? [A, B] : [A];
  const rows = [['group', 'colour', 'receptor', 'residue', 'number', 'generic_number', 'segment', 'plddt']];
  const add = (group, color, rec, set) => rec.residues.filter(r => set(r)).forEach(r =>
    rows.push([group, color, short(rec), `${r.aa}${r.seq}`, r.seq, r.display || '', r.segment, r.plddt ?? '']));
  view.custom.forEach(g => recs.forEach((rec, i) => {
    const seqs = resolveSelection(g.tokens, rec, i === 0).seqs;
    add(g.text, g.color, rec, r => seqs.has(r.seq));
  }));
  activeMotifs().forEach(m => recs.forEach(rec => add(m.name, m.color, rec, r => m.generic.includes(r.generic))));
  if (rows.length === 1) return 'nothing to list: add a highlight or turn on a motif first';
  return saveFile(exportName('csv').replace('_snakeplot', '_selection'), toCsv(rows));
}));
motifGroup.addEventListener('click', e => {
  const btn = e.target.closest('button');
  if (!btn) return;
  const id = btn.dataset.motif;
  view.motifs.has(id) ? view.motifs.delete(id) : view.motifs.add(id);
  btn.classList.toggle('on', view.motifs.has(id));
  btn.setAttribute('aria-pressed', String(view.motifs.has(id)));
  plotA.set('motifs', activeMotifs()); plotB?.set('motifs', activeMotifs());
  renderLegend();
  viewer?.mark(motifTargets());
});

// ------------------------------------------------- partner (compare with)
// A dropdown of every receptor of the same class, grouped by family.
// Default partner: the next receptor in A's family, else the class reference.
const compareSelect = document.getElementById('compare-with');
function fillPartnerSelect() {
  const same = index.filter(r => r.class === A.class_template && r.entry !== A.entry_name);
  const groups = d3.groups(same, r => r.family);
  compareSelect.innerHTML = groups.map(([fam, rs]) =>
    `<optgroup label="${fam}">${rs.map(r => `<option value="${r.entry}">${r.gene} · ${r.name}</option>`).join('')}</optgroup>`).join('');
}
async function choosePartner(preferred) {
  const meta = byEntry.get(A.entry_name);
  const fallback = index.find(r => r.family === meta?.family && r.entry !== A.entry_name)?.entry
    || (template.reference !== A.entry_name ? template.reference : index.find(r => r.class === meta?.class && r.entry !== A.entry_name).entry);
  const entry = preferred && preferred !== A.entry_name && byEntry.has(preferred) ? preferred : fallback;
  compareSelect.value = entry;
  B = await loadReceptor(entry);
}
compareSelect.addEventListener('change', async () => {
  B = await loadReceptor(compareSelect.value);
  await renderPlots(); await render3D();
});

// ------------------------------------------------------------------ routing
async function showReceptor(entry) {
  document.documentElement.dataset.view = 'receptor';
  A = await loadReceptor(entry);
  template = await loadTemplate(A.class_template);
  fillPartnerSelect();
  await choosePartner(null);   // new receptor: start from a family neighbour
  renderMotifButtons();
  await renderPlots();
  renderChips();
  viewer?.resize();   // the viewer was hidden while the index was showing
  await render3D();
  scrollTo(0, 0);
}

function showIndex() {
  document.documentElement.dataset.view = 'index';
  document.title = 'GPCR Snake Plot';
  tip.hidden = true;
}

async function route() {
  const entry = location.hash.slice(1);
  if (entry && byEntry.has(entry)) await showReceptor(entry);
  else showIndex();
}
addEventListener('hashchange', route);

// Toolbar buttons
const setOn = (group, value) =>
  group.querySelectorAll('button').forEach(b => b.classList.toggle('on', b.dataset.value === value));
document.querySelectorAll('.toolbar .group[data-control]').forEach(group => {
  group.addEventListener('click', async e => {
    const btn = e.target.closest('button');
    if (!btn) return;
    const control = group.dataset.control, value = btn.dataset.value;
    if (view[control] === value) return;
    const was = compareMode();
    view[control] = value;
    setOn(group, value);
    if (control === 'mode') {
      // Compare colouring only makes sense with two receptors; switch it on/off with the mode.
      if (compareMode() && !was) view.color = 'compare';
      if (!compareMode() && view.color === 'compare') view.color = 'segment';
      setOn(document.querySelector('.group[data-control="color"]'), view.color);
      await renderPlots();
      if (compareMode() !== was) await render3D();
    } else {
      plotA.set(control, value); plotB?.set(control, value);
      if (control === 'color') renderLegend();
    }
  });
});

renderIndex(index, {
  list: document.getElementById('index-list'),
  search: document.getElementById('search'),
  summary: document.getElementById('index-summary'),
});
await route();
