// Export the snake plot as SVG or PNG, and lists as CSV.
//
// The page styles the plot with CSS (classes and colour variables), which a
// saved SVG file would lose. So the export copies each element's computed
// style onto the copy as inline attributes first.

const STYLE_PROPS = [
  'fill', 'fill-opacity', 'stroke', 'stroke-width', 'stroke-opacity', 'stroke-dasharray',
  'opacity', 'font-family', 'font-size', 'font-weight', 'text-anchor', 'dominant-baseline',
  'letter-spacing', 'text-transform',
];

function inlineCopy(svgEl) {
  const copy = svgEl.cloneNode(true);
  copy.querySelectorAll('.hover').forEach(n => n.classList.remove('hover'));
  const src = [svgEl, ...svgEl.querySelectorAll('*')];
  const dst = [copy, ...copy.querySelectorAll('*')];
  src.forEach((el, i) => {
    const cs = getComputedStyle(el);
    const style = STYLE_PROPS.map(p => `${p}:${cs.getPropertyValue(p)}`).join(';');
    dst[i].setAttribute('style', style);
    dst[i].removeAttribute('class');
  });
  // text-transform is not an SVG property: apply it to the text itself.
  dst.forEach((el, i) => {
    if (el.tagName === 'text' && getComputedStyle(src[i]).textTransform === 'uppercase') el.textContent = el.textContent.toUpperCase();
  });
  return copy;
}

const esc = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

// One SVG document with the given plots side by side, a title line on top
// and a credit line at the bottom. Returns {text, width, height}.
export function buildSvg(svgEls, { title, subtitle }) {
  const TOP = 46, BOTTOM = 22, GAP = 20;
  const parts = svgEls.map(el => {
    const [x, y, w, h] = el.getAttribute('viewBox').split(/\s+/).map(Number);
    return { el, x, y, w, h };
  });
  const width = parts.reduce((s, p) => s + p.w, 0) + GAP * (parts.length - 1);
  const height = Math.max(...parts.map(p => p.h)) + TOP + BOTTOM;
  let offset = 0;
  const inner = parts.map(p => {
    const copy = inlineCopy(p.el);
    copy.setAttribute('x', offset); copy.setAttribute('y', TOP);
    copy.setAttribute('width', p.w); copy.setAttribute('height', p.h);
    copy.removeAttribute('id'); copy.removeAttribute('role');
    offset += p.w + GAP;
    return new XMLSerializer().serializeToString(copy);
  }).join('\n');
  const font = 'font-family:Inter,system-ui,-apple-system,Segoe UI,Roboto,sans-serif';
  const text = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">
<rect width="100%" height="100%" fill="#ffffff"/>
<text x="16" y="22" style="${font};font-size:16px;font-weight:600;fill:#1c2330">${esc(title)}</text>
<text x="16" y="38" style="${font};font-size:11px;fill:#5b6576">${esc(subtitle)}</text>
${inner}
<text x="16" y="${height - 8}" style="${font};font-size:9px;fill:#8a93a3">Generic numbers and segments: GPCRdb. Structure and pLDDT: AlphaFold DB (CC BY 4.0).</text>
</svg>`;
  return { text, width, height };
}

// Render an SVG string to a PNG Blob at `scale` times its size.
export function svgToPng({ text, width, height }, scale = 3) {
  return new Promise((resolve, reject) => {
    const img = new Image();
    const url = URL.createObjectURL(new Blob([text], { type: 'image/svg+xml' }));
    img.onload = () => {
      const canvas = document.createElement('canvas');
      canvas.width = Math.round(width * scale); canvas.height = Math.round(height * scale);
      const ctx = canvas.getContext('2d');
      ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
      URL.revokeObjectURL(url);
      canvas.toBlob(b => (b ? resolve(b) : reject(new Error('PNG encoding failed'))), 'image/png');
    };
    img.onerror = () => { URL.revokeObjectURL(url); reject(new Error('Could not render the SVG')); };
    img.src = url;
  });
}

export function toCsv(rows) {
  const cell = v => (v == null ? '' : /[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : String(v));
  return rows.map(r => r.map(cell).join(',')).join('\n') + '\n';
}

// Save a file. Inside a claude.ai artifact the page asks the viewer through
// the "downloads" capability; anywhere else (local server, GitHub Pages) it
// uses a normal download link.
let downloadsCap = null;
const getCap = () => {
  if (!downloadsCap && window.claude?.use) downloadsCap = window.claude.use('downloads').catch(() => null);
  return downloadsCap || Promise.resolve(null);
};
getCap();   // start asking early; the answer can take a moment
export async function saveFile(filename, data) {
  const cap = await getCap();
  if (cap) {
    try { await cap.save({ filename, data }); return 'saved'; }
    catch (e) { return e?.code === 'declined' ? 'declined' : `failed: ${e?.message || e?.code}`; }
  }
  const blob = data instanceof Blob ? data : new Blob([data]);
  const a = Object.assign(document.createElement('a'), { href: URL.createObjectURL(blob), download: filename });
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  return 'saved';
}
