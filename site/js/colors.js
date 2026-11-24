// Fill colours for the residue beads.

// AlphaFold's standard pLDDT bands, the same ones Mol* uses in 3D.
export const PLDDT_BANDS = [
  { min: 90, color: '#0053D6', label: 'Very high (≥90)' },
  { min: 70, color: '#65CBF3', label: 'Confident (70–90)' },
  { min: 50, color: '#FFDB13', label: 'Low (50–70)' },
  { min: -1, color: '#FF7D45', label: 'Very low (<50)' },
];

export function plddtColor(v) {
  if (v == null) return '#999';
  return PLDDT_BANDS.find(b => v >= b.min).color;
}

// The two receptors in a comparison: A is the selected tab, B the partner.
// The same colours are used for the two models in 3D.
export const RECEPTOR_COLORS = { a: '#2f6fd6', b: '#e0702b' };

// Comparing two residues at the same generic number.
// Groups of amino acids with similar side chains count as "similar".
const GROUPS = ['AVLIM', 'FWY', 'ST', 'NQ', 'DE', 'KRH', 'G', 'P', 'C'];
export const COMPARE_CLASSES = [
  { key: 'identical', color: '#bfe3c6', edge: '#4f9a5c', label: 'Identical' },
  { key: 'similar', color: '#fbe3a6', edge: '#c49a2c', label: 'Similar side chain' },
  { key: 'different', color: '#f5c4c4', edge: '#c25555', label: 'Different' },
  { key: 'unmatched', color: 'var(--loop)', edge: 'var(--loop-edge)', label: 'No counterpart' },
];

export function compareClass(aa, bb) {
  if (!aa || !bb) return 'unmatched';
  if (aa === bb) return 'identical';
  return GROUPS.some(g => g.includes(aa) && g.includes(bb)) ? 'similar' : 'different';
}

// Segment colours come from CSS variables so they follow the page tokens.
export function segmentKind(segment) {
  if (segment.startsWith('TM')) return 'tm';
  if (segment === 'H8') return 'h8';
  return 'loop';
}

// ---------------------------------------------------- residue properties
// Hydrophobicity: Kyte & Doolittle (1982) scale, from -4.5 (R, most
// hydrophilic) to +4.5 (I, most hydrophobic). Drawn blue -> white -> orange.
export const KYTE_DOOLITTLE = {
  I: 4.5, V: 4.2, L: 3.8, F: 2.8, C: 2.5, M: 1.9, A: 1.8, G: -0.4, T: -0.7, S: -0.8,
  W: -0.9, Y: -1.3, P: -1.6, H: -3.2, E: -3.5, Q: -3.5, D: -3.5, N: -3.5, K: -3.9, R: -4.5,
};
export const HYDRO_ENDS = { philic: '#2f6fd6', phobic: '#e0702b' };
export function hydroColor(aa) {
  const v = KYTE_DOOLITTLE[aa];
  return v == null ? '#cccccc' : hydroScale(v);
}
export function hydroScale(v) {
  return v < 0 ? d3.interpolateRgb('#ffffff', HYDRO_ENDS.philic)(-v / 4.5 * 0.85)
               : d3.interpolateRgb('#ffffff', HYDRO_ENDS.phobic)(v / 4.5 * 0.85);
}

// Charge at pH 7. His is mostly neutral at pH 7 but often protonated, so it
// gets a pale blue.
export const CHARGE_CLASSES = [
  { key: 'pos', aas: 'KR', color: '#3b6fd8', label: 'Positive (K, R)' },
  { key: 'his', aas: 'H', color: '#b7cdf3', label: 'His (partly positive)' },
  { key: 'neg', aas: 'DE', color: '#d64545', label: 'Negative (D, E)' },
  { key: 'neutral', aas: '', color: '#f2f3f5', label: 'Uncharged' },
];

// Residue type: side-chain chemistry, the groups used by most viewers.
export const TYPE_CLASSES = [
  { aas: 'AVLIM', color: '#f3d27c', label: 'Aliphatic (A V L I M)' },
  { aas: 'FWY', color: '#c7a3df', label: 'Aromatic (F W Y)' },
  { aas: 'STNQ', color: '#9fd8b0', label: 'Polar (S T N Q)' },
  { aas: 'KRH', color: '#8fb3ee', label: 'Basic (K R H)' },
  { aas: 'DE', color: '#f09a9a', label: 'Acidic (D E)' },
  { aas: 'G', color: '#e2e4e8', label: 'Glycine' },
  { aas: 'P', color: '#e9c19b', label: 'Proline' },
  { aas: 'C', color: '#f2ec86', label: 'Cysteine' },
];

// Dark fills get white letters (relative luminance below about 0.35).
export function isDark(color) {
  const c = d3.rgb(color);
  const l = (0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b) / 255;
  return l < 0.55;
}
const solid = color => ({ fill: color, stroke: d3.rgb(color).darker(0.9).formatHex(), whiteText: isDark(color) });

// r: the residue; partnerAa: the other receptor's amino acid at the same
// generic number (only used by the "compare" colouring).
export function fillFor(r, mode, partnerAa) {
  if (mode === 'plddt') return { fill: plddtColor(r.plddt), stroke: 'rgba(0,0,0,.35)', whiteText: r.plddt >= 90 };
  if (mode === 'compare') {
    const c = COMPARE_CLASSES.find(k => k.key === compareClass(r.aa, partnerAa));
    return { fill: c.color, stroke: c.edge, whiteText: false };
  }
  if (mode === 'hydro') return solid(hydroColor(r.aa));
  if (mode === 'charge') return solid((CHARGE_CLASSES.find(c => c.aas.includes(r.aa)) || CHARGE_CLASSES[3]).color);
  if (mode === 'type') return solid((TYPE_CLASSES.find(c => c.aas.includes(r.aa)) || { color: '#cccccc' }).color);
  const k = segmentKind(r.segment);
  return { fill: `var(--${k})`, stroke: `var(--${k}-edge)`, whiteText: false };
}
