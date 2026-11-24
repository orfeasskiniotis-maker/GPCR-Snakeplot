// Custom selections: turn what the user types into residues of a receptor.
//
// Accepted, separated by commas or spaces:
//   113  D113          residue number (letter optional)
//   110-120           a range of residue numbers
//   3x50  3.50  3.50x50  R3x50   GPCRdb generic number (helix x position)
//   3x48-3x52  3x48-52    a range of generic numbers on one helix
//   TM3  ECL2  H8  C-term   a whole segment
//
// Residue numbers only make sense on the receptor they were typed for, so
// they apply to the page's receptor. Generic numbers and segments apply to
// every receptor, so in compare modes they also mark the partner.

export function parseSelection(text) {
  const tokens = text.replace(/\s*-\s*/g, '-').split(/[\s,;]+/).filter(Boolean);
  return tokens.map(raw => {
    const t = raw.trim();
    let m;
    if ((m = t.match(/^[A-Za-z]?(\d+)$/))) return { raw, kind: 'seq', from: +m[1], to: +m[1] };
    if ((m = t.match(/^[A-Za-z]?(\d+)-[A-Za-z]?(\d+)$/))) return { raw, kind: 'seq', from: Math.min(+m[1], +m[2]), to: Math.max(+m[1], +m[2]) };
    const g = s => {   // "3x50", "3.50", "3.50x50", "R3x50" -> {helix: "3", pos: 50}
      let q;
      s = s.replace(/^[A-Za-z](?=\d)/, '');
      if ((q = s.match(/^(\d{1,2})\.\d{2,3}x(\d{2,3})$/))) return { helix: q[1], pos: +q[2] };
      if ((q = s.match(/^(\d{1,2})[x.](\d{2,3})$/))) return { helix: q[1], pos: +q[2] };
      return null;
    };
    if (g(t)) { const a = g(t); return { raw, kind: 'generic', helix: a.helix, from: a.pos, to: a.pos }; }
    if ((m = t.match(/^(.+)-(.+)$/))) {
      const a = g(m[1]), b = g(m[2]) || (/^\d{2,3}$/.test(m[2]) && a ? { helix: a.helix, pos: +m[2] } : null);
      if (a && b && a.helix === b.helix)
        return { raw, kind: 'generic', helix: a.helix, from: Math.min(a.pos, b.pos), to: Math.max(a.pos, b.pos) };
    }
    if (/^(TM[1-7]|ICL[1-3]|ECL[1-3]|H8|N-?term|C-?term)$/i.test(t))
      return { raw, kind: 'segment', segment: t.toUpperCase().replace(/^([NC])-?TERM$/, '$1-term') };
    return { raw, kind: 'bad' };
  });
}

// Generic key "3x50" or bulge "3x501" -> position 50 / 50.1 on helix "3".
function genericPos(key) {
  const [h, p] = key.split('x');
  return { helix: h, pos: p.length === 3 ? +p.slice(0, 2) + +p[2] / 10 : +p };
}

// Residues of `rec` matched by the parsed tokens. primary: the page's own
// receptor (residue numbers apply only there). Returns the set of residue
// numbers and the tokens that matched nothing.
export function resolveSelection(tokens, rec, primary) {
  const seqs = new Set(), missing = [];
  for (const tk of tokens) {
    let hits = [];
    if (tk.kind === 'seq') {
      if (!primary) continue;
      hits = rec.residues.filter(r => r.seq >= tk.from && r.seq <= tk.to);
    } else if (tk.kind === 'generic') {
      hits = rec.residues.filter(r => {
        if (!r.generic) return false;
        const g = genericPos(r.generic);
        return g.helix === tk.helix && g.pos >= tk.from && g.pos <= tk.to;
      });
    } else if (tk.kind === 'segment') {
      hits = rec.residues.filter(r => r.segment.toUpperCase() === tk.segment.toUpperCase());
    }
    hits.forEach(r => seqs.add(r.seq));
    if (primary && !hits.length) missing.push(tk.raw);
  }
  return { seqs, missing };
}
