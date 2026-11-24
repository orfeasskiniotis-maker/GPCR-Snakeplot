// Work out an (x, y) position for every residue of one receptor.
//
// Helix residues (TM1-7, H8) go straight into their fixed template slot,
// looked up by generic number. That is what makes plots comparable across
// receptors. Loops and termini have no fixed slots, so they are drawn here
// as curves between the helix ends they connect.

export const BEAD_STEP = 19;   // distance between neighbouring loop beads (px)
const MAX_LOOP = 20;           // loops longer than this are collapsed
const KEEP_EACH_SIDE = 5;      // residues kept visible at each end of a collapsed loop
const TAIL_KEEP = 10;          // residues kept visible on a collapsed terminus
const TAIL_ROW = 7;            // beads per row when a terminus folds back and forth

const isHelix = (r, slots) =>
  r.generic && slots[r.generic] && (r.segment.startsWith('TM') || r.segment === 'H8');

// Which side of the membrane a loop sits on: -1 outside (up), +1 inside (down).
const sideOf = segment =>
  segment.startsWith('ECL') || segment === 'N-term' ? -1 : 1;

export function layoutReceptor(receptor, template, expanded) {
  const slots = template.slots;
  const items = [];   // {kind:'residue', r, x, y, front} or {kind:'collapsed', id, count, x, y}

  // Split the sequence into runs: helix residues are placed, loops are queued.
  const runs = [];
  let current = null;
  receptor.residues.forEach(r => {
    if (isHelix(r, slots)) {
      current = null;
      const s = slots[r.generic];
      runs.push({ helix: true, r, x: s.x, y: s.y, front: s.front });
    } else {
      if (!current) { current = { helix: false, residues: [] }; runs.push(current); }
      current.residues.push(r);
    }
  });

  runs.forEach((run, i) => {
    if (run.helix) { items.push({ kind: 'residue', r: run.r, x: run.x, y: run.y, front: run.front }); return; }
    const prev = runs[i - 1], next = runs[i + 1];
    const id = `${run.residues[0].segment}-${run.residues[0].seq}`;
    const open = expanded.has(id);
    if (prev && next) items.push(...loopBetween(prev, next, run.residues, id, open));
    // N-terminus folds to the right above TM1; C-terminus folds back left under helix 8.
    else if (next) items.push(...tail(next, run.residues.slice().reverse(), +1, id, open).reverse());
    else if (prev) items.push(...tail(prev, run.residues, -1, id, open));
  });
  return items;
}

// A loop from helix end A to helix end B, drawn as a rounded "balloon"
// (a cubic Bézier curve whose control points lean outwards), sized so the
// beads come out evenly spaced. A balloon is wider and less tall than a U.
function loopBetween(a, b, residues, id, open) {
  const side = sideOf(residues[0].segment);
  let shown = residues, hidden = 0;
  if (!open && residues.length > MAX_LOOP) {
    hidden = residues.length - 2 * KEEP_EACH_SIDE;
    shown = [...residues.slice(0, KEEP_EACH_SIDE), null, ...residues.slice(-KEEP_EACH_SIDE)];
  }
  const n = shown.length;
  const target = (n + 1) * BEAD_STEP;
  const left = a.x <= b.x ? -1 : 1;  // which way is "away from B" for A
  const curve = d => bezier(a, { x: a.x + left * 0.8 * d, y: a.y + side * d },
                               { x: b.x - left * 0.8 * d, y: b.y + side * d }, b);

  // Binary search the depth d that gives the right curve length.
  let lo = 0, hi = 2000;
  if (curveLength(curve(0)) < target) {
    for (let k = 0; k < 40; k++) {
      const mid = (lo + hi) / 2;
      if (curveLength(curve(mid)) < target) lo = mid; else hi = mid;
    }
  }
  const points = evenlySpaced(curve(lo), n);
  return shown.map((r, k) => r
    ? { kind: 'residue', r, x: points[k].x, y: points[k].y, front: true }
    : { kind: 'collapsed', id, count: hidden, x: points[k].x, y: points[k].y });
}

// A free terminus: leave the anchor vertically, then fold back and forth in rows.
// residues are ordered starting next to the anchor; dir is the first row's
// horizontal direction (+1 right, -1 left).
function tail(anchor, residues, dir, id, open) {
  const side = sideOf(residues[0].segment);
  let shown = residues, hidden = 0;
  if (!open && residues.length > TAIL_KEEP + 2) {
    hidden = residues.length - TAIL_KEEP;
    shown = [...residues.slice(0, TAIL_KEEP), null];
  }
  const out = [];
  let x = anchor.x, y = anchor.y + side * BEAD_STEP * 1.2;
  let row = 0, col = 0;
  shown.forEach(r => {
    out.push(r
      ? { kind: 'residue', r, x, y, front: true }
      : { kind: 'collapsed', id, count: hidden, x, y });
    col++;
    if (col >= TAIL_ROW) { col = 0; row++; y += side * BEAD_STEP; }
    else x += (row % 2 === 0 ? dir : -dir) * BEAD_STEP;
  });
  return out;
}

// --- small curve helpers ---------------------------------------------------
function bezier(p0, p1, p2, p3) {
  return t => {
    const u = 1 - t;
    return {
      x: u * u * u * p0.x + 3 * u * u * t * p1.x + 3 * u * t * t * p2.x + t * t * t * p3.x,
      y: u * u * u * p0.y + 3 * u * u * t * p1.y + 3 * u * t * t * p2.y + t * t * t * p3.y,
    };
  };
}
function sample(f, steps = 400) {
  const pts = [];
  for (let i = 0; i <= steps; i++) pts.push(f(i / steps));
  return pts;
}
function curveLength(f) {
  const pts = sample(f);
  let len = 0;
  for (let i = 1; i < pts.length; i++) len += Math.hypot(pts[i].x - pts[i - 1].x, pts[i].y - pts[i - 1].y);
  return len;
}
// n points spread evenly along the curve, not counting its two ends.
function evenlySpaced(f, n) {
  const pts = sample(f, 800);
  const cum = [0];
  for (let i = 1; i < pts.length; i++) cum.push(cum[i - 1] + Math.hypot(pts[i].x - pts[i - 1].x, pts[i].y - pts[i - 1].y));
  const total = cum[cum.length - 1];
  const out = [];
  let j = 0;
  for (let k = 1; k <= n; k++) {
    const want = (total * k) / (n + 1);
    while (cum[j + 1] < want) j++;
    const t = (want - cum[j]) / (cum[j + 1] - cum[j] || 1);
    out.push({ x: pts[j].x + t * (pts[j + 1].x - pts[j].x), y: pts[j].y + t * (pts[j + 1].y - pts[j].y) });
  }
  return out;
}
