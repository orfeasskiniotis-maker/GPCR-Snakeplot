// Draw the snake plot as SVG with D3 and handle hover / click on residues.
//
// Three uses:
//   single        one receptor
//   side by side  two of these plots next to each other, each with `partner`
//                 set to the other receptor (for the compare colouring)
//   overlay       one plot with `overlay: true`: the partner's residues are
//                 drawn into the same template slots. Where the two amino
//                 acids differ the bead is split, A's letter left, B's right.
import { layoutReceptor } from './layout.js';
import { fillFor, isDark, RECEPTOR_COLORS } from './colors.js';

const R = 8.5;   // bead radius (px)
const PAD = 24;  // margin around the drawing

// frame: optional {x0, x1, y0, y1} the drawing must cover. Passing the same
// frame for every receptor keeps the helices at the same place on screen when
// switching tabs, which shows that the template slots really are shared.
export function createSnakePlot(svgEl, receptor, template, {
  onHover, onClick, frame = null, color = 'segment', label = 'aa', partner = null, overlay = false, motifs = [],
  custom = [],
}) {
  const svg = d3.select(svgEl);
  // custom: the user's highlight groups, [{color, seqs: Map(entry_name -> Set of seq)}].
  // A later group wins over an earlier one.
  const state = { color, label, motifs, custom, selA: null, selB: null, expanded: new Set() };
  const half = template.membrane_half_px;
  // Partner residue by generic number, for colouring and the overlay.
  const partnerByGeneric = new Map((partner?.residues || []).filter(r => r.generic).map(r => [r.generic, r]));
  const partnerAa = r => (r.generic ? partnerByGeneric.get(r.generic)?.aa : undefined);

  function buildItems() {
    const items = layoutReceptor(receptor, template, state.expanded);
    if (!overlay || !partner) return { items, ghosts: [] };
    // Attach B's residue to every A helix bead in the same slot.
    const slotted = new Set();
    items.forEach(d => {
      if (d.kind === 'residue' && d.r.generic && template.slots[d.r.generic]) {
        d.rb = partnerByGeneric.get(d.r.generic) || null;
        slotted.add(d.r.generic);
      }
    });
    // B's own layout gives its loops and any helix slots A does not use.
    const ghosts = layoutReceptor(partner, template, state.expanded)
      .filter(d => d.kind !== 'residue' || !slotted.has(d.r.generic))
      .map(d => ({ ...d, kind: d.kind === 'residue' ? 'ghost' : 'ghost-collapsed', rb: d.r, r: null }));
    return { items, ghosts };
  }

  function render() {
    const { items, ghosts } = buildItems();
    svg.selectAll('*').remove();

    const { x0, x1, y0, y1 } = bounds([...items, ...ghosts], half, frame);
    svg.attr('viewBox', `${x0} ${y0} ${x1 - x0} ${y1 - y0}`);

    // Membrane band with its two surfaces.
    svg.append('rect').attr('class', 'membrane')
      .attr('x', x0).attr('y', -half).attr('width', x1 - x0).attr('height', 2 * half);
    [-half, half].forEach(y => svg.append('line').attr('class', 'membrane-edge')
      .attr('x1', x0).attr('x2', x1).attr('y1', y).attr('y2', y));
    svg.append('text').attr('class', 'side-label').attr('text-anchor', 'end')
      .attr('x', x1 - 8).attr('y', -half - 10).text('Extracellular');
    svg.append('text').attr('class', 'side-label').attr('text-anchor', 'end')
      .attr('x', x1 - 8).attr('y', half + 18).text('Intracellular');

    // Helix names in a row along the top, above their columns.
    Object.entries(template.helices).forEach(([h, p]) => {
      svg.append('text').attr('class', 'helix-label').attr('x', p.X).attr('y', y0 + 18).text(`TM${h}`);
    });

    const line = d3.line().x(d => d.x).y(d => d.y).curve(d3.curveCatmullRom.alpha(0.5));

    // Overlay: B's chain as a dashed line, and B-only residues as small ghost beads.
    if (ghosts.length) {
      const bChain = layoutReceptor(partner, template, state.expanded);
      svg.append('path').attr('class', 'chain ghost-chain').attr('stroke', RECEPTOR_COLORS.b).attr('d', line(bChain));
      const g = svg.append('g').selectAll('g.ghost').data(ghosts).join('g')
        .attr('class', d => `ghost${d.kind === 'ghost-collapsed' ? ' collapsed' : ''}`)
        .attr('transform', d => `translate(${d.x},${d.y})`);
      g.append('circle').attr('r', d => (d.kind === 'ghost' && d.rb.generic ? R : 5))
        .attr('stroke', RECEPTOR_COLORS.b);
      g.filter(d => d.kind === 'ghost' && d.rb.generic).append('text').attr('fill', RECEPTOR_COLORS.b).text(d => d.rb.aa);
      g.filter(d => d.kind === 'ghost-collapsed').append('text').attr('fill', RECEPTOR_COLORS.b).text(d => `+${d.count}`);
      g.on('mouseenter', (event, d) => d.rb && onHover(null, event, d.rb))
        .on('mousemove', (event, d) => d.rb && onHover(null, event, d.rb))
        .on('mouseleave', () => onHover(null))
        .on('click', (event, d) => {
          if (d.kind === 'ghost-collapsed') { state.expanded.add(d.id); render(); return; }
          onClick(null, d.rb);
        });
    }

    // A's chain: one line through every bead in sequence order.
    svg.append('path').attr('class', 'chain').attr('d', line(items))
      .attr('stroke', overlay ? RECEPTOR_COLORS.a : null);

    // Back-facing helix residues first, so front-facing ones sit on top.
    const ordered = [...items].sort((a, b) => (a.front === b.front ? 0 : a.front ? 1 : -1));
    const beads = svg.append('g').selectAll('g.bead').data(ordered).join('g')
      .attr('class', d => [
        'bead',
        d.kind === 'collapsed' ? 'collapsed' : '',
        d.front === false ? 'back' : '',
        state.label === 'seq' ? 'small-text' : '',
      ].join(' '))
      .attr('transform', d => `translate(${d.x},${d.y})`);

    beads.append('circle').attr('r', R).each(function (d) {
      if (d.kind !== 'residue') return;
      const c = paint(d.r);
      d3.select(this).attr('fill', c.fill).attr('stroke', c.stroke);
    });

    beads.classed('in-custom', d => d.kind === 'residue' && !!paint(d.r).custom);

    const textColor = d => (d.kind === 'residue' && paint(d.r).whiteText ? '#fff' : null);
    // Split beads (overlay, amino acids differ): A's letter left, B's right.
    const split = d => overlay && d.kind === 'residue' && d.rb && d.rb.aa !== d.r.aa && state.label === 'aa';
    beads.filter(d => !split(d)).append('text')
      .text(d => d.kind === 'collapsed' ? `+${d.count}` : labelFor(d.r))
      .attr('style', d => (textColor(d) ? `fill:${textColor(d)}` : null));
    const halves = beads.filter(split);
    halves.append('line').attr('class', 'split').attr('y1', -R + 1.5).attr('y2', R - 1.5);
    halves.append('text').attr('class', 'half').attr('x', -3.4).attr('fill', RECEPTOR_COLORS.a).text(d => d.r.aa);
    halves.append('text').attr('class', 'half').attr('x', 3.4).attr('fill', RECEPTOR_COLORS.b).text(d => d.rb.aa);

    // Motif rings: one coloured ring per active motif the residue belongs to
    // (a residue in two motifs, like N7x49, gets two rings).
    beads.filter(d => d.kind === 'residue' && d.r.generic).each(function (d) {
      const hits = state.motifs.filter(m => m.generic.includes(d.r.generic));
      hits.forEach((m, i) => d3.select(this).insert('circle', ':first-child')
        .attr('class', 'motif-ring').attr('r', R + 3 + i * 2.6).attr('stroke', m.color));
      if (hits.length) d3.select(this).classed('in-motif', true);
    });

    // Mark the x50 reference residue of each helix with its generic number.
    items.filter(d => d.kind === 'residue' && /^[1-8]x50$/.test(d.r.generic)).forEach(d => {
      svg.append('text').attr('class', 'x50-label')
        .attr('x', d.x - R - 3).attr('y', d.y + 3).attr('text-anchor', 'end').text(d.r.generic);
    });

    beads
      .on('mouseenter', function (event, d) {
        d3.select(this).classed('hover', true);
        if (d.kind === 'residue') onHover(d.r, event, d.rb);
      })
      .on('mousemove', (event, d) => d.kind === 'residue' && onHover(d.r, event, d.rb))
      .on('mouseleave', function () { d3.select(this).classed('hover', false); onHover(null); })
      .on('click', (event, d) => {
        if (d.kind === 'collapsed') { state.expanded.add(d.id); render(); return; }
        onClick(d.r, d.rb);
      });

    markSelected();
  }

  // Bead colour: the user's highlight if the residue is in a group, else the colour mode.
  function paint(r) {
    const group = [...state.custom].reverse().find(g => g.seqs.get(receptor.entry_name)?.has(r.seq));
    if (group) return { fill: group.color, stroke: d3.rgb(group.color).darker(1).formatHex(), whiteText: isDark(group.color), custom: true };
    return fillFor(r, state.color, partnerAa(r));
  }

  function labelFor(r) {
    if (state.label === 'seq') return r.seq;
    if (state.label === 'generic') return r.generic && r.generic.includes('x') ? r.generic.split('x')[1] : '';
    return r.aa;
  }

  function markSelected() {
    svg.selectAll('g.bead, g.ghost').classed('selected', d =>
      (d.r && d.r.seq === state.selA) || (d.rb && d.rb.seq === state.selB));
    svg.selectAll('.selected').raise();
  }

  // Expand a collapsed loop when the residue to select is hidden inside it.
  function reveal(rec, seq) {
    if (seq == null || !rec) return false;
    const shown = layoutReceptor(rec, template, state.expanded).some(d => d.kind === 'residue' && d.r.seq === seq);
    if (shown) return false;
    const r = rec.residues.find(x => x.seq === seq);
    if (!r) return false;
    const first = rec.residues.find(x => x.segment === r.segment);
    state.expanded.add(`${r.segment}-${first.seq}`);
    return true;
  }

  // Make sure every residue of the active motifs is visible (e.g. C45x50 in ECL2).
  function revealMotifs() {
    const generics = new Set(state.motifs.flatMap(m => m.generic));
    [receptor, overlay ? partner : null].filter(Boolean).forEach(rec =>
      rec.residues.filter(r => generics.has(r.generic)).forEach(r => reveal(rec, r.seq)));
  }

  // ...and every residue of the user's highlight groups.
  function revealCustom() {
    state.custom.forEach(g => (g.seqs.get(receptor.entry_name) || []).forEach(seq => reveal(receptor, seq)));
  }

  revealMotifs();
  revealCustom();
  render();
  return {
    set(key, value) {
      state[key] = value;
      if (key === 'motifs') revealMotifs();
      if (key === 'custom') revealCustom();
      render();
    },
    // seqA: residue of this plot's receptor; seqB: of the partner (overlay only).
    select(seqA, seqB = null) {
      state.selA = seqA; state.selB = seqB;
      const changed = reveal(receptor, seqA) | (overlay ? reveal(partner, seqB) : false);
      if (changed) render(); else markSelected();
    },
  };
}

// The box (in template px) that a layout needs, grown to cover `frame` if given.
export function bounds(items, half, frame = null) {
  const xs = items.map(d => d.x), ys = items.map(d => d.y);
  const b = {
    x0: Math.min(...xs) - PAD - 40, x1: Math.max(...xs) + PAD,
    y0: Math.min(...ys, -half) - PAD - 26, y1: Math.max(...ys, half) + PAD,
  };
  if (frame) {
    b.x0 = Math.min(b.x0, frame.x0); b.x1 = Math.max(b.x1, frame.x1);
    b.y0 = Math.min(b.y0, frame.y0); b.y1 = Math.max(b.y1, frame.y1);
  }
  return b;
}
