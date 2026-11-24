// The 3D viewer (Mol*). Loads one or two AlphaFold models and lets the snake
// plot focus, highlight and pick residues by (model index, sequence number).
// Models come from the pipeline already superposed on the class reference,
// so two receptors loaded together sit on top of each other.
//
// The Mol* "viewer" bundle does not export its query helpers, so residues
// are found by walking the loaded structure directly (residueLoci below).

export async function createViewer3D(el, models, { onPick }) {
  let model = models[0];
  const viewer = await molstar.Viewer.create(el, {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowRemoteState: false,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowLeftPanel: false,
    viewportShowExpand: false,
    viewportShowSelectionMode: false,
    viewportShowAnimation: false,
    viewportShowTrajectoryControls: false,
    collapseLeftPanel: true,
  });
  const plugin = viewer.plugin;

  // Turn the model so it matches the snake plot: membrane horizontal,
  // extracellular side up, TM1 on the left. The vectors come from the pipeline.
  function orient(durationMs = 0) {
    const v = model?.view;
    if (!v) return plugin.managers.camera.reset();
    const camera = plugin.canvas3d.camera;
    const distance = 95;  // Å from the bundle centre, enough to fit the TM bundle
    const position = v.target.map((t, i) => t + v.toward_camera[i] * distance);
    camera.setState({ ...camera.getSnapshot(), target: v.target, position, up: v.up }, durationMs);
  }

  // White background (Mol* colours are 0xRRGGBB numbers).
  plugin.canvas3d?.setProps({ renderer: { backgroundColor: 0xffffff } });

  // Replace what is shown with these models. colors: one 0xRRGGBB per model
  // to paint each receptor in a single colour; leave it out for pLDDT colours.
  async function load(list, colors = null) {
    model = list[0];
    await plugin.clear();
    for (const m of list) await viewer.loadStructureFromUrl(m.file, 'mmcif', false);
    if (colors) {
      const structures = plugin.managers.structure.hierarchy.current.structures;
      for (let i = 0; i < structures.length; i++) {
        await plugin.managers.structure.component.updateRepresentationsTheme(
          structures[i].components, { color: 'uniform', colorParams: { value: colors[i] } });
      }
    }
    orient();
  }
  await load(models);

  const structures = () => plugin.managers.structure.hierarchy.current.structures.map(s => s.cell.obj?.data);

  // A Mol* "loci" (selection) covering every atom of one residue of model `index`.
  function residueLoci(index, seq) {
    const s = structures()[index];
    if (!s) return null;
    const elements = [];
    for (const unit of s.units) {
      const h = unit.model.atomicHierarchy;
      const idx = [];
      for (let i = 0; i < unit.elements.length; i++) {
        const res = h.residueAtomSegments.index[unit.elements[i]];
        if (h.residues.label_seq_id.value(res) === seq) idx.push(i);
      }
      if (idx.length) elements.push({ unit, indices: Int32Array.from(idx) });
    }
    return elements.length ? { kind: 'element-loci', structure: s, elements } : null;
  }

  // Read the residue number back out of a loci that Mol* hands us on click.
  function seqFromLoci(loci) {
    if (!loci || loci.kind !== 'element-loci' || !loci.elements.length) return null;
    const { unit, indices } = loci.elements[0];
    const first = typeof indices === 'number' ? intervalStart(indices) : indices[0];
    const res = unit.model.atomicHierarchy.residueAtomSegments.index[unit.elements[first]];
    return unit.model.atomicHierarchy.residues.label_seq_id.value(res);
  }

  plugin.behaviors.interaction.click.subscribe(({ current }) => {
    const seq = seqFromLoci(current.loci);
    if (seq == null) return;
    const s = current.loci.structure;
    const index = structures().findIndex(x => x === s || x === s.root);
    onPick(Math.max(index, 0), seq);
  });

  // targets: [{index, seq}, ...]; entries with seq == null are skipped.
  const lociFor = targets => targets.filter(t => t.seq != null)
    .map(t => residueLoci(t.index, t.seq)).filter(Boolean);

  return {
    load,
    resize() { viewer.handleResize(); },
    focus(targets) {
      const all = lociFor(targets);
      if (!all.length) return;
      plugin.managers.interactivity.lociSelects.deselectAll();
      all.forEach(loci => plugin.managers.interactivity.lociSelects.select({ loci }));
      plugin.managers.structure.focus.setFromLoci(all[0]);   // side chain + neighbours as sticks
      plugin.managers.camera.focusLoci(all[0], { extraRadius: 11, durationMs: 400 });
    },
    // Show a set of residues (e.g. a motif) as selected, without moving the camera.
    mark(targets) {
      plugin.managers.interactivity.lociSelects.deselectAll();
      lociFor(targets).forEach(loci => plugin.managers.interactivity.lociSelects.select({ loci }));
    },
    highlight(targets) {
      plugin.managers.interactivity.lociHighlights.clearHighlights();
      lociFor(targets || []).forEach(loci => plugin.managers.interactivity.lociHighlights.highlight({ loci }));
    },
    reset() {
      plugin.managers.interactivity.lociSelects.deselectAll();
      plugin.managers.structure.focus.clear();
      orient(400);
    },
  };
}

// Mol* stores a contiguous index range as two 32-bit ints packed into one
// float64 (its "Interval" type). This unpacks the first index.
const f64 = new Float64Array(1);
const i32 = new Int32Array(f64.buffer);
function intervalStart(packed) { f64[0] = packed; return i32[0]; }
