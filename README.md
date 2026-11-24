# GPCR Snake Plot Portal (prototype)

387 human non-olfactory GPCRs: Class A 286 (NPY6R, a pseudogene without a
TM7, is skipped), B1 15, B2 28, C 22, F 11, T2 25. Five adhesion receptors
(CELSR1-3, ADGRG4, ADGRV1) have no full-length AlphaFold DB model and are
missing. Each class has its own template keyed by GPCRdb generic numbers.
For A, B1, C and F the helix heights and phases are averaged over
inactive-state experimental structures placed in the membrane by OPM (one per
receptor: A 64, B1 5, C 8, F 4). B2 and T2 have no usable inactive
structures, so their templates are measured on one AlphaFold model (ADGRL3,
TAS2R46) and the page says so under the plot. Every class has motif presets; the exceptions found on
individual receptors are listed in site/js/motifs.js. The home page lists them by ligand type and family
with a search box; each receptor page is at #<entry_name>, e.g. #adrb2_human.
Every page uses the same plot frame, so a residue with the same generic number
sits on the same spot on every receptor.

View modes: **Single**, **Side by side** (two plots, linked selection) and
**Overlay** (both receptors in the same template slots; a split bead means the
amino acids differ). "Compare" colouring marks identical / similar / different
residues. In both compare modes the 3D viewer shows the two AlphaFold models
superposed on their shared TM C-alphas.

Motif presets (E/DRY, ionic lock, CWxP, PIF, Na⁺ site, NPxxY, disulfide)
ring their residues on the plot and select them in 3D; several can be on at
once. Each preset is a list of generic numbers in site/js/motifs.js.

Colour modes: segment, pLDDT, hydrophobicity (Kyte-Doolittle), charge,
residue type, and compare. **Highlight**: type residues (113, D113, 110-120),
generic numbers (3x50, 3.50, 6x48-6x52) or segments (TM3, ECL2), pick a
colour and add it as a group; groups also show in 3D and survive switching
receptor (generic numbers and segments carry over, residue numbers apply to
the receptor they were typed on). **Export**: SVG, PNG (3x) and a CSV list of
the highlighted and motif residues.

A template-based snake plot (D3) sits next to the AlphaFold model (Mol*).
Hover a residue for its details; click it to focus it in 3D; click a residue
in 3D to highlight it on the plot.

## Layout

```
pipeline/                Python: download data and write static JSON
  sources.py             GPCRdb / AlphaFold DB downloads, cached in pipeline/cache/
  generic_numbers.py     "3.50x50" -> key 3x50, helix 3, position 50 (bulges: 2x551 -> 55.1)
  build_template.py      Class template: generic number -> fixed (x, y) slot
  experimental.py        Averages helix heights/phases over OPM-oriented structures
  build_receptor.py      One receptor's residues, pLDDT, model file and 3D camera
  superpose.py           Moves each model onto the class reference (Kabsch on TM C-alphas)
  batch_build.py         Every receptor of a class + site/data/index.json
site/                    Static website, no build step
  index.html, style.css
  js/main.js             routing (index / receptor), data loading, plot <-> 3D, toolbar
  js/index-page.js       searchable receptor list
  js/layout.js           helix residues from template slots; loops/termini as curves
  js/snakeplot.js        SVG drawing and interaction (D3)
  js/viewer3d.js         Mol* viewer: load, focus, highlight, pick
  js/colors.js           segment, pLDDT and compare colours
  js/motifs.js           motif presets (generic numbers per class)
  js/selection.js        parse the Highlight box into residues
  js/export.js           SVG / PNG / CSV export (inlines CSS so files look the same)
  data/                  output of the pipeline
  vendor/                d3 7.9, Mol* 4.18
```

## Run it

```
cd pipeline
pip install -r requirements.txt
python build_template.py adrb2_human A   # writes site/data/templates/classA.json
python build_receptor.py adrb2_human     # writes site/data/receptors/ and models/

cd ../site
python3 -m http.server 8000              # open http://localhost:8000
```

All Class A receptors at once (about 20 minutes and ~100 MB the first time,
cached after that; also writes site/data/index.json):

```
python build_template.py adrb2_human A
python build_template.py glp1r_human B1   # needs gpcrdb.org and opm-assets.storage.googleapis.com
python batch_build.py A
python batch_build.py B1      # adds B1 to the same index
# same for B2 (agrl3_human), C (grm5_human), F (smo_human), T2 (t2r46_human)
```

One receptor: `python build_receptor.py <entry_name>` (then re-run
batch_build.py, or add it to site/data/index.json by hand). Set GPCR_CACHE to
keep the download cache somewhere else. A new class needs its own template and an
entry in REFERENCE (build_receptor.py) and CLASS_SLUG (batch_build.py).

## Hosting (GitHub Pages)

Live site: https://orfeasskiniotis-maker.github.io/GPCR-Snakeplot/

`.github/workflows/pages.yml` uploads the `site/` folder to GitHub Pages on
every push to `main` (no build step). One-time setup: Settings > Pages >
Source: **GitHub Actions**. All paths in the site are relative, so it works
under the `/GPCR-Snakeplot/` sub-path.

The AlphaFold models in `site/data/models/` (~150 MB, 387 files, largest
1.3 MB) are committed rather than fetched from AlphaFold DB at runtime: the
pipeline superposes every model onto its class reference, so the files differ
from the AlphaFold DB originals, and serving them from the same site avoids a
second server on first load. This sits well inside GitHub's limits (100 MB
per file, 1 GB per Pages site). Pipeline download caches (`pipeline/cache/`)
are not committed.

## Known simplifications

- Templates average only inactive-state structures that OPM has (its file
  store lacks many recent cryo-EM entries: about 100 Class A entries were
  skipped). B1 and F rest on 5 and 4 structures. B2 and T2 use one AlphaFold
  model each, with a membrane estimated from the model and assumed 30 Å thick.
- Helix tilt is ignored: every residue rises 1.5 Å, as in the roadmap.
- The angle formula carries the helix direction in θ as well as in y
  (θ = φ + 100° · d · (p − 50)); see build_template.py for why.
- Loops longer than 20 residues and termini longer than 12 are collapsed
  into a grey bead; click it to expand.

## Credits and licences

- Structures: [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk)
  (Jumper et al. 2021, Nature 596:583; Varadi et al. 2024, NAR 52:D368),
  licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The
  model files here are modified: their coordinates are superposed onto a
  class reference model.
- Generic residue numbers, receptor families and segments:
  [GPCRdb](https://gpcrdb.org) (Pándy-Szekeres et al. 2023, NAR 51:D395).
- Membrane orientation of experimental structures: [OPM](https://opm.phar.umich.edu).
- Sequences and accessions: [UniProt](https://www.uniprot.org).
- Libraries: [D3](https://d3js.org) (ISC) and [Mol*](https://molstar.org) (MIT).
