// Conserved motif presets for the toolbar (from the internship roadmap).
//
// Each preset is just a list of generic numbers, so it works on every
// receptor of its class without extra code: the template puts e.g. 3x50 on
// the same spot for all of them.
//
// Class A positions were checked against GPCRdb residue data for β2AR and
// D2R (2026-11-24): every position exists and carries the expected amino
// acid (D/E-R-Y, E6x30, C-W-x-P, P-I-F, D2x50/S3x39/N7x45/S7x46/N7x49,
// N-P-x-x-Y, C3x25/C45x50). Class B1 positions were checked on all 15 human
// B1 receptors: HETX is H2x50-E3x50-T6x42-Y7x57 (CALCRL has F7x57), the polar
// network R2x60-N3x43-H6x52-Q7x49 holds in the glucagon/PTH/VIP families and
// varies elsewhere (CRF receptors R-N-T-Q, calcitonin N-N-Q/E-Q), and PxxG is
// P6x47...G6x50 in all except PTH2R (L6x47).
//
// Classes B2, C, F and T2 (added 2026-11-24): motifs described in the
// literature, mapped to GPCRdb numbers and checked on every human receptor of
// the class in this portal. Exceptions:
//   B2 (28)  Disulfide C3x29-C45x50 in all. HETX: H2x50-E3x50 only in 14
//            (ADGRL/E/D/C/V-like); ADGRA, ADGRB, ADGRF and some ADGRG differ.
//            TM6 kink G6x50 in 22 (not ADGRA1-3, ADGRB1-3). Toggle W6x53 in 26
//            (ADGRD2 V, ADGRF2 T). TM7 Q7x49-G7x50 in 22 (not ADGRA1-3,
//            ADGRB1, ADGRF2, ADGRG5).
//   C (22)   Ionic lock K3x50-E6x35 in 11: GABA-B has K-D, taste receptors
//            R-E, the GPRC5 orphans and GPR156 lack it, GRM4 has no 6x35.
//            Disulfide C3x29-C45x50 in 18 (not GPR156, GPRC5A/D, TAS1R3).
//            W6x50 in 16 (GABA-B C, GPR156/158/179 and TAS1R2 differ).
//   F (11)   Switch R6x32-W7x55 in 9 (FZD4 and FZD9 have K6x32). KTxxxW
//            K8x49-T8x50-W8x54 in 10 (SMO has A8x49). Disulfide C3x25-C45x50 in all.
//   T2 (25)  LxxxR L2x46-R2x50 in 24 (TAS2R42 T2x50). FYxxK F3x49-Y3x50-K3x53
//            in 17 (TAS2R3/4/38 L, R30/R40 R, R7 Y, R19 C, R41 F). LxxSL
//            L5x61-S5x64-L5x65 in 24 (TAS2R14 M). HSxxL H7x49-S7x50-L7x53 in
//            15; seven have P7x50 (TAS2R30/31/41/43/45/46/60).

export const MOTIFS = [
  { id: 'dry', name: 'E/DRY', cls: 'A', color: '#d62f6c', generic: ['3x49', '3x50', '3x51'] },
  { id: 'ionic', name: 'Ionic lock', cls: 'A', color: '#7b3fc4', generic: ['3x50', '6x30'] },
  { id: 'cwxp', name: 'CWxP toggle', cls: 'A', color: '#0e8f7e', generic: ['6x47', '6x48', '6x49', '6x50'] },
  { id: 'pif', name: 'PIF', cls: 'A', color: '#e07b00', generic: ['3x40', '5x50', '6x44'] },
  { id: 'na', name: 'Na⁺ site', cls: 'A', color: '#2266cc', generic: ['2x50', '3x39', '7x45', '7x46', '7x49'] },
  { id: 'npxxy', name: 'NPxxY', cls: 'A', color: '#8a6d00', generic: ['7x49', '7x50', '7x51', '7x52', '7x53'] },
  { id: 'ss', name: 'Disulfide', cls: 'A', color: '#b8860b', generic: ['3x25', '45x50'] },
  // Class B1
  { id: 'hetx', name: 'HETX', cls: 'B1', color: '#d62f6c', generic: ['2x50', '3x50', '6x42', '7x57'] },
  { id: 'polar', name: 'Central polar network', cls: 'B1', color: '#2266cc', generic: ['2x60', '3x43', '6x52', '7x49'] },
  { id: 'pxxg', name: 'PxxG kink', cls: 'B1', color: '#0e8f7e', generic: ['6x47', '6x50'] },
  // Class B2 (adhesion)
  { id: 'hetx', name: 'HETX', cls: 'B2', color: '#d62f6c', generic: ['2x50', '3x50', '6x42', '7x57'] },
  { id: 'kink6', name: 'TM6 kink (xxG)', cls: 'B2', color: '#0e8f7e', generic: ['6x47', '6x50'] },
  { id: 'toggle', name: 'W6x53 toggle', cls: 'B2', color: '#7b3fc4', generic: ['6x53'] },
  { id: 'qg', name: 'TM7 QG', cls: 'B2', color: '#e07b00', generic: ['7x49', '7x50'] },
  { id: 'ss', name: 'Disulfide', cls: 'B2', color: '#b8860b', generic: ['3x29', '45x50'] },
  // Class C (glutamate)
  { id: 'ionic', name: 'Ionic lock', cls: 'C', color: '#7b3fc4', generic: ['3x50', '6x35'] },
  { id: 'toggle', name: 'W6x50 toggle', cls: 'C', color: '#0e8f7e', generic: ['6x50'] },
  { id: 'ss', name: 'Disulfide', cls: 'C', color: '#b8860b', generic: ['3x29', '45x50'] },
  // Class F (frizzled)
  { id: 'switch', name: 'R6x32–W7x55 switch', cls: 'F', color: '#d62f6c', generic: ['6x32', '7x55'] },
  { id: 'ktxxxw', name: 'KTxxxW', cls: 'F', color: '#2266cc', generic: ['8x49', '8x50', '8x54'] },
  { id: 'ss', name: 'Disulfide', cls: 'F', color: '#b8860b', generic: ['3x25', '45x50'] },
  // Class T2 (taste 2)
  { id: 'lxxxr', name: 'LxxxR (TM2)', cls: 'T2', color: '#2266cc', generic: ['2x46', '2x50'] },
  { id: 'fyxxk', name: 'FYxxK (TM3)', cls: 'T2', color: '#d62f6c', generic: ['3x49', '3x50', '3x53'] },
  { id: 'lxxsl', name: 'LxxSL (TM5)', cls: 'T2', color: '#e07b00', generic: ['5x61', '5x64', '5x65'] },
  { id: 'hsxxl', name: 'HSxxL (TM7)', cls: 'T2', color: '#8a6d00', generic: ['7x49', '7x50', '7x53'] },
];

export const motifsForClass = cls => MOTIFS.filter(m => m.cls === cls);
