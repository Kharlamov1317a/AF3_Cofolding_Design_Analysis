AF3 TARGET-GROUP ANALYSIS v2
============================

Input directory:
  /mnt/gs21/scratch/dolgikhb/af3/outputs/ppi_ids_szymanskiGH18s_run

Reference structure:
  /mnt/gs21/scratch/dolgikhb/af3/outputs/itr3_pd/ppi_cichi1_itr3_run/binder9_target1/binder9_target1/binder9_target1_model.cif

Predicted chains:
  binder = A
  target = B

Reference chains:
  binder = A
  target = B

Model mode:
  top

PRIMARY RMSD: target_orientation_rmsd_A
---------------------------------------
This metric is designed to answer:

  Does the predicted complex adopt a similar target placement/orientation
  relative to the binder as the reference complex?

For every predicted complex:

1. Predicted/reference binder residues are matched by global sequence alignment.
   This is appropriate when binder sequences are identical or highly similar.
2. The predicted binder is least-squares superposed onto the reference binder.
3. Predicted/reference target chains are structurally aligned with Biopython's
   sequence-independent Combinatorial Extension (CE) algorithm.
4. The CE alignment is used ONLY to identify structurally equivalent target
   C-alpha residue pairs. This makes correspondence robust when target sequence
   identity is low but structural similarity is moderate/high.
5. The CE target superposition is discarded.
6. The BINDER-derived rotation/translation from step 2 is applied to the
   predicted target.
7. RMSD is calculated between the binder-anchored predicted target and the
   reference target over the CE-derived structurally equivalent residue pairs.
   NO additional target fit is performed.

Consequently, target_orientation_rmsd_A remains sensitive to complex placement
and orientation rather than simply measuring whether the isolated target folds
can be superposed.

DIAGNOSTIC RMSD: target_structure_fit_rmsd_A
---------------------------------------------
This is the target-only best-fit RMSD over the CE-derived target residue pairs.
It is a diagnostic for intrinsic target structural similarity.

Useful interpretation:
  low target_structure_fit_rmsd_A + low target_orientation_rmsd_A
      -> similar target fold and similar complex orientation

  low target_structure_fit_rmsd_A + high target_orientation_rmsd_A
      -> similar target fold but different target placement/orientation

CE correspondence settings
---------------------------
Preferred CE window:
  8

CE max gap:
  30

Target structural-fit distance cutoff:
  5.0 Å

Minimum retained target CA pairs:
  8

Main outputs
------------
all_model_metrics.tsv
  One row per analyzed AF3 model.

target_group_summary.tsv
  Per-target descriptive statistics.

top10_highest_iptm_all_pairs.tsv
  Global top 10 highest-ipTM models across every binder-target pairing.

top10_lowest_orientation_rmsd_all_pairs.tsv
  Global top 10 lowest target-orientation-RMSD models across every
  binder-target pairing.

plots/iptm_distribution_by_target.png
plots/target_orientation_rmsd_distribution_by_target.png
plots/orientation_rmsd_vs_iptm_by_target.png

The scatter plot uses:
  green  = target orientation RMSD < 2 Å and ipTM > 0.8
  orange = not green, but RMSD < 5 Å and ipTM > 0.6
  red    = all other predictions

analysis_failures.tsv
  Any skipped runs/models and the corresponding error.

Run summary
-----------
Successfully analyzed models: 208
Target groups: 13
Targets: target1, target2, target3, target4, target5, target6, target7, target8, target9, target10, target11, target12, target13
Failures/skips: 0

Statistical-unit note
---------------------
With --model-mode top (default), each binder-target run contributes one AF3
top-ranked model, avoiding seed/sample predictions being treated as independent
binder observations.
