#!/usr/bin/env python3
"""
AF3 binder-target group analysis

For AlphaFold 3 co-folding output directories named like:
    binder1_target1/
    binder2_target1/
    binder3_target2/
    ...

the script:
  1. extracts ipTM from *_summary_confidences.json
  2. aligns the predicted binder onto the binder in a user-supplied reference complex
  3. structurally aligns predicted/reference targets with CE solely to identify
     structurally equivalent C-alpha residue pairs (sequence identity is not required)
  4. applies the binder-derived rigid-body transform to the predicted target
  5. calculates binder-anchored target orientation RMSD over the CE-derived
     target residue correspondences WITHOUT target-specific refitting
  6. generates target-group ipTM and orientation-RMSD distributions
  7. generates faceted RMSD-vs-ipTM plots with threshold-based point coloring
  8. writes two global top-N tables: highest ipTM and lowest orientation RMSD

Default interpretation:
  - predicted chain A = binder
  - predicted chain B = target
  - reference chain A = binder
  - reference chain B = target
  - one top-ranked AF3 model per binder-target run (--model-mode top)

Dependencies:
    numpy pandas matplotlib biopython

RMSD interpretation:
    The target chain is structurally aligned with Biopython's Combinatorial
    Extension (CE) algorithm ONLY to determine which target residues are
    structurally equivalent. The CE target superposition is not used for the
    reported orientation RMSD. Instead, the predicted binder is aligned onto
    the reference binder, that binder-derived transform is applied to the
    predicted target, and RMSD is measured over the CE-derived target residue
    pairs. This makes the primary RMSD sensitive to complex orientation while
    remaining robust to low target sequence identity.

Example:
    python af3_target_group_analysis.py \
      --input-dir /path/to/af3_outputs \
      --reference /path/to/reference_complex.cif \
      --outdir af3_target_analysis

If the reference uses different chain IDs:
    python af3_target_group_analysis.py \
      --input-dir /path/to/af3_outputs \
      --reference /path/to/reference_complex.pdb \
      --ref-binder-chain C \
      --ref-target-chain D \
      --outdir af3_target_analysis
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.cealign import run_cealign
from Bio.SVDSuperimposer import SVDSuperimposer
from matplotlib.lines import Line2D


SUMMARY_SUFFIX = "_summary_confidences.json"
MODEL_SUFFIX = "_model.cif"
SAMPLE_DIR_RE = re.compile(r"^seed-.+_sample-\d+$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Analyze AF3 binder-target co-folding runs by target group: ipTM, "
            "binder-anchored target orientation C-alpha RMSD, distributions, "
            "threshold-colored RMSD-vs-ipTM plots, and global top-model tables."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("--input-dir", required=True, type=Path,
                   help="Parent directory containing binder#_target# run directories.")
    p.add_argument("--reference", required=True, type=Path,
                   help="Reference complex structure (.cif/.mmcif/.pdb).")
    p.add_argument("--outdir", required=True, type=Path,
                   help="Output analysis directory.")

    p.add_argument("--pred-binder-chain", default="A",
                   help="Binder chain ID in predicted AF3 models.")
    p.add_argument("--pred-target-chain", default="B",
                   help="Target chain ID in predicted AF3 models.")
    p.add_argument("--ref-binder-chain", default="A",
                   help="Binder chain ID in reference complex.")
    p.add_argument("--ref-target-chain", default="B",
                   help="Target chain ID in reference complex.")

    p.add_argument(
        "--group-regex",
        default=r"^(?P<binder>binder\d+)_(?P<target>target\d+)$",
        help=(
            "Regex applied to run-directory basename. Must contain named groups "
            "'binder' and 'target'."
        ),
    )

    p.add_argument(
        "--model-mode",
        choices=["top", "samples", "all"],
        default="top",
        help=(
            "top = one top-ranked/root AF3 model per binder-target run; "
            "samples = all seed/sample models; all = both."
        ),
    )
    p.add_argument("--top-n", type=int, default=10,
                   help="Number of models retained in each global ranking table.")
    p.add_argument(
        "--distribution-kind",
        choices=["violin", "box"],
        default="violin",
        help="Distribution plot style for ipTM and orientation RMSD.",
    )

    p.add_argument("--min-binder-ca", type=int, default=3,
                   help="Minimum aligned binder C-alpha pairs needed for superposition.")
    p.add_argument("--min-target-ca", type=int, default=8,
                   help="Minimum CE-derived target C-alpha pairs needed for RMSD.")
    p.add_argument(
        "--ce-window-size",
        type=int,
        default=8,
        help=(
            "Preferred CE structural-alignment fragment window. It is automatically "
            "reduced for short target chains."
        ),
    )
    p.add_argument(
        "--ce-max-gap",
        type=int,
        default=30,
        help="Maximum CE structural-alignment gap."
    )
    p.add_argument(
        "--target-fit-distance-cutoff",
        type=float,
        default=5.0,
        help=(
            "After target-only CE superposition, retain structurally corresponding "
            "C-alpha pairs within this distance (Angstrom) for orientation-RMSD "
            "evaluation. Set <= 0 to disable distance filtering."
        ),
    )
    p.add_argument(
        "--figure-dpi",
        type=int,
        default=300,
        help="PNG resolution.",
    )

    args = p.parse_args()

    if args.top_n < 1:
        p.error("--top-n must be >= 1")
    if args.min_binder_ca < 3:
        p.error("--min-binder-ca must be >= 3")
    if args.min_target_ca < 3:
        p.error("--min-target-ca must be >= 3")
    if args.ce_window_size < 2:
        p.error("--ce-window-size must be >= 2")
    if args.ce_max_gap < 0:
        p.error("--ce-max-gap must be >= 0")

    return args


def compile_group_regex(pattern: str) -> re.Pattern:
    rx = re.compile(pattern)
    required = {"binder", "target"}
    missing = required - set(rx.groupindex)
    if missing:
        raise ValueError(
            f"--group-regex must contain named capture groups {sorted(required)}; "
            f"missing {sorted(missing)}"
        )
    return rx


def load_structure(path: Path):
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif suffix in {".pdb", ".ent"}:
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(
            f"Unsupported structure format for {path}. Use .cif/.mmcif/.pdb."
        )
    return parser.get_structure(path.stem, str(path))


def get_model0(structure):
    models = list(structure.get_models())
    if not models:
        raise ValueError("Structure contains no coordinate models.")
    return models[0]


def residue_to_one_letter(resname: str) -> str:
    key = resname.strip().upper()
    return protein_letters_3to1_extended.get(key, "X")


def extract_chain_ca(structure, chain_id: str) -> Tuple[str, np.ndarray, List[Tuple]]:
    """
    Return:
      sequence
      C-alpha coordinates (N x 3)
      residue identifiers
    Only residues containing a CA atom are retained.
    """
    model = get_model0(structure)

    if chain_id not in model:
        available = [c.id for c in model.get_chains()]
        raise KeyError(
            f"Chain {chain_id!r} not found. Available chains: {available}"
        )

    chain = model[chain_id]
    sequence = []
    coords = []
    residue_ids = []

    for residue in chain.get_residues():
        if "CA" not in residue:
            continue
        sequence.append(residue_to_one_letter(residue.get_resname()))
        coords.append(residue["CA"].get_coord().astype(float))
        residue_ids.append(residue.id)

    if not coords:
        raise ValueError(f"No C-alpha atoms found in chain {chain_id!r}")

    return "".join(sequence), np.asarray(coords, dtype=float), residue_ids


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


ALIGNER = make_aligner()


def sequence_correspondence(
    reference_seq: str,
    moving_seq: str,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Global sequence alignment, returning corresponding ungapped residue indices.

    The aligned coordinate blocks returned by Biopython contain positions that are
    directly paired. Mismatched amino acids are kept as structural correspondences.
    """
    if not reference_seq or not moving_seq:
        raise ValueError("Cannot align empty sequences.")

    alignment = ALIGNER.align(reference_seq, moving_seq)[0]
    ref_blocks, mov_blocks = alignment.aligned

    ref_idx: List[int] = []
    mov_idx: List[int] = []

    for (r0, r1), (m0, m1) in zip(ref_blocks, mov_blocks):
        r_len = int(r1 - r0)
        m_len = int(m1 - m0)
        if r_len != m_len:
            # This should not occur for ungapped aligned blocks, but guard anyway.
            n = min(r_len, m_len)
        else:
            n = r_len
        ref_idx.extend(range(int(r0), int(r0) + n))
        mov_idx.extend(range(int(m0), int(m0) + n))

    if not ref_idx:
        raise ValueError("Sequence alignment produced no paired residues.")

    ref_arr = np.asarray(ref_idx, dtype=int)
    mov_arr = np.asarray(mov_idx, dtype=int)

    matches = sum(
        reference_seq[i] == moving_seq[j]
        for i, j in zip(ref_arr.tolist(), mov_arr.tolist())
    )
    identity = matches / len(ref_arr)
    return ref_arr, mov_arr, identity


def structural_correspondence_ce(
    reference_seq: str,
    reference_ca: np.ndarray,
    moving_seq: str,
    moving_ca: np.ndarray,
    preferred_window_size: int,
    max_gap: int,
    distance_cutoff: float,
    min_pairs: int,
) -> Dict[str, object]:
    """
    Identify target residue correspondence by sequence-independent CE alignment.

    Important:
      The target-only CE superposition is used ONLY to establish structurally
      equivalent C-alpha residue pairs and to report a diagnostic target-fold
      fit RMSD. Its rotation/translation is NOT used for the primary complex
      orientation RMSD.
    """
    min_len = min(len(reference_ca), len(moving_ca))
    if min_len < 4:
        raise ValueError(
            f"Target chains are too short for structural correspondence: {min_len} CA atoms."
        )

    # CE needs enough residues for at least two windows. Adapt for short chains.
    window_size = min(preferred_window_size, max(2, min_len // 2))
    while window_size >= 2 and min_len < 2 * window_size:
        window_size -= 1
    if window_size < 2:
        raise ValueError(
            f"Could not choose a valid CE window for target length {min_len}."
        )

    alignments = run_cealign(
        reference_ca.tolist(),
        moving_ca.tolist(),
        window_size,
        max_gap,
    )
    if not alignments:
        raise ValueError("CE structural alignment produced no target alignment.")

    # Match CEAligner's own selection strategy: among longest CE paths, retain
    # the path with the lowest target-only best-fit RMSD.
    max_length = max(aln.length for aln in alignments)
    longest = [aln for aln in alignments if aln.length == max_length]

    best = None
    best_fit_rmsd = float("inf")
    best_rot = None
    best_tran = None
    best_ref_idx = None
    best_mov_idx = None

    for aln in longest:
        ref_idx = np.asarray(aln.path[0], dtype=int)
        mov_idx = np.asarray(aln.path[1], dtype=int)
        if len(ref_idx) < 3:
            continue

        sup = SVDSuperimposer()
        sup.set(reference_ca[ref_idx], moving_ca[mov_idx])
        sup.run()
        rms = float(sup.get_rms())

        if rms < best_fit_rmsd:
            best = aln
            best_fit_rmsd = rms
            best_rot, best_tran = sup.get_rotran()
            best_ref_idx = ref_idx
            best_mov_idx = mov_idx

    if best is None or best_ref_idx is None or best_mov_idx is None:
        raise ValueError("CE did not yield a usable target structural correspondence.")

    # Identify the structurally conserved portion after target-only best fit.
    moving_fit = moving_ca[best_mov_idx] @ best_rot + best_tran
    pair_distances = np.linalg.norm(
        moving_fit - reference_ca[best_ref_idx],
        axis=1,
    )

    if distance_cutoff > 0:
        keep = pair_distances <= distance_cutoff
    else:
        keep = np.ones(len(pair_distances), dtype=bool)

    ref_idx = best_ref_idx[keep]
    mov_idx = best_mov_idx[keep]

    if len(ref_idx) < min_pairs:
        raise ValueError(
            f"Only {len(ref_idx)} structurally corresponding target CA pairs remain "
            f"after CE/cutoff filtering; need >= {min_pairs}. "
            f"Try increasing --target-fit-distance-cutoff if appropriate."
        )

    # Refit only for a diagnostic measure of intrinsic target fold similarity.
    diagnostic_sup = SVDSuperimposer()
    diagnostic_sup.set(reference_ca[ref_idx], moving_ca[mov_idx])
    diagnostic_sup.run()
    structural_fit_rmsd = float(diagnostic_sup.get_rms())

    matches = sum(
        reference_seq[i] == moving_seq[j]
        for i, j in zip(ref_idx.tolist(), mov_idx.tolist())
    )
    structural_pair_identity = matches / len(ref_idx)

    return {
        "ref_idx": ref_idx,
        "mov_idx": mov_idx,
        "target_structure_fit_rmsd_A": structural_fit_rmsd,
        "target_correspondence_ca": int(len(ref_idx)),
        "target_ce_alignment_length": int(best.length),
        "target_ce_z_score": float(best.z_score),
        "target_structural_pair_seq_identity": float(structural_pair_identity),
        "target_ce_window_size_used": int(window_size),
    }


def binder_anchored_target_orientation_rmsd(
    pred_structure,
    ref_binder_seq: str,
    ref_binder_ca: np.ndarray,
    ref_target_seq: str,
    ref_target_ca: np.ndarray,
    pred_binder_chain: str,
    pred_target_chain: str,
    min_binder_ca: int,
    min_target_ca: int,
    ce_window_size: int,
    ce_max_gap: int,
    target_fit_distance_cutoff: float,
) -> Dict[str, float]:
    """
    Calculate target RMSD in a binder-anchored coordinate frame.

    Procedure
    ---------
    1. Match predicted/reference binder residues by global sequence alignment.
       This is appropriate because binders are expected to be identical or highly
       sequence-similar.
    2. Superpose the predicted binder C-alpha atoms onto the reference binder.
    3. Independently use CE target structural alignment ONLY to identify which
       target residues are structurally equivalent despite low target sequence
       identity.
    4. Apply ONLY the binder-derived transform to the predicted target.
    5. Compute C-alpha RMSD over the CE-derived target residue pairs with NO
       target-specific refitting.

    Therefore target_orientation_rmsd_A is sensitive to whether the predicted
    target occupies a similar position/orientation relative to the binder.
    """
    pred_binder_seq, pred_binder_ca, _ = extract_chain_ca(
        pred_structure, pred_binder_chain
    )
    pred_target_seq, pred_target_ca, _ = extract_chain_ca(
        pred_structure, pred_target_chain
    )

    # Binder anchoring: sequence correspondence is reliable for the expected
    # identical/highly similar binder sequences.
    rb_idx, pb_idx, binder_identity = sequence_correspondence(
        ref_binder_seq, pred_binder_seq
    )
    if len(rb_idx) < min_binder_ca:
        raise ValueError(
            f"Only {len(rb_idx)} aligned binder C-alpha pairs; "
            f"need >= {min_binder_ca}."
        )

    binder_sup = SVDSuperimposer()
    binder_sup.set(ref_binder_ca[rb_idx], pred_binder_ca[pb_idx])
    binder_sup.run()
    binder_rotation, binder_translation = binder_sup.get_rotran()
    binder_fit_rmsd = float(binder_sup.get_rms())

    # Target correspondence: sequence-independent structural matching.
    corr = structural_correspondence_ce(
        reference_seq=ref_target_seq,
        reference_ca=ref_target_ca,
        moving_seq=pred_target_seq,
        moving_ca=pred_target_ca,
        preferred_window_size=ce_window_size,
        max_gap=ce_max_gap,
        distance_cutoff=target_fit_distance_cutoff,
        min_pairs=min_target_ca,
    )

    rt_idx = corr["ref_idx"]
    pt_idx = corr["mov_idx"]

    # Critical step: use the BINDER transform, not the CE target transform.
    pred_target_binder_anchored = (
        pred_target_ca[pt_idx] @ binder_rotation + binder_translation
    )
    ref_target_eval = ref_target_ca[rt_idx]

    delta = pred_target_binder_anchored - ref_target_eval
    target_orientation_rmsd = float(
        np.sqrt(np.mean(np.sum(delta * delta, axis=1)))
    )

    return {
        "target_orientation_rmsd_A": target_orientation_rmsd,
        "binder_fit_rmsd_A": binder_fit_rmsd,
        "binder_aligned_ca": int(len(rb_idx)),
        "binder_seq_identity": float(binder_identity),
        "pred_binder_length_ca": int(len(pred_binder_ca)),
        "pred_target_length_ca": int(len(pred_target_ca)),
        "ref_target_length_ca": int(len(ref_target_ca)),
        "target_structure_fit_rmsd_A": float(corr["target_structure_fit_rmsd_A"]),
        "target_correspondence_ca": int(corr["target_correspondence_ca"]),
        "target_ce_alignment_length": int(corr["target_ce_alignment_length"]),
        "target_ce_z_score": float(corr["target_ce_z_score"]),
        "target_structural_pair_seq_identity": float(
            corr["target_structural_pair_seq_identity"]
        ),
        "target_ce_window_size_used": int(corr["target_ce_window_size_used"]),
    }


def matching_model_path(summary_path: Path) -> Optional[Path]:
    name = summary_path.name
    if not name.endswith(SUMMARY_SUFFIX):
        return None

    prefix = name[: -len(SUMMARY_SUFFIX)]
    exact = summary_path.with_name(prefix + MODEL_SUFFIX)
    if exact.exists():
        return exact

    # Conservative fallback within same directory.
    candidates = sorted(summary_path.parent.glob(prefix + "*_model.cif"))
    if len(candidates) == 1:
        return candidates[0]
    return None


def is_sample_summary(path: Path) -> bool:
    return bool(SAMPLE_DIR_RE.match(path.parent.name))


def choose_summary_files(run_dir: Path, mode: str) -> List[Path]:
    all_summaries = sorted(run_dir.rglob(f"*{SUMMARY_SUFFIX}"))
    if not all_summaries:
        return []

    sample_summaries = [p for p in all_summaries if is_sample_summary(p)]
    top_summaries = [p for p in all_summaries if not is_sample_summary(p)]

    if mode == "samples":
        # Fall back to top if a legacy/single-model run has no sample directories.
        return sample_summaries if sample_summaries else top_summaries

    if mode == "all":
        return all_summaries

    # mode == "top"
    if not top_summaries:
        # If only sample outputs exist, choose the first sample as a fallback.
        warnings.warn(
            f"No root/top summary found under {run_dir}; using first sample summary."
        )
        return sample_summaries[:1]

    # A wrapper directory may contain more than one non-sample summary. Prefer:
    #   1. summaries with a matching model CIF
    #   2. the shallowest relative path under the run directory
    valid = [p for p in top_summaries if matching_model_path(p) is not None]
    candidates = valid if valid else top_summaries
    candidates = sorted(
        candidates,
        key=lambda p: (len(p.relative_to(run_dir).parts), str(p)),
    )

    if len(candidates) > 1:
        warnings.warn(
            f"Multiple possible top-level summaries found under {run_dir}; "
            f"using {candidates[0]}"
        )

    return candidates[:1]


def discover_run_dirs(input_dir: Path, group_rx: re.Pattern) -> List[Tuple[Path, re.Match]]:
    hits: List[Tuple[Path, re.Match]] = []

    # Include input-dir itself in case the user points directly at a single run.
    candidates = [input_dir] + [p for p in input_dir.rglob("*") if p.is_dir()]

    for path in candidates:
        m = group_rx.fullmatch(path.name)
        if m:
            hits.append((path, m))

    # Avoid duplicate nested discovery of identically named wrappers by exact path only.
    hits = sorted({p.resolve(): (p, m) for p, m in hits}.values(), key=lambda x: str(x[0]))
    return hits


def read_summary(summary_path: Path) -> Dict:
    with summary_path.open("r") as fh:
        data = json.load(fh)

    if "iptm" not in data:
        raise KeyError(f"'iptm' not found in {summary_path}")

    return data


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def natural_key(value: str):
    parts = re.split(r"(\d+)", str(value))
    return [int(p) if p.isdigit() else p.lower() for p in parts]


def analyze(args: argparse.Namespace) -> Tuple[pd.DataFrame, pd.DataFrame]:
    group_rx = compile_group_regex(args.group_regex)

    if not args.input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {args.input_dir}")
    if not args.reference.exists():
        raise FileNotFoundError(f"Reference structure does not exist: {args.reference}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    plots_dir = args.outdir / "plots"
    plots_dir.mkdir(exist_ok=True)

    print(f"[INFO] Loading reference: {args.reference}")
    ref_structure = load_structure(args.reference)
    ref_binder_seq, ref_binder_ca, _ = extract_chain_ca(
        ref_structure, args.ref_binder_chain
    )
    ref_target_seq, ref_target_ca, _ = extract_chain_ca(
        ref_structure, args.ref_target_chain
    )

    print(
        f"[INFO] Reference binder chain {args.ref_binder_chain}: "
        f"{len(ref_binder_ca)} CA atoms"
    )
    print(
        f"[INFO] Reference target chain {args.ref_target_chain}: "
        f"{len(ref_target_ca)} CA atoms"
    )

    run_dirs = discover_run_dirs(args.input_dir, group_rx)
    if not run_dirs:
        raise RuntimeError(
            f"No run directories matched --group-regex under {args.input_dir}\n"
            f"Regex: {args.group_regex}"
        )

    print(f"[INFO] Found {len(run_dirs)} binder-target run directories.")

    rows = []
    failures = []

    for i, (run_dir, match) in enumerate(run_dirs, start=1):
        binder_group = match.group("binder")
        target_group = match.group("target")
        summaries = choose_summary_files(run_dir, args.model_mode)

        if not summaries:
            failures.append({
                "run_dir": str(run_dir),
                "summary_path": "",
                "model_path": "",
                "stage": "discovery",
                "error": "No *_summary_confidences.json found.",
            })
            print(f"[WARN] [{i}/{len(run_dirs)}] No summary JSON: {run_dir}")
            continue

        for summary_path in summaries:
            model_path = matching_model_path(summary_path)
            if model_path is None:
                failures.append({
                    "run_dir": str(run_dir),
                    "summary_path": str(summary_path),
                    "model_path": "",
                    "stage": "model_pairing",
                    "error": "Could not find matching *_model.cif.",
                })
                print(f"[WARN] Missing model CIF for {summary_path}")
                continue

            try:
                summary = read_summary(summary_path)
                iptm = safe_float(summary.get("iptm"))
                ranking_score = safe_float(summary.get("ranking_score"))
                ptm = safe_float(summary.get("ptm"))

                pred_structure = load_structure(model_path)
                rmsd_info = binder_anchored_target_orientation_rmsd(
                    pred_structure=pred_structure,
                    ref_binder_seq=ref_binder_seq,
                    ref_binder_ca=ref_binder_ca,
                    ref_target_seq=ref_target_seq,
                    ref_target_ca=ref_target_ca,
                    pred_binder_chain=args.pred_binder_chain,
                    pred_target_chain=args.pred_target_chain,
                    min_binder_ca=args.min_binder_ca,
                    min_target_ca=args.min_target_ca,
                    ce_window_size=args.ce_window_size,
                    ce_max_gap=args.ce_max_gap,
                    target_fit_distance_cutoff=args.target_fit_distance_cutoff,
                )

                try:
                    rel_model = model_path.relative_to(args.input_dir)
                except ValueError:
                    rel_model = model_path

                model_id = str(rel_model)
                seed_sample = (
                    summary_path.parent.name
                    if is_sample_summary(summary_path)
                    else "top_ranked"
                )

                row = {
                    "target_group": target_group,
                    "binder_group": binder_group,
                    "run_name": run_dir.name,
                    "model_id": model_id,
                    "seed_sample": seed_sample,
                    "iptm": iptm,
                    "ptm": ptm,
                    "ranking_score": ranking_score,
                    **rmsd_info,
                    "summary_path": str(summary_path.resolve()),
                    "model_path": str(model_path.resolve()),
                }
                rows.append(row)

            except Exception as exc:
                failures.append({
                    "run_dir": str(run_dir),
                    "summary_path": str(summary_path),
                    "model_path": str(model_path),
                    "stage": "analysis",
                    "error": f"{type(exc).__name__}: {exc}",
                })
                print(f"[WARN] Failed {model_path}: {exc}", file=sys.stderr)

    if not rows:
        raise RuntimeError(
            "No models were successfully analyzed. Check analysis_failures.tsv "
            "and chain IDs / directory layout."
        )

    df = pd.DataFrame(rows)
    fail_df = pd.DataFrame(failures)

    # Natural target/binder ordering.
    target_order = sorted(df["target_group"].unique(), key=natural_key)
    target_rank = {v: i for i, v in enumerate(target_order)}
    df["_target_order"] = df["target_group"].map(target_rank)
    # Numeric-aware binder ordering without requiring an external natsort package.
    df["_binder_order"] = df["binder_group"].map(
        lambda x: re.sub(r"(\d+)", lambda m: f"{int(m.group(1)):012d}", str(x).lower())
    )
    df = df.sort_values(
        ["_target_order", "_binder_order", "iptm"],
        ascending=[True, True, False],
    ).drop(columns=["_target_order", "_binder_order"])

    all_metrics_path = args.outdir / "all_model_metrics.tsv"
    df.to_csv(all_metrics_path, sep="\t", index=False)

    if not fail_df.empty:
        fail_df.to_csv(args.outdir / "analysis_failures.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(
            columns=["run_dir", "summary_path", "model_path", "stage", "error"]
        ).to_csv(args.outdir / "analysis_failures.tsv", sep="\t", index=False)

    # Per-target descriptive summaries.
    summary_rows = []

    for target in target_order:
        sub = df[df["target_group"] == target].copy()

        summary_rows.append({
            "target_group": target,
            "n_models": len(sub),
            "n_unique_binders": sub["binder_group"].nunique(),
            "iptm_mean": sub["iptm"].mean(),
            "iptm_median": sub["iptm"].median(),
            "iptm_std": sub["iptm"].std(ddof=1),
            "iptm_min": sub["iptm"].min(),
            "iptm_max": sub["iptm"].max(),
            "target_orientation_rmsd_A_mean": sub["target_orientation_rmsd_A"].mean(),
            "target_orientation_rmsd_A_median": sub["target_orientation_rmsd_A"].median(),
            "target_orientation_rmsd_A_std": sub["target_orientation_rmsd_A"].std(ddof=1),
            "target_orientation_rmsd_A_min": sub["target_orientation_rmsd_A"].min(),
            "target_orientation_rmsd_A_max": sub["target_orientation_rmsd_A"].max(),
            "target_structure_fit_rmsd_A_mean": sub["target_structure_fit_rmsd_A"].mean(),
            "target_structure_fit_rmsd_A_median": sub["target_structure_fit_rmsd_A"].median(),
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(
        args.outdir / "target_group_summary.tsv",
        sep="\t",
        index=False,
    )

    # Exactly two GLOBAL top-N ranking tables across every binder-target pairing.
    top_iptm = (
        df.dropna(subset=["iptm"])
        .sort_values(
            ["iptm", "target_orientation_rmsd_A"],
            ascending=[False, True],
        )
        .head(args.top_n)
        .copy()
    )
    top_iptm.insert(0, "rank", range(1, len(top_iptm) + 1))
    top_iptm.to_csv(
        args.outdir / f"top{args.top_n}_highest_iptm_all_pairs.tsv",
        sep="\t",
        index=False,
    )

    top_rmsd = (
        df.dropna(subset=["target_orientation_rmsd_A"])
        .sort_values(
            ["target_orientation_rmsd_A", "iptm"],
            ascending=[True, False],
        )
        .head(args.top_n)
        .copy()
    )
    top_rmsd.insert(0, "rank", range(1, len(top_rmsd) + 1))
    top_rmsd.to_csv(
        args.outdir / f"top{args.top_n}_lowest_orientation_rmsd_all_pairs.tsv",
        sep="\t",
        index=False,
    )

    return df, fail_df


def format_target_axis(ax, labels: Sequence[str]):
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(axis="y", alpha=0.25)


def distribution_plot(
    df: pd.DataFrame,
    value_col: str,
    ylabel: str,
    title: str,
    output_path: Path,
    kind: str,
    dpi: int,
):
    labels = sorted(df["target_group"].dropna().unique(), key=natural_key)
    grouped = [
        df.loc[df["target_group"] == label, value_col].dropna().to_numpy(dtype=float)
        for label in labels
    ]

    usable = [(lab, vals) for lab, vals in zip(labels, grouped) if len(vals) > 0]
    if not usable:
        print(f"[WARN] No finite data for plot: {value_col}")
        return

    labels = [x[0] for x in usable]
    grouped = [x[1] for x in usable]

    width = max(7.0, 0.8 * len(labels) + 2.0)
    fig, ax = plt.subplots(figsize=(width, 6.0))

    positions = np.arange(1, len(labels) + 1)

    if kind == "violin":
        vp = ax.violinplot(
            grouped,
            positions=positions,
            showmeans=False,
            showmedians=True,
            showextrema=True,
        )
        # Keep matplotlib defaults; just make bodies slightly transparent.
        for body in vp["bodies"]:
            body.set_alpha(0.65)
    else:
        ax.boxplot(grouped, positions=positions, showfliers=True)

    format_target_axis(ax, labels)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Target group")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def scatter_facets(
    df: pd.DataFrame,
    output_path: Path,
    dpi: int,
):
    """
    Plot orientation RMSD on x and ipTM on y, faceted by target.

    Point classes are applied in priority order:
      green  : RMSD < 2 A AND ipTM > 0.8
      orange : otherwise, RMSD < 5 A AND ipTM > 0.6
      red    : all remaining points
    """
    labels = sorted(df["target_group"].dropna().unique(), key=natural_key)
    n = len(labels)
    if n == 0:
        return

    ncols = min(3, n)
    nrows = math.ceil(n / ncols)
    fig_w = 5.0 * ncols
    fig_h = 4.2 * nrows

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(fig_w, fig_h),
        squeeze=False,
        sharex=True,
        sharey=True,
    )

    all_axes = axes.ravel()

    for ax, target in zip(all_axes, labels):
        sub = df[df["target_group"] == target].dropna(
            subset=["iptm", "target_orientation_rmsd_A"]
        ).copy()

        green = (
            (sub["target_orientation_rmsd_A"] < 2.0)
            & (sub["iptm"] > 0.8)
        )
        orange = (
            ~green
            & (sub["target_orientation_rmsd_A"] < 5.0)
            & (sub["iptm"] > 0.6)
        )

        point_colors = np.full(len(sub), "red", dtype=object)
        point_colors[orange.to_numpy()] = "orange"
        point_colors[green.to_numpy()] = "green"

        ax.scatter(
            sub["target_orientation_rmsd_A"],
            sub["iptm"],
            c=point_colors,
            alpha=0.8,
            s=34,
        )
        ax.set_title(f"{target} (n={len(sub)})")
        ax.set_xlabel("Binder-anchored target orientation Cα RMSD (Å)")
        ax.set_ylabel("ipTM")
        ax.grid(alpha=0.25)

    for ax in all_axes[n:]:
        ax.axis("off")

    legend_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="None", color="green",
            label="RMSD < 2 Å and ipTM > 0.8",
        ),
        Line2D(
            [0], [0], marker="o", linestyle="None", color="orange",
            label="RMSD < 5 Å and ipTM > 0.6",
        ),
        Line2D(
            [0], [0], marker="o", linestyle="None", color="red",
            label="All other predictions",
        ),
    ]

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=3,
        frameon=False,
    )
    fig.suptitle(
        "AF3 binder-anchored target orientation RMSD vs ipTM by target",
        y=1.075,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_readme(args: argparse.Namespace, df: pd.DataFrame, fail_df: pd.DataFrame):
    targets = sorted(df["target_group"].unique(), key=natural_key)

    text = f"""AF3 TARGET-GROUP ANALYSIS v2
============================

Input directory:
  {args.input_dir.resolve()}

Reference structure:
  {args.reference.resolve()}

Predicted chains:
  binder = {args.pred_binder_chain}
  target = {args.pred_target_chain}

Reference chains:
  binder = {args.ref_binder_chain}
  target = {args.ref_target_chain}

Model mode:
  {args.model_mode}

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
  {args.ce_window_size}

CE max gap:
  {args.ce_max_gap}

Target structural-fit distance cutoff:
  {args.target_fit_distance_cutoff} Å

Minimum retained target CA pairs:
  {args.min_target_ca}

Main outputs
------------
all_model_metrics.tsv
  One row per analyzed AF3 model.

target_group_summary.tsv
  Per-target descriptive statistics.

top{args.top_n}_highest_iptm_all_pairs.tsv
  Global top {args.top_n} highest-ipTM models across every binder-target pairing.

top{args.top_n}_lowest_orientation_rmsd_all_pairs.tsv
  Global top {args.top_n} lowest target-orientation-RMSD models across every
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
Successfully analyzed models: {len(df)}
Target groups: {len(targets)}
Targets: {", ".join(targets)}
Failures/skips: {len(fail_df)}

Statistical-unit note
---------------------
With --model-mode top (default), each binder-target run contributes one AF3
top-ranked model, avoiding seed/sample predictions being treated as independent
binder observations.
"""
    (args.outdir / "README_analysis.txt").write_text(text)


def main():
    args = parse_args()

    try:
        df, fail_df = analyze(args)

        plots_dir = args.outdir / "plots"
        distribution_plot(
            df=df,
            value_col="iptm",
            ylabel="ipTM",
            title="AF3 ipTM distribution by target group",
            output_path=plots_dir / "iptm_distribution_by_target.png",
            kind=args.distribution_kind,
            dpi=args.figure_dpi,
        )
        distribution_plot(
            df=df,
            value_col="target_orientation_rmsd_A",
            ylabel="Binder-anchored target orientation Cα RMSD (Å)",
            title="Target orientation Cα RMSD distribution by target group",
            output_path=plots_dir / "target_orientation_rmsd_distribution_by_target.png",
            kind=args.distribution_kind,
            dpi=args.figure_dpi,
        )
        scatter_facets(
            df=df,
            output_path=plots_dir / "orientation_rmsd_vs_iptm_by_target.png",
            dpi=args.figure_dpi,
        )
        write_readme(args, df, fail_df)

        print()
        print("[DONE] Analysis complete.")
        print(f"[DONE] Output directory: {args.outdir.resolve()}")
        print(f"[DONE] Successful models: {len(df)}")
        print(f"[DONE] Failures/skips: {len(fail_df)}")

    except Exception as exc:
        print(f"[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
