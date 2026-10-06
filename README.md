# AF3_Cofolding_Design_Analysis
Useful scripts for analyzing AF3 Co-folding outputs for Protein Design

## Sections

1. Using af3_target_group_analysis_v2.py - ipTM extraction, RMSD calculation, and plotting


## 1 - af3_target_group_analysis_v2.py

The following script is used to calculate and obtain information about AF3 Co-folding outputs, focusing on evaluating differences between target groups (e.g. select target groups of interest are co-folded with a variety of binder/inhibitors) compared to a user defined reference complex. Specifically,

- Extract ipTM values from complete af3 co-folding outputs and plot distributions between target groups
- Calculate binder-aligned target orientated C-alpha RMSD and plot distributions between target groups
- Plot ipTM values (X-axis) and RMSD (Y-axis) scatter plot for individual target groups
- Tabulate top 10 highest ipTM values from outputs
- Tabulate top 10 lowest RMSD values from outputs

### Setup

HPCC setup - Prepare a conda environment in HPCC:

```Bash
module purge

module load Miniforge3/25.11.0-1

conda create --name <env_name> conda-forge numpy pandas matplotlib biopython

conda activate <env_name>
```

Where <env_name> is your unique conda environment name.

### Basic Command and Input Structure

```Bash
python af3_target_group_analysis_v2.py \
    --input-dir /path/to/af3_outputs \ 
    --reference /path/to/reference_complex.cif \ # Only one reference complex is supported per analysis
    --outdir /path/to/af3_target_analysis
```

For `--input-dir`, the expected input directory layout should follow something like:

```
/path/to/af3_outputs/
        |---/binder1_target1
            |---/binder1_target1
                |---binder1_target1_summary_confidences.json
                |---binder1_target1_model.cif
        |---/binder1_target2
        |---/binder2_target1
        |---/binder2_target2
```

For `--reference`, choose one reference complex that encompasses a desired binding pose.

### Script Details

  1. Extracts ipTM from *_summary_confidences.json
  2. Aligns the predicted binder onto the binder in a user-supplied reference complex
  3. Structurally aligns predicted/reference targets with CE solely to identify
     structurally equivalent C-alpha residue pairs (sequence identity is not required)
  4. Applies the binder-derived rigid-body transform to the predicted target
  5. Calculates binder-anchored target orientation RMSD over the CE-derived
     target residue correspondences WITHOUT target-specific refitting
  6. Generates target-group ipTM and orientation-RMSD distributions
  7. Generates faceted RMSD-vs-ipTM plots with threshold-based point coloring
  8. Writes two global top-N tables: highest ipTM and lowest orientation RMSD

It is expected that your predicted model and reference chains are organized as such:

  - predicted chain A = binder
  - predicted chain B = target
  - reference chain A = binder
  - reference chain B = target

RMSD interpretation: target_orientation_rmsd_A
    
The target chain is structurally aligned with Biopython's Combinatorial
Extension (CE) algorithm ONLY to determine which target residues are
structurally equivalent. The CE target superposition is not used for the
reported orientation RMSD. Instead, the predicted binder is aligned onto
the reference binder, that binder-derived transform is applied to the
predicted target, and RMSD is measured over the CE-derived target residue
pairs. This makes the primary RMSD sensitive to complex orientation while
remaining robust to low target sequence identity.

Diagnostic RMSD: target_stucture_fit_rmsd_A
  
This is the target-only best-fit RMSD over the CE-derived target residue pairs.
It is a diagnostic for intrinsic target structural similarity

Useful interpretations for RMSD results:
  
Low target_structure_fit_rmsd_A + low target_orientation_rmsd_A
    -> similar target fold and similar complex orientation

Low target_structure_fit_rmsd_A + high target_orientation_rmsd_A
    -> similar target fold but different target placement/orientation

Both target_orientation_rmsd_A and target_structure_fit_rmsd_A are reported in `target_group_summary.tsv`. target_orientation_rmsd_A is the RMSD value used for plotting and recording in the global ranking table.

### Additional Arguments

| Argument | Options | Description | Required? |
|---|---|---|---|
| `--pred-binder-chain` | alphabetical | Default= `A` , Chain letter associated with prediction binder | OPTIONAL |
| `--pred-target-chain` | alphabetical | Default= `B` , Chain letter associated with prediction target | OPTIONAL |
| `--ref-binder-chain` | alphabetical | Default= `A` , Chain letter associated with reference binder | OPTIONAL |
| `--ref-target-chain` | alphabetical | Default= `B` , Chain letter associated with reference target | OPTIONAL |
| `--group-regex` | "^(?P<binder>binder\d+)_(?P<target>target\d+)$" | change `binder` and/or `target` outside angle brackets to the actual name of your output files, e.g. "^(?P<binder>inhibitor\d+)_(?P<target>VDAG_0000\d+)$" | OPTIONAL |
| `--model-mode` | `top`, `samples`, `all` | Default= `top` , top = one top-ranked/root AF3 model per binder-target run, samples = all seed/sample models, all = both | OPTIONAL |
| `--top-n` | number | Default=10, Number of models retained in each global ranking table | OPTIONAL |
| `--distribution-kind` | `violin`, `box` | Default= `violin` , Distribution plot style for ipTM and orientation RMSD | OPTIONAL |
| `--min-binder-ca` | number | Default=3, Minimum aligned binder C-alpha pairs needed for superposition | OPTIONAL |
| `--min-target-ca` | number | Default=8, Minimum CE-derived target C-alpha pairs needed for RMSD | OPTIONAL |
| `--ce-window-size` | number | Default=8, Preferred CE structural-alignment fragment window. It is automatically reduced for short target chains | OPTIONAL |
| `--ce-max-gap` | number | Default=30, Maximum CE structural-alignment gap | OPTIONAL |
| `--target-fit-distance-cutoff` | number w/ decimal | Default=5.0, After target-only CE superposition, retain structurally corresponding C-alpha pairs within this distance (Angstrom) for orientation-RMSD evaluation, Set <= 0 to disable distance filtering | OPTIONAL |
| `--figure-dpi` | number | Default=300, PNG resolution | OPTIONAL |

### Outputs

```
analysis/
├── all_model_metrics.tsv
├── target_group_summary.tsv
├── top10_highest_iptm_all_pairs.tsv
├── top10_lowest_orientation_rmsd_all_pairs.tsv
├── analysis_failures.tsv
├── README_analysis.txt
│
└── plots/
    ├── iptm_distribution_by_target.png
    ├── target_orientation_rmsd_distribution_by_target.png
    └── orientation_rmsd_vs_iptm_by_target.png
```

See `comparison_test` for example outputs.
