# Binning systematic scans

This folder contains shell scripts for running binning systematic scans with `Unfolding/Machine.C`.

Run all scripts from the main project directory, i.e. from the directory containing:

```bash
config.h
config_hist.h
Unfolding/Machine.C
Data/
```

Example:

```bash
cd D0JetAngularities200GeVUnfolding
bash Systematics/Binning/01_BinningJetPt.sh
```

Each script creates its own scan directory inside `Systematics/Binning/`, for example `scanJetPt`, `scanZ`, or `scanL11`.

Each scan directory contains:

```bash
overrides/       # generated override macros for individual runs
runs/            # ROOT logs for individual runs
summary.tsv      # list of tested binning configurations
stability.tsv    # stability output written by Machine.C
Output/          # ROOT output files written by Machine.C
```

The first real `Machine.C` call is run with `++` to force ACLiC recompilation. All following runs use `+`.

## Scripts

| Script | Observable scanned | What is changed |
|---|---|---|
| `01_BinningJetPt.sh` | jet transverse momentum, `p_{T,jet}` | Varies reco-level and true-level jet-pT binning. Other observable binnings stay unchanged. |
| `02_BinningZ.sh` | momentum sharing, `z` | Varies only `zRecoBinsVec` and `zMcBinsVecCustom`. Includes controlled coherent edge nudging by `0.05`. |
| `03_BinningL11.sh` | `lambda_1^1` / `l11` | Varies only angularity binning with index `0`: `angRecoBinsVec[*][0]` and `angMcBinsVecCustom[*][0]`. |
| `04_BinningL11p5.sh` | `lambda_{1.5}^1` / `l11p5` | Varies only angularity binning with index `1`: `angRecoBinsVec[*][1]` and `angMcBinsVecCustom[*][1]`. |
| `05_BinningL12.sh` | `lambda_2^1` / `l12` | Varies only angularity binning with index `2`: `angRecoBinsVec[*][2]` and `angMcBinsVecCustom[*][2]`. |
| `06_BinningL13.sh` | `lambda_3^1` / `l13` | Varies only angularity binning with index `3`: `angRecoBinsVec[*][3]` and `angMcBinsVecCustom[*][3]`. |
| `07_BinningL10p5.sh` | `lambda_{0.5}^1` / `l10p5` | Varies only angularity binning with index `4`: `angRecoBinsVec[*][4]` and `angMcBinsVecCustom[*][4]`. |
| `08_BinningpTD.sh` | momentum dispersion, `p_T^D` | Varies only angularity/binning variable with the pTD index. Other binnings should stay unchanged. |

## Important notes

- The override macros should change only the observable scanned by the given script.
- All other binnings should remain as defined in `config.h` / `config_hist.h`.
- The scripts are intended for nominal binning systematics, so they use:

```bash
USE_PRIOR_SHAPE_WEIGHTING=0
SYSTEMATIC_SPLOT=0
```

- The current `Machine.C` signature used by these scripts is:

```cpp
Machine(fonllJet, CutOfNegative, minJetPtRecoCut, savedIter,
        InputFile, OutputFile, minPtD0Cut, maxPtD0Cut,
        OverrideMacro, ScanDir, usePriorShapeWeighting, systematicSPlot)
```

