# TDIS-TPCproto_analysis

TPC cosmic-test inspection and analysis with two maintained implementations:
Python/PyROOT and native C++/ROOT. Both use the same pad mapping and manual-cut database.
The framework provides 1D/2D histograms, an event-display prototype, run-specific cleaning,
and an independent time-versus-voltage study. Event visualization and spatial reconstruction
remain under development; the current display is not a completed track reconstruction.

## Start with the included example

A compact **80 KB real-data ROOT file** is included at `data/root/ntuple1616.root`.
It retains complete selected events from run 1616, including event 5584 and the
run's known automatically vetoed event IDs. It is a subset, not a full-exposure run.
Selection, counts and SHA-256 hashes are in [data/example_provenance.json](data/example_provenance.json).

From the repository root, with ROOT installed:

```bash
root -l 'analysis_root/TPC.C("raw",1616)'            # raw 1D and 2D, windows + exports
root -l 'analysis_root/TPC.C("raw --1d",1616)'       # only 1D
root -l 'analysis_root/TPC.C("raw --2d",1616)'       # only 2D
root -l 'analysis_root/TPC.C("cuts",1616)'           # list manual rules, no plots
root -l 'analysis_root/TPC.C("raw-clean",1616)'      # apply DB + existing event vetoes
root -l 'analysis_root/TPC.C("event",1616,5584)'     # current event-display prototype
root -l -b -q 'analysis_root/TPC.C("save",1616,5584)' # file-only event export
```

For Python, use a Python environment with PyROOT from the ROOT installation,
then install `requirements.txt`:

```bash
python -m pip install -r requirements.txt
python analysis/tpc.py raw 1616
python analysis/tpc.py raw 1616 --2d --batch
python analysis/tpc.py raw 1616 --apply-cuts
python analysis/tpc.py event 1616 --event 5584 --save
```

**Raw plots deliberately include noisy events.** Mapping always applies, but raw
inspection bypasses manual cuts and automatic event vetoes. Cleaned inspection and
event displays apply manual DB rules and IDs from an existing `(run,event)` veto
ledger. They do not rerun noise detection. The independent `voltage` workflow
classifies the available listed runs afresh and writes the ledger. See the manual
for the exact thresholds, ledger lookup, and limitations for new datasets.

## Documentation and configuration

- [Illustrated PDF manual](docs/TPC_manual.pdf): command/option reference, cleaning,
  noise-veto implementation, current settings, examples and plots.
- [LaTeX source](docs/TPC_manual.tex).
- [Developer and documentation workflow](CONTRIBUTING.md).
- `config/run_cuts.ini`: global and per-run detector cleaning.
- `config/voltage_study_cuts.ini`: additional restrictions for the independent voltage study only.
- `config/raw_binning.csv`: 1D axes.
- `config/raw_2d_plots.csv`: 2D pairs, orientations and axis overrides.
- `analysis/settings.csv`, `analysis_root/settings.csv`: independent noise thresholds and histogram settings.

Printed fADC is one-based; ROOT Channel is zero-based; Pads are numbered 0–409.
A Channel cut removes the corresponding hits from mapped Pad plots too.
The six main 1D quantities, a separate electronics Channel view, and configured 2D
pairs are shown together by default. The current CSVs determine the precise binning
and pair orientation. Results are grouped under `results/python/` or `results/root/`.

Additional ROOT files may be placed under `data/root/` or selected with an input-directory override.
They are ignored by Git. The single-run sample can exercise the voltage command,
but the documented multi-voltage comparison requires additional runs.
The manual's full-run illustrations are reference figures, not plots of the reduced example file.
