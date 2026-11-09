# Development

The active analysis consists of three Python files and two native ROOT files.
`analysis/detector.py` and `analysis_root/TPCDetector.h` define input, mapping,
selection and noise classification. Drawing, commands and exports are in `tpc.py`
and `TPC.C`; `voltage.py` contains the independent Python voltage study.
The `PersonalTPC` C++ namespace is retained for API compatibility.

Preserve full pulse records for future reconstruction. The event display's first
passing hit per channel is a visualization policy, not the default analysis reduction.
Keep cuts on electronics values distinct from mapped coordinates.

## Checks

```bash
python -m unittest discover -s analysis/tests -v
python analysis/tpc.py check
root -l -b -q 'analysis_root/TPC.C("check")'
root -l -b -q 'analysis_root/tests/RunCutsCheck.C()'
root -l -b -q 'analysis_root/tests/Raw2DCheck.C()'
```

The Python unit checks use NumPy and fixed fixtures and do not require ROOT.
Native integration checks create temporary synthetic files. Validate changed plot
counts against both implementations when modifying binning, mapping, or selection.
Generated `results/` and large local datasets must not be committed.

## Manual

Update `docs/TPC_manual.tex` for user-visible command, selection, or behaviour changes.
The PDF is a user/developer reference; avoid conversational history and personal support answers.
Current settings, command help and palette source are generated excerpts.

```bash
python docs/build_manual.py
python docs/build_manual.py --check
```

Building requires `pdflatex`. Reference figures in `docs/figures/` remain available
when the full original datasets are absent. Regenerate the relevant plots and refresh
reference figures when their appearance changes. `--check` detects changes to tracked
sources, configurations, mapping and reference figures; generated excerpts do not
replace explanations of new semantics.
