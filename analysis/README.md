# TPC analysis — one program

See the [illustrated PDF manual](../docs/TPC_manual.pdf) and its
[LaTeX source](../docs/TPC_manual.tex). Rebuild with
`python docs/build_manual.py` after script/configuration changes.

From the data directory, run:

```bash
python analysis/tpc.py
```

Choose an item from the menu to browse events, save an event, recreate the
examples, or run the voltage analysis. No script selection is required.
The ROOT inputs are in `data/root/`; both input and output paths are resolved automatically.

For repeatable commands:

```bash
python analysis/tpc.py event 1616                         # browse numbered pads
python analysis/tpc.py event 1616 --event 5584 --save      # export PNG/PDF
python analysis/tpc.py event 1616 --event 5584 --mode 3d --save
python analysis/tpc.py examples                           # 1615/1616, 2D and 3D
python analysis/tpc.py voltage                            # voltage scan and audits
python analysis/tpc.py info                               # show current cuts/paths
python analysis/tpc.py check                              # regression checks
python analysis/tpc.py raw 1616                           # 1D + 2D together, UNCUT
python analysis/tpc.py raw 1616 --batch                   # save without a window
python analysis/tpc.py raw 1616 --1d                      # only 1D
python analysis/tpc.py raw 1616 --2d                      # only 2D
python analysis/tpc.py raw 1616 --apply-cuts              # both groups, cleaned
python analysis/tpc.py cuts 1616                          # show active DB rules
```

Raw inspection shows 1D and 2D together by default (13 windows), using one
input/mapping/selection and one window loop. `--1d` / `--2d` select a group;
`--dimensions both|1d|2d` is also available. The 1D group opens a 3-by-2
summary and six individual windows (Event,
Pad, Hit, Time, Charge, Pedestal). Pad always uses the established mapping,
0–409; unconnected channels are counted as Pad=-1 in underflow. Mapping
does not remove hit rows from the other distributions. Blockhead is omitted. All seven views
are accompanied by a separate original Channel diagnostic (electronics
0–575, including unconnected channels). All eight views are saved as PNG/PDF
and native canvases in `raw_histograms.root`, alongside their seven TH1
histograms. Close the windows to finish, or use `--batch` for unattended
exports. `raw 1616 --branches Channel` opens only the Channel diagnostic,
saving in `branches_Channel/` to preserve the full run's outputs.

The 2D group has a 2-by-2 summary and four individual canvases: Event vs
Pad, Hit vs Pad, Charge vs Pad, and the currently configured Charge vs Time.
Pads are horizontal in the first three; Time is horizontal in the last.
Axes and names are editable in the pair CSV. Outputs are in the run's
`2d/` or `cleaned/2d/` subfolder, including TH2D histograms and all five
canvases in `raw_2d_histograms.root`. Edit `config/raw_2d_plots.csv` for
pairs, axis binning and logarithmic colour scales. See the
[2D guide](../config/README.md#two-dimensional-inspection) for resolution,
flow accounting, focused views and file-only exports.

Saved events go to `results/python/events/runNNNN/`. Examples include `3d/` subfolders,
the comparison, and a CSV of the exact displayed hits. Voltage output goes
to `results/python/voltage/`. Results are kept separate from the input data.

Use `python analysis/tpc.py event --help`, `examples --help`, or
`voltage --help` for the options of that command. From elsewhere, use an
absolute script path and point at your data:

```bash
python /path/to/analysis/tpc.py event 1616 --event 5584 --save --input-dir /path/to/data --output-dir /path/to/plots
```

## What to edit

**Manual cuts:** edit `config/run_cuts.ini`; it has common `[global]` and
per-run blocks shared with ROOT. **Automatic noise thresholds and voltage
histogram limits:** edit `settings.csv`. See [cut/inspection guide](../config/README.md).
The independent voltage comparison keeps its additional restrictions in
`config/voltage_study_cuts.ini`; only the voltage command reads that file.

**Mapping:** `auxiliary/pad_map.csv` is the established mapping; no builder
needs to run. Pads 0–249 are medium, below large pads 250–409. The printed
one-based fADC label equals `Channel + 1`. Unconnected channels have no pad.
The top-row wiring corrections are already included.

Only three active source files remain:

| File | Edit it for |
| --- | --- |
| `tpc.py` | Menu/commands, event displays, plot styling, example export |
| `detector.py` | ROOT reading, selection, veto loading, pad geometry and noise detection |
| `voltage.py` | Voltage histograms, normalization, run choice and audit plots |

Comments inside each file describe the important stages. Tests are under
`tests/`; the archived mapping builder, ROOT macros, and earlier modules
are under `archive/previous_layout/` as historical reference. They are not
part of routine analysis and are not a second supported interface.

## Selections and future work

Event displays include **both** pad sections. Their current selection is
`Time < 228 && Channel != 271 && Channel != 149`, plus the selected run's manual DB rules and automatic event vetoes. The independent voltage scan additionally reads `config/voltage_study_cuts.ini`, which retains `Channel < 285` for that study.

Noise detection runs before these hit cuts. Its thresholds are in
`settings.csv`. Manual Time/Channel/Charge/Event rules are defined only in
`config/run_cuts.ini`. Vetoes use `(run,event)` pairs. The program prefers
`results/python/voltage/rejected_events.csv` in this workspace;
otherwise it uses the bundled snapshot. For different data, generate a
new ledger with the voltage command or pass `--veto-file` to event/examples.

The voltage scan selects the file with the most distinct events per
voltage. Run 1603 remains excluded. Run 1605 still has conflicting 8250/
8500 V entries; resolve it explicitly with `voltage --voltage-override
RUN=VOLTS`. `data/datafiles-voltage.txt` is the working run/voltage table. `--voltage-file` can select another table.

Physical displays retain the first passing hit per channel; voltage
histograms count all accepted raw hits. Missing pad values are drawn as
zero, which does not identify dead pads. Negative charge is clipped by
the display's zero-based color scale, not removed from the raw arrays.

For later reconstruction, reuse `detector.load_hit_arrays()` and
`detector.select_hit_mask()` to keep every passing pulse, including its
Event/Channel/Hit/Time/Charge. Pad centers follow from the mapping bounds.
Time remains in native units; drift calibration, t0, fitting, and track
output are future stages. No reconstruction is implemented here.

## Sharing

Copy `analysis/` together with its sibling `config/` folder, or share `bundles/analysis_bundle.zip`. Supply ROOT data separately.
Python needs ROOT/PyROOT, NumPy and Matplotlib. Tk is needed only for
interactive windows; exports work headlessly. `ENVIRONMENT.txt` records
the tested versions. Historical cut notes under `auxiliary/archive/` are
not loaded; the active settings are only in the top-level `settings.csv`.

Raw inspection uses `raw_binning.csv` in the shared config folder.
Use `raw --apply-cuts` to compare cleaning, and `--cuts-db` to select an
alternative complete manual-cut database. DB snapshots accompany exports.
