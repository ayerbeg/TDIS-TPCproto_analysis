# Manual cuts and raw-data inspection

`run_cuts.ini` is the **shared detector-cleaning database** for Python and
C++/ROOT, containing only global and per-run rules. The independent voltage
study has its own `voltage_study_cuts.ini`. `raw_binning.csv` controls the raw histograms in both programs.
The original ROOT inputs are always read-only.

The [illustrated PDF manual](../docs/TPC_manual.pdf) includes current
configuration snapshots and plot examples. Plot CSV headers/values may
have surrounding spaces or tabs for alignment in both engines; keep commas
as separators, including consecutive commas for empty fields.

## Cut database

Write a condition describing the values to **remove**:

```ini
[global]
reject = Time >= 228
reject = Channel in 271,149

[run 1616]
# Illustrations only: inspect data before enabling a new exclusion.
# reject = Channel == 12
# reject = Event in 1000..1020
# reject = Channel in 40..45
# reject = Channel == 44 && Time >= 90 && Time < 95
# reject = Charge < 0
# reject = Charge > 3500
```

- `[global]` applies to **every** run and cleaning workflow.
- `[run NNNN]` adds rules only for that exact run.
- Separate `reject` lines are **OR**: matching any line rejects the row.
- Conditions joined with `&&` on the same line are **AND**: all must match.
- `#` starts a comment, including at the end of a rule. Commenting out a
  line disables it. Each block name may appear once.

The detector database contains the existing common cuts and per-run blocks.
No additional run-specific exclusions have been enabled without inspection.
The old manual values have been removed from both `settings.csv` files so
there is no second hidden copy of the global cuts. The remaining settings
govern automatic-noise thresholds, excluded scan runs, and histogram limits.
If you extend the accepted Time range for a voltage scan, also extend its
`time_max` histogram limit in the relevant program's `settings.csv`.

Fields: `Event`, `Channel`, `Hit`, `Time`, `Charge`, `Pedestal`, `Blockhead`,
plus derived `FADC` and `Pad`. `Channel` is the zero-based ROOT value;
`FADC = Channel + 1`; `Pad` uses the active mapping and is `-1` for unmapped
electronics. Branch names are case-insensitive. Numeric values may be
negative, fractional or scientific notation; they remain in native units.

A Channel exclusion also removes its mapped pad's hits: all rules select
raw hit rows before any cleaned histogram or event display is populated.
For the current mapping, `Channel == 152`, `FADC == 153`, and `Pad == 0`
are equivalent ways to reject the same pad. One rule is sufficient. The
numbers themselves are not interchangeable: ROOT Channel152 is Pad0,
not Pad152. Uncut raw views intentionally bypass all these DB rules.

Operators: `==`, `!=`, `<`, `<=`, `>`, `>=`, `in`, `not in`.
`in 40..45` includes **both** endpoints. `in 12,44,271` is a value list.
For a histogram bin with edges 90 and 95, use `Time >= 90 && Time < 95`.
There are no arbitrary functions, arithmetic, parentheses, or `||`; express
OR by adding another rejection line. Invalid syntax, unknown fields, and
missing input branches produce errors with the DB filename and line.

An Event rule removes all rows with that event ID. A channel/time rule
removes matching **hit rows**, not the rest of the event. Rules run before
the display keeps the first accepted hit per channel, so a later good pulse
can still be shown if the earlier one was removed.

Automatic event-noise vetoes remain a separate generated ledger with their
diagnostic evidence. Clean data must pass both that ledger and the manual
DB. Uncut inspection bypasses **all** manual rules and automatic vetoes.

## Inspect and then clean

From the workspace root:

```bash
# Default: both 1D and 2D, with summary and individual windows.
python analysis/tpc.py raw 1616
python analysis/tpc.py raw 1616 --1d
python analysis/tpc.py raw 1616 --2d

# Focus on channels, change binning, or zoom onto an event interval.
python analysis/tpc.py raw 1616 --branches Pad --bins 410 --range -0.5 409.5
python analysis/tpc.py raw 1616 --branches Channel
python analysis/tpc.py raw 1616 --branches Event --range 8000 9000 --bin-width 1 --log-y

# List the detector-cleaning rules active for this run.
python analysis/tpc.py cuts 1616

# After editing the DB, compare cleaned raw histograms and rerun a display.
python analysis/tpc.py raw 1616 --apply-cuts
python analysis/tpc.py event 1616 --event 5584 --save
```

The six plotted quantities are Event, **Pad**, Hit, Time, Charge, and Pedestal;
Blockhead is omitted. Pad is derived from Channel using the established
mapping, always, even without cleaning: ROOT Channel -> printed FADC
(Channel+1) -> zero-based Pad. Medium pads occupy 0–249, large pads 250–409.
The original electronics **Channel** branch is also plotted in an extra
standalone diagnostic canvas, with all 576 channels (0–575), including
unconnected channels. Its axis explicitly states that printed fADC is
Channel+1. `--branches Channel` / ROOT `raw:Channel` inspect electronics;
`--branches Pad` / `raw:Pad` inspect mapped pads.
The 1D group has one 3-by-2 summary canvas, six individual mapped-analysis
canvases, and the separate Channel canvas (eight windows total).
Each view is exported as PNG/PDF. Both programs also
save the seven TH1 histograms and all eight native TCanvas objects inside
`raw_histograms.root`, so they can be reopened and edited in ROOT.

Default inspection opens both groups together: **13 windows** (eight 1D
plus five 2D). Python reads/maps/selects the data once and prepares both
groups before entering one shared window loop. Close the windows to finish.
Add `--1d` or `--2d` to select a group, or use
`--dimensions both|1d|2d` (`--dim` is an alias). Add `--batch` to save without windows,
for example `python analysis/tpc.py raw 1616 --batch`. If a graphical backend
is unavailable, the program explains why and still saves the outputs.

Native ROOT equivalents:

```bash
root -l 'analysis_root/TPC.C("raw",1616)'
root -l 'analysis_root/TPC.C("raw --1d",1616)'
root -l 'analysis_root/TPC.C("raw --2d",1616)'
root -l 'analysis_root/TPC.C("raw:Pad",1616)'
root -l 'analysis_root/TPC.C("raw:Channel",1616)'
root -l 'analysis_root/TPC.C("cuts",1616)'
root -l -b -q 'analysis_root/TPC.C("raw-clean",1616)'
root -l -b -q 'analysis_root/TPC.C("save",1616,5584)'
```

ROOT opens the summary and individual windows without `-b`; batch mode
saves the same plots without windows. Raw outputs go to
`results/python/raw/runNNNN/` or `results/root/raw/runNNNN/`, with cleaned
views in `cleaned/`. They contain PNG/PDF overviews, native ROOT histograms,
per-branch bin/count CSVs and a summary. `raw_overview` is the summary;
`raw_Event`, `raw_Pad`, etc. are the individual exports. To browse the
native canvases interactively, use:

```bash
root -l results/root/raw/run1616/raw_histograms.root
```

Then run `new TBrowser;` at the ROOT prompt and double-click a canvas.

Focused `--branches` / `raw:NAME` commands open only the requested
individual views. They save under `branches_NAME/` (or
`branches_NAME1_NAME2/` in Python), so they do not overwrite the full run's
summary or ROOT file. Cleaned focused views are under `cleaned/branches_NAME/`.
Likewise, `--plots NAME` alone selects focused 2D inspection. Explicit
`--dimensions both` can combine focused 1D and 2D selections. Python
`--bins`, `--range`, `--bin-width` and `--log-y` control the 1D group;
2D axis overrides remain in `raw_2d_plots.csv`.

Edit `raw_binning.csv` for persistent per-branch choices. `BinWidth` is in
native units; a nonempty `Bins` takes precedence. `Min`/`Max` are histogram
**edges**, not rejection cuts. Empty limits use the full raw extent; log
scales are controlled by `LogY`. Integer branches default to one bin/value;
Pad includes the full 0–409 detector, Channel includes electronics 0–575,
and Charge includes negative data.
Counts below/above a plotting range are reported as under/overflow. All
bins, including the final one, use ROOT's `[low, high)` convention.
Automatic raw/clean axes are identical so removed tails remain visible.

Mapping is not an additional hit cut. Unconnected channels have Pad=-1,
recorded in the Pad histogram's underflow and explicitly reported as
unmapped. The other raw distributions still include those rows. The summary
records total and selected unmapped counts; `pad_map_used.csv` carries the
actual wiring table. The original Channel array and native units stay
unchanged. DB `Channel` conditions still refer to ROOT electronics values;
use `Pad` conditions when selecting by the pad numbers shown in the plot.
Plot labels use **Raw Hits**, **Underflow**, **Overflow**, and
**Hits on unconnected channels**. Counts have thousands separators.

Python options `--bins`, `--bin-width`, `--range`, `--log-y`, `--linear-y`
override those choices for one invocation. The same per-axis binning
helpers are used by 2D inspection.

## Two-dimensional inspection

The 2D group opens one 2-by-2 summary and four individual canvases:

| Plot name | Horizontal X axis | Vertical Y axis |
| --- | --- | --- |
| `Event_vs_Pad` | Pad | Event ID |
| `Hit_vs_Pad` | Pad | Hit (stored pulse index) |
| `Charge_vs_Pad` | Pad | Charge |
| `Charge_vs_Time` | Time | Charge |

Every raw pulse fills one histogram cell; the colour scale is **Raw Hits**.
Multiple pulses on the same pad/event remain separate contributions.
`Hit` is the original Hit branch, not an event's total number of pads.
Colour scales are logarithmic by default; empty cells are white.

```bash
# Both groups are the default; flags select just 2D when wanted.
python analysis/tpc.py raw 1616 --2d
python analysis/tpc.py raw 1616 --2d --apply-cuts
root -l 'analysis_root/TPC.C("raw --2d",1616)'
root -l 'analysis_root/TPC.C("raw-clean --2d",1616)'

# Focus on one individual view without replacing the full 2D output.
python analysis/tpc.py raw 1616 --plots Charge_vs_Pad
root -l 'analysis_root/TPC.C("raw2d:Charge_vs_Pad",1616)'
```

Outputs live in `results/python/raw/runNNNN/2d/` and
`results/root/raw/runNNNN/2d/`, or `cleaned/2d/` for cleaned views. Each
contains summary/individual PNG/PDF files, four native TH2D histograms and
five editable TCanvas objects in `raw_2d_histograms.root`, nonzero-bin CSV
tables, a summary (Python JSON / ROOT CSV), and configuration/mapping
snapshots. Focused exports use a `plots_NAME/` subfolder.

Mapping and cleaning follow exactly the 1D convention: Pad is always
0–409, with unconnected hits at Pad=-1 in underflow. Time vs Charge includes
those hits because it has no Pad axis. Raw views bypass the DB and vetoes;
cleaned views apply the detector global/run DB and automatic event vetoes.
The independent voltage-study restrictions are never applied here.
The previous Python `raw2d` and ROOT `raw2d` / `raw2d-clean` commands
remain 2D-only shortcuts. The normal menu now has one combined inspection
entry; ROOT `raw-clean` also shows both groups by default.

Edit **`raw_2d_plots.csv`** to change or add pairs. `X` and `Y` are field
names. `XBinWidth`, `XBins`, `XMin`, `XMax` and their Y counterparts override
the corresponding 1D axis; empty cells inherit `raw_binning.csv`.
An explicit bin count takes precedence over width. `LogZ` controls the
colour scale; X and Y remain linear. Python `--plots-file` selects an
alternate pair/axis file.

Current resolution is **64 event IDs per bin**; the remaining axes inherit
their one-unit bins from `raw_binning.csv`. Bin widths are
printed on the axes. To inspect individual event IDs in a smaller region,
set `YBinWidth=1`, `YMin=8000`, `YMax=9000` on the Event-vs-Pad row.
These are plotting edges, not cuts. Raw/cleaned axes always use the same
uncut data extent. A limit of four million cells per histogram prevents
accidentally allocating an enormous dense Event/Pad array.

ROOT retains all underflow/overflow cells, including both-axis corners.
The sparse CSV contains every nonzero cell: bin 0 is underflow and N+1
is overflow, with infinite open-end bounds. **Outside plotted ranges**
counts each hit once, even when it overflows both axes. Per-axis flow
counts in the summary can therefore overlap. Native ROOT histograms retain
true entries and moments from the original accepted pulse coordinates.

## Reproducibility and alternatives

### Independent voltage study

The voltage comparison reads `run_cuts.ini` for detector cleaning, then adds
restrictions from **`voltage_study_cuts.ini`**. Only the voltage command
opens this second file. Event displays and raw/cleaned inspection use the
detector database alone, so both pad sections remain available.

The study file preserves the original comparison's electronics limit:

```ini
# config/voltage_study_cuts.ini -- used ONLY by the voltage study
[global]
reject = Channel >= 285
```

Here `[global]` means all runs within that study. Both files use the same
rule syntax; study-specific per-run rules can also go in the study file.
General noisy-channel/event cleaning belongs in `run_cuts.ini`.

```bash
python analysis/tpc.py voltage --study-cuts-db config/voltage_study_cuts.ini
root -l -b -q 'analysis_root/TPC.C("voltage")'
```

ROOT uses that study file by default; its optional eighth `studyCutsDB`
argument overrides it. Voltage exports save both `cuts_used.ini` (detector
cleaning) and `voltage_study_cuts_used.ini` (study restrictions). Removing
the study-file rule changes only the independent voltage study.

### Detector database overrides

Each invocation reads the DB; reopen an interactive browser after edits.
Event/example/cleaned-raw exports and voltage analysis save the parsed detector DB
text alongside their output. The snapshot is the version actually used.

An alternate file is a **complete replacement**, including its global
block, not an append to the default file:

```bash
python analysis/tpc.py event 1616 --event 5584 --save --cuts-db /path/cuts.ini
root -l -b -q 'analysis_root/TPC.C("save",1616,5584,"","","","/path/cuts.ini")'
```

Both shareable ZIPs include this `config/` folder. Keep it beside `analysis/`
or `analysis_root/` when copying the programs.
