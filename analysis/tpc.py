#!/usr/bin/env python3
"""One command for routine TPC analysis.

python analysis/tpc.py                         # simple text menu
python analysis/tpc.py event 1616               # interactive numbered pads
python analysis/tpc.py event 1616 --event 5584 --save
python analysis/tpc.py examples                # both established examples, 2D/3D
python analysis/tpc.py voltage                 # original voltage comparison
python analysis/tpc.py info                    # paths, mapping and current cuts

Edit ../config/run_cuts.ini for manual cuts and settings.csv for automatic
noise thresholds. This file controls displays and the
command/menu interface. detector.py supplies reusable raw-hit/geometry
functions; voltage.py contains the longer voltage-scan plotting workflow.
No track reconstruction is performed. Raw Time remains in native units.
"""
import argparse
import csv
import json
import math
from array import array
import select
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # interactive event/raw commands switch to TkAgg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import numpy as np
from detector import (DEFAULT_INPUT_DIR, DEFAULT_PAD_MAP, SETTINGS_FILE,
    BASE_CUT, VOLTAGE_CUT, TIME_CUT_MAX, BOARD_WIDTH_MM, BOARD_HEIGHT_MM,
    default_veto_file, load_pad_map, load_pad_records, channel_to_pad,
    pad_rect_mm, pad_number, passes_hit_cuts, select_event_hits,
    rejected_events_for_run, load_hit_arrays, results_directory, DEFAULT_CUTS_DB, DEFAULT_BINNING,
    CutDatabase, cut_database, select_hit_mask, raw_histogram_spec, DEFAULT_VOLTAGE_STUDY_DB,
    DEFAULT_2D_PLOTS, load_2d_plots, raw_2d_axes)

# ---- Drawing only: change styles here without changing hit selection ----
def render_2d(ax, pad_map, values_by_channel, cmap='plasma', title="", show_pad_numbers=True):
    """values_by_channel: dict {channel: value}. Every pad in pad_map is
    drawn at its true rect; pads with no entry in values_by_channel are
    drawn as 0 (baseline), matching the ROOT version's convention."""
    # A missing pad value is drawn at 0 for visualization. That does NOT
    # mean a measured Time of zero or a known dead channel. Scales start
    # at zero, so negative Charge values are visually clipped at baseline.
    vals = list(values_by_channel.values())
    vmax = max(vals) if vals else 1.0
    patches, colors = [], []
    for ch, (size, row, col) in pad_map.items():
        xlo, ylo, w, h = pad_rect_mm(size, row, col)
        val = values_by_channel.get(ch, 0)
        patches.append(Rectangle((xlo, ylo), w, h))
        colors.append(val)
        if show_pad_numbers:
            ax.text(xlo + w / 2, ylo + h / 2, str(pad_number(size, row, col)),
                    ha='center', va='center', fontsize=4.5, color='white' if val < vmax / 2 else 'black', zorder=3)
    pc = PatchCollection(patches, cmap=cmap, edgecolor='black', linewidth=0.2)
    pc.set_array(np.array(colors))
    pc.set_clim(0, vmax if vmax > 0 else 1)
    ax.add_collection(pc)
    ax.set_xlim(0, BOARD_WIDTH_MM)
    ax.set_ylim(0, BOARD_HEIGHT_MM)
    ax.set_aspect('equal')
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    ax.set_title(title)
    ax.axhline(40, color='white', linewidth=1, linestyle='--')  # medium/large boundary
    ax.figure.colorbar(pc, ax=ax, fraction=0.04, pad=0.02)
    return pc


def render_3d(ax, pad_map, values_by_channel, cmap='plasma', title=""):
    """Each pad drawn as one 3D bar at its true position/size; height =
    value. Only ~410 bars total (one per pad), not one per mm -- much
    lighter than a rasterized fine grid would be."""
    xs, ys, dxs, dys, dzs = [], [], [], [], []
    for ch, (size, row, col) in pad_map.items():
        xlo, ylo, w, h = pad_rect_mm(size, row, col)
        val = values_by_channel.get(ch, 0)
        xs.append(xlo); ys.append(ylo); dxs.append(w); dys.append(h); dzs.append(val)
    dzs = np.array(dzs)
    vmax = dzs.max() if dzs.max() > 0 else 1.0
    colors = plt.get_cmap(cmap)(dzs / vmax)
    # Each bar footprint follows the pad bounds; height is Charge or Time.
    # Do not use this display's vertical axis as a reconstructed z axis.
    ax.bar3d(xs, ys, np.zeros_like(dzs), dxs, dys, dzs, color=colors, shade=True)
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_zlabel("value")
    ax.set_title(title)
    ax.set_xlim(0, BOARD_WIDTH_MM); ax.set_ylim(0, BOARD_HEIGHT_MM)
    ax.set_box_aspect((BOARD_WIDTH_MM, BOARD_HEIGHT_MM, 50))


# ---- Single-event navigation/export ----
def stdin_has_input(timeout_sec):
    ready, _, _ = select.select([sys.stdin], [], [], timeout_sec)
    return bool(ready)


def read_line_keeping_plots_alive(fig_charge, fig_time):
    """Polls for stdin input while repeatedly pausing on the matplotlib
    figures, so their GUI event loops (including 3D rotation) keep
    running instead of freezing while the prompt waits for you to type."""
    while not stdin_has_input(0.05):
        plt.pause(0.05)  # lets matplotlib process GUI events (rotation, resize, etc.)
    return sys.stdin.readline().strip()


def load_tree_arrays(root_filename, tree_name="ntuple"):
    """Tuple for plotting; detector.load_hit_arrays() also exposes all raw branches."""
    arrays, f = load_hit_arrays(root_filename, tree_name)
    return arrays['Event'], arrays['Channel'], arrays['Time'], arrays['Charge'], f


def run_event(args):
    if args.save and args.output_dir is None:
        args.output_dir = results_directory(args.input_dir) / 'events' / f'run{args.run_number}'
    if not args.output_dir:
        plt.switch_backend('TkAgg')

    pad_map = load_pad_map(args.pad_map)
    print(f"Loaded {len(pad_map)} channel mappings from {args.pad_map}")

    # Paths are independent of the current shell directory. --output-dir
    # exports a single event; omitting it opens the interactive navigator.
    filename = str(args.input_dir / f"ntuple{args.run_number:04d}.root")
    print(f"Opening {filename} ...")
    arrays, _f = load_hit_arrays(filename)
    events, channels, times, charges = (arrays[name] for name in ('Event', 'Channel', 'Time', 'Charge'))
    database = cut_database(args.cuts_db)
    print(f"Loaded {len(events)} hit entries.")

    vetoes = rejected_events_for_run(args.run_number, args.veto_file or default_veto_file(args.input_dir))
    accepted = select_hit_mask(events, channels, times, vetoes, run=args.run_number,
                               cuts_db=database, arrays=arrays, pad_map=args.pad_map)
    mapped = channel_to_pad(channels, args.pad_map) >= 0
    event_set = sorted(set(events[accepted & mapped]))
    print(f"DB: {database.path}; {len(database.active(args.run_number))} active rules")
    print(f"Cuts: {database.describe(args.run_number)}; {len(vetoes)} saved event vetoes")
    if not event_set:
        print(f"ERROR: no events found with any hit passing Time<{TIME_CUT_MAX}.")
        return
    if args.event is not None and args.event not in event_set:
        raise ValueError(f'Event {args.event} has no mapped hits after DB cuts and automatic vetoes')
    print(f"Found {len(event_set)} distinct events with hits passing Time<{TIME_CUT_MAX} "
          f"(range {event_set[0]} to {event_set[-1]}).")

    render_fn = render_3d if args.mode == "3d" else render_2d

    if not args.output_dir:
        plt.ion()
    fig_charge = plt.figure(figsize=(10, 9))
    fig_time = plt.figure(figsize=(10, 9))

    def draw_event(target_event):
        # The first passing hit per Channel is a display convention, not a
        # physical charge sum. Later reconstruction should use raw arrays
        # and select_hit_mask() from detector.py to preserve all hits.
        hits = select_event_hits(events, channels, times, charges, target_event, vetoes, accepted_mask=accepted)
        hits = [(ch, q, t) for ch, q, t in hits if ch in pad_map]
        # Values remain keyed by zero-based ROOT Channel. The geometry
        # lookup decides positions and the pad-number labels in the plots.
        charge_by_channel = {ch: q for ch, q, t in hits}
        time_by_channel = {ch: t for ch, q, t in hits}

        fig_charge.clear()
        ax_c = fig_charge.add_subplot(111, projection='3d' if args.mode == '3d' else None)
        render_fn(ax_c, pad_map, charge_by_channel,
                  title=f"Charge - Run {args.run_number} - Event {target_event}")
        fig_charge.canvas.draw_idle()

        fig_time.clear()
        ax_t = fig_time.add_subplot(111, projection='3d' if args.mode == '3d' else None)
        render_fn(ax_t, pad_map, time_by_channel, cmap='viridis',
                  title=f"Time - Run {args.run_number} - Event {target_event}")
        fig_time.canvas.draw_idle()

        if not args.output_dir:
            plt.pause(0.01)
        print(f"Event {target_event}: {len(hits)} hit(s) drawn.")

    current_idx = event_set.index(args.event) if args.event in event_set else 0
    if args.output_dir:
        # Exporting never changes the original ntuple. PNG is convenient
        # for inspection; PDF keeps labels suitable for sharing/printing.
        target = args.event if args.event is not None else event_set[0]
        if target not in event_set:
            raise ValueError(f"Event {target} is vetoed or has no mapped hits passing the cuts")
        args.output_dir.mkdir(parents=True, exist_ok=True)
        draw_event(target)
        for name, fig in [('charge', fig_charge), ('time', fig_time)]:
            for suffix in ('png', 'pdf'):
                fig.savefig(args.output_dir / f"event_{args.run_number}_{target}_{name}.{suffix}", dpi=180)
        database.snapshot(args.output_dir / f"event_{args.run_number}_{target}_cuts.ini")
        print(f"Saved numbered event display to {args.output_dir}")
        return
    draw_event(event_set[current_idx])

    while True:
        print(f"Event {event_set[current_idx]} ({current_idx+1} of {len(event_set)}). "
              f"[Enter]=next  p=previous  <number>=jump  q=quit: ", end="", flush=True)
        line = read_line_keeping_plots_alive(fig_charge, fig_time)

        if line == "":
            if current_idx + 1 < len(event_set):
                current_idx += 1
            else:
                print("Already at the last event with hits.")
        elif line == "p":
            if current_idx > 0:
                current_idx -= 1
            else:
                print("Already at the first event with hits.")
        elif line == "q":
            break
        else:
            try:
                target = int(line)
                if target in event_set:
                    current_idx = event_set.index(target)
                draw_event(target)
                continue
            except ValueError:
                print(f"Not a recognized command ('{line}'). Use Enter, p, q, or an event number.")
                continue

        draw_event(event_set[current_idx])


# ---- Reproduce examples and the exact displayed-hit table ----
def run_examples(args):
    examples = [tuple(map(int, s.split(':'))) for s in (args.example or ['1615:42535', '1616:5584'])]
    out = args.output_dir or results_directory(args.input_dir) / 'events'
    comparison_dir = out / 'comparisons'
    comparison_dir.mkdir(parents=True, exist_ok=True)
    pad_map = load_pad_map(args.pad_map)
    database = cut_database(args.cuts_db)
    comparison, axes = plt.subplots(len(examples), 2, figsize=(18, 7.5 * len(examples)), squeeze=False)
    rows = []
    for i, (run, event) in enumerate(examples):
        arrays, root_file = load_hit_arrays(args.input_dir / f'ntuple{run}.root')
        events, channels, times, charges = (arrays[name] for name in ('Event', 'Channel', 'Time', 'Charge'))
        vetoes = rejected_events_for_run(run, args.veto_file or default_veto_file(args.input_dir))
        if event in vetoes:
            raise ValueError(f'Run {run}, event {event} is noise-vetoed')
        hits = [(ch, q, t) for ch, q, t in select_event_hits(events, channels, times, charges, event, vetoes,
                run=run, cuts_db=database, arrays=arrays, pad_map=args.pad_map)
                if ch in pad_map]
        if not hits:
            raise ValueError(f'Run {run}, event {event} has no selected mapped hits')
        # First-hit-per-channel reduction is only for current displays.
        # The CSV records precisely these plotted hits, not all raw pulses.
        pads = channel_to_pad([h[0] for h in hits], args.pad_map)
        for pad, (channel, charge, time) in zip(pads, hits):
            rows.append(dict(run=run, event=event, pad=int(pad), fadc=int(channel)+1,
                             root_channel=int(channel), charge=int(charge), time=int(time)))
        for mode, renderer in [('2d', render_2d), ('3d', render_3d)]:
            dest = out / f'run{run}'
            if mode == '3d':
                dest /= '3d'
            dest.mkdir(parents=True, exist_ok=True)
            for j, (name, values, cmap) in enumerate([
                ('charge', {ch: q for ch, q, t in hits}, 'plasma'),
                ('time', {ch: t for ch, q, t in hits}, 'viridis')]):
                fig = plt.figure(figsize=(10, 9))
                ax = fig.add_subplot(111, projection='3d' if mode == '3d' else None)
                renderer(ax, pad_map, values, cmap=cmap, title=f'{name.title()} - Run {run} - Event {event}')
                for ext in ('png', 'pdf'):
                    fig.savefig(dest / f'event_{run}_{event}_{name}.{ext}', dpi=180)
                plt.close(fig)
                if mode == '2d':
                    render_2d(axes[i, j], pad_map, values, cmap=cmap,
                              title=f'Run {run} | Event {event} | {name.title()} | {len(hits)} pads')
        root_file.Close()
        print(f'Run {run}, event {event}: {(pads < 250).sum()} medium, {(pads >= 250).sum()} large pads', flush=True)
    # Comparison filenames encode the chosen runs for easy sharing.
    tag = '_'.join(str(run) for run, event in examples)
    comparison.suptitle('Numbered event displays — medium pads 0–249 below large pads 250–409', fontsize=17)
    comparison.text(0.5, 0.025, database.describe() + '; per-run DB rules' + '; saved event vetoes; first passing hit per channel. Both sections included.', ha='center', fontsize=10)
    comparison.tight_layout(rect=(0, 0.05, 1, 0.96))
    for ext in ('png', 'pdf'):
        comparison.savefig(comparison_dir / f'examples_{tag}.{ext}', dpi=180)
    plt.close(comparison)
    with (comparison_dir / f'example_hits_{tag}.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    database.snapshot(comparison_dir / f'examples_{tag}_cuts.ini')
    print(f'Saved examples to {out}')


# ---- Common inspection input: one reading and selection for both groups ----
def inspection_data(args):
    """Keep mapping, detector cuts and event vetoes identical in 1D and 2D."""
    arrays, source = load_hit_arrays(args.input_dir / f'ntuple{args.run_number}.root')
    try:
        arrays['Pad'] = channel_to_pad(arrays['Channel'], args.pad_map)
        arrays['FADC'] = arrays['Channel'] + 1
        database = cut_database(args.cuts_db) if args.apply_cuts else None
        keep = np.ones(len(arrays['Event']), dtype=bool)
        if args.apply_cuts:
            vetoes = rejected_events_for_run(args.run_number, args.veto_file or default_veto_file(args.input_dir))
            keep = select_hit_mask(arrays['Event'], arrays['Channel'], arrays['Time'], vetoes,
                                   run=args.run_number, cuts_db=database, arrays=arrays, pad_map=args.pad_map)
        return dict(arrays=arrays, source=source, database=database, keep=keep)
    except Exception:
        source.Close()
        raise


# ---- 1D drawing/export; common entry point can defer display until 2D is ready ----
def run_raw(args, *, inspection=None, defer_show=False):
    # Raw inspection normally opens a zoomable window AND saves its outputs.
    # --batch keeps unattended exports independent of a graphical desktop.
    show_window = not args.batch
    if show_window and not defer_show:
        try:
            plt.switch_backend('TkAgg')
        except (ImportError, RuntimeError) as error:
            print(f'Cannot open a plot window: {error}. Saving outputs only; '
                  'use --batch for unattended exports.', flush=True)
            show_window = False
    owns_input = inspection is None
    data = inspection if inspection is not None else inspection_data(args)
    arrays, root_file = data['arrays'], data['source']
    # Mapping is a coordinate conversion, independent of every cut. Keep
    # electronics Channel unchanged for DB rules; use derived Pad in plots.
    # Unconnected electronics map to -1, reported in the Pad underflow.
    preferred = ['Event', 'Pad', 'Hit', 'Time', 'Charge', 'Pedestal']
    # Channel is a separate electronics diagnostic; Pad stays in the six
    # summary panels. Both representations are useful and retain all rows.
    available = preferred + ['Channel']
    branches = list(dict.fromkeys(args.branches)) if args.branches else [b for b in available if b in arrays]
    for name in branches:
        if name not in available or name not in arrays:
            raise ValueError(f'Branch {name!r} is unavailable; choose from {[b for b in available if b in arrays]}')
    summary_branches = [b for b in preferred if b in arrays] if not args.branches else []
    database, keep = data['database'], data['keep']
    profile = 'clean' if args.apply_cuts else 'raw'
    out = args.output_dir or results_directory(args.input_dir) / 'raw' / f'run{args.run_number}'
    if args.apply_cuts:
        out = out / 'cleaned'
    if args.branches:
        # Focused inspection must not overwrite the full run's ROOT file.
        out = out / ('branches_' + '_'.join(branches))
    out.mkdir(parents=True, exist_ok=True)
    fig = None
    figures = []
    if summary_branches:
        columns = min(3, len(summary_branches)); rows = math.ceil(len(summary_branches) / columns)
        fig, axes = plt.subplots(rows, columns, figsize=(5 * columns, 3.7 * rows), squeeze=False)
        figures.append(fig)
    import ROOT
    # The Python windows use Matplotlib. Also persist native ROOT canvases,
    # created in batch mode to avoid opening duplicate ROOT GUI windows.
    previous_root_batch = ROOT.gROOT.IsBatch()
    ROOT.gROOT.SetBatch(True)
    root_out = ROOT.TFile(str(out / 'raw_histograms.root'), 'RECREATE')
    canvas_stem = f'TPC_raw_{args.run_number}' + ('_clean' if args.apply_cuts else '')
    root_summary = None
    if summary_branches:
        root_summary = ROOT.TCanvas(canvas_stem, f'Run {args.run_number} | {profile}', columns * 500, rows * 360)
        root_summary.Divide(columns, rows)
    ROOT.TNamed('selection', 'UNCUT: all raw hit rows' if profile == 'raw' else database.describe(args.run_number)).Write()
    ROOT.TNamed('mapping', f'{args.pad_map}: Pad from ROOT Channel; FADC=Channel+1; Pad=-1 is unmapped (underflow)').Write()
    metadata = dict(run=args.run_number, selection=profile, total_rows=len(keep), selected_rows=int(keep.sum()),
                    cuts_applied=bool(args.apply_cuts), cuts_db=str(database.path) if args.apply_cuts else None,
                    pad_map=str(args.pad_map.resolve()), unmapped_rows=int((arrays['Pad'] < 0).sum()),
                    selected_unmapped_rows=int(((arrays['Pad'] < 0) & keep).sum()),
                    summary_branches=summary_branches,
                    histograms=[])
    # Carry the actual wiring table alongside the outputs for reproducibility.
    (out / 'pad_map_used.csv').write_bytes(args.pad_map.read_bytes())
    for branch in branches:
        values = arrays[branch][keep]
        # Base axes on the uncut column so raw/clean comparisons have
        # identical automatic edges, including the bins removed by cuts.
        edges, log_y = raw_histogram_spec(arrays[branch], branch, args.bins, args.bin_width, args.range, args.log_y, args.binning_file)
        finite = values[np.isfinite(values)]
        # Match ROOT's [low, high) convention, including the LAST bin.
        inside = (finite >= edges[0]) & (finite < edges[-1])
        counts, _ = np.histogram(finite[inside], bins=edges)
        under = int((finite < edges[0]).sum()); over = int((finite >= edges[-1]).sum())
        title = f'Run {args.run_number} | {profile.upper()} | {branch}'
        axis_label = 'Pad number (0-based)' if branch == 'Pad' else f'{branch} (native values)'
        if branch == 'Channel':
            axis_label = 'ROOT Channel (0-based, fADC = Channel + 1)'
        annotation_lines = [f'Underflow: {under:,}   Overflow: {over:,}']
        if branch in ('Pad', 'Channel'):
            annotation_lines.append(f"Hits on unconnected channels: {metadata['selected_unmapped_rows']:,}")

        def draw_axis(target):
            """Use the same counts/edges in the summary and individual view."""
            target.stairs(counts, edges, color='#1559D4', linewidth=1.2)
            target.set(xlabel=axis_label, ylabel='Raw Hits', title=title)
            if log_y and counts.max(initial=0):
                target.set_yscale('log'); target.set_ylim(bottom=.7)
            target.grid(alpha=.2)
            target.text(.02, .96, '\n'.join(annotation_lines), transform=target.transAxes, va='top', fontsize=8)

        if branch in summary_branches:
            draw_axis(axes.flat[summary_branches.index(branch)])
        individual, single_axis = plt.subplots(figsize=(10, 7))
        figures.append(individual)
        draw_axis(single_axis)
        individual.tight_layout()
        for ext in ('png', 'pdf'):
            individual.savefig(out / f'raw_{branch}.{ext}', dpi=160)

        hist = ROOT.TH1D(f'raw_{branch}', f'{title};{axis_label};Raw Hits', len(counts), array('d', edges))
        hist.SetDirectory(0)
        for i, n in enumerate(counts, 1):
            hist.SetBinContent(i, int(n)); hist.SetBinError(i, np.sqrt(n))
        hist.SetBinContent(0, under); hist.SetBinContent(len(counts)+1, over); hist.SetEntries(len(finite))
        root_out.cd(); hist.Write()

        def draw_root(pad):
            pad.cd()
            pad.SetBottomMargin(.14)
            pad.SetLogy(bool(log_y and counts.max(initial=0)))
            if log_y and counts.max(initial=0): hist.SetMinimum(.7)
            hist.SetLineColor(ROOT.kBlue + 1)
            hist.SetStats(False)
            # Copies belong to their canvas and survive file closure.
            hist.DrawCopy('HIST')
            note = ROOT.TLatex()
            note.SetNDC(); note.SetTextSize(.028)
            for index, line in enumerate(annotation_lines):
                note.DrawLatex(.13, .84 - .045 * index, line)

        if branch in summary_branches:
            draw_root(root_summary.cd(summary_branches.index(branch) + 1))
        root_individual = ROOT.TCanvas(f'{canvas_stem}_{branch}', title, 1000, 700)
        draw_root(root_individual)
        root_out.cd(); root_individual.Write(); root_individual.Close()
        with (out / f'{branch}_bins.csv').open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(['low_inclusive', 'high_exclusive', 'count'])
            writer.writerows(zip(edges[:-1], edges[1:], counts))
        metadata['histograms'].append(dict(branch=branch, bins=len(counts), low=float(edges[0]), high=float(edges[-1]),
                                            underflow=under, overflow=over, nonfinite=int(len(values)-len(finite)), log_y=log_y))
    if fig is not None:
        for ax in list(axes.flat)[len(summary_branches):]: ax.set_visible(False)
        fig.suptitle(f'Run {args.run_number} — ' + ('UNCUT raw data (all rows, no vetoes)' if profile == 'raw' else 'After global/run cuts and automatic event vetoes'))
        fig.tight_layout(rect=(0,0,1,.96))
        for ext in ('png', 'pdf'): fig.savefig(out / f'raw_overview.{ext}', dpi=160)
        root_out.cd(); root_summary.Write(); root_summary.Close()
    root_out.Close()
    if owns_input: root_file.Close()
    ROOT.gROOT.SetBatch(previous_root_batch)
    (out / 'raw_summary.json').write_text(json.dumps(metadata, indent=2) + '\n')
    if args.apply_cuts: database.snapshot(out / 'cuts_used.ini')
    print(f'Run {args.run_number}: {len(keep):,} raw rows; {keep.sum():,} plotted rows; selection={profile}; saved to {out}', flush=True)
    if defer_show: return figures
    if show_window:
        # Block until all inspection windows close; files are already saved.
        # Each figure stays alive for independent toolbar zoom/pan.
        plt.show(block=True)
    for figure in figures:
        plt.close(figure)


# ---- Two-dimensional inspection; one raw pulse fills one TH2 cell ----
def run_raw2d(args, *, inspection=None, defer_show=False):
    """Mapped coordinates, optional detector cleaning, native TH2s and canvases.

    Pair definitions and axis overrides are data in raw_2d_plots.csv. Adding
    another pair does not require another script or a new drawing function.
    """
    show_window = not args.batch
    if show_window and not defer_show:
        try:
            plt.switch_backend('TkAgg')
        except (ImportError, RuntimeError) as error:
            print(f'Cannot open plot windows: {error}. Saving outputs only.', flush=True)
            show_window = False
    owns_input = inspection is None
    data = inspection if inspection is not None else inspection_data(args)
    arrays, source = data['arrays'], data['source']
    definitions = load_2d_plots(args.plots_file)
    if args.plots:
        unknown = set(args.plots) - {p['Name'] for p in definitions}
        if unknown: raise ValueError(f'Unknown 2D plot names: {sorted(unknown)}')
        definitions = [p for p in definitions if p['Name'] in args.plots]
    database, keep = data['database'], data['keep']
    out = args.output_dir or results_directory(args.input_dir) / 'raw' / f'run{args.run_number}'
    if args.apply_cuts: out /= 'cleaned'
    out /= '2d'
    if args.plots: out /= 'plots_' + '_'.join(p['Name'] for p in definitions)
    out.mkdir(parents=True, exist_ok=True)
    profile = 'CLEANED' if args.apply_cuts else 'UNCUT'
    unmapped = int(((arrays['Pad'] < 0) & keep).sum())
    metadata = dict(run=args.run_number, selection=profile, total_rows=len(keep),
                    selected_rows=int(keep.sum()), selected_unmapped_hits=unmapped,
                    cuts_db=str(database.path) if database else None, histograms=[])
    for path, name in [(args.pad_map, 'pad_map_used.csv'), (args.plots_file, 'raw_2d_plots_used.csv'),
                       (args.binning_file, 'raw_binning_used.csv')]:
        (out / name).write_bytes(path.read_bytes())

    import ROOT
    from matplotlib.colors import LogNorm, Normalize
    previous_batch = ROOT.gROOT.IsBatch()
    ROOT.gROOT.SetBatch(True)
    root_file = ROOT.TFile(str(out / 'raw_2d_histograms.root'), 'RECREATE')
    figures, root_canvases = [], []
    summary = root_summary = None
    stem = f'TPC_raw2d_{args.run_number}' + ('_clean' if args.apply_cuts else '')
    try:
        if not args.plots:
            columns = min(2, len(definitions)); rows = math.ceil(len(definitions) / columns)
            summary, axes = plt.subplots(rows, columns, figsize=(7 * columns, 5.5 * rows), squeeze=False)
            figures.append(summary)
            root_summary = ROOT.TCanvas(stem, f'Run {args.run_number} | {profile} | 2D summary', columns * 650, rows * 500)
            root_summary.Divide(columns, rows); root_canvases.append(root_summary)
        ROOT.TNamed('selection', database.describe(args.run_number) if database else 'UNCUT: all raw hits; no vetoes').Write()
        ROOT.TNamed('mapping', 'Pad from established mapping; Pad=-1 for unconnected channels, retained in flow cells').Write()
        for index, plot in enumerate(definitions):
            name, xfield, yfield = plot['Name'], plot['X'], plot['Y']
            xedges, yedges = raw_2d_axes(arrays, plot, args.binning_file)
            nx, ny = len(xedges) - 1, len(yedges) - 1
            def label(field, edges):
                base = {'Pad': 'Pad number (0-based)', 'Event': 'Event ID', 'Hit': 'Hit (pulse index)'}.get(field, f'{field} (native values)')
                return base + f' [bin width {edges[1]-edges[0]:g}]'
            xlabel, ylabel = label(xfield, xedges), label(yfield, yedges)
            title = f'Run {args.run_number} | {profile} | {yfield} vs {xfield}'
            histogram = ROOT.TH2D(f'raw2d_{name}', f'{title};{xlabel};{ylabel};Raw Hits',
                                  nx, float(xedges[0]), float(xedges[-1]), ny, float(yedges[0]), float(yedges[-1]))
            histogram.SetDirectory(0); histogram.SetStats(False); histogram.SetOption('COLZ')
            x, y = arrays[xfield][keep], arrays[yfield][keep]
            finite = np.isfinite(x) & np.isfinite(y)
            # ROOT FillN defines the bin/flow conventions for BOTH versions.
            # It also retains true moments, entries and corner-flow cells.
            x = np.ascontiguousarray(x[finite], dtype=np.float64)
            y = np.ascontiguousarray(y[finite], dtype=np.float64)
            if len(x): histogram.FillN(len(x), x, y, np.ones(len(x), dtype=np.float64))
            full = np.frombuffer(histogram.GetArray(), dtype=np.float64,
                                 count=histogram.GetNcells()).reshape(ny + 2, nx + 2).copy()
            visible = full[1:-1, 1:-1]
            inside = int(visible.sum()); outside = len(x) - inside
            notes = [f'Outside plotted ranges: {outside:,}']
            if 'Pad' in (xfield, yfield): notes.append(f'Hits on unconnected channels: {unmapped:,}')
            maximum = float(visible.max(initial=0))
            log_z = plot['LogZ'] and maximum > 0
            norm = LogNorm(vmin=1, vmax=max(2., maximum)) if log_z else Normalize(vmin=0, vmax=max(1., maximum))
            def draw_axis(ax, figure):
                mesh = ax.pcolormesh(xedges, yedges, np.ma.masked_equal(visible, 0),
                                     cmap='viridis', norm=norm, shading='flat', rasterized=True)
                ax.set(xlabel=xlabel, ylabel=ylabel)
                ax.set_title(title + '\n' + '\n'.join(notes), fontsize=10)
                figure.colorbar(mesh, ax=ax, label='Raw Hits')
            if summary is not None: draw_axis(axes.flat[index], summary)
            figure, ax = plt.subplots(figsize=(10, 8)); figures.append(figure)
            draw_axis(ax, figure); figure.tight_layout()
            for ext in ('png', 'pdf'): figure.savefig(out / f'raw2d_{name}.{ext}', dpi=160)

            def draw_root(pad):
                pad.cd(); pad.SetLeftMargin(.13); pad.SetRightMargin(.17)
                pad.SetBottomMargin(.15); pad.SetTopMargin(.19); pad.SetLogz(log_z)
                if log_z:
                    histogram.SetMinimum(1); histogram.SetMaximum(max(2., maximum))
                histogram.DrawCopy('COLZ')
                text = ROOT.TLatex(); text.SetNDC(); text.SetTextSize(.028)
                for line, note in enumerate(notes): text.DrawLatex(.13, .91 - line * .045, note)
            if root_summary: draw_root(root_summary.cd(index + 1))
            individual = ROOT.TCanvas(f'{stem}_{name}', title, 1000, 800)
            root_canvases.append(individual); draw_root(individual)
            root_file.cd(); histogram.Write(); individual.Write()
            # Sparse CSV includes nonzero visible AND under/overflow cells.
            # Bin 0 is underflow; N+1 is overflow, including both-axis corners.
            def bounds(edges, number):
                if number == 0: return (-math.inf, edges[0])
                if number == len(edges): return (edges[-1], math.inf)
                return (edges[number-1], edges[number])
            with (out / f'{name}_bins.csv').open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(['x_bin', 'y_bin', 'x_low_inclusive', 'x_high_exclusive', 'y_low_inclusive', 'y_high_exclusive', 'count'])
                for iy, ix in np.argwhere(full > 0):
                    writer.writerow([ix, iy, *bounds(xedges, ix), *bounds(yedges, iy), int(full[iy, ix])])
            metadata['histograms'].append(dict(name=name, x=xfield, y=yfield, x_bins=nx, y_bins=ny,
                x_low=float(xedges[0]), x_high=float(xedges[-1]), y_low=float(yedges[0]), y_high=float(yedges[-1]),
                entries=len(x), visible_hits=inside, outside_hits=outside, nonfinite_hits=int((~finite).sum()),
                x_underflow=int(full[:, 0].sum()), x_overflow=int(full[:, -1].sum()),
                y_underflow=int(full[0, :].sum()), y_overflow=int(full[-1, :].sum()), log_z=plot['LogZ']))
        if summary is not None:
            for ax in list(axes.flat)[len(definitions):]: ax.set_visible(False)
            summary.tight_layout()
            for ext in ('png', 'pdf'): summary.savefig(out / f'raw2d_overview.{ext}', dpi=160)
            root_file.cd(); root_summary.Write()
        (out / 'raw_2d_summary.json').write_text(json.dumps(metadata, indent=2) + '\n')
        if database: database.snapshot(out / 'cuts_used.ini')
    except Exception:
        for figure in figures: plt.close(figure)
        raise
    finally:
        for canvas in root_canvases: canvas.Close()
        root_file.Close()
        if owns_input: source.Close()
        ROOT.gROOT.SetBatch(previous_batch)
    print(f'Run {args.run_number}: {profile}, {keep.sum():,} hits; {len(definitions)} 2D views saved to {out}', flush=True)
    if defer_show: return figures
    if show_window: plt.show(block=True)
    for figure in figures: plt.close(figure)


def run_inspection(args):
    """One analysis session: show both groups together, or select dimensions.

    The drawing functions remain separate edit points. Input/mapping/cuts
    are evaluated once, then all requested windows share one show() loop.
    """
    dimensions = args.dimensions
    if dimensions is None:
        dimensions = '1d' if args.branches and not args.plots else '2d' if args.plots and not args.branches else 'both'
    if dimensions == '1d' and args.plots:
        raise ValueError('--plots needs --2d or --dimensions both')
    if dimensions == '2d' and (args.branches or args.bins is not None or args.bin_width is not None or args.range or args.log_y is not None):
        raise ValueError('1D branch/bin options need --1d or --dimensions both; use --plots-file for 2D axis overrides')
    show_window = not args.batch
    if show_window:
        try:
            plt.switch_backend('TkAgg')
        except (ImportError, RuntimeError) as error:
            print(f'Cannot open plot windows: {error}. Saving outputs only.', flush=True)
            show_window = False
    data = inspection_data(args)
    figures = []
    try:
        if dimensions in ('both', '1d'):
            figures.extend(run_raw(args, inspection=data, defer_show=True))
        if dimensions in ('both', '2d'):
            figures.extend(run_raw2d(args, inspection=data, defer_show=True))
        if show_window: plt.show(block=True)
    finally:
        for figure in figures: plt.close(figure)
        data['source'].Close()


# ---- The only entry point: commands or a small menu ----
def common_arguments(parser):
    """Keep the same path options for event and example workflows."""
    parser.add_argument('--input-dir', type=Path, default=DEFAULT_INPUT_DIR,
                        help='Directory containing ntupleNNNN.root; defaults to data/root/')
    parser.add_argument('--output-dir', type=Path,
                        help='Where to save plots; default under results/python/ for the chosen workflow')
    parser.add_argument('--pad-map', type=Path, default=DEFAULT_PAD_MAP)
    parser.add_argument('--cuts-db', type=Path, default=DEFAULT_CUTS_DB, help='Shared detector global/run rejection database')
    parser.add_argument('--veto-file', type=Path,
                        help='Run/event veto CSV; otherwise use live output or bundled auxiliary snapshot')


def make_parser():
    parser = argparse.ArgumentParser(description='TPC analysis: events, raw inspection, examples and voltage scans.',
                                     epilog='Run without arguments for a menu. Edit config/run_cuts.ini for manual cuts.')
    commands = parser.add_subparsers(dest='command')
    event = commands.add_parser('event', help='Browse events, or save one with --save')
    event.add_argument('run_number', type=int, help='Run number, for example 1616')
    event.add_argument('--event', type=int, help='Event to start with or export')
    event.add_argument('--mode', choices=['2d', '3d'], default='2d')
    event.add_argument('--save', action='store_true', help='Export without opening windows')
    common_arguments(event)
    examples = commands.add_parser('examples', help='Save the established 1615/1616 examples and comparison')
    examples.add_argument('--example', action='append', metavar='RUN:EVENT',
                          help='Repeat to choose your own examples; default is 1615:42535 and 1616:5584')
    common_arguments(examples)
    raw = commands.add_parser('raw', help='Inspect 1D and 2D together; no cuts/vetoes unless --apply-cuts')
    raw.add_argument('run_number', type=int)
    raw.add_argument('--batch', action='store_true', help='Save outputs without opening a plot window')
    dimensions = raw.add_mutually_exclusive_group()
    dimensions.add_argument('--dimensions', '--dim', choices=['both', '1d', '2d'],
                            help='Default: both; --branches or --plots alone focuses its own group')
    dimensions.add_argument('--1d', dest='dimensions', action='store_const', const='1d', help='Show only 1D')
    dimensions.add_argument('--2d', dest='dimensions', action='store_const', const='2d', help='Show only 2D')
    dimensions.add_argument('--both', dest='dimensions', action='store_const', const='both', help='Show both groups (default)')
    raw.add_argument('--plots', nargs='+', help='Focus on selected 2D plot names')
    raw.add_argument('--plots-file', type=Path, default=DEFAULT_2D_PLOTS)
    raw.add_argument('--branches', nargs='+', help='Focus on individual quantities; default: six-panel mapped summary plus separate Channel diagnostic')
    bin_choice = raw.add_mutually_exclusive_group()
    bin_choice.add_argument('--bins', type=int, help='Override 1D bin count')
    bin_choice.add_argument('--bin-width', type=float, help='Override 1D native bin width')
    raw.add_argument('--range', type=float, nargs=2, metavar=('LOW', 'HIGH'), help='1D histogram edges, not a data cut')
    scale = raw.add_mutually_exclusive_group()
    scale.add_argument('--log-y', dest='log_y', action='store_true')
    scale.add_argument('--linear-y', dest='log_y', action='store_false')
    raw.set_defaults(log_y=None)
    raw.add_argument('--binning-file', type=Path, default=DEFAULT_BINNING)
    raw.add_argument('--apply-cuts', action='store_true', help='Compare cleaned data using DB plus automatic vetoes')
    common_arguments(raw)
    raw2d = commands.add_parser('raw2d', help='Shortcut for raw --dimensions 2d')
    raw2d.add_argument('run_number', type=int)
    raw2d.add_argument('--batch', action='store_true', help='Save without opening windows')
    raw2d.add_argument('--apply-cuts', action='store_true', help='Apply detector global/run DB and automatic vetoes')
    raw2d.add_argument('--plots', nargs='+', help='Focus on names from raw_2d_plots.csv')
    raw2d.add_argument('--plots-file', type=Path, default=DEFAULT_2D_PLOTS)
    raw2d.add_argument('--binning-file', type=Path, default=DEFAULT_BINNING)
    common_arguments(raw2d)
    cuts = commands.add_parser('cuts', help='List the DB rejection rules active for a run')
    cuts.add_argument('run_number', type=int)
    cuts.add_argument('--cuts-db', type=Path, default=DEFAULT_CUTS_DB)
    scan = commands.add_parser('voltage' , help='Independent voltage study and noise audits')
    # ROOT is needed only by this workflow or by the lazy input reader.
    # Keeping its lengthy plotting code separate keeps routine editing simple.
    import voltage
    voltage.add_arguments(scan)
    info = commands.add_parser('info', help='Print current cuts, mapping size, and input locations')
    info.add_argument('--input-dir', type=Path, default=DEFAULT_INPUT_DIR)
    info.add_argument('--cuts-db', type=Path, default=DEFAULT_CUTS_DB)
    commands.add_parser('check', help='Run the detector/mapping regression checks')
    return parser


def menu_arguments():
    """No command memorization needed: turn menu choices into the same CLI."""
    print('\nTPC analysis')
    print('  1  Browse events')
    print('  2  Save one event')
    print('  3  Recreate the 1615/1616 examples (2D and 3D)')
    print('  4  Run the independent voltage study')
    print('  5  Show current cuts and locations')
    print('  6  Check the mapping and selection')
    print('  7  Inspect raw data (1D + 2D)')
    print('  0  Quit')
    choice = input('Choose [0]: ').strip() or '0'
    if choice == '0':
        return None
    if choice in ('1', '2'):
        run = input('Run number [1616]: ').strip() or '1616'
        # The known example is a convenience, not an event-selection cut.
        suggested = {'1615': '42535', '1616': '5584', '1610': '687'}.get(run, '')
        event = input(f'Event [{suggested or "first accepted event"}]: ').strip() or suggested
        args = ['event', run]
        if event:
            args += ['--event', event]
        if choice == '2':
            args += ['--save']
        return args
    if choice == '7':
        run = input('Run number [1616]: ').strip() or '1616'
        return ['raw', run]
    commands = {'3': 'examples', '4': 'voltage', '5': 'info', '6': 'check'}
    if choice not in commands:
        raise ValueError('Please choose a number from 0 to 7')
    return [commands[choice]]


def main(argv=None):
    parser = make_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        # Piped/non-interactive execution should print help instead of hanging.
        if not sys.stdin.isatty():
            parser.print_help()
            return
        argv = menu_arguments()
        if argv is None:
            return
    args = parser.parse_args(argv)
    if args.command == 'event':
        run_event(args)
    elif args.command == 'examples':
        run_examples(args)
    elif args.command == 'raw':
        run_inspection(args)
    elif args.command == 'raw2d':
        run_raw2d(args)
    elif args.command == 'cuts':
        database = cut_database(args.cuts_db)
        print(database.path)
        for rule in database.active(args.run_number):
            print(f'[{rule.section}] line {rule.line}: REJECT {rule.expression}' + (f' # {rule.note}' if rule.note else ''))
        print(f'\nRules listed only; no plots opened. To apply cuts and plot: python analysis/tpc.py raw {args.run_number} --apply-cuts')
        print('Keep the same --cuts-db override if you supplied one.')
    elif args.command == 'voltage':
        import voltage
        voltage.run(args)
    elif args.command == 'info':
        print(f'Settings: {SETTINGS_FILE}')
        print(f'Mapping: {DEFAULT_PAD_MAP} ({len(load_pad_records())} pads)')
        print(f'ROOT input: {args.input_dir.resolve()}')
        print(f'Veto ledger: {default_veto_file(args.input_dir)}')
        database = cut_database(args.cuts_db)
        print(f'Cut DB: {database.path}')
        print(f'Event displays: {database.describe()}')
        print(f'Independent voltage-study cuts: {DEFAULT_VOLTAGE_STUDY_DB}')
        print('Event displays and cleaned raw plots use only the detector global/run DB and automatic vetoes.')
    elif args.command == 'check':
        import unittest
        tests = unittest.defaultTestLoader.discover(str(Path(__file__).parent / 'tests'))
        result = unittest.TextTestRunner(verbosity=2).run(tests)
        if not result.wasSuccessful():
            raise SystemExit(1)


if __name__ == '__main__':
    main()
