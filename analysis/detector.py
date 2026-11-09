"""Detector data and selection, independent of the user interface.

Routine use: python analysis/tpc.py (menu), or tpc.py --help.
Edit ../config/run_cuts.ini for manual cuts and settings.csv for automatic
noise thresholds. The established auxiliary/pad_map.csv is the
mapping source; no builder is needed during analysis.

For later reconstruction, load_hit_arrays() retains Event/Channel/Hit/
Time/Charge and select_hit_mask() keeps ALL passing pulses. The first-hit
reduction in select_event_hits() is ONLY for current pad displays. Pad
coordinates are transverse positions; native Time is not calibrated z.
"""
import csv
import re
import math
from dataclasses import dataclass
from pathlib import Path
import numpy as np

# Resolve paths from this file, never from the shell's working directory.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
AUXILIARY_DIR = SCRIPT_DIR / 'auxiliary'
DEFAULT_INPUT_DIR = PROJECT_DIR / 'data' / 'root'
# A shared source folder can still be copied beside a flat ROOT dataset.
if not DEFAULT_INPUT_DIR.is_dir():
    DEFAULT_INPUT_DIR = PROJECT_DIR
DEFAULT_PAD_MAP = AUXILIARY_DIR / 'pad_map.csv'
SETTINGS_FILE = SCRIPT_DIR / 'settings.csv'
DEFAULT_CUTS_DB = PROJECT_DIR / 'config' / 'run_cuts.ini'
DEFAULT_VOLTAGE_STUDY_DB = PROJECT_DIR / 'config' / 'voltage_study_cuts.ini'
DEFAULT_BINNING = PROJECT_DIR / 'config' / 'raw_binning.csv'
DEFAULT_2D_PLOTS = PROJECT_DIR / 'config' / 'raw_2d_plots.csv'
with SETTINGS_FILE.open(newline='') as stream:
    SETTINGS = {r['Key']: r['Value'] for r in csv.DictReader(stream)}

# Histogram limits and automatic-noise settings remain here. Manual
# Time/Channel/etc. rejections come exclusively from the shared cut DB.
TIME_CUT_MAX = int(SETTINGS['time_max'])
NOISE_CHANNELS = int(SETTINGS['noise_channels'])
MODULE_CHANNELS = int(SETTINGS['module_channels'])
MODULE_SYNCHRONOUS = int(SETTINGS['module_synchronous'])
MODULE_HITS_PER_CHANNEL = int(SETTINGS['module_hits_per_channel'])
CHANNELS_PER_MODULE = int(SETTINGS['channels_per_module'])
READOUT_MODULES = int(SETTINGS['readout_modules'])
EXCLUDED_RUNS = [int(r) for r in SETTINGS['excluded_runs'].split(';') if r]
# ---- Input and run-specific noise vetoes ----

def results_directory(input_dir=DEFAULT_INPUT_DIR):
    """Keep outputs outside data/root/ in the organized workspace.

    For a different organized dataset, infer its project from data/root/.
    A flat external dataset gets its own results/python/ directory.
    Every export command also accepts an explicit --output-dir override.
    """
    source = Path(input_dir).resolve()
    if source == DEFAULT_INPUT_DIR.resolve():
        project = PROJECT_DIR
    elif source.name == 'root' and source.parent.name == 'data':
        project = source.parent.parent
    else:
        project = source
    return project / 'results' / 'python'


def default_veto_file(input_dir=DEFAULT_INPUT_DIR):
    """Use a freshly generated run ledger when present, else the bundled copy.

    The CSV is indexed by (run,event), never by event number alone. Use an
    explicit --veto-file when working with unrelated data having reused run
    numbers. auxiliary/veto_provenance.json documents the bundled snapshot.
    """
    generated = results_directory(input_dir) / 'voltage' / 'rejected_events.csv'
    if generated.exists():
        return generated
    # Accept the former layout when using a separate, older dataset.
    legacy = Path(input_dir) / 'time_voltage_plots' / 'rejected_events.csv'
    return legacy if legacy.exists() else AUXILIARY_DIR / 'rejected_events.csv'

def load_hit_arrays(root_filename, tree_name='ntuple',
                    branches=None):
    import ROOT
    root_file = ROOT.TFile.Open(str(root_filename), 'READ')
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f'Could not open {root_filename}')
    tree = root_file.Get(tree_name)
    if not tree:
        root_file.Close()
        raise ValueError(f'Tree {tree_name!r} is missing from {root_filename}')
    if branches is None:
        numeric = {'Int_t', 'UInt_t', 'Long64_t', 'ULong64_t', 'Short_t', 'UShort_t',
                   'Char_t', 'UChar_t', 'Bool_t', 'Float_t', 'Double_t', 'Long_t', 'ULong_t'}
        branches = [b.GetName() for b in tree.GetListOfBranches()
                    if b.GetLeaf(b.GetName()) and b.GetLeaf(b.GetName()).GetLenStatic() == 1
                    and not b.GetLeaf(b.GetName()).GetLeafCount()
                    and b.GetLeaf(b.GetName()).GetTypeName() in numeric]
    missing = [name for name in branches if not tree.GetBranch(name)]
    if missing:
        root_file.Close()
        raise ValueError(f'Missing branches in {root_filename}: {missing}')
    arrays = ROOT.RDataFrame(tree).AsNumpy(list(branches))
    return arrays, root_file

def rejected_events_for_run(run, path=None):
    path = Path(path) if path else default_veto_file()
    if not path.exists():
        raise FileNotFoundError(f"Missing veto ledger: {path}. Generate it with tpc.py voltage or pass --veto-file.")
    with path.open(newline='') as stream:
        return {int(row['event']) for row in csv.DictReader(stream) if int(row['run']) == run}

# ---- Manual cut database: no eval(), no arbitrary ROOT expressions ----
_NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
_FIELDS = {name.lower(): name for name in
           ('Event', 'Channel', 'Hit', 'Time', 'Charge', 'Pedestal', 'Blockhead', 'FADC', 'Pad')}


@dataclass(frozen=True)
class CutCondition:
    field: str
    operator: str
    values: tuple
    is_range: bool = False

    def matches(self, values):
        if self.operator in ('in', 'not in'):
            answer = ((values >= self.values[0]) & (values <= self.values[1])
                      if self.is_range else np.isin(values, self.values))
            return ~answer if self.operator == 'not in' else answer
        value = self.values[0]
        return {'==': np.equal, '!=': np.not_equal, '<': np.less, '<=': np.less_equal,
                '>': np.greater, '>=': np.greater_equal}[self.operator](values, value)


@dataclass(frozen=True)
class CutRule:
    section: str
    expression: str
    conditions: tuple
    line: int
    note: str = ''


class CutDatabase:
    """Human-editable rejection rules, shared with the native ROOT program.

    [global] always applies. [run N] applies only to that run. Independent
    studies load their own additional file, outside this detector database.
    ANY matching reject line removes a row; conditions joined by && on
    one line must ALL match. Ranges are inclusive. No expression is eval'd.
    """
    def __init__(self, path=DEFAULT_CUTS_DB):
        self.path = Path(path).resolve()
        self.text = self.path.read_text()
        self.rules = []
        section = None
        sections = set()
        for line_number, original in enumerate(self.text.splitlines(), 1):
            code, _, note = original.partition('#')
            code = code.strip()
            if not code:
                continue
            try:
                if code.startswith('['):
                    match = re.fullmatch(r'\[\s*(global|run\s+\d+)\s*\]', code, re.I)
                    if not match:
                        raise ValueError('Expected [global] or [run NNNN]; study cuts belong in a separate file')
                    section = ' '.join(match[1].lower().split())
                    if section.startswith('run '): section = f'run {int(section[4:])}'
                    if section in sections:
                        raise ValueError(f'Duplicate block [{section}]')
                    sections.add(section)
                    continue
                if section is None:
                    raise ValueError('Put rules inside a named block')
                match = re.fullmatch(r'reject\s*=\s*(.+)', code, re.I)
                if not match:
                    raise ValueError('Expected reject = CONDITION')
                conditions = tuple(self._condition(term.strip()) for term in match[1].split('&&'))
                self.rules.append(CutRule(section, match[1].strip(), conditions, line_number, note.strip()))
            except ValueError as error:
                raise ValueError(f'{self.path}:{line_number}: {error}') from error

    @staticmethod
    def _condition(term):
        match = re.fullmatch(r'([A-Za-z_][A-Za-z_0-9]*)\s*(not\s+in|in|==|!=|<=|>=|<|>)\s*(.+)', term)
        if not match or match[1].lower() not in _FIELDS:
            raise ValueError(f'Invalid field or condition: {term!r}')
        field, operator, rhs = _FIELDS[match[1].lower()], ' '.join(match[2].split()), match[3].strip()
        def number(text):
            if not re.fullmatch(_NUMBER, text.strip()):
                raise ValueError(f'Expected a finite number, got {text!r}')
            value = float(text)
            if not math.isfinite(value):
                raise ValueError('Non-finite cut value')
            return value
        if operator in ('in', 'not in'):
            limits = re.fullmatch(rf'({_NUMBER})\s*\.\.\s*({_NUMBER})', rhs)
            if limits:
                values = (number(limits[1]), number(limits[2]))
                if values[0] > values[1]:
                    raise ValueError('Range lower endpoint exceeds upper endpoint')
                return CutCondition(field, operator, values, True)
            return CutCondition(field, operator, tuple(number(x) for x in rhs.split(',')))
        return CutCondition(field, operator, (number(rhs),))

    def active(self, run=None):
        sections = {'global'}
        if run is not None:
            sections.add(f'run {int(run)}')
        return [rule for rule in self.rules if rule.section in sections]

    def describe(self, run=None):
        return ' && '.join(f'NOT({r.expression})' for r in self.active(run)) or 'all rows'

    def keep(self, arrays, run=None, pad_map=DEFAULT_PAD_MAP):
        """Return one bool/raw row; never aggregate pulses or expand event IDs."""
        fields = {name: np.asarray(values) for name, values in arrays.items()}
        size = len(next(iter(fields.values())))
        rejected = np.zeros(size, dtype=bool)
        for rule in self.active(run):
            bad = np.ones(size, dtype=bool)
            for condition in rule.conditions:
                name = condition.field
                if name not in fields and name == 'FADC' and 'Channel' in fields:
                    fields[name] = fields['Channel'] + 1
                if name not in fields and name == 'Pad' and 'Channel' in fields:
                    fields[name] = channel_to_pad(fields['Channel'], pad_map)
                if name not in fields:
                    raise ValueError(f'{self.path}:{rule.line}: branch {name} is missing from this input')
                bad &= condition.matches(fields[name])
            rejected |= bad
        return ~rejected

    def snapshot(self, destination):
        """Save the text actually parsed, even if the on-disk DB changes later."""
        Path(destination).write_text(self.text)


_DEFAULT_DATABASE = None

def cut_database(database=None):
    global _DEFAULT_DATABASE
    if isinstance(database, CutDatabase):
        return database
    if database is not None:
        return CutDatabase(database)
    if _DEFAULT_DATABASE is None:
        _DEFAULT_DATABASE = CutDatabase()
    return _DEFAULT_DATABASE


# Common cut expressions now come ONLY from the DB, not settings.csv.
BASE_CUT = 'Global and run-specific cuts from run_cuts.ini'
EVENT_CUT = BASE_CUT
VOLTAGE_CUT = 'Detector global/run DB + independent voltage-study DB'


def passes_hit_cuts(channel, time):
    """Legacy scalar helper for global Channel/Time checks; use full masks for analysis."""
    return bool(cut_database().keep({'Channel': [channel], 'Time': [time]})[0])


def select_hit_mask(events, channels, times, rejected_events=(), voltage=False,
                    run=None, cuts_db=None, arrays=None, pad_map=DEFAULT_PAD_MAP,
                    study_cuts_db=None):
    """Preserve ALL passing raw pulses; combine detector DB and run auto vetoes.

    voltage=True opts into the separate study_cuts_db, used only by the
    independent voltage analysis. General detector cleaning leaves it off.
    """
    fields = dict(arrays or {})
    fields.update(Event=np.asarray(events), Channel=np.asarray(channels), Time=np.asarray(times))
    mask = cut_database(cuts_db).keep(fields, run, pad_map=pad_map)
    # Only the independent voltage workflow opts into this second DB.
    # Ordinary cleaning never opens the study file, even if it is missing.
    if voltage:
        study = cut_database(study_cuts_db if study_cuts_db is not None else DEFAULT_VOLTAGE_STUDY_DB)
        mask &= study.keep(fields, run, pad_map=pad_map)
    return mask & ~np.isin(events, list(rejected_events))


def select_event_hits(events, channels, times, charges, target_event, rejected_events=(),
                      run=None, cuts_db=None, arrays=None, pad_map=DEFAULT_PAD_MAP,
                      accepted_mask=None):
    """First passing hit/channel for DISPLAY only, after DB cuts and event vetoes.

    Applying the DB before seen-channel reduction matters: removing a first
    pulse must not hide a later acceptable pulse on that same channel.
    """
    if accepted_mask is None:
        fields = dict(arrays or {})
        fields['Charge'] = np.asarray(charges)
        accepted_mask = select_hit_mask(events, channels, times, rejected_events,
                                        run=run, cuts_db=cuts_db, arrays=fields, pad_map=pad_map)
    seen, result = set(), []
    for event, channel, time, charge, keep in zip(events, channels, times, charges, accepted_mask):
        if not keep or event != target_event or channel in seen:
            continue
        seen.add(channel)
        result.append((channel, charge, time))
    return result

# ---- Established pad mapping and geometry ----
MED_ROWS, MED_COLS, MED_PAD_MM = 10, 25, 4.0
LG_ROWS, LG_COLS, LG_PAD_MM = 8, 20, 5.0
BOARD_WIDTH_MM, BOARD_HEIGHT_MM = 100.0, 80.0

def load_pad_records(csv_path=DEFAULT_PAD_MAP):
    """Full CSV records include pad ID, Channel, row/column, and mm bounds.

    A later track fit can derive pad centers as X_mm + Width_mm/2 and
    Y_mm + Height_mm/2, while retaining pad widths for uncertainties. Time
    is still in native units; z requires drift calibration and t0, neither
    of which this geometry module assumes or calculates.
    """
    with open(csv_path, newline="") as stream:
        records = list(csv.DictReader(stream))
    for r in records:
        for key in ("Pad", "FADC", "Channel", "Row", "Col"):
            r[key] = int(r[key])
        for key in ("X_mm", "Y_mm", "Width_mm", "Height_mm"):
            r[key] = float(r[key])
    if sorted(r['Pad'] for r in records) != list(range(410)):
        raise ValueError("Expected exactly 410 continuous pad numbers, 0–409")
    if len({r['Channel'] for r in records}) != len(records):
        raise ValueError("Each electronics channel must map to only one pad")
    return sorted(records, key=lambda r: r['Pad'])

def load_pad_map(csv_path=DEFAULT_PAD_MAP):
    """ROOT Channel -> (section, bottom-based row, column)."""
    return {r['Channel']: (r['PadSize'], r['Row'], r['Col']) for r in load_pad_records(csv_path)}

def channel_to_pad(channels, csv_path=DEFAULT_PAD_MAP):
    """Vectorized lookup of zero-based ROOT channels; -1 for unmapped values."""
    records = load_pad_records(csv_path)
    lookup = np.full(max(r['Channel'] for r in records) + 1, -1, dtype=np.int32)
    for r in records:
        lookup[r['Channel']] = r['Pad']
    values = np.asarray(channels)
    pads = np.full(values.shape, -1, dtype=np.int32)
    valid = (values >= 0) & (values < len(lookup)) & (values == np.floor(values))
    pads[valid] = lookup[values[valid].astype(int)]
    return pads

def pad_number(pad_size, row, col):
    if pad_size == "Medium":
        return row * MED_COLS + col
    if pad_size == "Large":
        return 250 + row * LG_COLS + col
    raise ValueError(f"Unknown pad size: {pad_size}")

def pad_rect_mm(pad_size, row, col):
    # Detector coordinates: origin at lower-left medium pad, x to the
    # right, y upward. This is a pad footprint, NOT a 3D track coordinate.
    if pad_size == "Medium":
        return col * MED_PAD_MM, row * MED_PAD_MM, MED_PAD_MM, MED_PAD_MM
    if pad_size == "Large":
        return col * LG_PAD_MM, 40 + row * LG_PAD_MM, LG_PAD_MM, LG_PAD_MM
    raise ValueError(f"Unknown pad size: {pad_size}")

# ---- Whole-event noise analysis, before any hit cuts ----
def module_floods(event, channel, time, pairs, min_channels, min_synchronous):
    """Detect module flooding BEFORE time cuts, channel cuts, or geometry.

    A module needs >=min_channels distinct channels AND either enough
    simultaneous channels or >=MODULE_HITS_PER_CHANNEL hits per channel.
    The returned dictionary retains evidence for each (event,module).
    Electronics modules are not the physical pad rows in the wiring map.
    """
    groups, occupancy = np.unique(pairs[:, 0].astype(np.int64) * READOUT_MODULES + pairs[:, 1] // CHANNELS_PER_MODULE,
                                  return_counts=True)
    row_groups = event.astype(np.int64) * READOUT_MODULES + channel // CHANNELS_PER_MODULE
    keys, rows = np.unique(row_groups, return_counts=True)
    details = {}
    for key, nch in zip(groups[occupancy >= min_channels], occupancy[occupancy >= min_channels]):
        mask = row_groups == key
        # Distinct channels per exact time, so repeated hits cannot inflate synchrony.
        tc = np.unique(np.column_stack((time[mask], channel[mask])), axis=0)
        _, synchronous = np.unique(tc[:, 0], return_counts=True)
        nrows = int(rows[np.searchsorted(keys, key)])
        nsync = int(synchronous.max())
        if nsync >= min_synchronous or nrows >= MODULE_HITS_PER_CHANNEL * nch:
            ev, module = divmod(int(key), READOUT_MODULES)
            details.setdefault(ev, []).append(dict(module=module, channels=int(nch),
                                                   hits=nrows, synchronous_channels=nsync))
    return details


# ---- Shared 1D/2D histogram axes, independent of hit selection ----
def plot_config_rows(path):
    """Read comma-separated plot settings with optional spaces/tabs for alignment.

    Commas, including consecutive commas for empty fields, are mandatory.
    Whitespace surrounding a header/value is formatting, not part of its name.
    """
    with Path(path).open(newline='') as stream:
        for row in csv.DictReader(stream, skipinitialspace=True):
            if None in row and any(value.strip() for value in row[None]):
                raise ValueError(f'{path}: too many CSV fields; keep empty columns and commas aligned')
            yield {key.strip(): (value or '').strip() for key, value in row.items() if key is not None}


def raw_histogram_spec(values, branch, bins=None, bin_width=None, limits=None, log_y=None,
                       binning_file=DEFAULT_BINNING, overrides=None):
    """Return explicit bin EDGES and scale. All histogram upper edges are exclusive.

    Default integer fields use one native value per bin. Channel includes
    the configured full readout span; Time/Event/Hit start at zero where
    possible. Charge/Pedestal retain their observed negative values.
    A plotting range changes only the view: under/overflow are reported.
    """
    row = {}
    for item in plot_config_rows(binning_file):
        if item['Branch'] == branch:
            row = item
            break
    if overrides:
        # Per-axis 2D choices override only nonempty fields. A new width
        # replaces an inherited bin count; an explicit count takes precedence.
        if overrides.get('BinWidth') and not overrides.get('Bins'): row['Bins'] = ''
        row.update({key: value for key, value in overrides.items() if value != ''})
    finite = np.asarray(values)[np.isfinite(values)]
    low = float(finite.min()) if finite.size else 0.
    high = float(finite.max()) if finite.size else 0.
    if branch in ('Event', 'Channel', 'Hit', 'Time', 'Blockhead', 'FADC', 'Pad'):
        low = min(0., low)
    if branch == 'Channel':
        high = max(high, CHANNELS_PER_MODULE * READOUT_MODULES - 1)
    low, high = math.floor(low) - .5, math.ceil(high) + .5
    if row.get('Min', '').strip(): low = float(row['Min'])
    if row.get('Max', '').strip(): high = float(row['Max'])
    if limits is not None: low, high = map(float, limits)
    if not math.isfinite(low) or not math.isfinite(high) or high <= low:
        raise ValueError(f'{branch}: histogram edges must be finite and increasing')
    if bins is None and row.get('Bins', '').strip(): bins = int(row['Bins'])
    if bin_width is not None: bins = None
    if bins is None:
        width = float(bin_width if bin_width is not None else row.get('BinWidth') or 1)
        if not math.isfinite(width) or width <= 0: raise ValueError('Bin width must be positive')
        bins = math.ceil((high - low) / width)
        high = low + bins * width
    if int(bins) != bins or not 1 <= int(bins) <= 1000000:
        raise ValueError(f'{branch}: bin count must be 1..1000000; choose --bins, --bin-width, or --range')
    if log_y is None:
        scale = row.get('LogY', 'false').strip().lower()
        if scale not in ('true', 'false', '1', '0'): raise ValueError('LogY must be true or false')
        log_y = scale in ('true', '1')
    return np.linspace(low, high, int(bins) + 1), bool(log_y)


# ---- Declarative 2D plot definitions; shared with native ROOT ----
def load_2d_plots(path=DEFAULT_2D_PLOTS):
    """Each CSV row defines Name, X, Y, per-axis overrides, and LogZ.

    Empty axis options inherit raw_binning.csv. Coarser Event/Charge defaults
    keep the 2D canvases practical; the 1D binning is never changed here.
    """
    plots, names = [], set()
    for line, original in enumerate(plot_config_rows(path), 2):
        row = {k: (v or '').strip() for k, v in original.items() if k}
        if not any(row.values()):
            continue
        try:
            name = row.get('Name', '')
            if not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', name) or name in names:
                raise ValueError('Plot names must be unique filename-safe identifiers')
            for axis in ('X', 'Y'):
                field = row.get(axis, '')
                if field.lower() not in _FIELDS:
                    raise ValueError(f'Unknown {axis} field: {field}')
                row[axis] = _FIELDS[field.lower()]
            scale = row.get('LogZ', 'true').lower()
            if scale not in ('true', 'false', '1', '0'):
                raise ValueError('LogZ must be true or false')
            row['LogZ'] = scale in ('true', '1')
            names.add(name); plots.append(row)
        except ValueError as error:
            raise ValueError(f'{path}:{line}: {error}') from error
    if not plots:
        raise ValueError(f'No 2D plot definitions in {path}')
    return plots


def raw_2d_axes(arrays, plot, binning_file=DEFAULT_BINNING):
    """Uniform ROOT-compatible edges, always derived from the uncut input.

    Min/Max change plotting limits only. The ROOT histograms keep all flow
    cells, including corners. Guard allocation before creating a dense TH2.
    """
    axes = []
    for axis in ('X', 'Y'):
        field = plot[axis]
        if field not in arrays:
            raise ValueError(f'{plot["Name"]}: missing branch {field}')
        options = {key: plot.get(axis + key, '') for key in ('BinWidth', 'Bins', 'Min', 'Max')}
        edges, _ = raw_histogram_spec(arrays[field], field, binning_file=binning_file, overrides=options)
        axes.append(edges)
    if (len(axes[0]) + 1) * (len(axes[1]) + 1) > 4_000_000:
        raise ValueError(f'{plot["Name"]}: more than 4 million 2D cells; increase BinWidth or narrow limits in raw_2d_plots.csv')
    return tuple(axes)
