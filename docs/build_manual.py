#!/usr/bin/env python3
"""Refresh current configuration/help/examples and compile the maintained manual.
Run from any folder: python docs/build_manual.py [--check].
Compilation uses a temporary directory; only the PDF and intentional sources remain.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
CONFIGS = ['run_cuts.ini', 'voltage_study_cuts.ini', 'raw_binning.csv', 'raw_2d_plots.csv']
FIGURES = {
    'voltage_overlay.png': 'results/root/voltage/time_overlay_counts.png',
    'raw_1d.png': 'results/root/raw/run1616/raw_overview.png',
    'raw_2d.png': 'results/root/raw/run1616/2d/raw2d_overview.png',
    'event_charge.png': 'results/root/events/run1616/event_1616_5584_charge.png',
    'event_charge_3d.png': 'results/root/events/run1616/3d/event_1616_5584_charge.png',
}

def inputs():
    paths = [HERE/'TPC_manual.tex', Path(__file__).resolve()]
    paths += [PROJECT/'analysis'/name for name in ('tpc.py','detector.py','voltage.py','settings.csv')]
    paths += [PROJECT/'analysis_root'/name for name in ('TPC.C','TPCDetector.h','settings.csv')]
    paths += [PROJECT/'config'/name for name in CONFIGS]
    paths += [PROJECT/engine/'auxiliary/pad_map.csv' for engine in ('analysis','analysis_root')]
    result={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.exists()}
    for name,source in FIGURES.items():
        path=PROJECT/source if (PROJECT/source).exists() else HERE/'figures'/name
        result[source]=hashlib.sha256(path.read_bytes()).hexdigest()
    return result

def tex_escape(value):
    return ''.join({'_':r'\_', '&':r'\&', '%':r'\%', '#':r'\#', '$':r'\$', '{':r'\{', '}':r'\}'}.get(c,c) for c in str(value))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Fail if source/config/example changes have made the PDF stale')
    args=parser.parse_args()
    manifest=HERE/'generated/manifest.json'
    current=inputs()
    if args.check:
        old=json.loads(manifest.read_text()) if manifest.exists() else {}
        # A source bundle intentionally carries only its own engine.
        for engine in ('analysis','analysis_root'):
            if not (PROJECT/engine).exists():old={k:v for k,v in old.items() if not k.startswith(engine+'/')}
        changed=sorted(k for k in set(current)|set(old) if current.get(k)!=old.get(k))
        if changed or not (HERE/'TPC_manual.pdf').exists():
            print('Manual needs rebuilding: '+', '.join(changed or ['missing PDF']))
            return 1
        print('Manual PDF is current with scripts, settings, mapping and example figures.')
        return 0
    generated=HERE/'generated';generated.mkdir(exist_ok=True)
    figures=HERE/'figures';figures.mkdir(exist_ok=True)
    for name in CONFIGS:shutil.copy2(PROJECT/'config'/name, generated/name)
    for name,source in FIGURES.items():
        if (PROJECT/source).exists():shutil.copy2(PROJECT/source, figures/name)
    for engine in ('analysis','analysis_root'):
        if (PROJECT/engine/'settings.csv').exists():shutil.copy2(PROJECT/engine/'settings.csv',generated/(engine+'_settings.csv'))
    environment=dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLCONFIGDIR='/tmp/cosmic-matplotlib')
    for command in ('raw','event','voltage'):
        if not (PROJECT/'analysis/tpc.py').exists():continue  # Reuse bundled help in ROOT-only distributions.
        result=subprocess.run([sys.executable,str(PROJECT/'analysis/tpc.py'),command,'--help'],text=True,capture_output=True,env=environment,check=True)
        (generated/(command+'_help.txt')).write_text(result.stdout)
    # Include the user's actual palette definition, without duplicating RGB values in prose.
    if (PROJECT/'analysis_root/TPC.C').exists():
        source=(PROJECT/'analysis_root/TPC.C').read_text()
        start=source.index('void set_color_env(');opening=source.index('{',start);depth=1;end=opening+1
        while depth:
            depth += (source[end]=='{')-(source[end]=='}');end+=1
        (generated/'palette.C').write_text(source[start:end]+'\n')
    with (PROJECT/'config/raw_2d_plots.csv').open(newline='') as stream:
        rows=[{k.strip():(v or '').strip() for k,v in r.items()} for r in csv.DictReader(stream)]
    table=[r'\begin{longtable}{p{.27\linewidth}p{.12\linewidth}p{.12\linewidth}p{.33\linewidth}}',r'\textbf{Name} & \textbf{X} & \textbf{Y} & \textbf{Overrides} \\\hline']
    for row in rows:
        options=', '.join(f'{k}={v}' for k,v in row.items() if v and k not in ('Name','X','Y'))
        table.append(' & '.join(tex_escape(v) for v in (row['Name'],row['X'],row['Y'],options))+r' \\')
    table.append(r'\end{longtable}')
    (generated/'pairs_table.tex').write_text('\n'.join(table)+'\n')
    revision=hashlib.sha256(json.dumps(current,sort_keys=True).encode()).hexdigest()[:12]
    (generated/'revision.tex').write_text(r'\newcommand{\manualrevision}{'+revision+'}\n')
    with tempfile.TemporaryDirectory(prefix='tpc-manual-build-') as directory:
        command=['pdflatex','-interaction=nonstopmode','-halt-on-error','-file-line-error','-output-directory='+directory,'TPC_manual.tex']
        for _ in range(2):
            result=subprocess.run(command,cwd=HERE,text=True,capture_output=True)
            if result.returncode:
                print(result.stdout[-8000:]);return result.returncode
        log=(Path(directory)/'TPC_manual.log').read_text()
        warnings=[line for line in log.splitlines() if 'Overfull' in line or 'undefined' in line]
        if warnings:print('\n'.join(warnings))
        shutil.copy2(Path(directory)/'TPC_manual.pdf',HERE/'TPC_manual.pdf')
    manifest.write_text(json.dumps(current,indent=2,sort_keys=True)+'\n')
    print(f'Built {HERE/"TPC_manual.pdf"}; revision {revision}')
    return 0

if __name__=='__main__':raise SystemExit(main())
