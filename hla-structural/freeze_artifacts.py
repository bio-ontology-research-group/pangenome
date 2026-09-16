#!/usr/bin/env python3
"""Freeze local code, results, measurement archives and dependency fingerprints."""
import datetime,hashlib,json,platform
from pathlib import Path
import numpy,matplotlib

ROOT=Path(__file__).resolve().parent
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    dest=ROOT/'results/artifact_manifest.json'
    files=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p!=dest]
    dependencies=[ROOT.parent/'hla-analysis/results/artifact_manifest.json',ROOT.parent/'hla-analysis/results/frozen_panel.tsv',ROOT.parent/'hla-analysis/results/rccx_eligibility.tsv',ROOT.parent/'hla-analysis/results/haplotype_signatures.tsv',ROOT.parent/'hla-analysis/source/rccx_depth_targets.bed']
    obj=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),software=dict(python=platform.python_version(),numpy=numpy.__version__,matplotlib=matplotlib.__version__),dependencies=[dict(path=str(p),sha256=sha(p)) for p in dependencies],files=files)
    dest.write_text(json.dumps(obj,indent=2)+'\n');print('Frozen',len(files),'artifacts')

if __name__=='__main__':main()
