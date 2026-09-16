#!/usr/bin/env python3
"""Hash the remote assembly and read files actually used by this analysis."""
from pathlib import Path
import hashlib
import json
import datetime


def main():
    work=Path('/home/leechuck/hla/codex-analysis/rccx');paths=[]
    for entry in (work/'rccx_panel_ids.txt').read_text().splitlines():
        paths.append(Path('/home/leechuck/hla/mhc_all')/(entry+'.mhc.fa'))
    for donor in ['HG01252','HG01358','HG01943','HG02976','NA18608','NA18620','NA18943','NA19007','NA19159','HG03130','HG03195','HG03139','NA21144']:
        paths.append(Path('/home/asianhla/data/upload/1000G_MHC/cram')/(donor+'.mhc.cram'))
    paths.append(Path('/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa'))
    rows=[]
    for p in paths:
        if not p.exists():rows.append(dict(path=str(p),status='absent'));continue
        digest=hashlib.sha256()
        with p.open('rb') as f:
            for chunk in iter(lambda:f.read(1024*1024),b''):digest.update(chunk)
        rows.append(dict(path=str(p),status='present',bytes=p.stat().st_size,sha256=digest.hexdigest()))
    (work/'remote_source_manifest.json').write_text(json.dumps(dict(recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=rows),indent=2)+'\n')
    print(len(rows),'remote sources recorded')


if __name__=='__main__':main()
