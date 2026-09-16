#!/usr/bin/env python3
"""Compute-node provenance capture without exporting source reads."""
import hashlib,json,subprocess,datetime
from pathlib import Path

def main():
    paths=[Path('/home/asianhla/data/upload/HLA/mhc_graph/MHC.full.gfa.gz'),Path('/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa')]
    paths += [Path('/home/leechuck/hla/mhc_all')/(x+'.mhc.fa') for x in Path('assembly_entries.txt').read_text().splitlines()]
    paths += [Path('/home/asianhla/data/upload/1000G_MHC/cram')/(x+'.mhc.cram') for x in Path('evaluation_donors.txt').read_text().splitlines()]
    records=[]
    for p in paths:
        before=p.stat();h=hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
        after=p.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
        records.append(dict(path=str(p),bytes=after.st_size,mtime_ns=after.st_mtime_ns,sha256=h.hexdigest()))
    versions={cmd:subprocess.check_output([cmd,'--version'],text=True).splitlines()[0] for cmd in ['/home/leechuck/hla/mm/envs/hla/bin/samtools','/home/leechuck/hla/mm/envs/hla/bin/python3']}
    Path('input_manifest.json').write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),versions=versions,files=records),indent=2)+'\n')
    print('COMPLETE',len(records),'input hashes',flush=True)

if __name__=='__main__':main()
