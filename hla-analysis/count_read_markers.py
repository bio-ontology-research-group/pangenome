#!/usr/bin/env python3
"""Count exact genomic 31-mers in paired/singleton primary reads, BQ>=20.

One count per read, with read-name fragment deduplication across mates/files.
No genomic uniqueness is assumed. Duplicate-marked reads are excluded upstream.
"""
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path


def main():
    markers,donor,out,*fastqs=sys.argv[1:]
    rows=[r for r in csv.DictReader(open(markers),delimiter='\t') if r['donor']==donor]
    probes=defaultdict(list)
    complement=str.maketrans('ACGT','TGCA')
    for r in rows:
        for allele in ['assembly','alternative']:
            k=r[allele+'_kmer']; probes[k].append((r['marker_id'],allele))
            rc=k.translate(complement)[::-1]
            if rc!=k: probes[rc].append((r['marker_id'],allele))
    support=defaultdict(set); n=0
    for path in fastqs:
        with gzip.open(path,'rt') as f:
            while True:
                name=f.readline().strip()
                if not name:break
                seq=f.readline().strip(); plus=f.readline(); qual=f.readline().strip()
                assert name.startswith('@') and plus.startswith('+') and len(seq)==len(qual)
                n+=1; fragment=name[1:].split()[0].removesuffix('/1').removesuffix('/2')
                for i in range(len(seq)-30):
                    hit=probes.get(seq[i:i+31])
                    if hit and min(map(ord,qual[i:i+31]))>=53:
                        for key in hit:support[key].add(fragment)
    result=[]
    for r in rows:
        a=support[r['marker_id'],'assembly']; b=support[r['marker_id'],'alternative']
        result.append(dict(r,assembly_fragments=len(a),alternative_fragments=len(b),fragments_supporting_both=len(a&b),primary_reads_scanned=n))
    with open(out,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
    print(donor,n,'reads',len(result),'markers',flush=True)


if __name__=='__main__':main()
