#!/usr/bin/env python3
"""Training-independent deterministic marker sketch; fold-specific filtering is later."""
import csv, gzip, hashlib, json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
MASK=(1<<62)-1
U64=(1<<64)-1

def mix(v):
    v=((v^(v>>30))*0xbf58476d1ce4e5b9)&U64
    v=((v^(v>>27))*0x94d049bb133111eb)&U64
    return v^(v>>31)

def kmers(s):
    f=r=n=0
    for c in s:
        b={'A':0,'C':1,'G':2,'T':3}.get(c)
        if b is None: f=r=n=0;continue
        f=((f<<2)|b)&MASK;r=(r>>2)|((3-b)<<60);n+=1
        if n>=31:
            k=min(f,r)
            if mix(k)&63==0:yield k

def main():
    rows=list(csv.DictReader(open(ROOT/'source/catalogue.tsv'),delimiter='\t'))
    rows=[r for r in rows if r['eligible']=='1']
    seq={}
    with gzip.open(ROOT/'source/locus_sequences.fa.gz','rt') as f:
        for line in f:
            if line.startswith('>'):name=line.strip()[1:];seq[name]=''
            else:seq[name]+=line.strip()
    profiles=[Counter(kmers(seq[r['sequence_id']])) for r in rows]
    # Global union is only a measurement vocabulary, never the fold eligibility mask.
    keys=sorted(set().union(*(p.keys() for p in profiles)))
    lookup={k:i for i,k in enumerate(keys)}
    matrix=np.zeros((len(rows),len(keys)),dtype=np.uint16)
    for i,p in enumerate(profiles):
        for k,v in p.items():matrix[i,lookup[k]]=v
    np.save(ROOT/'source/path_marker_counts.npy',matrix)
    (ROOT/'source/markers.txt').write_text(''.join(str(k)+'\n' for k in keys))
    (ROOT/'source/matrix_rows.json').write_text(json.dumps(rows,indent=2)+'\n')
    bydonor=defaultdict(dict)
    for r in rows:bydonor[r['donor']].setdefault(r['locus'],[]).append(r)
    candidates=[]
    for donor,loci in bydonor.items():
        if not all(len(loci.get(l,[]))==2 and all(r['reads_available']=='1' for r in loci[l]) for l in ['RCCX','DRB']):continue
        group=tuple(sorted(int(r['copy_number']) for r in loci['RCCX']))
        candidates.append((group,hashlib.sha256(('20260916:'+donor).encode()).hexdigest(),donor))
    # Three donors per observed RCCX haploid-count pair, selected without read results.
    strata=defaultdict(list)
    for group,h,donor in sorted(candidates):strata[group].append(donor)
    selected=set(d for v in strata.values() for d in v[:3])
    # Positive control and prior assembly-matched structural comparison donor.
    for donor in ['HG03139','HG03195','NA21144']:
        if donor in {x[2] for x in candidates}:selected.add(donor)
    (ROOT/'source/donors.txt').write_text(''.join(d+'\n' for d in sorted(selected)))
    all_donors=sorted(x[2] for x in candidates)
    (ROOT/'source/evaluation_donors.txt').write_text(''.join(d+'\n' for d in all_donors))
    (ROOT/'source/validation_donors.txt').write_text(''.join(d+'\n' for d in all_donors if d not in selected))
    (ROOT/'source/assembly_entries.txt').write_text(''.join(h.replace('#','_')+'\n' for h in sorted({r['hap_id'] for r in rows})))
    summary=dict(k=31,selection='splitmix64(canonical_2bit_kmer)&63 == 0',markers=len(keys),paths=len(rows),eligible_diploid_donors=len(candidates),selected_donors=sorted(selected),strata={str(k):v[:3] for k,v in strata.items()},per_fold_minimum_independent_families=3)
    (ROOT/'source/design.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
