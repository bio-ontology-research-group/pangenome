#!/usr/bin/env python3
"""Check actual assembled RCCX intervals for gaps and contig-edge proximity."""
from pathlib import Path
import csv
import hashlib
from marker_specificity import fasta


def main():
    rows=list(csv.DictReader(open('rccx_screen.tsv'),delimiter='\t'));out=[]
    for r in rows:
        if not r['span_contig']:continue
        sample,hap=r['hap_id'].split('#')
        records=dict(fasta(Path(f'/home/leechuck/hla/mhc_all/{sample}_{hap}.mhc.fa')))
        seq=records[r['span_contig']];lo=int(r['span_start0']);hi=int(r['span_end0'])
        region=seq[lo:hi]
        assert len(region)==hi-lo
        out.append(dict(hap_id=r['hap_id'],contig=r['span_contig'],start0=lo,end0=hi,span_bp=len(region),N_bases=region.count('N'),other_ambiguous_bases=sum(c not in 'ACGTN' for c in region),contig_edge_distance=min(lo,len(seq)-hi),span_sha256=hashlib.sha256(region.encode()).hexdigest()))
    with open('rccx_span_qc.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
    print(len(out),'assembled RCCX spans checked')


if __name__=='__main__':main()
