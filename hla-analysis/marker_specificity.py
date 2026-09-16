#!/usr/bin/env python3
"""Locate probes in both donor MHC assemblies; not a genome-wide uniqueness test."""
import csv
import json
from pathlib import Path


def fasta(path):
    name=None; pieces=[]
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                if name: yield name,''.join(pieces).upper()
                name=line[1:].split()[0];pieces=[]
            else:pieces.append(line.strip())
    if name:yield name,''.join(pieces).upper()


def main():
    rows=list(csv.DictReader(open('read_markers.tsv'),delimiter='\t'))
    complement=str.maketrans('ACGT','TGCA'); cache={}; out=[]
    for r in rows:
        donor=r['donor']
        if donor not in cache:
            cache[donor]=[record for hap in [1,2] for record in fasta(Path(f'/home/leechuck/hla/mhc_all/{donor}_{hap}.mhc.fa'))]
        result=dict(marker_id=r['marker_id'],donor=donor)
        for allele in ['assembly','alternative']:
            probe=r[allele+'_kmer']; hits=[]
            for name,seq in cache[donor]:
                for strand,k in [('+',probe),('-',probe.translate(complement)[::-1])]:
                    pos=seq.find(k)
                    while pos!=-1:
                        hits.append(dict(contig=name,start_1based=pos+1,strand=strand))
                        pos=seq.find(k,pos+1)
            result[allele+'_mhc_hits']=len(hits)
            result[allele+'_locations']=json.dumps(hits)
        out.append(result)
    with open('marker_specificity.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
    print(len(out),'probes checked')


if __name__=='__main__':main()
