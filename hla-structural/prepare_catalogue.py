#!/usr/bin/env python3
"""Freeze coordinate requests and provisional structural endpoints before reading WGS."""
import csv, gzip, hashlib, json, re
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'hla-analysis'

def read(path):
    return list(csv.DictReader(open(path), delimiter='\t'))

def write(path, rows):
    with open(path, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)

def main():
    panel = {r['hap_id']: r for r in read(OLD/'results/frozen_panel.tsv')}
    names = {}
    for line in open(ROOT/'source/contig_names.tsv'):
        sample, hap, short, original = line.rstrip().split('\t')
        names[(sample+'#'+hap, original)] = sample+'#'+hap+'#'+short
    rows = []
    for r in read(OLD/'results/rccx_eligibility.tsv'):
        if r['strict_gapfree_rccx_eligible'] != '1': continue
        p = panel[r['hap_id']]
        rows.append(dict(locus='RCCX', hap_id=r['hap_id'], donor=p['donor_id'], family=p['validation_family_group'], cohort=p['cohort'], contig=r['span_contig'], graph_path=names[(r['hap_id'],r['span_contig'])], start0=int(r['span_start0'])-1000, end0=int(r['span_end0'])+1000, structural_signature=r['screen_signature'], copy_number=int(r['c4_annotated_copies']), truth_status='assembly_homology_and_annotation'))
    for hap, p in panel.items():
        if p['strict_structural_screen_eligible'] != '1': continue
        sample, h = hap.split('#')
        genes = []
        with gzip.open(OLD/f'source/gtf_all/{sample}_{h}.gtf.gz','rt') as f:
            for line in f:
                if line.startswith('#'): continue
                x = line.rstrip().split('\t')
                if len(x)<9 or x[2]!='gene': continue
                name = re.search(r'gene_name "([^"]+)"',x[8])
                if name and name[1] in ['HLA-DRA','HLA-DRB1','HLA-DRB3','HLA-DRB4','HLA-DRB5']:
                    genes.append((x[0],int(x[3])-1,int(x[4]),name[1],x[6]))
        counts=Counter(g[3] for g in genes)
        if counts['HLA-DRA']!=1 or counts['HLA-DRB1']!=1 or len({g[0] for g in genes})!=1: continue
        genes.sort(key=lambda g:g[1]); lo=min(g[1] for g in genes)-1000; hi=max(g[2] for g in genes)+1000
        if lo<0: continue
        signature='|'.join(g[3]+g[4] for g in genes)
        rows.append(dict(locus='DRB',hap_id=hap,donor=p['donor_id'],family=p['validation_family_group'],cohort=p['cohort'],contig=genes[0][0],graph_path=names[(hap,genes[0][0])],start0=lo,end0=hi,structural_signature=signature,copy_number=sum(counts[g] for g in ['HLA-DRB1','HLA-DRB3','HLA-DRB4','HLA-DRB5']),truth_status='provisional_annotation_absence_not_confirmed_deletion'))
    for r in rows:
        r['structure_id']=r['locus']+'-STRUCT-'+hashlib.sha256(r['structural_signature'].encode()).hexdigest()[:16]
    rows.sort(key=lambda r:(r['locus'],r['hap_id']))
    write(ROOT/'source/intervals.tsv',rows)
    print(json.dumps({'intervals':len(rows),'loci':dict(Counter(r['locus'] for r in rows))}))

if __name__=='__main__': main()
