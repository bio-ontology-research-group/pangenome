#!/usr/bin/env python3
"""Coordinate, provenance and benchmark invariants; no accuracy thresholds."""
import csv
import hashlib
import json
import gzip
import subprocess
import tempfile
import sys
from pathlib import Path
from collections import defaultdict
from Bio.Seq import Seq
from sequence_catalogue import BASE, OUT, oriented_slice, merged_intervals


def table(path):
    return list(csv.DictReader(path.open(),delimiter='\t'))


def main():
    # Independent forward-genome example, including reverse-strand intervals.
    genome='AACCGGTTAACG'
    for strand,seq in [('+',genome),('-',str(Seq(genome).reverse_complement()))]:
        observed=oriented_slice(seq,103,108,101,112,strand)
        expected=genome[2:8]
        if strand=='-': expected=str(Seq(expected).reverse_complement())
        assert observed==expected
    assert merged_intervals([(1,6),(4,9),(12,14)])==[(1,9),(12,14)]
    assert merged_intervals([(1,6),(7,9)])==[(1,9)]
    manifest=json.loads((BASE/'source/imgt/manifest.json').read_text())
    for r in manifest['files']:
        assert hashlib.sha256((BASE/'source/imgt'/r['path']).read_bytes()).hexdigest()==r['sha256']
    seq=table(OUT/'sequence_catalogue.tsv')
    assert len(seq)==6670 and all(r['annotation_matches']=='1' for r in seq)
    index={r['name']:r for r in seq}
    assert index['apr003#1#HLA-A']['cds_complete']=='0'
    assert 'partial_CDS' in index['apr003#1#HLA-A']['qc_flags']
    assert 'DQB1*02:180' in index['HG02717#1#HLA-DQB1']['exact_cds_alleles']
    candidates=table(OUT/'coding_candidates.tsv')
    for r in candidates:
        assert hashlib.sha256(r['protein_sequence'].encode()).hexdigest()==r['protein_sha256']
        assert int(r['nearest_reference_protein_edit_distance'])>0
        assert r['assessment']=='unresolved_database_absent_protein'
    if (OUT/'read_markers.tsv').exists():
        # Exercise read-quality, strand and paired-fragment behavior independently
        # of real counts; counting mates twice would exaggerate candidate support.
        marker=table(OUT/'read_markers.tsv')[0]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            with (p/'markers.tsv').open('w') as f:
                w=csv.DictWriter(f,fieldnames=marker,delimiter='\t');w.writeheader();w.writerow(marker)
            k=marker['assembly_kmer']; alt=marker['alternative_kmer']
            rc=str(Seq(k).reverse_complement())
            with gzip.open(p/'reads.fq.gz','wt') as f:
                for name,s,q in [('pair/1',k,'I'*31),('pair/2',k,'I'*31),('reverse',rc,'I'*31),('alt',alt,'I'*31),('low',alt,'!'*31)]:
                    f.write('@'+name+'\n'+s+'\n+\n'+q+'\n')
            subprocess.run([sys.executable,str(BASE/'count_read_markers.py'),str(p/'markers.tsv'),marker['donor'],str(p/'out.tsv'),str(p/'reads.fq.gz')],check=True,capture_output=True)
            result=table(p/'out.tsv')[0]
            assert result['assembly_fragments']=='2' and result['alternative_fragments']=='1'
    if (OUT/'typing_predictions.tsv').exists():
        predictions=table(OUT/'typing_predictions.tsv')
        groups={r['hap_id']:r for r in table(BASE.parent/'hla-pilot/results/validation_groups.tsv')}
        seen=defaultdict(set)
        methods=defaultdict(set)
        for r in predictions:
            assert r['hap_id'] in groups
            seen[r['gene'],r['split'],r['method']].add(r['hap_id'])
            methods[r['gene'],r['split']].add(r['method'])
            assert int(r['correct'])<=int(r['unambiguous'])
        for (gene,split),ms in methods.items():
            assert len({frozenset(seen[gene,split,m]) for m in ms})==1
        for r in table(OUT/'typing_metrics.tsv'):
            assert int(r['correct'])<=int(r['called'])<=int(r['n'])
    if (OUT/'rccx_eligibility.tsv').exists():
        rccx={r['hap_id']:r for r in table(OUT/'rccx_eligibility.tsv')}
        assert len(rccx)==754
        reference=rccx['GRCh38#0']
        assert all(reference[k]=='2' for k in ['c4_annotated_copies','cyp21_like_loci','tnx_like_loci','whr1_like_loci'])
        for h in ['ksa006#1','ksa008#1','HG02735#2','NA21144#2']:
            assert rccx[h]['strict_gapfree_rccx_eligible']=='1'
        for r in rccx.values():
            if r['strict_gapfree_rccx_eligible']=='1':
                assert r['source_type']=='assembly_or_reference'
                assert r['N_bases']=='0' and int(r['contig_edge_distance'])>=1000
        assert rccx['ksa004#2']['strict_gapfree_rccx_eligible']=='0'
    print('PASS: strand coordinates, interval union, frozen reference hashes, partial-CDS control, current allele control, candidate identity, and available benchmark invariants')


if __name__=='__main__':
    main()
