#!/usr/bin/env python3
"""Focused scientific consistency checks on the completed pilot artifacts."""
from pathlib import Path
from collections import defaultdict
import csv
import importlib.util
import tempfile
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

P = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('pilot', P/'run_pilot.py')
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)

assert pilot.norm('HLA-DRB1*15:new') is None
assert pilot.norm('HLA-DRB1*15:00:new') is None
assert pilot.norm('HLA-DRB1*15:03:01:new') == 'HLA-DRB1*15:03'
assert pilot.norm('DRB1*15:03N') == 'HLA-DRB1*15:03'

# Independent reverse-complement control: same fixed features for either strand.
import random
rng=random.Random(734)
sequence=''.join(rng.choice('ACGT') for _ in range(2000))
with tempfile.TemporaryDirectory() as d:
    path=Path(d)/'test.fa'
    SeqIO.write([SeqRecord(Seq(sequence),id='X#1#GENE'),
                 SeqRecord(Seq(sequence).reverse_complement(),id='Y#1#GENE')],path,'fasta')
    features=pilot.kmers(path)
    assert features['X#1'] and features['X#1']==features['Y#1']
    sim=pilot.similarity(features,['X#1','Y#1'])
    assert sim[0,1]==1

def read(name):
    with (P/'results'/name).open() as f:
        return list(csv.DictReader(f,delimiter='\t'))

groups=read('validation_groups.tsv')
donors=defaultdict(set)
for r in groups:donors[r['donor']].add(r['group'])
assert all(len(v)==1 for v in donors.values())
assert not any('.JaSaPaGe' in r['hap_id'] for r in groups)
pred=read('predictions.tsv')
assert len(pred)==len({(r['hap_id'],r['gene'],r['method'],r['split']) for r in pred})
for r in pred:
    assert 0<=float(r['nearest_similarity'])<=1
    assert int(r['correct'])<=int(r['unambiguous'])
    if not int(r['allele_seen']):assert not int(r['correct'])
for m in read('metrics.tsv'):
    rr=[r for r in pred if all(r[k]==m[k] for k in ['gene','method','split']) and (m['cohort']=='ALL' or r['cohort']==m['cohort'])]
    assert len(rr)==int(m['n'])
    assert sum(int(r['correct']) for r in rr)==int(m['correct'])
for r in read('experimental_genotypes.tsv'):
    assert int(r['match'])<=int(r['called'])
    if r['called']=='1':
        truth=[set(':'.join(a.split(':')[:2]) for a in s.split('/')) for s in r['truth'].split(';')]
        calls=[s.split('*')[1] for s in r['prediction'].split(';')]
        match=(calls[0] in truth[0] and calls[1] in truth[1]) or (calls[1] in truth[0] and calls[0] in truth[1])
        assert int(match)==int(r['match'])

source=P.parent/'hla-audit/2026-09-16/source'
order=list(dict.fromkeys(line.split('\t')[0] for line in (source/'classII.bed').open() if not line.startswith('#')))
tree=Phylo.read(source/'classII.nwk','newick');tree.ladderize()
expected=['#'.join(order[int(t.name)].split('#')[:2]) for t in tree.get_terminals()]
actual=read('corrected_dendrogram/data/dendrogram_tip_order.tsv')
assert [r['hap_id'] for r in actual]==expected
print('PASS: allele parsing, strand invariance, duplicate exclusion, donor groups, metric denominators, experimental pair matching, and all 753 dendrogram rows.')
