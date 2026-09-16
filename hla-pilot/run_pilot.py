#!/usr/bin/env python3
"""Fixed-feature nearest-neighbour pilot; annotations are not experimental truth.

Class-II bundles are transductive (constructed on the complete panel).
Canonical 31-mer sampling is sequence-local and never fits on test labels.
Excludes target donor/family from neighbours and preserves tied-label ambiguity.
"""
from pathlib import Path
from collections import defaultdict, Counter
import csv
import json
import re
import zlib
import numpy as np
from scipy.sparse import csr_matrix
from Bio import SeqIO

BASE = Path(__file__).resolve().parent
AUDIT = BASE.parent / 'hla-audit/2026-09-16'
OUT = BASE / 'results'
GENES = ['HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1']


def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def save(name, rows):
    if not rows:
        (OUT / name).write_text('')
        return
    with (OUT / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def norm(value):
    # Preserve expression suffixes in source tables; this endpoint compares
    # numeric two-field identity only, not expression status or full sequence.
    match = re.match(r'(?:HLA-)?([A-Z0-9]+)\*(\d+):(\d+)', value)
    if match and int(match[3]) > 0:
        return f'HLA-{match[1]}*{match[2]}:{match[3]}'
    return None


def similarity(features, ids):
    vocabulary = {}
    rr, cc = [], []
    for i, h in enumerate(ids):
        for feature in features.get(h, set()):
            rr.append(i)
            cc.append(vocabulary.setdefault(feature, len(vocabulary)))
    # The vocabulary merely indexes a fixed feature universe; no weighting,
    # filtering, normalization or tuning is fitted on the held-out donors.
    X = csr_matrix((np.ones(len(rr), dtype=np.float32), (rr, cc)), shape=(len(ids), len(vocabulary)))
    overlap = (X @ X.T).toarray()
    size = np.asarray(X.sum(axis=1)).ravel()
    union = size[:, None] + size[None, :] - overlap
    return np.divide(overlap, union, out=np.zeros_like(overlap), where=union > 0)


def kmers(path):
    features = {}
    complement = str.maketrans('ACGT', 'TGCA')
    for rec in SeqIO.parse(path, 'fasta'):
        if rec.id.split('#')[-1].endswith('_2'):
            continue
        h = '#'.join(rec.id.split('#')[:2])
        seq = str(rec.seq).upper()
        rc = seq.translate(complement)[::-1]
        n = len(seq)
        selected = set()
        for i in range(n - 30):
            k = min(seq[i:i+31], rc[n-i-31:n-i])
            if any(b not in 'ACGT' for b in k):
                continue
            if zlib.crc32(k.encode()) % 32 == 0:
                selected.add(k)  # retain actual kmer, so no hash collisions
        if h in features:
            raise ValueError(f'Multiple gene sequences for {h}: {path}')
        features[h] = selected
    return features


def main():
    OUT.mkdir(exist_ok=True)
    panel = {r['hap_id']: r for r in table(AUDIT / 'panel_manifest.tsv')}
    calls = defaultdict(list)
    all_calls = defaultdict(list)
    for r in table(AUDIT / 'source/hla_calls.tsv'):
        calls[r['hap_id'], r['gene']].append(r)
        all_calls[r['hap_id']].append(r)
    # Select one assembly source per known donor; keep both haplotypes together.
    ids = sorted(h for h, r in panel.items() if r['cohort'] != 'REF' and '.JaSaPaGe' not in r['sample'])
    parent = {}
    def root(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = root(parent[x])
        return parent[x]
    def join(a, b):
        a, b = root(a), root(b)
        parent[max(a, b)] = min(a, b)
    pedigree = {}
    for r in csv.DictReader((BASE / 'source/1000G.ped').open(), delimiter=' '):
        pedigree[r['SampleID']] = r
        join(r['SampleID'], 'family:' + r['FamilyID'])
        for k in ['FatherID', 'MotherID']:
            if r[k] not in {'0', '-9', ''}:
                join(r['SampleID'], r[k])
    # Conservative name-based grouping, explicitly not verified pedigree.
    for donor in ['KPanRef-076', 'KPanRef-076-M', 'KPanRef-076-P']:
        join(donor, 'provisional:KPanRef-076')
    for donor in ['apr-f', 'apr-m', 'apr-s']:
        join(donor, 'provisional:APR_trio')
    group = {h: root(panel[h]['donor_id']) for h in ids}
    save('validation_groups.tsv', [dict(hap_id=h, donor=panel[h]['donor_id'], group=group[h],
        pedigree_status='1000G_record' if panel[h]['donor_id'] in pedigree else 'unverified_outside_1000G',
        cohort=panel[h]['cohort']) for h in ids])
    label = {}
    for h in ids:
        for gene in GENES:
            r = calls[h, gene]
            label[h, gene] = norm(r[0]['consensus']) if len(r) == 1 else None
    bundle = {}
    for line in (BASE / 'source/classII.ord').read_text().splitlines():
        name, vec = line.split('\t')[:2]
        bundle['#'.join(name.split('#')[:2])] = {i for i, v in enumerate(vec.split(',')) if int(v)}
    matrices = {'classII_bundle_presence_transductive': similarity(bundle, ids)}
    features_by_gene = {}
    for gene in GENES:
        print('Extracting fixed 31-mer features:', gene, flush=True)
        features_by_gene[gene] = kmers(BASE / 'source' / (gene + '.fa'))
    predictions = []
    for gene in GENES:
        methods = dict(matrices)
        methods['gene_flanks_31mer_fixed'] = similarity(features_by_gene[gene], ids)
        for method, matrix in methods.items():
            feature = bundle if method.startswith('classII') else features_by_gene[gene]
            for mode in ['leave_family_out', 'leave_cohort_out']:
                for i, h in enumerate(ids):
                    truth = label[h, gene]
                    if truth is None or not feature.get(h):
                        continue
                    train = [j for j, v in enumerate(ids) if group[v] != group[h] and label[v, gene]
                             and feature.get(v) and (mode != 'leave_cohort_out' or panel[v]['cohort'] != panel[h]['cohort'])]
                    assert all(group[ids[j]] != group[h] for j in train)
                    if not train:
                        continue
                    scores = matrix[i, train]
                    best = float(scores.max())
                    ties = np.asarray(train)[np.isclose(scores, best, rtol=0, atol=1e-7)].tolist()
                    labels = sorted({label[ids[j], gene] for j in ties}) if best > 0 else []
                    majority = Counter(label[ids[j], gene] for j in train)
                    top = max(majority.values())
                    majority_labels = {k for k, v in majority.items() if v == top}
                    predictions.append(dict(hap_id=h, donor=panel[h]['donor_id'], cohort=panel[h]['cohort'],
                        gene=gene, method=method, split=mode, annotation=truth,
                        prediction='|'.join(labels), unambiguous=int(len(labels)==1),
                        correct=int(labels == [truth]), contains_annotation=int(truth in labels),
                        nearest_similarity=best, tied_neighbours=len(ties), reference_haplotypes=len(train),
                        allele_seen=int(truth in majority), majority_correct=int(majority_labels == {truth})))
    save('predictions.tsv', predictions)
    summaries = []
    for method, gene, split, cohort in sorted({(r['method'], r['gene'], r['split'], c) for r in predictions for c in ['ALL', r['cohort']]}):
        rows = [r for r in predictions if (r['method'], r['gene'], r['split']) == (method,gene,split) and (cohort=='ALL' or r['cohort']==cohort)]
        called = sum(r['unambiguous'] for r in rows)
        correct = sum(r['correct'] for r in rows)
        summaries.append(dict(method=method,gene=gene,split=split,cohort=cohort,n=len(rows),called=called,
            correct=correct,accuracy_all=correct/len(rows),call_rate=called/len(rows),
            accuracy_called=correct/called if called else 0,
            unseen_alleles=sum(not r['allele_seen'] for r in rows),
            majority_accuracy=sum(r['majority_correct'] for r in rows)/len(rows)))
    save('metrics.tsv', summaries)
    # Independent diploid comparison: retain ambiguity lists and both pairings.
    truth_rows = list(csv.DictReader((BASE/'source/1000G_HLA.txt').open(), delimiter=' '))
    pred_index = {(r['donor'],r['hap_id'].split('#')[1],r['gene'],r['method'],r['split']):r for r in predictions}
    experiment = []
    for t in truth_rows:
        for gene in ['HLA-DRB1','HLA-DQB1']:
            g = gene.removeprefix('HLA-')
            allowed = [{norm(g+'*'+a) for a in t[k].split('/')} for k in [g,g+'.1']]
            if any(None in s for s in allowed):
                continue
            for method in ['classII_bundle_presence_transductive','gene_flanks_31mer_fixed','assembly_annotation']:
                for split in (['leave_family_out','leave_cohort_out'] if method!='assembly_annotation' else ['none']):
                    ps=[]
                    for hap in ['1','2']:
                        if method=='assembly_annotation':
                            ps.append(label.get((t['id']+'#'+hap,gene)))
                        else:
                            r=pred_index.get((t['id'],hap,gene,method,split))
                            ps.append(r['prediction'] if r and r['unambiguous'] else None)
                    # Only evaluate donors with two numeric assembly labels.
                    if any(label.get((t['id']+'#'+hap,gene)) is None for hap in ['1','2']):
                        continue
                    called=all(ps)
                    matched=bool(called and ((ps[0] in allowed[0] and ps[1] in allowed[1]) or (ps[1] in allowed[0] and ps[0] in allowed[1])))
                    experiment.append(dict(donor=t['id'],gene=gene,method=method,split=split,
                        truth=t[g]+';'+t[g+'.1'],prediction=';'.join(x or 'NO_CALL' for x in ps),called=int(called),match=int(matched)))
    save('experimental_genotypes.tsv',experiment)
    # Structural candidate screen: signed gene order on a single contig,
    # oriented by DRA; do not interpret absent annotations as confirmed losses.
    structures=[]
    for h in ids:
        gene_rows=[r for r in all_calls[h] if r['gene'].startswith(('HLA-DR','HLA-DQ'))]
        anchor=calls[h,'HLA-DRA']
        if len(anchor)!=1 or not gene_rows or len({r['contig'] for r in gene_rows})!=1:
            continue
        reverse=anchor[0]['strand']=='-'
        ordered=sorted(gene_rows,key=lambda r:int(r['start']),reverse=reverse)
        signature='>'.join(r['gene']+('+' if (r['strand']=='-')==reverse else '-') for r in ordered)
        if any(not label[h,g] for g in GENES):
            continue
        structures.append(dict(hap_id=h,donor=panel[h]['donor_id'],cohort=panel[h]['cohort'],
            classical_combination='|'.join(label[h,g] for g in GENES),signature=signature,
            upstream_clipped=int(panel[h]['source_type']=='clipped_graph_paths'),
            span_bp=max(int(r['end']) for r in ordered)-min(int(r['start']) for r in ordered)+1))
    save('structural_haplotypes.tsv',structures)
    candidates=[]
    for combo in sorted({r['classical_combination'] for r in structures}):
        rows=[r for r in structures if r['classical_combination']==combo]
        signatures={r['signature'] for r in rows}
        if len(signatures)<2:
            continue
        for sig in sorted(signatures):
            rr=[r for r in rows if r['signature']==sig]
            candidates.append(dict(classical_combination=combo,signature=sig,
                donors=len({r['donor'] for r in rr}),original_assembly_donors=len({r['donor'] for r in rr if not r['upstream_clipped']}),
                haplotypes=';'.join(r['hap_id'] for r in rr),status='annotation_candidate_requires_sequence_validation'))
    save('within_type_structure_candidates.tsv',candidates)
    result=dict(haplotypes=len(ids),donors=len({panel[h]['donor_id'] for h in ids}),
        family_groups=len(set(group.values())),prediction_rows=len(predictions),
        structural_haplotypes=len(structures),structural_combinations_with_multiple_signatures=len({r['classical_combination'] for r in candidates}),
        bundle_status='transductive: bundles fitted on whole panel',
        kmer_status='fixed canonical k=31; crc32 selection modulo 32; exact retained strings; gene plus existing flanks',
        limitations=['Not a read-based typer','Pedigree incomplete outside 1000G','Old IPD annotation labels','Unresolved numeric two-field labels excluded per locus; noncoding novelty retained','No parameter tuning','Structural candidates unvalidated'])
    (OUT/'run_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
