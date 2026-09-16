#!/usr/bin/env python3
"""Fixed-panel assembly typing benchmark. No training labels from query families.

Bundles were constructed using the full frozen panel: this is transductive, not
unseen-graph validation. Current exact CDS matches supply benchmark labels;
assembly labels are not experimental truth. Numeric two-field endpoint only.
"""
from pathlib import Path
from collections import defaultdict, Counter
import csv
import sys
import zlib
import numpy as np
from Bio import SeqIO
from sequence_catalogue import GENES, write_table

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parent / 'hla-pilot'))
from run_pilot import similarity, norm, table

OUT = BASE / 'results'


def features(path, sampling=32):
    out = {}
    complement = str.maketrans('ACGT', 'TGCA')
    for rec in SeqIO.parse(path, 'fasta'):
        s = str(rec.seq).upper()
        selected = set()
        for i in range(len(s)-30):
            k = s[i:i+31]
            if set(k) - set('ACGT'):
                continue
            k = min(k, k.translate(complement)[::-1])
            if sampling == 1 or zlib.crc32(k.encode()) % sampling == 0:
                selected.add(k)
        out[rec.id] = selected
    return out


def bundles(path):
    out = {}
    for line in path.read_text().splitlines():
        rid, vector = line.split('\t')[:2]
        out[rid] = {i for i, v in enumerate(vector.split(',')) if int(v)}
    return out


def by_haplotype(records, gene):
    out = {}
    for rid, f in records.items():
        bits = rid.split('#')
        if bits[-1] == gene:
            out['#'.join(bits[:2])] = f
    return out


def main():
    groups = {r['hap_id']: r for r in table(BASE.parent/'hla-pilot/results/validation_groups.tsv')}
    ids = sorted(groups)
    catalogue = table(OUT/'sequence_catalogue.tsv')
    counts = Counter((r['hap_id'], r['gene']) for r in catalogue)
    labels, eligibility = {}, []
    for r in catalogue:
        h, g = r['hap_id'], r['gene']
        if h not in groups:
            continue
        alleles = {norm(x) for x in r['exact_cds_alleles'].split(';') if x}
        ok = counts[h,g] == 1 and r['cds_complete'] == '1' and len(alleles) == 1 and None not in alleles
        labels[h,g] = next(iter(alleles)) if ok else None
        eligibility.append(dict(hap_id=h,gene=g,eligible_current_cds=int(ok),current_twofield=labels[h,g] or '',old_twofield=norm(r['old_consensus']) or '',qc_flags=r['qc_flags'],assessment=r['assessment']))
    write_table(OUT/'typing_eligibility.tsv',eligibility)
    sequence = {k: features(OUT/p) for k,p in [('gene_sequence','gene_sequences.fa'),('coding_sequence','coding_sequences.fa'),('flanks_only','flanking_sequences.fa')]}
    sequence['coding_sequence_all_31mers'] = features(OUT/'coding_sequences.fa',sampling=1)
    regional = {'#'.join(k.split('#')[:2]):v for k,v in bundles(BASE.parent/'hla-pilot/source/classII.ord').items()}
    calls = defaultdict(Counter)
    for r in table(BASE.parent/'hla-audit/2026-09-16/source/hla_calls.tsv'):
        calls[r['hap_id']][r['gene']] += 1
    content = {h:{(g,i) for g,n in calls[h].items() for i in range(1,n+1)} for h in ids}
    predictions = []
    for suffix in GENES:
        g = 'HLA-'+suffix
        methods = {k:by_haplotype(v,g) for k,v in sequence.items()}
        methods['gene_bundles_fixed_panel'] = by_haplotype(bundles(BASE/'source'/f'{g}.ord'),g)
        methods['gene_content'] = content
        if suffix.startswith(('DR','DQ','DP')):
            methods['classII_bundles_fixed_panel'] = regional
        # Same query and reference eligibility for every method at this locus.
        eligible = [h for h in ids if labels.get((h,g)) and all(f.get(h) for f in methods.values())]
        group_arr = np.array([groups[h]['group'] for h in eligible])
        cohort_arr = np.array([groups[h]['cohort'] for h in eligible])
        truth = np.array([labels[h,g] for h in eligible])
        for method, feature in methods.items():
            mat = similarity(feature,eligible)
            for split in ['leave_family_out','leave_cohort_out']:
                for i,h in enumerate(eligible):
                    mask = group_arr != group_arr[i]
                    if split == 'leave_cohort_out':
                        mask &= cohort_arr != cohort_arr[i]
                    train = np.flatnonzero(mask)
                    assert all(groups[eligible[j]]['group'] != groups[h]['group'] for j in train)
                    distribution = Counter(truth[train])
                    best = float(mat[i,train].max()) if len(train) else 0
                    ties = train[np.isclose(mat[i,train],best,rtol=0,atol=1e-7)] if best>0 else []
                    predicted = sorted(set(truth[ties])) if len(ties) else []
                    top = max(distribution.values(),default=0)
                    majority = sorted(k for k,v in distribution.items() if v==top)
                    donor_support = len({groups[eligible[j]]['donor'] for j in train if truth[j]==truth[i]})
                    predictions.append(dict(hap_id=h,donor=groups[h]['donor'],cohort=groups[h]['cohort'],gene=g,method=method,split=split,annotation=truth[i],prediction='|'.join(predicted),unambiguous=int(len(predicted)==1),correct=int(predicted==[truth[i]]),contains_annotation=int(truth[i] in predicted),nearest_similarity=best,tied_neighbours=len(ties),reference_haplotypes=len(train),reference_allele_donors=donor_support,majority_correct=int(majority==[truth[i]])))
        print(g, 'common eligible',len(eligible),flush=True)
    write_table(OUT/'typing_predictions.tsv',predictions)
    bins = defaultdict(list)
    for r in predictions:
        rarity = 'unseen' if r['reference_allele_donors']==0 else 'rare_1_to_5_donors' if r['reference_allele_donors']<=5 else 'more_than_5_donors'
        for c in ['ALL',r['cohort']]:
            bins[r['gene'],r['method'],r['split'],c,'ALL'].append(r)
        bins[r['gene'],r['method'],r['split'],'ALL',rarity].append(r)
    metrics=[]
    for (g,m,s,c,rarity),rows in sorted(bins.items()):
        n=len(rows); called=sum(r['unambiguous'] for r in rows); correct=sum(r['correct'] for r in rows)
        metrics.append(dict(gene=g,method=m,split=s,cohort=c,rarity=rarity,n=n,called=called,correct=correct,accuracy_all=correct/n,call_rate=called/n,accuracy_called=correct/called if called else '',majority_accuracy=sum(r['majority_correct'] for r in rows)/n))
    write_table(OUT/'typing_metrics.tsv',metrics)
    # Experimental five-locus truth, independent of assembly-assigned labels.
    index={(r['hap_id'],r['gene'],r['method'],r['split']):r for r in predictions}
    experiments=[]
    for t in csv.DictReader((BASE.parent/'hla-pilot/source/1000G_HLA.txt').open(),delimiter=' '):
        for suffix in ['A','B','C','DRB1','DQB1']:
            g='HLA-'+suffix
            allowed=[{norm(suffix+'*'+a) for a in t[k].split('/')} for k in [suffix,suffix+'.1']]
            if any(None in x for x in allowed): continue
            if any(not labels.get((t['id']+'#'+hap,g)) for hap in ['1','2']): continue
            configs=sorted({(r['method'],r['split']) for r in predictions if r['gene']==g})
            for m,s in configs+[('current_exact_cds','none')]:
                ps=[]
                for hap in ['1','2']:
                    h=t['id']+'#'+hap
                    r=index.get((h,g,m,s))
                    ps.append(labels[h,g] if m=='current_exact_cds' else r['prediction'] if r and r['unambiguous'] else None)
                called=all(ps)
                match=bool(called and ((ps[0] in allowed[0] and ps[1] in allowed[1]) or (ps[1] in allowed[0] and ps[0] in allowed[1])))
                experiments.append(dict(donor=t['id'],gene=g,method=m,split=s,truth=t[suffix]+';'+t[suffix+'.1'],prediction=';'.join(p or 'NO_CALL' for p in ps),called=int(called),match=int(match)))
    write_table(OUT/'experimental_genotypes.tsv',experiments)


if __name__=='__main__':
    main()
