#!/usr/bin/env python3
"""Describe observed annotation/sequence differences conditional on HLA labels.

Zero annotations are not confirmed absence. Hash differences identify sequence
differences, not functional effects. C4 calls do not reconstruct full RCCX modules.
"""
from collections import defaultdict, Counter
import csv
import hashlib
import json
import sys
from sequence_catalogue import BASE, OUT, GENES, write_table
sys.path.insert(0,str(BASE.parent/'hla-pilot'))
from run_pilot import norm, table


def main():
    groups={r['hap_id']:r for r in table(BASE.parent/'hla-pilot/results/validation_groups.tsv')}
    seq=defaultdict(list)
    calls=defaultdict(list)
    for r in table(OUT/'sequence_catalogue.tsv'):
        seq[r['hap_id'],r['gene']].append(r)
    for r in table(BASE.parent/'hla-audit/2026-09-16/source/hla_calls.tsv'):
        calls[r['hap_id']].append(r)
    signatures=[]
    for h in sorted(groups):
        labels={}; reliable=True
        for suffix in GENES:
            g='HLA-'+suffix
            entries=seq[h,g]
            if not entries:
                labels[g]='NO_ANNOTATION'
                if suffix not in ['DRB3','DRB4','DRB5']: reliable=False
                continue
            values=[]
            for r in entries:
                exact={norm(a) for a in r['exact_cds_alleles'].split(';') if a}
                if r['cds_complete']!='1' or len(exact)!=1 or None in exact:
                    reliable=False
                    values.append('UNRESOLVED')
                else: values.append(next(iter(exact)))
            labels[g]='&'.join(sorted(values))
        core=['HLA-DRB1','HLA-DQA1','HLA-DQB1']
        drdq='|'.join(labels[g] for g in core)
        if any(len(seq[h,g])!=1 or 'UNRESOLVED' in labels[g] or labels[g]=='NO_ANNOTATION' for g in core): drdq=''
        full='|'.join(g+'='+labels[g] for g in sorted(labels)) if reliable else ''
        annotations=calls[h]
        counts=Counter(r['gene'] for r in annotations)
        samecontig=len({r['contig'] for r in annotations})==1
        ordered=sorted(annotations,key=lambda r:int(r['start']))
        anchor=[r for r in annotations if r['gene']=='HLA-A']
        oriented=samecontig and len(anchor)==1
        reverse=oriented and anchor[0]['strand']=='-'
        if reverse: ordered.reverse()
        order='|'.join(r['gene']+('+' if (r['strand']=='+')!=reverse else '-') for r in ordered) if oriented else ''
        c4=[r for r in ordered if r['gene'].startswith('C4')]
        c4contig=len({r['contig'] for r in c4})==1
        c4signature='|'.join(r['gene'] for r in c4) if c4contig and oriented else ''
        row=dict(hap_id=h,donor=groups[h]['donor'],cohort=groups[h]['cohort'],drdq_label=drdq,all_classical_label=full,all_classical_label_eligible=int(reliable),gene_content='|'.join(f'{g}:{n}' for g,n in sorted(counts.items())),gene_order=order,c4_annotation_order=c4signature,c4_annotation_count=len(c4),interpretation='Annotation signatures; zero calls are not confirmed absence; C4 is not complete RCCX')
        for g in ['HLA-'+x for x in GENES]:
            for kind in ['gene','cds','protein','flank']:
                row[g+'_'+kind]='|'.join(sorted(r[kind+'_sha256'] for r in seq[h,g]))
        signatures.append(row)
    write_table(OUT/'haplotype_signatures.tsv',signatures)
    contrasts=[]
    fields=['gene_content','gene_order','c4_annotation_order']+[g+'_'+k for g in ['HLA-'+x for x in GENES] for k in ['gene','cds','protein','flank']]
    for level in ['drdq_label','all_classical_label']:
        bins=defaultdict(list)
        for r in signatures:
            if r[level]: bins[r[level]].append(r)
        for label,rows in sorted(bins.items()):
            if len({r['donor'] for r in rows})<2: continue
            for field in fields:
                variants=defaultdict(list)
                for r in rows:
                    if r[field]: variants[r[field]].append(r)
                if len(variants)<2: continue
                stable=hashlib.sha256((level+'\n'+label+'\n'+field).encode()).hexdigest()[:16]
                contrasts.append(dict(contrast_id='HLA-WITHIN-'+stable,matching_level=level,matching_label=label,feature=field,distinct_signatures=len(variants),donors=len({r['donor'] for v in variants.values() for r in v}),haplotypes=sum(map(len,variants.values())),signature_members=json.dumps({k:[r['hap_id'] for r in v] for k,v in variants.items()},sort_keys=True),assessment='observed_difference_requires_quality_and_biological_assessment'))
    write_table(OUT/'within_type_contrasts.tsv',contrasts)
    summary=dict(haplotypes=len(signatures),all_classical_label_eligible=sum(r['all_classical_label_eligible'] for r in signatures),contrasts=dict(Counter(r['matching_level'] for r in contrasts)),limitations=['DRB3/4/5 no-annotation states are not validated biological absence','All available classical label matching uses numeric two-field calls, not expression suffix or full-gene identity','C4 annotations alone do not describe CYP21/TNX/STK19 module structure','Flanking sequence differences do not establish regulatory function','Multiple features from a single haplotype difference are not independent discoveries'])
    (OUT/'within_type_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
