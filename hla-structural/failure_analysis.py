#!/usr/bin/env python3
"""Describe remaining errors without fitting another model to validation outcomes."""
import csv,json
from collections import Counter
from pathlib import Path
from benchmark import write

ROOT=Path(__file__).resolve().parent
def read(path):return list(csv.DictReader(open(path),delimiter='\t'))

def main():
    cat={r['structure_id']:r['structural_signature'] for r in read(ROOT/'source/catalogue.tsv')}
    rows=read(ROOT/'results/calibrated_predictions.tsv');out=[]
    def forms(pair):return sorted(tuple(g.rstrip('+-') for g in cat[s].split('|') if g.startswith('C4')) for s in pair.split(';'))
    def dosage(pair):return Counter(g for path in pair for g in path)
    def marginal(count):return dict(A=sum(v for k,v in count.items() if k.startswith('C4A')),B=sum(v for k,v in count.items() if k.startswith('C4B')),long=sum(v for k,v in count.items() if k.endswith('L')),short=sum(v for k,v in count.items() if k.endswith('S')))
    for r in rows:
        if r['pair_correct']=='1':continue
        truth=forms(r['truth_pair']);pred=forms(r['predicted_pair']) if r['predicted_pair'] else []
        if not pred:category='no_call'
        elif truth==pred:category='non_C4_prototype_signature_only'
        elif dosage(truth)==dosage(pred):category='C4_order_or_chromosome_distribution'
        elif marginal(dosage(truth))==marginal(dosage(pred)):category='C4_joint_form_phase_with_matching_marginal_dosage'
        else:category='C4_marginal_form_dosage'
        out.append(dict(donor=r['donor'],category=category,represented_in_training=r['represented_in_training'],truth_C4_pair=json.dumps(truth),predicted_C4_pair=json.dumps(pred),truth_marginal_dosage=json.dumps(marginal(dosage(truth))),predicted_marginal_dosage=json.dumps(marginal(dosage(pred))),relative_structure_margin=r['relative_structure_margin'],normalized_rmse=r['normalized_rmse']))
    write(ROOT/'results/arrangement_errors.tsv',out)
    summary=dict(errors=len(out),categories=dict(Counter(r['category'] for r in out)),unrepresented_errors=sum(r['represented_in_training']=='0' for r in out),unique_optima_are_not_confidence=True)
    (ROOT/'results/arrangement_errors.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Remaining structural-typing errors','',f"The corrected model has {len(out)} full-signature errors among 88 additional donors.",'']
    lines += [f'- {k}: {v}' for k,v in summary['categories'].items()]
    lines += ['',f"{summary['unrepresented_errors']} errors involve a true structural signature absent from the family-excluded reference catalogue. The method still selected a unique optimum: a numerical optimum must not be called a confident genotype.",'',
        'The sampled marker vocabulary supplies conserved total-C4 dosage markers but no qualifying markers for the joint C4AL/C4AS/C4BL/C4BS copy counts. Those joint forms combine A/B sequence identity with long/short status at separated positions; failure to find a universal 31-mer for a joint form does not mean the underlying component variants are unmeasurable. The independent joint-form marker baseline is therefore an identifiability diagnostic, not a fair substitute for a dedicated C4 assay.','',
        'Next assay development should explicitly include separate A/B diagnostic sites, long/short insertion junctions, and module-boundary evidence, then infer their phase. Preserve ambiguous arrangement sets where short reads do not distinguish them. Validate with independent long reads or targeted assays, including candidate structures absent from the reference. Freeze any revised assay and test a new untouched sample set; these 88 donors have now been inspected and are not a fresh validation set for further tuning.','']
    (ROOT/'FAILURE_ANALYSIS.md').write_text('\n'.join(lines));print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
