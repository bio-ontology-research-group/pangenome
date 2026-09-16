#!/usr/bin/env python3
"""Untuned tie-breaking diagnostic, exploratory after the baseline comparison.

Only intersect a tied all-CDS-kmer label set with the bundle label set. Preserve
the original call if the intersection is empty; never replace a single-label call.
This rule's assessment is exploratory, not a preregistered independent validation.
"""
from collections import defaultdict
from sequence_catalogue import OUT, write_table
from check_results import table


def labels(r):
    return set(r['prediction'].split('|')) if r['prediction'] else set()


def main():
    index=defaultdict(dict)
    for r in table(OUT/'typing_predictions.tsv'):
        index[r['hap_id'],r['gene'],r['split']][r['method']]=r
    details=[]
    for (h,g,s),methods in sorted(index.items()):
        base=methods['coding_sequence_all_31mers']; seq=labels(base)
        for method in ['gene_bundles_fixed_panel','classII_bundles_fixed_panel']:
            if method not in methods: continue
            b=labels(methods[method]); combined=seq
            if len(seq)>1 and seq&b:combined=seq&b
            truth=base['annotation']; correct=int(combined=={truth})
            details.append(dict(hap_id=h,gene=g,split=s,bundle_method=method,sequence_prediction=base['prediction'],bundle_prediction=methods[method]['prediction'],combined_prediction='|'.join(sorted(combined)),sequence_correct=int(base['correct']),combined_correct=correct,sequence_called=int(base['unambiguous']),combined_called=int(len(combined)==1),rescued=int(correct and not int(base['correct'])),lost=int(int(base['correct']) and not correct),false_call_from_ambiguity=int(len(seq)>1 and len(combined)==1 and not correct)))
    write_table(OUT/'bundle_increment_predictions.tsv',details)
    bins=defaultdict(list)
    for r in details:bins[r['gene'],r['split'],r['bundle_method']].append(r)
    summary=[]
    for (g,s,m),rows in sorted(bins.items()):
        summary.append(dict(gene=g,split=s,bundle_method=m,n=len(rows),**{k:sum(r[k] for r in rows) for k in ['sequence_correct','combined_correct','sequence_called','combined_called','rescued','lost','false_call_from_ambiguity']}))
    write_table(OUT/'bundle_increment_summary.tsv',summary)
    for split in ['leave_family_out','leave_cohort_out']:
        for method in ['gene_bundles_fixed_panel','classII_bundles_fixed_panel']:
            rows=[r for r in details if r['split']==split and r['bundle_method']==method]
            print(split,method,'rescued',sum(r['rescued'] for r in rows),'false calls',sum(r['false_call_from_ambiguity'] for r in rows))


if __name__=='__main__':main()
