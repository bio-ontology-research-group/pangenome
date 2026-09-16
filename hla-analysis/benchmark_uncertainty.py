#!/usr/bin/env python3
"""Family-cluster bootstrap; shared draws for paired method comparisons."""
from collections import defaultdict
import numpy as np
from sequence_catalogue import BASE, OUT, write_table
from check_results import table


def main():
    groups={r['hap_id']:r['group'] for r in table(BASE.parent/'hla-pilot/results/validation_groups.tsv')}
    bins=defaultdict(lambda:defaultdict(list))
    for r in table(OUT/'typing_predictions.tsv'):
        bins[r['gene'],r['split']][r['method']].append(r)
    rng=np.random.default_rng(20260916)
    results=[]
    for (gene,split),methods in sorted(bins.items()):
        cluster=sorted({groups[r['hap_id']] for rows in methods.values() for r in rows})
        ci={g:i for i,g in enumerate(cluster)}
        counts=np.zeros(len(cluster))
        reference=methods['gene_sequence']
        for r in reference: counts[ci[groups[r['hap_id']]]]+=1
        draws=rng.multinomial(len(cluster),np.full(len(cluster),1/len(cluster)),size=2000)
        denominator=draws@counts
        baseline=np.zeros(len(cluster))
        for r in reference: baseline[ci[groups[r['hap_id']]]]+=int(r['correct'])
        for method,rows in sorted(methods.items()):
            assert {r['hap_id'] for r in rows}=={r['hap_id'] for r in reference}
            correct=np.zeros(len(cluster))
            for r in rows: correct[ci[groups[r['hap_id']]]]+=int(r['correct'])
            sampled=(draws@correct)/denominator
            difference=(draws@(correct-baseline))/denominator
            lo,hi=np.quantile(sampled,[.025,.975])
            dlo,dhi=np.quantile(difference,[.025,.975])
            results.append(dict(gene=gene,split=split,method=method,n=len(rows),family_clusters=len(cluster),accuracy=correct.sum()/counts.sum(),accuracy_ci_low=lo,accuracy_ci_high=hi,difference_vs_sampled_gene_sequence=(correct-baseline).sum()/counts.sum(),difference_ci_low=dlo,difference_ci_high=dhi,replicates=2000,seed=20260916,scope='Conditional on fixed predictions and panel; does not refit graph or capture unknown relatedness'))
    write_table(OUT/'typing_uncertainty.tsv',results)


if __name__=='__main__':main()
