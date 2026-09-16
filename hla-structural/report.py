#!/usr/bin/env python3
"""Catalogue description, baseline comparisons, uncertainties and honest scope."""
import csv,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from benchmark import dosage_labels,write

ROOT=Path(__file__).resolve().parent
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def integer(r,k):return int(r[k])

def main():
    cat=read(ROOT/'source/catalogue.tsv');pred=read(ROOT/'results/predictions.tsv');dose=read(ROOT/'results/dosage_baseline.tsv')
    summaries=[]
    for locus in ['DRB','RCCX']:
        grouped=defaultdict(list)
        for r in cat:
            if r['locus']==locus:grouped[r['structure_id']].append(r)
        for sid,rs in sorted(grouped.items()):
            summaries.append(dict(locus=locus,structure_id=sid,signature=rs[0]['structural_signature'],copy_number=rs[0]['copy_number'],haplotypes=len(rs),donors=len({r['donor'] for r in rs}),families=len({r['family'] for r in rs}),sequence_haplotypes=len({r['sequence_sha256'] for r in rs}),cohorts=';'.join(sorted({r['cohort'] for r in rs})),truth_status=rs[0]['truth_status']))
    write(ROOT/'results/structural_types.tsv',summaries)
    labels={r['hap_id']:r for r in read(ROOT.parent/'hla-analysis/results/haplotype_signatures.tsv') if r['all_classical_label_eligible']=='1'}
    all_labels={r['hap_id']:r for r in read(ROOT.parent/'hla-analysis/results/haplotype_signatures.tsv')}
    grouped_haps=defaultdict(dict)
    for r in cat:grouped_haps[r['hap_id']][r['locus']]=r
    extended=[]
    for hap,loci in sorted(grouped_haps.items()):
        r=loci['RCCX'];d=loci['DRB'];label=all_labels.get(hap,{})
        extended.append(dict(hap_id=hap,donor=r['donor'],cohort=r['cohort'],classical_HLA=label.get('all_classical_label',''),complete_classical_label=label.get('all_classical_label_eligible','0'),DRB_structure_id=d['structure_id'],DRB_gene_order=d['structural_signature'],RCCX_structure_id=r['structure_id'],RCCX_signature=r['structural_signature'],RCCX_component_count=r['copy_number'],C4_order='|'.join(g for g in r['structural_signature'].split('|') if g.startswith('C4')),DRB_sequence_sha256=d['sequence_sha256'],RCCX_sequence_sha256=r['sequence_sha256'],evidence='source_assembly_and_exact_graph_walk; not independent structural validation',phase_scope='assembly_haplotype; read-inferred structure pairs are not phased to classical HLA'))
    write(ROOT/'results/extended_HLA_panel.tsv',extended)
    contrasts=[]
    for locus in ['DRB','RCCX']:
        groups=defaultdict(list)
        for r in cat:
            if r['locus']==locus and r['hap_id'] in labels:groups[labels[r['hap_id']]['all_classical_label']].append(r)
        for label,rs in groups.items():
            if len({r['structure_id'] for r in rs})>1:
                contrasts.append(dict(locus=locus,classical_label=label,haplotypes=';'.join(r['hap_id'] for r in rs),structure_ids=';'.join(sorted({r['structure_id'] for r in rs})),copy_numbers=';'.join(sorted({r['copy_number'] for r in rs}))))
    if contrasts:write(ROOT/'results/same_hla_structural_contrasts.tsv',contrasts)
    # Reference-mapped read-depth baseline, not tuned against test labels.
    positions=defaultdict(set)
    for line in open(ROOT.parent/'hla-analysis/source/rccx_depth_targets.bed'):
        chrom,start,end,name=line.rstrip().split('\t');positions[name].update(range(int(start)+1,int(end)+1))
    depth=[]
    for p in sorted((ROOT/'source/read_counts').glob('*.depth.tsv.gz')):
        donor=p.name.split('.')[0]
        with gzip.open(p,'rt') as f:d={int(x[1]):int(x[2]) for x in (l.split('\t') for l in f)}
        assert all(ps<=d.keys() for ps in positions.values()),donor
        means={k:np.mean([d[i] for i in ps]) for k,ps in positions.items()}
        c4=(means['C4A']+means['C4B'])/2
        truth=sum(int(r['copy_number']) for r in cat if r['donor']==donor and r['locus']=='RCCX')
        for control in ['WHR1_unique_region','TNXB_distal_region']:
            value=4*c4/means[control] if means[control]>0 else None
            call=int(np.floor(value+0.5)) if value is not None else None
            depth.append(dict(donor=donor,locus='RCCX',control=control,truth=truth,estimate=value if value is not None else '',rounded_call=call if call is not None else '',correct=int(call==truth),control_depth=means[control]))
    write(ROOT/'results/depth_baseline.tsv',depth)
    # Fairness diagnostic: give the reference-depth baseline the same scalar
    # calibration recipe, using pilot labels only, and score the other 88 donors.
    pilot=set((ROOT/'source/donors.txt').read_text().splitlines())
    calibrated_depth=[]
    for control in ['WHR1_unique_region','TNXB_distal_region']:
        training=[r for r in depth if r['control']==control and r['donor'] in pilot]
        factor=float(np.median([r['estimate']/r['truth'] for r in training]))
        for r in depth:
            if r['control']!=control or r['donor'] in pilot:continue
            value=r['estimate']/factor;call=int(np.floor(value+.5))
            calibrated_depth.append(dict(**r,calibration_factor=factor,calibrated_estimate=value,calibrated_call=call,calibrated_correct=int(call==r['truth'])))
    write(ROOT/'results/calibrated_depth_baseline.tsv',calibrated_depth)
    # Dosage alone resolves only configurations uniquely compatible with its estimates.
    independent=[]
    for p in pred:
        if p['method']!='path_multiplicity':continue
        donor,locus=p['donor'],p['locus'];family=p['family']
        ds={r['feature']:r for r in dose if r['donor']==donor and r['locus']==locus and r['rounded_call']!=''}
        types={r['structure_id']:r for r in cat if r['locus']==locus and r['family']!=family}
        keys=sorted(types);allowed=[]
        for a,ka in enumerate(keys):
            for kb in keys[a:]:
                rs=[types[ka],types[kb]]
                features={k:sum(dosage_labels(r)[k] for r in rs) for k in dosage_labels(rs[0])}
                features['total']=sum(int(r['copy_number']) for r in rs)
                if ds and all(features[k]==int(v['rounded_call']) for k,v in ds.items()):allowed.append((ka,kb))
        call=';'.join(allowed[0]) if len(allowed)==1 else ''
        independent.append(dict(donor=donor,locus=locus,features=len(ds),compatible_pairs=len(allowed),predicted_pair=call,truth_pair=p['truth_pair'],pair_correct=int(call==p['truth_pair']),status='unique' if len(allowed)==1 else 'ambiguous' if allowed else 'inconsistent_or_missing_dosage'))
    write(ROOT/'results/independent_dosage_pairs.tsv',independent)
    summed=[]
    for donor in sorted({r['donor'] for r in dose}):
        ds=[r for r in dose if r['donor']==donor and r['locus']=='DRB' and r['feature']!='total']
        call=sum(int(r['rounded_call']) for r in ds) if len(ds)==4 and all(r['rounded_call']!='' for r in ds) else None
        truth=sum(int(r['truth']) for r in ds)
        summed.append(dict(donor=donor,locus='DRB',truth=truth,call=call if call is not None else '',correct=int(call==truth)))
    write(ROOT/'results/summed_dosage_baseline.tsv',summed)
    # Post-pilot fairness diagnostic for component dosage, using the same fixed
    # correction as the constrained path method. No validation-specific fitting.
    model=json.loads((ROOT/'source/dosage_calibration.json').read_text())
    calibrated_independent=[]
    for p in pred:
        if p['method']!='path_multiplicity' or p['locus']!='RCCX' or p['evaluation_set']!='additional_validation':continue
        records=[r for r in model['pilot_records'] if r['family']!=p['family']]
        factor=float(np.median([r['estimated_copies']/r['assembly_copies'] for r in records]))
        ds={r['feature']:int(np.floor(float(r['estimate'])/factor+.5)) for r in dose if r['donor']==p['donor'] and r['locus']=='RCCX' and r['estimate']!=''}
        types={r['structure_id']:r for r in cat if r['locus']=='RCCX' and r['family']!=p['family']};keys=sorted(types);allowed=[]
        for a,ka in enumerate(keys):
            for kb in keys[a:]:
                rs=[types[ka],types[kb]];features={k:sum(dosage_labels(r)[k] for r in rs) for k in dosage_labels(rs[0])};features['total']=sum(int(r['copy_number']) for r in rs)
                if ds and all(features[k]==v for k,v in ds.items()):allowed.append((ka,kb))
        call=';'.join(allowed[0]) if len(allowed)==1 else ''
        calibrated_independent.append(dict(donor=p['donor'],features=len(ds),compatible_pairs=len(allowed),predicted_pair=call,truth_pair=p['truth_pair'],pair_correct=int(call==p['truth_pair']),status='unique' if len(allowed)==1 else 'ambiguous' if allowed else 'inconsistent_or_missing_dosage'))
    write(ROOT/'results/calibrated_independent_pairs.tsv',calibrated_independent)
    # Uncertainty conditional on this enriched pilot: resample donor family clusters.
    rng=np.random.default_rng(20260916);ci=[]
    for locus in ['RCCX','DRB']:
        for method in ['path_multiplicity','binary_paths']:
            ps=[r for r in pred if r['locus']==locus and r['method']==method];fams=sorted({r['family'] for r in ps})
            for endpoint in ['cn_correct','pair_correct']:
                groups={fam:[int(r[endpoint]) for r in ps if r['family']==fam] for fam in fams}
                vals=[np.mean([v for fam in rng.choice(fams,len(fams),replace=True) for v in groups[fam]]) for _ in range(2000)]
                low,high=np.quantile(vals,[.025,.975]);ci.append(dict(locus=locus,method=method,endpoint=endpoint,n=len(ps),correct=sum(int(r[endpoint]) for r in ps),low=float(low),high=float(high)))
    write(ROOT/'results/uncertainty.tsv',ci)
    splits=[]
    for locus in ['RCCX','DRB']:
        for subset in ['pilot','additional_validation']:
            for method in ['path_multiplicity','binary_paths','assembly_profile_oracle']:
                ps=[r for r in pred if r['locus']==locus and r['method']==method and r['evaluation_set']==subset]
                splits.append(dict(locus=locus,evaluation_set=subset,method=method,n=len(ps),cn_correct=sum(int(r['cn_correct']) for r in ps),pair_correct=sum(int(r['pair_correct']) for r in ps),unrepresented_truth=sum(r['represented_in_training']=='0' for r in ps)))
    write(ROOT/'results/split_summary.tsv',splits)
    calibrated=read(ROOT/'results/calibrated_predictions.tsv') if (ROOT/'results/calibrated_predictions.tsv').exists() else []
    type_sequences={r['structure_id']:r['structural_signature'] for r in cat}
    c4=[]
    for p in pred+calibrated:
        if p['locus']!='RCCX':continue
        def c4pair(text):
            return ';'.join(sorted('|'.join(g for g in type_sequences[sid].split('|') if g.startswith('C4')) for sid in text.split(';'))) if text else ''
        truth=c4pair(p['truth_pair']);call=c4pair(p['predicted_pair'])
        c4.append(dict(donor=p['donor'],evaluation_set=p['evaluation_set'],method=p['method'],truth_C4_order_pair=truth,predicted_C4_order_pair=call,correct=int(bool(call) and call==truth)))
    write(ROOT/'results/c4_order_pairs.tsv',c4)
    if calibrated:
        intervals=[]
        for endpoint in ['cn_correct','pair_correct']:
            n=len(calibrated);k=sum(int(r[endpoint]) for r in calibrated);prop=k/n;z=1.959963984540054
            centre=(prop+z*z/(2*n))/(1+z*z/n);half=z*np.sqrt(prop*(1-prop)/n+z*z/(4*n*n))/(1+z*z/n)
            intervals.append(dict(method='pilot_calibrated_CN_constrained_paths',endpoint=endpoint,n=n,correct=k,wilson_low=float(centre-half),wilson_high=float(centre+half),assumption='88 distinct recorded families; approximate independent-donor binomial interval'))
        write(ROOT/'results/calibrated_uncertainty.tsv',intervals)
    lines=['# Structural MHC typing pilot','',
        'This is a technical pilot using held-out families and MHC-recruited public short reads. It does not establish clinical accuracy, complete SV reconstruction, a new structural class or disease relevance.','',
        '**Main result:** coarse DRB gene-content typing matches all 106 assembly labels. On 88 additional RCCX donors, pilot-calibrated total copy number matches 88/88, the full annotation-signature pair 73/88, and ordered C4 forms 74/88. Calibrated ordinary reference depth also recovers copy number in 88/88: no graph-specific copy-number advantage is established. Arrangement remains the harder, incompletely resolved endpoint.','',
        '## Catalogue','',
        'All 1,182 intervals (591 haplotypes at each locus) match the corresponding full-graph walk and assembly sequence exactly and contain no ambiguous bases. The frozen panel yields four coarse DRB gene-content types and 38 RCCX annotation signatures. RCCX signatures include nearest-prototype CYP21/TNX assignments, not validated functional alleles. The DRB interval omits a full pseudogene classification; annotation absence remains provisional.','',
        f'{len(contrasts)} classical-HLA label groups contain differing structural signatures in the strict subset; see `results/same_hla_structural_contrasts.tsv`. These are panel observations, not population frequencies or proof of disease effects.','',
        '`results/extended_HLA_panel.tsv` adds the structural layer to classical HLA labels for each of the 591 eligible assembled haplotypes: DRB content/order, RCCX component count, ordered C4 forms, stable structure IDs and sequence hashes. These are experimental descriptors, not new official HLA allele names. Assembly phase links these descriptors within a source haplotype; the short-read assay returns unordered structural pairs and does not phase them to classical HLA alleles.','',
        '## Prospective benchmark','',
        'Eighteen donors were selected deterministically before read results, enriching uncommon RCCX copy-number pairs and including known controls. Before computing predictions, the unchanged method was locked for the remaining 88 eligible donors. All 106 donors contribute two loci; per-set results are below. Candidate paths, marker family support and off-locus specificity exclude the query family. Unknown relatives outside available pedigree metadata may remain. The original graph alignment used all donors; this is a fixed-graph experiment, not a rebuilt-graph benchmark. Inference uses only retained training sequence profiles, not target-only nodes or target structural labels.','',
        'Exact canonical 31-mers are sampled deterministically at 1/64. Markers must occur in at least three training families, have no extra hits outside the locus in any retained training MHC assembly, and have at most eight occurrences per training path. This is MHC specificity, not whole-genome uniqueness. Reads exclude duplicate-marked, secondary, supplementary, unmapped and QC-failed records; all marker bases require Q20, and matching mates count once per fragment.','',
        'The prototype exhaustively fits unordered pairs of observed local paths by squared error against normalized marker dosage. Depth is estimated from training-defined near-universal single-copy markers. Scores are not calibrated genotype probabilities. Sequence-wide matching can impute arrangement through linked variation without a read spanning its defining junction; it is not independent physical phasing.','',
        '| Locus | Method | Copy number correct | Structural pair correct | Unique best pair |','|---|---|---:|---:|---:|']
    for r in read(ROOT/'results/summary.tsv'):
        if r['method']=='flat_sequence_catalogue':continue
        lines.append(f"| {r['locus']} | {r['method']} | {r['cn_correct']}/{r['n']} | {r['pair_correct']}/{r['n']} | {r['unique_calls']}/{r['n']} |")
    lines+=['','The binary ablation discards within-path marker multiplicity. The flat-sequence control is analytically identical to multiplicity-aware path inference because verified path sequences and markers are identical; its copied result rows record that equivalence, not an independently tested third algorithm. No graph-storage advantage is demonstrated. The assembly-profile oracle supplies exact held-out assembly marker counts to the same reference search: it diagnoses panel/representation limits under perfect measurements and is not a deployable read-typing method.','',
        '| Locus | Evaluation set | Method | Copy number correct | Structural pair correct |','|---|---|---|---:|---:|']
    for r in splits:lines.append(f"| {r['locus']} | {r['evaluation_set']} | {r['method']} | {r['cn_correct']}/{r['n']} | {r['pair_correct']}/{r['n']} |")
    if calibrated:
        model=json.loads((ROOT/'source/dosage_calibration.json').read_text())
        nc=sum(int(r['cn_correct']) for r in calibrated);npair=sum(int(r['pair_correct']) for r in calibrated)
        lines+=['','## Pilot-trained correction, independent additional evaluation','',
            f"The pilot showed systematic overestimation of C4 dosage. A single factor ({model['factor']:.4f}) was fitted from the 18 pilot donors as the median marker-dosage estimate divided by assembly copy number. The secondary method divides dosage by this factor and restricts path pairs to the rounded corrected copy number. It was frozen before any predictions on the remaining 88 donors were inspected; any overlapping query family is excluded from the calibration records. The original benchmark above is unchanged.", '',
            f'On the 88 additional donors, the corrected method recovers total RCCX copy number in **{nc}/88** and the complete annotation-signature pair in **{npair}/88**. These are assembly-label agreements, not independent breakpoint or functional validation. The copy-number result is supplied by calibrated marker dosage; graph path inference cannot claim credit for that improvement.', '',
            'C4-only ordered forms are separately tabulated in `results/c4_order_pairs.tsv`, avoiding reliance on uncertain CYP21/TNX prototype distinctions. The structural pair still derives from sequence-panel inference and is not direct long-range phasing.','']
        cn_dist=Counter(r['truth_cn'] for r in calibrated)
        lines += [f"Validation copy-number distribution: {dict(sorted(cn_dist.items()))}. The two- and six-copy donors occurred in the pilot, so this validation does not establish performance at those extremes. Each query family is excluded separately; other evaluation donors' assembled reference paths can remain in its reference panel.", '',
            'As a fairness diagnostic, the reference-depth baselines receive the same pilot-only median-ratio calibration recipe. This comparison was added after observing the secondary method result; no validation outcomes were used to fit or choose its factors.','']
        for control in ['WHR1_unique_region','TNXB_distal_region']:
            ds=[r for r in calibrated_depth if r['control']==control]
            lines.append(f"- Calibrated reference depth / {control}: {sum(r['calibrated_correct'] for r in ds)}/{len(ds)} total-copy calls correct.")
        lines += ['- The independent RCCX dosage sketch has no qualifying markers for the combined A/B-and-long/short forms. Total copy number alone does not uniquely resolve the signature pair. This is an identifiability diagnostic, not a competitive dedicated C4-typing baseline.','',
            f"For the corrected method, Wilson 95% intervals are {intervals[0]['wilson_low']:.1%}–{intervals[0]['wilson_high']:.1%} for total copy number and {intervals[1]['wilson_low']:.1%}–{intervals[1]['wilson_high']:.1%} for the full signature pair. These intervals do not include annotation, recruitment or future-cohort uncertainty.",'']
    lines+=['',
        '## Simpler baselines','', '| Locus | Baseline | Endpoint | Correct / attempted |','|---|---|---|---:|']
    for locus in ['RCCX','DRB']:
        ds=[r for r in dose if r['locus']==locus and r['feature']=='total']
        if locus=='RCCX':lines.append(f"| {locus} | Independent total-copy markers | Total copies | {sum(int(r['correct']) for r in ds)}/{len(ds)} |")
        else:lines.append(f"| DRB | Sum of independent paralog dosages | Total copies | {sum(r['correct'] for r in summed)}/{len(summed)} |")
        ps=[r for r in independent if r['locus']==locus]
        if locus=='DRB':lines.append(f"| {locus} | Independent component dosage | Structural pair (no forced ambiguity resolution) | {sum(r['pair_correct'] for r in ps)}/{len(ps)} |")
    for control in ['WHR1_unique_region','TNXB_distal_region']:
        ds=[r for r in depth if r['control']==control]
        lines.append(f"| RCCX | Reference depth / {control} | Total copies | {sum(r['correct'] for r in ds)}/{len(ds)} |")
    lines+=['','The depth baselines use original GRCh38 mappings, C4 exon means, two separately reported control regions and nearest-integer rounding, without tuning against these samples. Reference and paralog mapping bias remain. Independent dosage estimates require at least ten eligible component markers; missing/ambiguous calls count as failures in all-attempted accuracy.','',
        '## Limits and next decision','',
        'See `results/predictions.tsv` for represented/unrepresented types, margins and residuals, `results/uncertainty.tsv` for family-bootstrap intervals, and `results/dosage_baseline.tsv` for per-component estimates. A unique least-squares optimum is not a confident clinical call. Bootstrap intervals describe this enriched pilot only.','',
        'Neither an existing general-purpose SV caller nor graph read alignment is benchmarked here. Novel structures absent from training can be forced onto a known pair; residuals are recorded, but there is no validated novelty detector. Full chromosome phase, exact breakpoints, functional CYP21/TNX identity, gene conversion and UKB phenotype association are not established.','',
        'Before deployment, validate discordant structures with long reads or independent assemblies, extend beyond recruited reads, calibrate no-call thresholds on separate development donors, and test the locked assay on additional donors and ancestries. Any claim of graph-specific improvement requires a graph-alignment or topology-aware method against this flat-sequence baseline.','',
        'The [failure analysis](FAILURE_ANALYSIS.md) separates C4 marginal dosage errors, joint-form phase errors, chromosome-distribution errors and uncertain prototype distinctions. The 88 evaluation donors have now been inspected; further assay tuning requires a fresh validation set.','',
        'Methodological precedent: [PanGenie, Ebler et al. 2022](https://www.nature.com/articles/s41588-022-01043-w); direct MHC/RCCX precedent: [Logsdon et al. 2025](https://www.nature.com/articles/s41586-025-09140-6). This custom pilot does not implement PanGenie.','']
    (ROOT/'REPORT.md').write_text('\n'.join(lines))
    print('Report written; same-HLA contrasting groups',len(contrasts))

if __name__=='__main__':main()
