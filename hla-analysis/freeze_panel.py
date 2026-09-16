#!/usr/bin/env python3
"""Explicit full-panel endpoint eligibility; unresolved identities remain flagged."""
from collections import defaultdict, Counter
import json
from sequence_catalogue import BASE, OUT, GENES, write_table
from check_results import table


def main():
    panel=table(BASE.parent/'hla-audit/2026-09-16/panel_manifest.tsv')
    groups={r['hap_id']:r for r in table(BASE.parent/'hla-pilot/results/validation_groups.tsv')}
    signatures={r['hap_id']:r for r in table(OUT/'haplotype_signatures.tsv')}
    eligibility={(r['hap_id'],r['gene']):r for r in table(OUT/'typing_eligibility.tsv')}
    predictions=table(OUT/'typing_predictions.tsv')
    evaluated={(r['hap_id'],r['gene']) for r in predictions}
    seq=defaultdict(list)
    for r in table(OUT/'sequence_catalogue.tsv'):seq[r['hap_id'],r['gene']].append(r)
    frozen=[];endpoints=[]
    for r in panel:
        h=r['hap_id'];g=groups.get(h);s=signatures.get(h)
        structural=bool(g and r['source_type']=='assembly_or_reference' and r['graph_input_contigs']=='1' and 'previously_reported_duplicated_haplotype' not in r['flags'])
        frozen.append(dict(r,retained_benchmark_donor=int(bool(g)),validation_family_group=g['group'] if g else '',pedigree_audit=g['pedigree_status'] if g else 'not_an_independent_benchmark_donor',relationship_limitation='Unknown relatives outside available pedigree are not excluded',strict_structural_screen_eligible=int(structural),current_all_classical_labels=int(s['all_classical_label_eligible']) if s else 0))
        for gene in ['HLA-'+x for x in GENES]:
            e=eligibility.get((h,gene));rows=seq[h,gene];reason=[]
            if not g:reason.append('reference_or_duplicate_source')
            if not rows:reason.append('no_gene_sequence_annotation_not_proven_absence')
            elif len(rows)>1:reason.append('multiple_gene_sequences')
            else:
                if rows[0]['cds_complete']!='1':reason.append('cds_completeness_or_abnormality')
                if not rows[0]['exact_cds_alleles']:reason.append('no_current_exact_cds_match')
            if g and e and e['eligible_current_cds']=='1' and (h,gene) not in evaluated:reason.append('missing_feature_in_common_method_comparison')
            ok=(h,gene) in evaluated
            assert not (ok and reason),(h,gene,reason)
            if not ok and not reason:reason.append('ambiguous_current_twofield_label')
            endpoints.append(dict(hap_id=h,gene=gene,retained_donor=int(bool(g)),annotated_sequence_copies=len(rows),typing_common_methods_eligible=int(ok),typing_exclusion_reasons=';'.join(reason),strict_structural_screen_eligible=int(structural),all_classical_matched_comparison_eligible=int(bool(structural and s and s['all_classical_label_eligible']=='1')),source_type=r['source_type'],panel_flags=r['flags']))
    write_table(OUT/'frozen_panel.tsv',frozen);write_table(OUT/'endpoint_eligibility.tsv',endpoints)
    summary=dict(panel_haplotypes=len(frozen),retained_haplotypes=sum(r['retained_benchmark_donor'] for r in frozen),retained_donors=len({r['donor_id'] for r in frozen if r['retained_benchmark_donor']}),strict_structural_haplotypes=sum(r['strict_structural_screen_eligible'] for r in frozen),pedigree_status_haplotypes=dict(Counter(r['pedigree_audit'] for r in frozen)),typing_eligible_by_locus=dict(Counter(r['gene'] for r in endpoints if r['typing_common_methods_eligible'])),limitations=['Donor reconciliation uses names and public pedigree, not genome-wide identity/kinship estimation','Strict structural eligibility is a conservative source/contiguity screen, not per-base assembly validation','Typing uses supplied gene sequences and fixed whole-panel bundles; downstream graph clipping does not alter these input sequences, while upstream clipping may'])
    (OUT/'frozen_panel_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
