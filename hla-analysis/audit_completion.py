#!/usr/bin/env python3
"""Check deliverable coverage and internal consistency, not biological truth."""
from collections import Counter
import json
from sequence_catalogue import BASE, OUT
from check_results import table


def main():
    evidence={}
    panel=table(OUT/'frozen_panel.tsv');endpoints=table(OUT/'endpoint_eligibility.tsv')
    assert len(panel)==754 and len(endpoints)==754*11
    assert len({r['hap_id'] for r in panel})==754
    assert sum(r['retained_benchmark_donor']=='1' for r in panel)==742
    remote=json.loads((BASE/'source/remote_source_manifest.json').read_text())['files']
    assert len([r for r in remote if r['path'].endswith('.mhc.fa') and r['status']=='present'])==754
    evidence['1_frozen_panel']={'haplotypes':len(panel),'locus_endpoint_rows':len(endpoints),'remote_source_records':len(remote),'unknown_kinship':'explicitly unresolved'}
    pred=table(OUT/'typing_predictions.tsv');metrics=table(OUT/'typing_metrics.tsv')
    required={'gene_sequence','gene_bundles_fixed_panel','classII_bundles_fixed_panel','gene_content','coding_sequence_all_31mers'}
    assert required<={r['method'] for r in pred}
    assert len({r['gene'] for r in pred})==11
    assert {r['split'] for r in pred}=={'leave_family_out','leave_cohort_out'}
    assert {'unseen','rare_1_to_5_donors','more_than_5_donors'}<={r['rarity'] for r in metrics}
    assert len({r['gene'] for r in table(OUT/'experimental_genotypes.tsv')})==5
    assert len(table(OUT/'experimental_discordances.tsv'))==14
    evidence['2_benchmark']={'prediction_rows':len(pred),'methods':sorted({r['method'] for r in pred}),'loci':11,'truth_loci':5,'scope':'fixed panel, numeric two-field, eligible annotated genes'}
    contrast=table(OUT/'within_type_evidence.tsv')
    assert {r['matching_level'] for r in contrast}=={'drdq_label','all_classical_label'}
    assert any(r['feature'].endswith('_flank') for r in contrast)
    assert len(table(OUT/'rccx_eligibility.tsv'))==754
    evidence['3_beyond_labels']={'contrasts':len(contrast),'rccx_entries':754,'functional_module_and_regulatory_identity':'not established'}
    candidates=table(OUT/'candidate_evidence_ranking.tsv')
    assert len(candidates)==47 and len({r['candidate_id'] for r in candidates})==47
    assert len(table(OUT/'published_candidate_sequence_check.tsv'))==50
    assert len(table(OUT/'flank_validation_summary.tsv'))==4
    assert len(table(OUT/'rccx_depth_copy_check.tsv'))==2
    assert any(r['assessment']=='unresolved_read_assembly_discordance' for r in table(OUT/'flank_validation_summary.tsv'))
    evidence['4_candidate_assessment']={'protein_candidates':len(candidates),'categories':dict(Counter(r['category'] for r in candidates)),'validated_new_alleles':0,'validation_scope':'local SNPs, cross-source assemblies, known DRB5 control; structural depth exploratory'}
    for name in ['FINAL_REPORT.md','REPRODUCIBILITY.md','UKB_VALIDATION_PLAN.md','COMPLETION_AUDIT.md']:
        assert (BASE/name).is_file() and (BASE/name).stat().st_size>100
    log=json.loads((OUT/'reproduction/run.json').read_text())
    from run_analysis import STEPS
    assert [r['script'] for r in log]==STEPS
    assert all(r['exit_code']==0 for r in log)
    evidence['5_report_workflow']={'reproduced_stages':len(log),'ukb_execution':'prospective, not run','final_report':'FINAL_REPORT.md'}
    (OUT/'completion_evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(evidence,indent=2))


if __name__=='__main__':main()
