#!/usr/bin/env python3
"""Cross-artifact assertions: no dropped donors or family leakage, exact control."""
import csv,json,hashlib
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))

def main():
    cat=read(ROOT/'source/catalogue.tsv');pred=read(ROOT/'results/predictions.tsv')
    donors=set((ROOT/'source/evaluation_donors.txt').read_text().splitlines())
    assert len(cat)==1182 and all(r['eligible']=='1' and r['sequence_matches_assembly']=='1' for r in cat)
    assert len(donors)==106
    assert len(pred)==len(donors)*2*4
    assert {r['donor'] for r in pred}==donors
    assert len({(r['donor'],r['locus'],r['method']) for r in pred})==len(pred)
    by={(r['donor'],r['locus'],r['method']):r for r in pred}
    for donor in donors:
        for locus in ['RCCX','DRB']:
            a=by[donor,locus,'path_multiplicity'];b=by[donor,locus,'flat_sequence_catalogue']
            assert {k:v for k,v in a.items() if k!='method'}=={k:v for k,v in b.items() if k!='method'}
            reference=[r for r in cat if r['locus']==locus and r['family']!=a['family']]
            assert len(reference)==int(a['training_haplotypes'])
            assert all(r['donor']!=donor for r in reference)
            assert a['truth_cn']==str(sum(int(r['copy_number']) for r in cat if r['locus']==locus and r['donor']==donor))
    markers=len((ROOT/'source/markers.txt').read_text().splitlines())
    for donor in donors:assert (ROOT/'source/read_counts'/f'{donor}.bin').stat().st_size==4*markers
    assert len(read(ROOT/'results/depth_baseline.tsv'))==2*len(donors)
    calibration=json.loads((ROOT/'source/dosage_calibration.json').read_text())
    assert hashlib.sha256((ROOT/'calibrated_benchmark.py').read_bytes()).hexdigest()==calibration['script_sha256']
    cp=read(ROOT/'results/calibrated_predictions.tsv');cs=json.loads((ROOT/'results/calibrated_summary.json').read_text())
    pilot={r['donor'] for r in calibration['pilot_records']}
    assert len(cp)==88 and len(pilot)==18 and not pilot&{r['donor'] for r in cp}
    assert {r['donor'] for r in cp}==donors-pilot
    for r in cp:
        baseline=by[r['donor'],'RCCX','path_multiplicity']
        assert r['truth_pair']==baseline['truth_pair'] and r['truth_cn']==baseline['truth_cn']
        allowed=[x for x in calibration['pilot_records'] if x['family']!=r['family']]
        assert int(r['calibration_donors'])==len(allowed)
        assert abs(float(r['calibration_factor'])-np.median([x['estimated_copies']/x['assembly_copies'] for x in allowed]))<1e-10
    assert cs['cn_correct']==sum(int(r['cn_correct']) for r in cp)
    assert cs['pair_correct']==sum(int(r['pair_correct']) for r in cp)
    freeze=json.loads((ROOT/'source/pre_evaluation_freeze.json').read_text())
    assert hashlib.sha256((ROOT/'source/benchmark_frozen_before_predictions.py.txt').read_bytes()).hexdigest()==freeze['diagnostic_addendum']['benchmark_sha256']
    summary=dict(status='passed',catalogue_intervals=len(cat),donors=len(donors),primary_predictions=len(pred),calibrated_validation_predictions=len(cp),markers=markers,checks=['graph_sequence_concordance','no_dropped_donors','query_family_exclusion','flat_path_equivalence','copy_number_truth_consistency','complete_read_measurements','depth_baseline_denominators','pilot_validation_disjoint','calibration_family_exclusion','calibration_code_frozen_before_validation','original_model_snapshot_hash'])
    (ROOT/'results/checks.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
