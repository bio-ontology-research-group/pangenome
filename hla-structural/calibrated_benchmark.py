#!/usr/bin/env python3
"""Pilot-trained scalar dosage calibration; frozen before 88-donor prediction review."""
import csv,json,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
from pathlib import Path
import numpy as np
from benchmark import fit_pair,write

ROOT=Path(__file__).resolve().parent

def main():
    calibration=json.loads((ROOT/'source/dosage_calibration.json').read_text())
    rows=json.loads((ROOT/'source/matrix_rows.json').read_text())
    A=np.load(ROOT/'source/path_marker_counts.npy')
    ids=[i for i,r in enumerate(rows) if r['locus']=='RCCX'];rr=[rows[i] for i in ids];A=A[ids]
    B=np.stack([np.fromfile(ROOT/'source/background_counts'/f'{r["hap_id"].replace("#","_")}.bin',dtype='<u4') for r in rr])
    donors=(ROOT/'source/validation_donors.txt').read_text().splitlines();out=[]
    for donor in donors:
        test=[i for i,r in enumerate(rr) if r['donor']==donor];assert len(test)==2
        fam=rr[test[0]]['family'];train=[i for i,r in enumerate(rr) if r['family']!=fam]
        T=A[train];families=sorted({rr[i]['family'] for i in train});support=np.zeros(A.shape[1],dtype=np.uint16)
        for f in families:support+=(T[[j for j,i in enumerate(train) if rr[i]['family']==f]]>0).any(axis=0)
        keep=(support>=3)&(T==B[train]).all(axis=0)&(T.max(axis=0)<=8)
        anchors=keep&((T==1).mean(axis=0)>=0.98)
        raw=np.fromfile(ROOT/'source/read_counts'/f'{donor}.bin',dtype='<u4').astype(float)
        scale=float(np.median(raw[anchors])/2)
        cn=np.array([int(rr[i]['copy_number']) for i in train])
        markers=keep&((T==cn[:,None]).mean(axis=0)>=0.98)&((T[cn>0]==cn[cn>0,None]).mean(axis=0)>=0.95)
        # If a validation family overlaps a pilot family, exclude that calibration record.
        records=[r for r in calibration['pilot_records'] if r['family']!=fam]
        factor=float(np.median([r['estimated_copies']/r['assembly_copies'] for r in records]))
        estimate=float(np.median(raw[markers])/scale/factor)
        target=int(np.floor(estimate+0.5))
        scores=fit_pair(T[:,keep],raw[keep]/scale)
        scores[cn[:,None]+cn[None,:]!=target]=np.inf
        truth=sorted(rr[i]['structure_id'] for i in test);truthcn=sum(int(rr[i]['copy_number']) for i in test)
        pair=[];pairset=set();margin='';residual=''
        status='no_compatible_copy_number'
        if scale>=3 and markers.sum()>=10 and anchors.sum()>=30 and np.isfinite(scores).any():
            ij=np.unravel_index(np.argmin(scores),scores.shape);best=float(scores[ij])
            pair=sorted(rr[train[i]]['structure_id'] for i in ij)
            tied=np.argwhere(np.isclose(scores,best,rtol=1e-9,atol=1e-7))
            pairset={tuple(sorted(rr[train[int(j)]]['structure_id'] for j in ab)) for ab in tied}
            types=np.array([rr[i]['structure_id'] for i in train]);same=((types[:,None]==pair[0])&(types[None,:]==pair[1]))|((types[:,None]==pair[1])&(types[None,:]==pair[0]))
            alternate=float(np.min(np.where(same,np.inf,scores)));margin=(alternate-best)/max(best,1);residual=float(np.sqrt(max(best,0)/keep.sum()))
            status='unique_best_uncalibrated' if len(pairset)==1 else 'ambiguous'
        out.append(dict(donor=donor,family=fam,locus='RCCX',evaluation_set='additional_validation',method='pilot_calibrated_CN_constrained_paths',truth_pair=';'.join(truth),truth_cn=truthcn,represented_in_training=int(all(t in {rr[i]['structure_id'] for i in train} for t in truth)),calibration_factor=factor,calibration_donors=len(records),dosage_markers=int(markers.sum()),copy_estimate=estimate,predicted_cn=target if pairset else '',predicted_pair=';'.join(pair) if len(pairset)==1 else '',call_status=status,cn_correct=int(bool(pairset) and target==truthcn),pair_correct=int(len(pairset)==1 and pair==truth),compatible_structural_pairs=len(pairset),relative_structure_margin=margin,normalized_rmse=residual))
    write(ROOT/'results/calibrated_predictions.tsv',out)
    summary=dict(n=len(out),cn_correct=sum(r['cn_correct'] for r in out),pair_correct=sum(r['pair_correct'] for r in out),unique_calls=sum(r['call_status']=='unique_best_uncalibrated' for r in out),calibration='18 pilot donors, excluded overlapping query families',evaluation='88 additional donors; no evaluation-outcome tuning')
    (ROOT/'results/calibrated_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
