#!/usr/bin/env python3
"""Family-excluded diploid path inference, with independent dosage and binary controls.

Read marker counting is a sketch-based prototype, not PanGenie, a calibrated
likelihood, de novo SV calling or a clinical genotyper. Exact sequence and graph
path representations intentionally share an inference engine.
"""
import argparse, csv, json, os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent

def write(path, rows):
    with open(path,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def dosage_labels(row):
    s=row['structural_signature']
    if row['locus']=='RCCX':
        return {g:s.count(g+'+')+s.count(g+'-') for g in ['C4AL','C4AS','C4BL','C4BS']}
    return {g:s.count(g+'+')+s.count(g+'-') for g in ['HLA-DRB1','HLA-DRB3','HLA-DRB4','HLA-DRB5']}

def fit_pair(X,y):
    """All unordered path pairs, squared loss, including homozygotes."""
    X=X.astype(np.float64)
    gram=X@X.T; z=np.diag(gram)-2*(X@y)
    scores=z[:,None]+z[None,:]+2*gram+float(y@y)
    scores[np.tril_indices(len(X),-1)]=np.inf
    return scores

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--donors',type=Path,default=ROOT/'source/evaluation_donors.txt')
    parser.add_argument('--output',type=Path,default=ROOT/'results')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    rows=json.loads((ROOT/'source/matrix_rows.json').read_text())
    A=np.load(ROOT/'source/path_marker_counts.npy')
    donors=args.donors.read_text().splitlines()
    pilot=set((ROOT/'source/donors.txt').read_text().splitlines())
    # Entire assembled MHC background removes markers with off-locus hits.
    background={h:np.fromfile(ROOT/'source/background_counts'/f'{h.replace("#","_")}.bin',dtype='<u4') for h in {r['hap_id'] for r in rows}}
    assert all(len(a)==A.shape[1] for a in background.values())
    B=np.stack([background[r['hap_id']] for r in rows])
    predictions=[]; dosage=[]; folds=[]
    for donor in donors:
        raw=np.fromfile(ROOT/'source/read_counts'/f'{donor}.bin',dtype='<u4').astype(float)
        assert len(raw)==A.shape[1]
        for locus in ['RCCX','DRB']:
            test=[i for i,r in enumerate(rows) if r['donor']==donor and r['locus']==locus]
            assert len(test)==2
            family=rows[test[0]]['family']
            train=[i for i,r in enumerate(rows) if r['locus']==locus and r['family']!=family]
            assert all(rows[i]['donor']!=donor and rows[i]['family']!=family for i in train)
            T=A[train]; back=B[train]
            families=sorted({rows[i]['family'] for i in train})
            support=np.zeros(A.shape[1],dtype=np.uint16)
            for fam in families:
                support+=(T[[j for j,i in enumerate(train) if rows[i]['family']==fam]]>0).any(axis=0)
            # These conditions use training rows only; no target assembly features.
            specific=(T==back).all(axis=0)
            keep=(support>=3)&specific&(T.max(axis=0)<=8)
            anchors=keep & ((T==1).mean(axis=0)>=0.98)
            scale=float(np.median(raw[anchors])/2) if anchors.sum() else 0
            truth=sorted(rows[i]['structure_id'] for i in test)
            truth_cn=sum(int(rows[i]['copy_number']) for i in test)
            train_types={rows[i]['structure_id'] for i in train}
            represented=all(t in train_types for t in truth)
            common=dict(donor=donor,family=family,locus=locus,evaluation_set='pilot' if donor in pilot else 'additional_validation',truth_pair=';'.join(truth),truth_cn=truth_cn,represented_in_training=int(represented),training_haplotypes=len(train),markers=int(keep.sum()),anchors=int(anchors.sum()),haploid_marker_depth=scale)
            folds.append(common.copy())
            if scale<3 or keep.sum()<100 or anchors.sum()<30:
                for method in ['path_multiplicity','binary_paths','flat_sequence_catalogue','assembly_profile_oracle']:
                    predictions.append(dict(**common,method=method,call_status='insufficient_evidence',predicted_pair='',predicted_cn='',pair_correct=0,cn_correct=0,compatible_structural_pairs=0,relative_structure_margin='',normalized_rmse=''))
                continue
            y=raw[keep]/scale
            X=T[:,keep].astype(float)
            # No empirical tuning on the evaluation read sets.
            for method in ['path_multiplicity','binary_paths','assembly_profile_oracle']:
                Z=(X>0).astype(float) if method=='binary_paths' else X
                observed=A[test][:,keep].sum(axis=0).astype(float) if method=='assembly_profile_oracle' else y
                scores=fit_pair(Z,observed)
                ij=np.unravel_index(np.argmin(scores),scores.shape);best=float(scores[ij])
                pair=sorted(rows[train[j]]['structure_id'] for j in ij)
                cn=sum(int(rows[train[j]]['copy_number']) for j in ij)
                # Determine distinct structural pairs in all score ties; never break
                # a structural ambiguity by arbitrary donor ordering.
                tied=np.argwhere(np.isclose(scores,best,rtol=1e-9,atol=1e-7))
                pairset={tuple(sorted(rows[train[int(j)]]['structure_id'] for j in ab)) for ab in tied}
                cnset={sum(int(rows[train[int(j)]]['copy_number']) for j in ab) for ab in tied}
                typeids=np.array([rows[i]['structure_id'] for i in train])
                same=((typeids[:,None]==pair[0])&(typeids[None,:]==pair[1]))|((typeids[:,None]==pair[1])&(typeids[None,:]==pair[0]))
                alternate=float(np.min(np.where(same,np.inf,scores)))
                p=dict(**common,method=method,call_status='unique_best_uncalibrated' if len(pairset)==1 else 'ambiguous',predicted_pair=';'.join(pair) if len(pairset)==1 else '',predicted_cn=cn if len(cnset)==1 else '',pair_correct=int(len(pairset)==1 and pair==truth),cn_correct=int(len(cnset)==1 and cn==truth_cn),compatible_structural_pairs=len(pairset),relative_structure_margin=(alternate-best)/max(best,1),normalized_rmse=float(np.sqrt(max(best,0)/len(y))))
                predictions.append(p)
                if method=='path_multiplicity':
                    predictions.append(dict(p,method='flat_sequence_catalogue'))
            # Copy-number baseline: marker counts matching each component dosage
            # in >=98% of training haplotypes, with >=5 positive families.
            targets={'total':np.array([int(rows[i]['copy_number']) for i in train])}
            for key in dosage_labels(rows[test[0]]):
                targets[key]=np.array([dosage_labels(rows[i])[key] for i in train])
            for feature,target in targets.items():
                marker=keep & ((T==target[:,None]).mean(axis=0)>=0.98)
                positive_fams={rows[train[j]]['family'] for j in range(len(train)) if target[j]>0}
                if (target>0).any():marker &= ((T[target>0]==target[target>0,None]).mean(axis=0)>=0.95)
                truthdose=truth_cn if feature=='total' else sum(dosage_labels(rows[i])[feature] for i in test)
                n=int(marker.sum())
                estimate=float(np.median(raw[marker])/scale) if n>=10 and len(positive_fams)>=5 else None
                called=int(np.floor(estimate+0.5)) if estimate is not None else None
                dosage.append(dict(donor=donor,locus=locus,feature=feature,truth=truthdose,markers=n,estimate=estimate if estimate is not None else '',rounded_call=called if called is not None else '',correct=int(called==truthdose) if called is not None else 0,status='uncalibrated_dosage' if called is not None else 'insufficient_markers'))
            print(donor,locus,'markers',int(keep.sum()),'depth',round(scale,2),flush=True)
    write(args.output/'predictions.tsv',predictions)
    write(args.output/'dosage_baseline.tsv',dosage)
    write(args.output/'folds.tsv',folds)
    summary=[]
    for locus in ['RCCX','DRB']:
        for method in ['path_multiplicity','binary_paths','flat_sequence_catalogue','assembly_profile_oracle']:
            p=[r for r in predictions if r['locus']==locus and r['method']==method]
            summary.append(dict(locus=locus,method=method,n=len(p),pair_correct=sum(r['pair_correct'] for r in p),cn_correct=sum(r['cn_correct'] for r in p),unique_calls=sum(r['call_status']=='unique_best_uncalibrated' for r in p),unrepresented_truth=sum(not r['represented_in_training'] for r in p)))
    write(args.output/'summary.tsv',summary)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
