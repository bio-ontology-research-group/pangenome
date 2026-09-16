#!/usr/bin/env python3
"""Cross-source donor comparisons without assuming matching haplotype phase."""
from collections import defaultdict
import openpyxl
from sequence_catalogue import BASE, OUT, write_table
from check_results import table


def main():
    seq=table(OUT/'sequence_catalogue.tsv'); bysample=defaultdict(list)
    donors=sorted({r['hap_id'].split('#')[0].removesuffix('.JaSaPaGe') for r in seq if '.JaSaPaGe#' in r['hap_id']})
    for r in seq:bysample[r['hap_id'].split('#')[0],r['gene']].append(r)
    comparisons=[]
    for d in donors:
        for gene in sorted({r['gene'] for r in seq}):
            a=bysample[d,gene];b=bysample[d+'.JaSaPaGe',gene]
            if not a and not b:continue
            for feature in ['gene_sha256','cds_sha256']:
                comparisons.append(dict(donor=d,gene=gene,feature=feature,hprc_entries=len(a),jasapage_entries=len(b),hprc_signature=';'.join(sorted(r[feature] for r in a)),jasapage_signature=';'.join(sorted(r[feature] for r in b)),unordered_match=int(sorted(r[feature] for r in a)==sorted(r[feature] for r in b)),interpretation='Cross-source assemblies of one donor, not independent biological replication; shared raw-data provenance not excluded'))
    write_table(OUT/'cross_source_assembly_comparison.tsv',comparisons)
    path=BASE.parent/'literature/2026-09-16/logsdon2025_supplementary_tables.xlsx'
    workbook=openpyxl.load_workbook(path,read_only=True,data_only=True)
    rows=list(workbook['53'].values); hdr=rows[4]
    data=[dict(zip(hdr,r)) for r in rows[5:] if r and r[0] and '.hap' in str(r[0])]
    write_table(OUT/'published_rccx_logsdon2025.tsv',data)
    ours=defaultdict(list);published=defaultdict(list)
    for r in table(OUT/'haplotype_signatures.tsv'):
        if r['c4_annotation_order']:ours[r['donor']].append(r['c4_annotation_order'])
    for r in data:published[r['Haplotype'].split('.hap')[0]].append('|'.join(g for g in r['RCCX_Architecture'].split('_') if g.startswith('C4')))
    out=[]
    for d in sorted(ours.keys()&published.keys()):
        if len(ours[d])==len(published[d])==2:
            out.append(dict(donor=d,panel_unphased_c4=';'.join(sorted(ours[d])),published_unphased_c4=';'.join(sorted(published[d])),match=int(sorted(ours[d])==sorted(published[d])),interpretation='Unordered diploid C4 annotation comparison; assembly identity and full RCCX structure not established'))
    write_table(OUT/'published_c4_comparison.tsv',out)
    print('Cross-source comparisons',len(comparisons),'published C4 overlap',len(out),'matches',sum(r['match'] for r in out))


if __name__=='__main__':main()
