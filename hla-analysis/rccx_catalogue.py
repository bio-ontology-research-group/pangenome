#!/usr/bin/env python3
"""Conservative RCCX homology screen; not functional paralog or fusion typing.

Retain >=85% prototype span and >=95% alignment identity. Collapse competing
homologous prototypes at overlapping loci. Existing C4 calls retain long/short
and A/B assignments; these assignments are not independently inferred here.
"""
from collections import defaultdict, Counter
import json
import hashlib
from sequence_catalogue import BASE, OUT, write_table
from check_results import table


def main():
    raw=[]; grouped=defaultdict(list)
    for p in sorted((BASE/'source/rccx_paf').glob('*.paf')):
        for line in p.read_text().splitlines():
            c=line.split('\t')
            gene=c[0]
            if gene.startswith('C4'):continue
            qlen,qs,qe=int(c[1]),int(c[2]),int(c[3])
            coverage=(qe-qs)/qlen;identity=int(c[9])/int(c[10])
            family='CYP21' if gene.startswith('CYP21') else 'TNX' if gene.startswith('TNX') else 'WHR1'
            r=dict(entry=p.stem,hap_id=p.stem.rsplit('_',1)[0]+'#'+p.stem.rsplit('_',1)[1],prototype=gene,family=family,contig=c[5],start0=int(c[7]),end0=int(c[8]),strand=c[4],query_length=qlen,query_coverage=coverage,identity=identity,mapq=int(c[11]),eligible=int(coverage>=.85 and identity>=.95))
            raw.append(r)
            if r['eligible']:grouped[r['hap_id'],c[5],family].append(r)
    if not raw:raise ValueError('No RCCX alignments found')
    write_table(OUT/'rccx_homology_hits.tsv',raw)
    loci=[]
    for (h,contig,family),hits in grouped.items():
        # Greedy representative selection, preferring complete TNXB/WHR1 over
        # contained short pseudogene alignments; CYP21 identity chooses prototype.
        hits.sort(key=lambda r:(r['identity'],r['query_coverage']) if family=='CYP21' else (r['query_length'],r['identity']),reverse=True)
        selected=[]
        for r in hits:
            overlaps=[x for x in selected if r['strand']==x['strand'] and max(0,min(r['end0'],x['end0'])-max(r['start0'],x['start0']))>=.8*min(r['end0']-r['start0'],x['end0']-x['start0'])]
            if overlaps:
                overlaps[0]['competing_prototypes'].add(r['prototype'])
            else:selected.append(dict(r,competing_prototypes={r['prototype']}))
        for r in selected:
            r['competing_prototypes']=';'.join(sorted(r['competing_prototypes']))
            loci.append(r)
    write_table(OUT/'rccx_homology_loci.tsv',loci)
    perhap=defaultdict(list);calls=defaultdict(list)
    for r in loci:perhap[r['hap_id']].append(r)
    for r in table(BASE.parent/'hla-audit/2026-09-16/source/hla_calls.tsv'):
        if r['gene'].startswith('C4'):calls[r['hap_id']].append(r)
    panel=table(BASE.parent/'hla-audit/2026-09-16/panel_manifest.tsv')
    out=[]
    for p in panel:
        h=p['hap_id'];hits=perhap[h];c4=calls[h]; counts=Counter(r['family'] for r in hits)
        events=[dict(contig=r['contig'],start0=r['start0'],end0=r['end0'],strand=r['strand'],label=r['prototype']+'-like') for r in hits]
        events += [dict(contig=r['contig'],start0=int(r['start'])-1,end0=int(r['end']),strand=r['strand'],label=r['gene']) for r in c4]
        same=len({r['contig'] for r in events})==1
        forward=len({r['strand'] for r in c4})==1 and bool(c4)
        reverse=forward and c4[0]['strand']=='-'
        ordered=sorted(events,key=lambda r:r['start0'],reverse=reverse)
        signature='|'.join(r['label']+('+' if (r['strand']=='+')!=reverse else '-') for r in ordered) if same and forward else ''
        balanced=bool(c4) and same and counts['CYP21']==counts['TNX']==counts['WHR1']==len(c4)
        out.append(dict(hap_id=h,donor=p['donor_id'],cohort=p['cohort'],source_type=p['source_type'],panel_flags=p['flags'],c4_annotated_copies=len(c4),cyp21_like_loci=counts['CYP21'],tnx_like_loci=counts['TNX'],whr1_like_loci=counts['WHR1'],same_contig=int(same),span_contig=ordered[0]['contig'] if same and ordered else '',span_start0=min((x['start0'] for x in events),default=0) if same else '',span_end0=max((x['end0'] for x in events),default=0) if same else '',balanced_module_components=int(balanced),screen_signature=signature,signature_id='RCCX-SCREEN-'+hashlib.sha256(signature.encode()).hexdigest()[:16] if signature else '',classification='component_counts_consistent' if balanced else 'unresolved_component_or_contiguity',limitation='Nearest prototype is not functional gene identity; no fusion or gene-conversion classification; C4 annotation inherited'))
    write_table(OUT/'rccx_screen.tsv',out)
    summary=dict(panel_entries=len(out),balanced_component_entries=sum(r['balanced_module_components'] for r in out),classifications=dict(Counter(r['classification'] for r in out)),component_counts=dict(Counter(str((r['c4_annotated_copies'],r['cyp21_like_loci'],r['tnx_like_loci'],r['whr1_like_loci'])) for r in out)))
    (OUT/'rccx_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
