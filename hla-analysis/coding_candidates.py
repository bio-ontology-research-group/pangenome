#!/usr/bin/env python3
"""Conservative protein candidate triage; database absence is not novelty proof."""
from collections import defaultdict
from Bio import SeqIO, Align
from sequence_catalogue import BASE, OUT, reference_index, write_table
import csv
import json


def main():
    rows = list(csv.DictReader((OUT/'sequence_catalogue.tsv').open(),delimiter='\t'))
    cds = {r.id:r.seq for r in SeqIO.parse(OUT/'coding_sequences.fa','fasta')}
    panel = {r['hap_id']:r for r in csv.DictReader((BASE.parent/'hla-audit/2026-09-16/panel_manifest.tsv').open(),delimiter='\t')}
    groups = defaultdict(list)
    for r in rows:
        if r['assessment']=='unresolved_no_exact_match' and r['cds_complete']=='1':
            groups[r['gene'],r['protein_sha256']].append(r)
    aligner = Align.PairwiseAligner(mode='global',match_score=0,mismatch_score=-1,open_gap_score=-1,extend_gap_score=-1)
    cache={}; out=[]
    for (gene,sha),entries in sorted(groups.items()):
        if gene not in cache:
            cache[gene]=reference_index(gene.removeprefix('HLA-'),'prot')
        protein=str(cds[entries[0]['name']].translate()).removesuffix('*')
        best=None; closest=[]
        for ref,alleles in cache[gene].items():
            # Length difference bounds global unit-cost edit distance.
            if best is not None and abs(len(ref)-len(protein))>best:
                continue
            distance=int(-aligner.score(protein,ref))
            if best is None or distance<best:
                best=distance; closest=list(alleles)
            elif distance==best:
                closest.extend(alleles)
        donors=sorted({panel[r['hap_id']]['donor_id'] for r in entries if panel[r['hap_id']]['cohort']!='REF'})
        flags=sorted({f for r in entries for f in panel[r['hap_id']]['flags'].split(';') if f})
        out.append(dict(candidate_id='HLA-PROT-'+gene.removeprefix('HLA-')+'-'+sha[:16],gene=gene,protein_sha256=sha,protein_sequence=protein,haplotypes=';'.join(r['hap_id'] for r in entries),independent_name_reconciled_donors=len(donors),donors=';'.join(donors),cohorts=';'.join(sorted({r['cohort'] for r in entries})),nearest_reference_protein_edit_distance=best,nearest_reference_alleles=';'.join(sorted(closest)),panel_flags=';'.join(flags),assessment='unresolved_database_absent_protein',validation='not_yet_read_validated',interpretation='Exact absence from frozen IPD 3.65.0; annotation, assembly, literature and independent evidence still require assessment'))
        print(gene,sha[:12],best,len(donors),flush=True)
    write_table(OUT/'coding_candidates.tsv',out)
    (OUT/'coding_candidates_summary.json').write_text(json.dumps(dict(candidate_proteins=len(out),entries=sum(len(v) for v in groups.values()),read_validated=0),indent=2)+'\n')


if __name__=='__main__':
    main()
