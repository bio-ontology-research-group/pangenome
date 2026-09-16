#!/usr/bin/env python3
"""Apply conservative source QC and known-database checks to within-type contrasts."""
from collections import defaultdict, Counter
from Bio import SeqIO
import json
from sequence_catalogue import OUT, write_table
from check_results import table


def main():
    panel={r['hap_id']:r for r in table(OUT/'frozen_panel.tsv')}
    seq=defaultdict(list)
    for r in table(OUT/'sequence_catalogue.tsv'):seq[r['hap_id'],r['gene']].append(r)
    flanks={r.id:str(r.seq) for r in SeqIO.parse(OUT/'flanking_sequences.fa','fasta')}
    out=[]
    for r in table(OUT/'within_type_contrasts.tsv'):
        variants=json.loads(r['signature_members'])
        kept={k:[h for h in v if panel[h]['strict_structural_screen_eligible']=='1'] for k,v in variants.items()}
        kept={k:v for k,v in kept.items() if v};haps=[h for v in kept.values() for h in v]
        donors={panel[h]['donor_id'] for h in haps}
        status='unresolved';checks=[];gene=''
        if len(kept)<2 or len(donors)<2:status='excluded_after_strict_source_qc'
        elif r['feature'] in ['gene_content','gene_order','c4_annotation_order']:
            status='known_variation_class_structural_candidate'
            checks.append('Annotation-level structural distinction; candidate-level module/sequence validation required; not a novel structural class')
        else:
            gene,kind=r['feature'].rsplit('_',1)
            entries=[x for h in haps for x in seq[h,gene]]
            if kind=='gene' and all(x['exact_genomic_alleles'] for x in entries):
                status='known_database_genomic_sequence_variation'
            elif kind in ['cds','protein']:
                field='exact_cds_alleles' if kind=='cds' else 'exact_protein_alleles'
                status='known_database_sequence_variation' if all(x[field] for x in entries) else 'unresolved_sequence_difference'
            elif kind=='flank':
                status='unresolved_flanking_sequence_difference'
                pieces=[flanks[x['name']].split('N'*31) for x in entries]
                checks.append('all_flanks_2000bp_each='+str(all(len(v)==2 and all(len(s)==2000 for s in v) for v in pieces)))
                checks.append('all_flanks_unambiguous='+str(all(set(s)<=set('ACGT') for v in pieces for s in v)))
                checks.append('No regulatory function or sequence novelty inferred')
            else:status='unresolved_genomic_sequence_difference'
        recurrence=max((len({panel[h]['donor_id'] for h in v}) for v in kept.values()),default=0)
        out.append(dict(contrast_id=r['contrast_id'],matching_level=r['matching_level'],feature=r['feature'],category=status,retained_variants=len(kept),retained_donors=len(donors),retained_haplotypes=len(haps),maximum_variant_donor_support=recurrence,signature_members=json.dumps(kept,sort_keys=True),checks=';'.join(checks),matching_label=r['matching_label']))
    write_table(OUT/'within_type_evidence.tsv',out)
    summary={level:dict(Counter(r['category'] for r in out if r['matching_level']==level)) for level in ['drdq_label','all_classical_label']}
    (OUT/'within_type_evidence_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
