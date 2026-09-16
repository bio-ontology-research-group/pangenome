#!/usr/bin/env python3
"""Conservative evidence ranking from targeted read and assembly checks."""
from collections import defaultdict
import json
from sequence_catalogue import OUT, write_table
from check_results import table


def main():
    root=OUT/'read_support'
    specific={r['marker_id']:r for r in table(root/'marker_specificity.tsv')}
    support={r['marker_id']:r for p in root.glob('*/support.tsv') for r in table(p)}
    markers=table(OUT/'read_markers.tsv'); details=[]; bycase=defaultdict(list)
    for r in markers:
        s=support.get(r['marker_id']); t=specific[r['marker_id']]
        tested=bool(s); count=int(s['assembly_fragments']) if s else 0
        unique=int(t['assembly_mhc_hits'])==1
        # Exploratory evidence threshold, not a calibrated clinical call rule.
        supported=tested and count>=5 and unique
        item=dict(marker_id=r['marker_id'],case_id=r['case_id'],donor=r['donor'],hap_id=r['hap_id'],gene=r['gene'],cds_position=r['cds_position'],reads_available=int(tested),assembly_fragments=count if tested else '',alternative_fragments=s['alternative_fragments'] if s else '',assembly_mhc_hits=t['assembly_mhc_hits'],alternative_mhc_hits=t['alternative_mhc_hits'],local_support_at_least_5_and_mhc_unique=int(supported),status='local_read_support' if supported else 'nonunique_mhc_probe' if tested and not unique else 'insufficient_reads' if tested else 'local_cram_unavailable')
        details.append(item);bycase[r['case_id']].append(item)
    write_table(OUT/'read_marker_evidence.tsv',details)
    audit=defaultdict(list)
    for r in table(OUT/'read_marker_alignment_audit.tsv'):audit[r['case_id']].append(r)
    ranked=[]
    for r in table(OUT/'coding_candidates.tsv'):
        hits=bycase[r['candidate_id']]
        validated_donors=[]
        for donor in sorted({x['donor'] for x in hits}):
            local=[x for x in hits if x['donor']==donor]
            complete=all(x['local_support_at_least_5_and_mhc_unique'] for x in local)
            audits=[x for x in audit[r['candidate_id']] if x['name'].startswith(donor+'#')]
            complete &= bool(audits) and all(x['indel_blocks']=='[]' and int(x['snps'])==len([y for y in local if y['hap_id']+'#'+y['gene']==x['name']]) for x in audits)
            if complete:validated_donors.append(donor)
        category='plausible_novel_allele_candidate' if validated_donors else 'unresolved'
        ranked.append(dict(candidate_id=r['candidate_id'],gene=r['gene'],category=category,complete_cds=1,exact_current_protein_match=0,assembly_donors=r['independent_name_reconciled_donors'],donors_with_all_distinguishing_snps_locally_supported=';'.join(validated_donors),markers_tested=sum(x['reads_available'] for x in hits),markers_mhc_unique_and_supported=sum(x['local_support_at_least_5_and_mhc_unique'] for x in hits),nearest_reference_protein_edit_distance=r['nearest_reference_protein_edit_distance'],limitations='Local short-read SNP support does not phase the whole allele; MHC-only recruitment and MHC-only uniqueness; database absence is not proof of unpublished novelty; no expression/function evidence'))
    ranked.sort(key=lambda r:(r['category']!='plausible_novel_allele_candidate',-int(r['assembly_donors']),r['candidate_id']))
    write_table(OUT/'candidate_evidence_ranking.tsv',ranked)
    report=['# Targeted read validation checkpoint','',
        'DDBJ Slurm array 20604953 examined primary MHC-recruited reads after excluding secondary, supplementary, duplicate-marked and QC-failed reads. Six tasks completed; three failed because the local CRAM files for HG01252, HG01358 and HG01943 do not exist. No missing data are treated as negative evidence. Specificity job 20605048 completed.', '',
        'Genomic 31-mer probes compare annotated assembly bases with the closest selected reference CDS. Alternative probes carry all alternative CDS SNPs within the window, but retain assembly-derived intronic context. Counts require all 31 bases at Phred ≥20 and count each read-name fragment once across mates. Probes are checked against both donor MHC assemblies. This is not a genome-wide mapping or allele-phasing assay. Reference-selection ties and alignments are retained in `read_marker_alignment_audit.tsv`.', '',
        '## Coding candidates','']
    for r in ranked:
        if r['category']=='plausible_novel_allele_candidate':report.append(f"- {r['candidate_id']}: all tested distinguishing SNPs have at least five supporting fragments and a unique match in the donor MHC assemblies ({r['donors_with_all_distinguishing_snps_locally_supported']}). The full allele and novelty remain unconfirmed.")
    report += ['', '## Experimental disagreements', '',
        'NA18943 HLA-A and DRB1 have identical coding sequences in the HPRC and JaSaPaGe assemblies; local distinguishing SNPs have 9–15 supporting fragments and unique donor-MHC matches. NA19007 HLA-A has 14 supporting fragments at its distinguishing SNP, with no alternative-probe fragments. This supports the assembly bases without proving why the historical experimental labels differ.', '',
        'NA18608 DRB1 has supporting reads at the tested bases, but several probes occur in two or three places in the donor MHC assemblies. Its diploid allele discrepancy remains unresolved; these markers cannot assign all relevant bases to the intended DRB1 copy. Alternative counts can originate from the other haplotype or paralogs and are not automatically contradictory evidence.', '',
        'Raw counts, coordinates, alleles, probe sequences and specificity locations are preserved. The ≥5-fragment rule is an exploratory ranking threshold; it has not been calibrated for sensitivity or specificity. Nearby markers can count the same fragments and are not independent confirmations.', '']
    (OUT/'read_validation.md').write_text('\n'.join(report))
    print('plausible candidates with local SNP support',sum(r['category']=='plausible_novel_allele_candidate' for r in ranked))


if __name__=='__main__':main()
