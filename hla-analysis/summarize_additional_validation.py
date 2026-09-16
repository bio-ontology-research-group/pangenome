#!/usr/bin/env python3
"""Summarise local flank evidence and exploratory, reference-based C4 depth."""
from collections import defaultdict
import json
import numpy as np
from sequence_catalogue import BASE, OUT, write_table
from check_results import table


def main():
    root=OUT/'flank_validation'
    specificity={r['marker_id']:r for r in table(root/'marker_specificity.tsv')}
    flank=[]
    for f in sorted(root.glob('*/support.tsv')):
        for r in table(f):
            s=specificity[r['marker_id']]
            flank.append(dict(contrast_id=r['case_id'],donor=r['donor'],gene=r['gene'],flank_offset_1based=r['flank_offset_1based'],assembly_fragments=int(r['assembly_fragments']),alternative_fragments=int(r['alternative_fragments']),assembly_mhc_hits=int(s['assembly_mhc_hits']),alternative_mhc_hits=int(s['alternative_mhc_hits']),assessment='unresolved_read_assembly_discordance' if int(r['assembly_fragments'])<5 and int(r['alternative_fragments'])>=5 else 'local_base_supported_phase_not_independently_validated'))
    write_table(OUT/'flank_validation_summary.tsv',flank)
    targets=defaultdict(set)
    for line in (BASE/'source/rccx_depth_targets.bed').read_text().splitlines():
        chrom,start,end,name=line.split('\t')
        targets[name].update(range(int(start)+1,int(end)+1))
    structural=table(OUT/'rccx_screen.tsv');depthrows=[];cn=[]
    for f in sorted((OUT/'rccx_read_depth').glob('*.depth.tsv')):
        donor=f.name.split('.')[0];depth={int(c[1]):int(c[2]) for c in (line.split('\t') for line in f.read_text().splitlines())}
        means={}
        for name,positions in targets.items():
            assert positions<=depth.keys(),(donor,name,len(positions-depth.keys()))
            values=np.array([depth[p] for p in positions]);means[name]=float(values.mean())
            depthrows.append(dict(donor=donor,region=name,bases=len(values),mean_depth=float(values.mean()),median_depth=float(np.median(values)),zero_bases=int((values==0).sum())))
        c4=(means['C4A']+means['C4B'])/2
        expected=sum(int(r['c4_annotated_copies']) for r in structural if r['donor']==donor)
        cn.append(dict(donor=donor,assembly_diploid_c4=expected,c4_mean_depth=c4,estimated_copy_using_whr1=4*c4/means['WHR1_unique_region'],estimated_copy_using_tnxb=4*c4/means['TNXB_distal_region'],assessment='exploratory_depth_not_validated_copy_number',limitation='GRCh38-recruited reads; reference/paralog mapping and coverage bias; two controls expose normalization sensitivity; no haplotype or junction validation'))
    write_table(OUT/'rccx_depth_summary.tsv',depthrows);write_table(OUT/'rccx_depth_copy_check.tsv',cn)
    report=['# Additional candidate validation','',
        'Flank read array 20605749 and specificity job 20605750 completed. Probes use exact 31-base genomic windows, Phred ≥20 across the window and deduplicated read names after removal of duplicate-marked/QC-failed/secondary/supplementary reads. Both donor MHC assemblies were searched for probe matches. Whole-allele phase and regulatory function are not tested.', '',
        '## Flanking candidates','',
        'HLA-C contrast HLA-WITHIN-456a27fd141204be: the downstream position 1564 allele has 13 supporting fragments in HG03195 and 17 in HG03130. HG03195 has 11 alternative-probe fragments, consistent with its two different assembled MHC haplotypes. HG03130 has no alternative-probe fragments and its observed probe occurs on both assembled haplotypes. This supports the local sequence difference; it does not independently phase it to the classical HLA allele.', '',
        'DRB1 contrast HLA-WITHIN-f7a4b2ab22ceb184: HG03130 has 12 fragments supporting its downstream position 1812 allele and zero alternative-probe fragments. HG03195 has only two assembly-probe fragments versus six alternative-probe fragments; the alternative probe is absent from both assembled MHC haplotypes. Retain this as unresolved read/assembly discordance, with a possible assembly or phasing artefact. Do not promote it as a validated regulatory variant.', '',
        '## RCCX depth','',
        'Job 20605817 measured original GRCh38-aligned primary read depth across C4A/B coding exons and two nearby control regions. Base quality ≥20; mapping quality ≥0 intentionally retains repetitive-region alignments. The rough diploid copy estimate is four times the mean of C4A/B exon depths divided by the control depth, because GRCh38 represents two C4 genes per haploid reference. This estimator is exploratory and uncalibrated.', '',
        '| Donor | Assembly diploid C4 copies | WHR1-normalised estimate | TNXB-normalised estimate |','|---|---:|---:|---:|']
    for r in cn:report.append(f"| {r['donor']} | {r['assembly_diploid_c4']} | {r['estimated_copy_using_whr1']:.2f} | {r['estimated_copy_using_tnxb']:.2f} |")
    report += ['', 'The control-dependent estimates are supporting context, not validated CN calls. No new structural classification is based on them. HG02735, the other donor in the one-versus-two matched-haplotype comparison, has no local CRAM; the Saudi pair also lacks a validated read-based structural assay. Module phase, duplication junctions and CYP21/TNX functional identity remain unresolved.', '']
    (OUT/'additional_validation.md').write_text('\n'.join(report))
    print(json.dumps(cn,indent=2))


if __name__=='__main__':main()
