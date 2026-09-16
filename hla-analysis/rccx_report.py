#!/usr/bin/env python3
from collections import Counter
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sequence_catalogue import BASE, OUT, write_table
from check_results import table


def main():
    panel={r['hap_id']:r for r in table(OUT/'frozen_panel.tsv')}
    qc={r['hap_id']:r for r in table(OUT/'rccx_span_qc.tsv')}
    rows=table(OUT/'rccx_screen.tsv');eligible=[]
    for r in rows:
        h=r['hap_id'];q=qc.get(h)
        ok=bool(panel[h]['strict_structural_screen_eligible']=='1' and q and int(q['N_bases'])==0 and int(q['other_ambiguous_bases'])==0 and int(q['contig_edge_distance'])>=1000 and r['balanced_module_components']=='1' and r['screen_signature'])
        eligible.append(dict(r,strict_gapfree_rccx_eligible=int(ok),N_bases=q['N_bases'] if q else '',contig_edge_distance=q['contig_edge_distance'] if q else '',span_sha256=q['span_sha256'] if q else ''))
    write_table(OUT/'rccx_eligibility.tsv',eligible)
    counts=Counter(int(r['c4_annotated_copies']) for r in eligible if r['strict_gapfree_rccx_eligible'])
    summary=dict(strict_gapfree_entries=sum(counts.values()),c4_module_component_counts=dict(sorted(counts.items())),screened_entries=len(rows),limitations=['Homology count consistency does not establish functional paralog identity or gene conversion','Copy number is assembly-derived; source-read depth/junction validation remains unavailable for the matched-type structural leads','Upstream-clipped and fragmented sources excluded from strict subset','Counts describe this panel, not population frequencies'])
    (OUT/'rccx_eligibility_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    selected=['ksa006#1','ksa008#1','HG02735#2','NA21144#2']
    loci=table(OUT/'rccx_homology_loci.tsv')
    calls=table(BASE.parent/'hla-audit/2026-09-16/source/hla_calls.tsv')
    fig,ax=plt.subplots(figsize=(12,5))
    colors={'WHR1':'#b88b20','C4':'#2471a3','CYP21':'#bc4b51','TNX':'#40916c'}
    for y,h in enumerate(selected):
        q=qc[h];origin=int(q['start0']);end=(int(q['end0'])-origin)/1000
        ax.plot([0,end],[y,y],color='#aaaaaa',lw=1,zorder=1)
        events=[(int(r['start0']),int(r['end0']),r['prototype'],r['family']) for r in loci if r['hap_id']==h]
        events += [(int(r['start'])-1,int(r['end']),r['gene'],'C4') for r in calls if r['hap_id']==h and r['gene'].startswith('C4')]
        for lo,hi,label,family in events:
            x=(lo-origin)/1000;width=(hi-lo)/1000
            ax.broken_barh([(x,width)],(y-.13,.26),facecolors=colors[family],zorder=2)
            if family=='C4':ax.text(x+width/2,y-.21,label,ha='center',va='top',fontsize=8)
    ax.set_yticks(range(4),selected);ax.invert_yaxis();ax.set_xlabel('Position relative to observed RCCX span start (kb)')
    ax.set_title('RCCX component differences within matching classical HLA types\nEach adjacent pair matches all available numeric two-field calls',fontsize=11)
    ax.legend(handles=[Patch(color=c,label=g+'-like' if g!='C4' else 'C4 annotation') for g,c in colors.items()],loc='upper left',bbox_to_anchor=(0,-.16),frameon=False,ncol=4,fontsize=8)
    ax.spines[['top','right','left']].set_visible(False);ax.grid(axis='x',alpha=.15)
    fig.tight_layout();fig.savefig(OUT/'rccx_matched_hla.png',dpi=180);fig.savefig(OUT/'rccx_matched_hla.pdf');plt.close(fig)
    report=['# C4/RCCX structural screen','',
        'Reference prototypes: GRCh38 UCSC RefSeq gene spans plus the Ensembl WHR1B pseudogene span. WHR1/WHR1B are the current names corresponding to historical STK19/STK19B. Prototype coordinates, accessions, sequence and source snapshots are frozen in `source/rccx_*`. DDBJ minimap2 job 20605384 and span-QC job 20605526 completed.', '',
        'This screen locates homologous components, retaining ≥85% prototype-span coverage and ≥95% alignment identity and collapsing overlapping alternative-prototype matches. C4 A/B and long/short labels are inherited from the original annotations. CYP21A1P-like versus CYP21A2-like is nearest-prototype similarity, not a functional or gene-conversion call. Shared TNX/WHR1 alignments contained in a longer homolog are not counted twice.', '',
        f"Across 754 entries, 750 have equal observed component counts on one contig. The conservative subset has {sum(counts.values())} retained, non-clipped, single-contig, gap-free entries with matching counts, a usable order signature and ≥1 kb distance from contig edges. Component-count distribution: {dict(sorted(counts.items()))}. These are panel counts, not population frequencies.", '',
        '![RCCX within matching HLA types](rccx_matched_hla.png)', '',
        'The ksa006#1 versus ksa008#1 contrast has three versus two C4/CYP21/TNX/WHR1 components. HG02735#2 versus NA21144#2 has one versus two. Both pairs match all available classical numeric two-field calls, including annotated DRB3/4/5 states. These distinctions therefore extend beyond the three-gene DR–DQ label comparison. The spans have no ambiguous bases. Their module-count classes are known biology; neither a new structural class nor a phenotype effect is established.', '',
        'Four entries need caution: HIFI032007D#2 and HIFI032462D#2 have unbalanced components from upstream-clipped graph inputs; HIFI032164D#1 has five of each component but lacks a consistent C4 orientation signature and is also from clipped inputs; ksa004#2 has no RCCX calls in a fragmented ten-contig MHC extraction. These are not validated deletions, inversions or five-module discoveries.', '',
        'Four four-component entries (HG02392#1, NA18948#2, NA18952#1, apr048#1) merit targeted structural validation. The reference-based screen does not resolve chimeric CYP21/TNX genes, functional status, pathogenicity, or exact duplication breakpoints. Short-read depth/junction or independent long-read evidence is still needed before genotype deployment.', '',
        'Logsdon et al. 2025 provides direct prior art and a separate published RCCX catalogue. Its C4 patterns agree with all four overlapping evaluable diploid donors; our broader analysis should not claim the method or module classes as novel. Stable screen-signature IDs and per-haplotype coordinates are in `rccx_eligibility.tsv` and `rccx_homology_loci.tsv`.','']
    (OUT/'rccx_report.md').write_text('\n'.join(report));print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
