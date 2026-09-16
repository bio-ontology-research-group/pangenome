#!/usr/bin/env python3
"""Isogenic 31-mer probes for SNPs distinguishing assembly CDS from alternatives.

Probes use genomic context, not artificial exon junctions. They test local bases,
not full allele phasing. Indels and adjacent differences require separate review.
"""
from collections import defaultdict
from Bio import SeqIO, Align
import json
from sequence_catalogue import BASE, OUT, read_gtf, merged_intervals, reference_index, write_table
from check_results import table


def main():
    catalogue={r['name']:r for r in table(OUT/'sequence_catalogue.tsv')}
    cds={r.id:str(r.seq) for r in SeqIO.parse(OUT/'coding_sequences.fa','fasta')}
    genomic={r.id:str(r.seq) for r in SeqIO.parse(OUT/'gene_sequences.fa','fasta')}
    tasks=[]
    for r in table(OUT/'coding_candidates.tsv'):
        for h in r['haplotypes'].split(';'):
            if h.startswith(('HG','NA')):
                tasks.append((r['candidate_id'],h+'#'+r['gene'],r['nearest_reference_alleles'].split(';'),'candidate'))
    for r in table(OUT/'experimental_discordances.tsv'):
        if r['assessment']!='unresolved_antigen_exon_discordance': continue
        for hap in ['1','2']:
            prefixes=[r['gene'].removeprefix('HLA-')+'*'+a for call in r['experimental'].split(';') for a in call.split('/')]
            tasks.append(('EXPERIMENT-'+r['donor']+'-'+r['gene'],r['donor']+'#'+hap+'#'+r['gene'],prefixes,'discordance'))
    # Prefer substitutions over equally costly compensating indels.
    aligner=Align.PairwiseAligner(mode='global',match_score=0,mismatch_score=-1,open_gap_score=-1.01,extend_gap_score=-1)
    cache={}; markers=[]; audits=[]
    for case,rid,alternatives,kind in tasks:
        row=catalogue[rid]; gene=row['gene']; h=row['hap_id']; sample,hap=h.split('#')
        if gene not in cache: cache[gene]=reference_index(gene.removeprefix('HLA-'),'nuc')
        refs={s:[a for a in aa if any(a==x or a.startswith(x+':') for x in alternatives)] for s,aa in cache[gene].items()}
        refs={s:aa for s,aa in refs.items() if aa}
        query=cds[rid]
        ranked=sorted((int(-aligner.score(query,s)),s,aa) for s,aa in refs.items())
        if not ranked: raise ValueError(rid)
        distance,ref,alleles=ranked[0]
        # Stable lexical choice among equally near reference CDS; all ties recorded.
        tied=[aa for d,s,aa in ranked if d==distance]
        ann=[g for g in read_gtf(BASE/'source/gtf_all'/f'{sample}_{hap}.gtf.gz') if g.get('gene_name')==gene]
        assert len(ann)==1,rid
        g=ann[0]; spans=merged_intervals(g['features']['CDS']+g['features']['stop_codon'])
        positions=[p for lo,hi in spans for p in range(lo,hi+1)]
        if g['strand']=='-': positions.reverse()
        offsets=[p-g['start'] if g['strand']=='+' else g['end']-p for p in positions]
        assert ''.join(genomic[rid][p] for p in offsets)==query,rid
        aln=aligner.align(query,ref)[0]
        coords=aln.coordinates
        snps=[]; indels=[]
        for i in range(coords.shape[1]-1):
            a,b=map(int,coords[:,i]); z,w=map(int,coords[:,i+1])
            if z-a==w-b:
                snps.extend((x,b+x-a) for x in range(a,z) if query[x]!=ref[b+x-a])
            else: indels.append([a,z,b,w])
        for qpos,rpos in snps:
            pos=offsets[qpos]; seq=genomic[rid]
            if pos<15 or pos+15>=len(seq):continue
            observed=seq[pos-15:pos+16]
            alt=list(observed)
            for nearby_q,nearby_r in snps:
                local=offsets[nearby_q]-(pos-15)
                if 0<=local<31: alt[local]=ref[nearby_r]
            alt=''.join(alt)
            assert observed[15]==query[qpos]
            other_variants=sum(abs(offsets[p]-pos)<=30 for p,_ in snps)-1
            markers.append(dict(marker_id=f'{case}-{h}-cds{qpos+1}',case_id=case,kind=kind,donor=sample,hap_id=h,gene=gene,cds_position=qpos+1,genomic_contig=g['contig'],genomic_position=positions[qpos],gene_strand=g['strand'],assembly_base=query[qpos],alternative_base=ref[rpos],assembly_kmer=observed,alternative_kmer=alt,nearest_reference_alleles=';'.join(alleles),cds_edit_distance=distance,other_snps_within_30bp=other_variants,indel_blocks=len(indels),interpretation='Local allele-specific evidence; isogenic alternative context; not full allele or haplotype validation'))
        audits.append(dict(case_id=case,name=rid,kind=kind,minimum_cds_edit_distance=distance,nearest_reference_alleles=';'.join(alleles),equally_near_reference_alleles=json.dumps(tied),snps=len(snps),indel_blocks=json.dumps(indels)))
    write_table(OUT/'read_markers.tsv',markers)
    write_table(OUT/'read_marker_alignment_audit.tsv',audits)
    (OUT/'read_marker_donors.txt').write_text('\n'.join(sorted({r['donor'] for r in markers}))+'\n')
    print(len(markers),'markers',len({r['donor'] for r in markers}),'donors')


if __name__=='__main__':main()
