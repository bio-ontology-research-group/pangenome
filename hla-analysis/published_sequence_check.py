#!/usr/bin/env python3
"""Positive exact-sequence checks against published HGSVC3 allele sequences.

Ordered exact CDS exon blocks establish sequence presence. No match is inconclusive:
synonymous changes, annotation boundaries and other publications are not exhausted.
"""
from collections import defaultdict
from Bio import SeqIO
import gzip
import json
from sequence_catalogue import BASE, OUT, read_gtf, merged_intervals, write_table
from check_results import table


def ordered_match(blocks,sequence):
    end=0
    for block in blocks:
        start=sequence.find(block,end)
        if start<0:return False
        end=start+len(block)
    return True


def main():
    published=defaultdict(list)
    with gzip.open(BASE/'source/HGSVCv3_HLA_alleles.fasta.gz','rt') as f:
        for r in SeqIO.parse(f,'fasta'):
            if '_HLA-' not in r.id:continue
            gene='HLA-'+r.id.split('_HLA-')[1].split('_')[0]
            published[gene].append((r.id,str(r.seq).upper(),r.description))
    cds={r.id:str(r.seq) for r in SeqIO.parse(OUT/'coding_sequences.fa','fasta')}
    genomic={r.id:str(r.seq) for r in SeqIO.parse(OUT/'gene_sequences.fa','fasta')}
    known=set(genomic.values())
    controls=[name for records in published.values() for name,s,desc in records if s[1000:-1000] in known]
    assert controls, 'Published archive produced no known sequence controls'
    (OUT/'published_sequence_positive_controls.json').write_text(json.dumps(dict(exact_genomic_matches=len(controls),published_records=controls),indent=2)+'\n')
    result=[]
    for r in table(OUT/'coding_candidates.tsv'):
        for h in r['haplotypes'].split(';'):
            rid=h+'#'+r['gene'];sample,hap=h.split('#')
            ann=[g for g in read_gtf(BASE/'source/gtf_all'/f'{sample}_{hap}.gtf.gz') if g.get('gene_name')==r['gene']]
            assert len(ann)==1
            g=ann[0];spans=merged_intervals(g['features']['CDS']+g['features']['stop_codon'])
            if g['strand']=='-':spans.reverse()
            lengths=[hi-lo+1 for lo,hi in spans];q=cds[rid];blocks=[];pos=0
            for n in lengths:blocks.append(q[pos:pos+n]);pos+=n
            assert pos==len(q)
            gene_matches=[];coding_matches=[]
            for name,s,desc in published[r['gene']]:
                if genomic[rid] in s:gene_matches.append(name)
                if ordered_match(blocks,s):coding_matches.append(name)
            result.append(dict(candidate_id=r['candidate_id'],hap_id=h,gene=r['gene'],published_exact_genomic=';'.join(gene_matches),published_ordered_exact_cds_blocks=';'.join(coding_matches),status='previously_published_coding_sequence' if coding_matches else 'not_resolved_by_exact_block_search',source='Logsdon2025_HGSVC3_allele_archive',limitation='Negative exact search does not establish novelty; other publications and synonymous variants are not exhausted'))
    write_table(OUT/'published_candidate_sequence_check.tsv',result)
    found={r['candidate_id'] for r in result if r['published_ordered_exact_cds_blocks']}
    print(len(result),'candidate haplotypes checked;',len(found),'candidate proteins with published exact coding-block matches')
    for r in result:
        if r['published_ordered_exact_cds_blocks']:print(r['candidate_id'],r['hap_id'],r['published_ordered_exact_cds_blocks'])


if __name__=='__main__':main()
