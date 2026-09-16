#!/usr/bin/env python3
"""Reproducible single-SNP flank probes for two all-classical-matched donors."""
import re
from Bio import SeqIO
from sequence_catalogue import BASE, OUT, write_table
from check_results import table


def main():
    outdir=OUT/'flank_validation';outdir.mkdir(exist_ok=True)
    seq={r.id:str(r.seq) for r in SeqIO.parse(OUT/'flanking_sequences.fa','fasta')}
    out=[]
    for gene,case in [('HLA-C','HLA-WITHIN-456a27fd141204be'),('HLA-DRB1','HLA-WITHIN-f7a4b2ab22ceb184')]:
        regions={r['name']:r for r in table(BASE/'source'/f'{gene}.regions.tsv')}
        a=seq['HG03195#1#'+gene].split('N'*31)[1];b=seq['HG03130#1#'+gene].split('N'*31)[1]
        differences=[i for i,(x,y) in enumerate(zip(a,b)) if x!=y]
        assert len(a)==len(b)==2000 and len(differences)==1
        i=differences[0];assert 15<=i<1985
        for donor,own,other in [('HG03195',a,b),('HG03130',b,a)]:
            r=regions[donor+'#1#'+gene];m=re.fullmatch(r'(.*):(\d+)-(\d+)',r['region'])
            lo,hi=int(m[2]),int(m[3]);pos=hi-1999+i if r['strand']=='+' else lo+1999-i
            out.append(dict(marker_id=case+'-'+donor,case_id=case,kind='flank',donor=donor,hap_id=donor+'#1',gene=gene,cds_position='',flank_side='downstream',flank_offset_1based=i+1,genomic_contig=m[1],genomic_position=pos,gene_strand=r['strand'],assembly_base=own[i],alternative_base=other[i],assembly_kmer=own[i-15:i+16],alternative_kmer=other[i-15:i+16],interpretation='Single SNP in 2kb gene-oriented downstream flank; no regulatory effect inferred'))
    write_table(outdir/'read_markers.tsv',out)
    (outdir/'read_marker_donors.txt').write_text('HG03130\nHG03195\n')


if __name__=='__main__':main()
