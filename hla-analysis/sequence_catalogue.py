#!/usr/bin/env python3
"""Audit annotated sequence against a pinned IPD release; no de novo annotation.

Exact unmatched sequences are unresolved, not automatically novel. GTF coordinates
are 1-based closed; supplied region FASTAs are already gene-strand oriented.
"""
from pathlib import Path
from collections import defaultdict, Counter
import csv
import gzip
import hashlib
import json
import re
from Bio import SeqIO
from Bio.Seq import Seq

BASE = Path(__file__).resolve().parent
SRC = BASE / 'source'
OUT = BASE / 'results'
GENES = ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DRB3', 'DRB4', 'DRB5']


def attrs(s):
    return dict(re.findall(r'(\w+)\s+"([^"]*)"', s))


def read_gtf(path):
    genes = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            col = line.rstrip().split('\t')
            a = attrs(col[8])
            gid = a.get('gene_id')
            if not gid:
                continue
            g = genes.setdefault(gid, {'features': defaultdict(list)})
            g['features'][col[2]].append((int(col[3]), int(col[4])))
            if col[2] in ('gene', 'transcript'):
                g.update(a)
                g.update(contig=col[0], start=int(col[3]), end=int(col[4]), strand=col[6])
    return list(genes.values())


def oriented_slice(seq, lo, hi, region_lo, region_hi, strand):
    assert region_lo <= lo <= hi <= region_hi
    if strand == '+':
        return seq[lo-region_lo:hi-region_lo+1]
    return seq[region_hi-hi:region_hi-lo+1]


def merged_intervals(intervals):
    out = []
    for lo, hi in sorted(intervals):
        if out and lo <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(hi, out[-1][1]))
        else:
            out.append((lo, hi))
    return out


def reference_index(gene, kind):
    stem = 'DRB345' if gene in ('DRB3', 'DRB4', 'DRB5') and kind != 'gen' else gene
    index = defaultdict(list)
    for rec in SeqIO.parse(SRC / 'imgt/fasta' / f'{stem}_{kind}.fasta', 'fasta'):
        allele = rec.description.split()[1]
        if allele.split('*')[0] == gene:
            index[str(rec.seq).upper()].append(allele)
    return index


def write_table(path, rows):
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True)
    annotations = {p.name.removesuffix('.gtf.gz'): read_gtf(p)
                   for p in sorted((SRC / 'gtf_all').glob('*.gtf.gz'))}
    allrows, summary = [], {}
    with (OUT / 'gene_sequences.fa').open('w') as genomic_out, (OUT / 'coding_sequences.fa').open('w') as cds_out, (OUT / 'flanking_sequences.fa').open('w') as flank_out:
        for gene in GENES:
            name = 'HLA-' + gene
            refs = {k: reference_index(gene, k) for k in ('gen', 'nuc', 'prot')}
            fasta = {r.id: str(r.seq).upper() for r in SeqIO.parse(SRC / f'{name}.fa', 'fasta')}
            with (SRC / f'{name}.regions.tsv').open() as f:
                regions = list(csv.DictReader(f, delimiter='\t'))
            rows = []
            for region in regions:
                rid = region['name']
                m = re.fullmatch(r'(.*):(\d+)-(\d+)', region['region'])
                contig, lo, hi = m[1], int(m[2]), int(m[3])
                seq = fasta[rid]
                assert len(seq) == hi-lo+1, rid
                ann = annotations[f"{region['sample']}_{region['haplotype']}"]
                matches = [g for g in ann if g.get('gene_name') == name and g['contig'] == contig and g['strand'] == region['strand'] and lo <= g['start'] <= g['end'] <= hi]
                row = dict(name=rid, hap_id=f"{region['sample']}#{region['haplotype']}", gene=name, cohort=region['cohort'], old_consensus=region['consensus'], region=region['region'], strand=region['strand'], annotation_matches=len(matches), template_warning='', gene_bp='', cds_bp='', protein_aa='', cds_complete=0, qc_flags='', gene_sha256='', cds_sha256='', protein_sha256='', flank_sha256='', exact_genomic_alleles='', exact_cds_alleles='', exact_protein_alleles='', assessment='unresolved_annotation_mapping')
                if len(matches) != 1:
                    rows.append(row)
                    continue
                g = matches[0]
                extract = lambda a, b: oriented_slice(seq, a, b, lo, hi, g['strand'])
                genomic = extract(g['start'], g['end'])
                # Stops may be separate GTF features: union prevents counting twice.
                spans = merged_intervals(g['features']['CDS'] + g['features']['stop_codon'])
                if g['strand'] == '-':
                    spans.reverse()
                cds = ''.join(extract(a, b) for a, b in spans)
                warning = g.get('template_warning', '')
                flags = []
                if warning not in ('', 'NA'):
                    flags.append('annotation:' + warning)
                if not cds.startswith('ATG'):
                    flags.append('no_start')
                if len(cds) % 3:
                    flags.append('incomplete_codon')
                if cds[-3:] not in ('TAA', 'TAG', 'TGA'):
                    flags.append('no_terminal_stop')
                if set(cds) - set('ACGT'):
                    flags.append('ambiguous_cds_bases')
                protein = str(Seq(cds[:len(cds)//3*3]).translate())
                if '*' in protein[:-1]:
                    flags.append('internal_stop')
                protein = protein.removesuffix('*')
                # Separator keeps the two flanks distinct and prevents junction k-mers.
                flank = extract(lo, g['start']-1) if lo < g['start'] else ''
                other = extract(g['end']+1, hi) if g['end'] < hi else ''
                if g['strand'] == '-':
                    flank, other = other, flank
                flank += 'N'*31 + other
                exact = {k: refs[k].get(s, []) for k, s in [('gen', genomic), ('nuc', cds), ('prot', protein)]}
                state = 'known_exact_genomic' if exact['gen'] else 'known_exact_cds' if exact['nuc'] else 'known_exact_protein' if exact['prot'] else 'unresolved_no_exact_match'
                if flags and not any(exact.values()):
                    state = 'unresolved_incomplete_or_abnormal_cds'
                row.update(template_warning=warning, gene_bp=len(genomic), cds_bp=len(cds), protein_aa=len(protein), cds_complete=int(not flags), qc_flags=';'.join(flags), assessment=state)
                for label, sequence in [('gene', genomic), ('cds', cds), ('protein', protein), ('flank', flank)]:
                    row[label+'_sha256'] = hashlib.sha256(sequence.encode()).hexdigest()
                for label, kind in [('genomic', 'gen'), ('cds', 'nuc'), ('protein', 'prot')]:
                    row['exact_'+label+'_alleles'] = ';'.join(sorted(exact[kind]))
                genomic_out.write(f'>{rid}\n{genomic}\n')
                cds_out.write(f'>{rid}\n{cds}\n')
                flank_out.write(f'>{rid}\n{flank}\n')
                rows.append(row)
            allrows.extend(rows)
            summary[name] = dict(n=len(rows), complete=sum(r['cds_complete'] for r in rows), assessments=dict(Counter(r['assessment'] for r in rows)), exact_genomic=sum(bool(r['exact_genomic_alleles']) for r in rows), exact_cds=sum(bool(r['exact_cds_alleles']) for r in rows), exact_protein=sum(bool(r['exact_protein_alleles']) for r in rows))
            print(name, summary[name], flush=True)
    write_table(OUT / 'sequence_catalogue.tsv', allrows)
    (OUT / 'sequence_summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
