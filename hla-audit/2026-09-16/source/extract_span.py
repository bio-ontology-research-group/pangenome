#!/usr/bin/env python3
"""Cut a multi-gene span (e.g. the MHC class II region from HLA-DRA to HLA-DMA)
out of every MHC haplotype using the Immuannot gene coordinates in hla_calls.tsv.

For each haplotype the span runs from the outermost coordinate of the first gene
to the outermost coordinate of the last gene (+/- flank), and is written in the
orientation of the first gene, so all records read in the same direction.
Haplotypes missing either gene, or where the two genes sit on different contigs,
are skipped (listed on stderr). Records are named sample#haplotype#<name>.
"""
import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--from-gene", required=True)
    ap.add_argument("--to-gene", required=True)
    ap.add_argument("--name", required=True, help="record suffix and output prefix, e.g. classII")
    ap.add_argument("--flank", type=int, default=20000)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    fai = a.fasta + ".fai"
    if not os.path.exists(fai):
        subprocess.run(["samtools", "faidx", a.fasta], check=True)
    ctg_len = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}

    genes = defaultdict(dict)  # hap -> gene -> (contig, start, end, strand)
    with open(a.calls) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] in (a.from_gene, a.to_gene) and r["contig"] in ctg_len:
                hap = (r["cohort"], r["sample"], r["haplotype"])
                if r["gene"] not in genes[hap]:  # first copy only
                    genes[hap][r["gene"]] = (r["contig"], int(r["start"]), int(r["end"]), r["strand"])

    plan = []
    for hap in sorted(genes):
        g = genes[hap]
        if a.from_gene not in g or a.to_gene not in g:
            print(f"{hap[1]}#{hap[2]}: missing {a.from_gene if a.from_gene not in g else a.to_gene}", file=sys.stderr)
            continue
        c1, s1, e1, st1 = g[a.from_gene]
        c2, s2, e2, st2 = g[a.to_gene]
        if c1 != c2:
            print(f"{hap[1]}#{hap[2]}: genes on different contigs", file=sys.stderr)
            continue
        s = max(1, min(s1, s2) - a.flank)
        e = min(ctg_len[c1], max(e1, e2) + a.flank)
        plan.append((f"{hap[1]}#{hap[2]}#{a.name}", c1, s, e, st1, hap))

    os.makedirs(a.outdir, exist_ok=True)
    out_fa = os.path.join(a.outdir, f"{a.name}.fa")
    out_tsv = os.path.join(a.outdir, f"{a.name}.regions.tsv")
    with open(out_fa, "w") as fa, open(out_tsv, "w") as tsv:
        tsv.write("name\tcohort\tsample\thaplotype\tcontig\tstart\tend\tstrand\tlength\n")
        for name, c, s, e, st, hap in plan:
            cmd = ["samtools", "faidx"] + (["-i"] if st == "-" else []) + [a.fasta, f"{{{c}}}:{s}-{e}"]
            seq = "".join(l.strip() for l in subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.splitlines()[1:])
            fa.write(f">{name} {c}:{s}-{e}({st})\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i + 80] + "\n")
            tsv.write(f"{name}\t{hap[0]}\t{hap[1]}\t{hap[2]}\t{c}\t{s}\t{e}\t{st}\t{e - s + 1}\n")
    print(f"{a.name}: {len(plan)} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
