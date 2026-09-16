#!/usr/bin/env python3
"""Reconcile saved panel metadata; no sequence inference or cluster execution.

Run: python3 hla-audit/2026-09-16/audit_panel.py
Input provenance and checksums: source/manifest.json.
"""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent
SRC = BASE / "source"


def read(name):
    with (SRC / name).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(name, rows):
    with (BASE / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    calls = read("hla_calls.tsv")
    by_hap = defaultdict(list)
    for row in calls:
        by_hap[row["hap_id"]].append(row)
    graph = defaultdict(list)
    for line in (SRC / "graph_contig_names.tsv").read_text().splitlines():
        sample, hap, renamed, original = line.split("\t")
        graph[f"{sample}#{hap}"].append((renamed, original))
    graph_paths = {}
    for kind in ("full", "clipped"):
        paths = defaultdict(list)
        for line in (SRC / f"graph_{kind}_path_lengths.tsv").read_text().splitlines():
            name, length = line.split("\t")
            paths["#".join(name.split("#")[:2])].append((name, int(length)))
        assert set(paths) == set(graph), f"{kind} graph path IDs differ from input IDs"
        graph_paths[kind] = paths
    regions = {r["sample"] + "#" + r["haplotype"]: r for r in read("classII.regions.tsv")}
    bundles = {r["#ctg"].rsplit("#", 1)[0]: r for r in read("classII.ctg.summary.tsv")}
    core = ("HLA-DRB1", "HLA-DQA1", "HLA-DQB1")
    known_duplicate_samples = {"NA18940", "NA18943", "NA18945", "NA18952", "NA18970"}
    manifest = []
    for hap_id, contigs in sorted(graph.items()):
        records = by_hap.get(hap_id, [])
        if not records:
            raise ValueError(f"No annotations for {hap_id}")
        first = records[0]
        sample = first["sample"]
        donor = sample
        if sample.endswith(".JaSaPaGe") and sample.split(".")[0] in known_duplicate_samples:
            donor = sample.split(".")[0]
        genes = defaultdict(list)
        for record in records:
            genes[record["gene"]].append(record)
        core_counts = [len(genes[g]) for g in core]
        unique_core = all(n == 1 for n in core_counts)
        one_contig = unique_core and len({genes[g][0]["contig"] for g in core}) == 1
        flags = []
        reference = first["cohort"] == "REF"
        clipped_source = first["cohort"] in {"CPC-Chinese", "KPanRef-Korean"}
        if reference:
            flags.append("reference_not_independent_donor")
        if clipped_source:
            flags.append("upstream_clipped_graph_source")
        if donor in known_duplicate_samples:
            flags.append("cross_cohort_duplicate_donor")
        if sample == "NA18952.JaSaPaGe":
            flags.append("previously_reported_duplicated_haplotype")
        if not unique_core:
            flags.append("DR_DQ_core_missing_or_multiple_annotations")
        elif not one_contig:
            flags.append("DR_DQ_core_split_contigs")
        if hap_id not in regions:
            flags.append("absent_from_existing_DRA_DMA_span")
        novel = any(r["novel"].lower() == "true" for g in core for r in genes[g])
        if novel:
            flags.append("DR_DQ_candidate_novel_annotation")
        # Syntactic screening only: no sequence or nomenclature validation here.
        available_twofield = unique_core and all(
            genes[g][0]["allele_2field"] not in {"", "NA"}
            and ":" in genes[g][0]["allele_2field"]
            and all(v.isdigit() and int(v) > 0 for v in genes[g][0]["allele_2field"].split("*")[-1].split(":"))
            for g in core
        )
        full_bp = sum(n for _, n in graph_paths["full"][hap_id])
        clipped_bp = sum(n for _, n in graph_paths["clipped"][hap_id])
        if clipped_bp < full_bp:
            flags.append("downstream_graph_clipping_loss")
        row = dict(
            hap_id=hap_id, sample=sample, haplotype=first["haplotype"],
            cohort=first["cohort"], donor_id=donor,
            donor_identity_status="name_based_pedigree_not_verified",
            graph_input_contigs=len(contigs),
            full_graph_paths=len(graph_paths["full"][hap_id]),
            clipped_graph_paths=len(graph_paths["clipped"][hap_id]),
            full_graph_path_bp=full_bp, clipped_graph_path_bp=clipped_bp,
            downstream_clipped_bp=full_bp - clipped_bp,
            source_type="clipped_graph_paths" if clipped_source else "assembly_or_reference",
            classII_span_present=int(hap_id in regions),
            bundle_summary_present=int(hap_id in bundles),
            bundle_coverage_pct=bundles.get(hap_id, {}).get("total_bundle_coverage_percentage", "NA"),
            DR_DQ_one_copy_each=int(unique_core), DR_DQ_same_contig=int(one_contig),
            DR_DQ_numeric_twofield_labels=int(available_twofield),
            independent_truth_status="not_audited", flags=";".join(flags),
        )
        for gene in core:
            row[gene + "_copies"] = len(genes[gene])
            row[gene + "_calls"] = "|".join(r["consensus"] for r in genes[gene]) or "NA"
        manifest.append(row)
    assert set(by_hap) == set(graph), "Annotation and graph input IDs differ"
    assert len(graph) == len((SRC / "graph_seqfile.txt").read_text().splitlines()), "Input count mismatch"
    write("panel_manifest.tsv", manifest)

    donors = defaultdict(set)
    for row in manifest:
        if row["cohort"] != "REF":
            donors[row["donor_id"]].add(row["sample"])
    summary = dict(
        graph_input_haplotype_entries=len(graph),
        graph_input_contigs=sum(map(len, graph.values())),
        annotation_haplotype_entries=len(by_hap),
        cohort_counts=dict(Counter(r["cohort"] for r in manifest)),
        biological_sample_labels=len({r["sample"] for r in manifest if r["cohort"] != "REF"}),
        provisional_unique_donor_ids=len(donors),
        cross_cohort_duplicate_donors={k: sorted(v) for k, v in donors.items() if len(v) > 1},
        classII_span_entries=len(regions),
        missing_classII_span=sorted(set(graph) - set(regions)),
        DR_DQ_one_copy_same_contig=sum(r["DR_DQ_same_contig"] for r in manifest),
        DR_DQ_numeric_twofield_labels=sum(r["DR_DQ_numeric_twofield_labels"] for r in manifest),
        full_graph_paths=sum(len(p) for p in graph_paths["full"].values()),
        clipped_graph_paths=sum(len(p) for p in graph_paths["clipped"].values()),
        downstream_clipped_bp=sum(r["downstream_clipped_bp"] for r in manifest),
        flag_counts=dict(Counter(flag for r in manifest for flag in r["flags"].split(";") if flag)),
        status="Metadata screening only; no held-out accuracy, novel-allele validation or pedigree verification",
    )
    try:
        from Bio import Phylo
    except ImportError:
        summary["tree_annotation_check"] = "not run: Biopython unavailable"
    else:
        tree = Phylo.read(SRC / "classII.nwk", "newick")
        before = [t.name for t in tree.get_terminals()]
        tree.ladderize()
        after = [t.name for t in tree.get_terminals()]
        summary["tree_tip_count"] = len(before)
        summary["tip_positions_changed_by_ladderize"] = sum(a != b for a, b in zip(before, after))
    (BASE / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
