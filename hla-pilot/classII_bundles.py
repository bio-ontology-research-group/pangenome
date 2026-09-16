#!/usr/bin/env python3
"""Class II region (HLA-DRA .. HLA-DMA) pgr-tk principal-bundle analysis, after Chin (ASHI 2023):
 - MAP-graph (pmapg.gfa) drawn as a graph, nodes = principal bundles
 - bundle-distance dendrogram (classII.nwk) with leaves coloured by DRB1 allele group / secondary DRB gene
 - association between sequence-level clusters and gene-level haplotype strings
Run inside hla-viz: data/classII/classII.{nwk,pmapg.gfa,ctg.summary.tsv,bed}, data/classII_haplotype_strings.tsv
"""
import re, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import networkx as nx
from Bio import Phylo
from io import StringIO
from cohorts import SHORT, COLOR, present

H = pd.read_csv("data/classII_haplotype_strings.tsv", sep="\t", index_col=0)
H["DRB1grp"] = H["HLA-DRB1"].str.split(":").str[0]
H["DRhap"] = H["DRB345"].map(lambda v: "DRB3" if v.startswith("DRB3") else "DRB4" if v.startswith("DRB4") else "DRB5" if v.startswith("DRB5") else "none")

# ---------- dendrogram
tree = Phylo.read("data/classII/classII.nwk", "newick")
# newick leaves are integer indices in order of first appearance of the contigs in the bed / dist files
ctg_order = []
for line in open("data/classII/classII.bed"):
    if line.startswith("#"): continue
    c = line.split("\t")[0]
    if not ctg_order or ctg_order[-1] != c:
        if c not in ctg_order: ctg_order.append(c)
def hap_of(name):
    if name is not None and name.isdigit(): name = ctg_order[int(name)]
    m = re.match(r"([^#]+)#(\d)#", name)
    return f"{m.group(1)}#{m.group(2)}" if m else name
grp_col = {}
cmap = plt.get_cmap("tab20")
groups = sorted(H.DRB1grp.unique())
for i, g in enumerate(groups):
    grp_col[g] = cmap(i % 20)
drcol = {"DRB3": "#1f77b4", "DRB4": "#ff7f0e", "DRB5": "#2ca02c", "none": "#7f7f7f"}
tree.ladderize()
tips = list(tree.get_terminals())  # labels must follow the displayed leaf order
fig = plt.figure(figsize=(14, 42))
ax = fig.add_axes([0.02, 0.01, 0.6, 0.98])
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: "", show_confidence=False)
ax.set_xlabel("bundle distance"); ax.set_ylabel("")
# tip order (y positions as drawn by Phylo.draw: 1..n top-down in ladderized order)
order = [hap_of(t.name) for t in tips]
ymax = ax.get_ylim()
for i, h in enumerate(order):
    y = i + 1
    r = H.loc[h] if h in H.index else None
    if r is None: continue
    ax.add_patch(plt.Rectangle((ax.get_xlim()[1] * 0.92, y - 0.5), ax.get_xlim()[1] * 0.03, 1, color=grp_col[r.DRB1grp], lw=0))
    ax.add_patch(plt.Rectangle((ax.get_xlim()[1] * 0.955, y - 0.5), ax.get_xlim()[1] * 0.03, 1, color=drcol[r.DRhap], lw=0))
    ax.add_patch(plt.Rectangle((ax.get_xlim()[1] * 0.99, y - 0.5), ax.get_xlim()[1] * 0.03, 1, color=COLOR.get(r.cohort, "#000000"), lw=0))
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=grp_col[g]) for g in groups] + [plt.Rectangle((0, 0), 1, 1, color=drcol[k]) for k in drcol],
          labels=groups + [f"secondary {k}" for k in drcol], fontsize=7, loc="upper left", ncol=2, title="DRB1 group | secondary DRB gene")
ax.add_artist(ax.get_legend())
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=COLOR[c]) for c in present(H.cohort)], labels=[SHORT[c] for c in present(H.cohort)], fontsize=7, loc="lower left", title="cohort (third strip)")
ax.set_title(f"pgr-tk bundle-distance dendrogram of the class II region (DRA..DMA), {len(order)} haplotypes; colour strips: DRB1 group, secondary DRB gene, cohort", fontsize=10, loc="left")
fig.savefig("figures/fig9_classII_dendrogram.png", dpi=110)
pd.DataFrame({"y": range(1, len(order)+1), "hap_id": order}).to_csv("data/dendrogram_tip_order.tsv", sep="\t", index=False)
if "--dendrogram-only" in sys.argv:
    sys.exit(0)

# ---------- cluster / gene-level association: cut the tree into clusters at a distance threshold
def clusters_at(tree, k):
    """cut into k clusters by greedily splitting the root-most clades"""
    clades = [tree.root]
    while len(clades) < k:
        # split the clade with the largest branch-depth spread among its children
        best = max((c for c in clades if not c.is_terminal()), key=lambda c: max((ch.branch_length or 0) for ch in c.clades), default=None)
        if best is None: break
        clades.remove(best); clades.extend(best.clades)
    return clades
rows = []
for k in (10, 20, 40, 80):
    cl = clusters_at(tree, k)
    pur_d, pur_s, sizes = [], [], []
    for c in cl:
        members = [hap_of(t.name) for t in c.get_terminals()]
        members = [m for m in members if m in H.index]
        if len(members) < 2: continue
        vc = H.loc[members, "HLA-DRB1"].value_counts(); pur_d.append(vc.iloc[0] / len(members))
        vs = H.loc[members, ["DRB345", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1"]].agg("-".join, axis=1).value_counts(); pur_s.append(vs.iloc[0] / len(members))
        sizes.append(len(members))
    rows.append((k, len(sizes), np.average(pur_d, weights=sizes), np.average(pur_s, weights=sizes)))
pur = pd.DataFrame(rows, columns=["clusters_requested", "clusters_with_2plus", "purity_DRB1_2field", "purity_DRB345-DRB1-DQA1-DQB1"])
print(pur.to_string())
pur.to_csv("data/classII_cluster_purity.tsv", sep="\t", index=False)

# ---------- principal-bundle graph: nodes = principal bundles (bed column 4 "bundle:len:dir:..."), edges = consecutive
# bundles along each haplotype, weighted by the number of haplotypes traversing them (bundle-level MAP-graph view)
G = nx.DiGraph(); blen = {}; nhap = {}
prev = None; cur = None
for line in open("data/classII/classII.bed"):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t"); b_id = f[3].split(":")[0]
    blen[b_id] = max(blen.get(b_id, 0), int(f[2]) - int(f[1])); nhap.setdefault(b_id, set()).add(f[0])
    if f[0] != cur: cur, prev = f[0], None
    if prev is not None and prev != b_id:
        if G.has_edge(prev, b_id): G[prev][b_id]["w"] += 1
        else: G.add_edge(prev, b_id, w=1)
    prev = b_id
print("bundle graph nodes", G.number_of_nodes(), "edges", G.number_of_edges())
pos = nx.kamada_kawai_layout(G.to_undirected())
fig, ax = plt.subplots(figsize=(16, 13))
nodes = list(G.nodes()); sizes = np.array([blen[n] for n in nodes]); nh = np.array([len(nhap[n]) for n in nodes])
w = np.array([G[u][v]["w"] for u, v in G.edges()])
nx.draw_networkx_edges(G, pos, ax=ax, arrows=False, alpha=0.35, width=0.3 + 2.5 * w / w.max(), edge_color="#555555")
sc = nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=nodes, node_size=8 + 120 * np.log10(sizes + 1) / np.log10(sizes.max() + 1), node_color=nh, cmap="viridis", linewidths=0)
plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.01, label="haplotypes carrying the bundle")
ax.set_axis_off(); ax.set_title(f"pgr-tk principal-bundle graph of the class II region (DRA..DMA), {len(order)} haplotypes: {G.number_of_nodes()} bundles, {G.number_of_edges()} adjacencies\nnode size = bundle length, colour = number of haplotypes; edge width = haplotypes sharing the adjacency", fontsize=10, loc="left")
fig.savefig("figures/fig10_classII_bundle_graph.png", dpi=140, bbox_inches="tight")

# ---------- PCA on principal-bundle presence (sequence-level structure), from the .ord file
vec = {}
for line in open("data/classII/classII.ord"):
    if line.startswith("#") or not line.strip(): continue
    n, v = line.rstrip("\n").split("\t")[:2]
    vec[hap_of(n)] = np.array([int(x) for x in v.split(",")])
B = pd.DataFrame(vec).T
B = B.loc[[h for h in B.index if h in H.index]]
X = B.values - B.values.mean(0); U, S, Vt = np.linalg.svd(X, full_matrices=False); pc = U[:, :2] * S[:2]; ev = S ** 2 / np.sum(S ** 2)
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(19, 8.5))
samp = pd.Series([h.split("#")[0] for h in B.index], index=B.index)
for smp, idx in samp.groupby(samp).groups.items():
    if len(idx) == 2: ax.plot(pc[[B.index.get_loc(i) for i in idx], 0], pc[[B.index.get_loc(i) for i in idx], 1], color="#bbbbbb", lw=0.5, zorder=1)
for k, col in drcol.items():
    m = (H.loc[B.index, "DRhap"] == k).values
    ax.scatter(pc[m, 0], pc[m, 1], s=14, color=col, alpha=0.8, zorder=2, label=f"secondary {k}" if k != "none" else "no DRB3/4/5")
for grp, g in H.loc[B.index].groupby("DRB1grp"):
    if len(g) >= 8:
        ii = [B.index.get_loc(i) for i in g.index]; ax.text(np.median(pc[ii, 0]), np.median(pc[ii, 1]), grp, fontsize=8, ha="center", va="center", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
ax.set_xlabel(f"PC1 ({ev[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({ev[1]*100:.1f}%)"); ax.legend(fontsize=8)
coh = H.loc[B.index, "cohort"].values
for c in present(coh)[::-1]:
    m = coh == c; ax2.scatter(pc[m, 0], pc[m, 1], s=16 if c == "HPRC-Rest" else 24, color=COLOR[c], alpha=0.55 if c == "HPRC-Rest" else 0.9, zorder=1 if c == "HPRC-Rest" else 2, edgecolors="none", label=f"{SHORT[c]} (n={m.sum()})")
ax2.set_xlabel(f"PC1 ({ev[0]*100:.1f}%)"); ax2.set_ylabel(f"PC2 ({ev[1]*100:.1f}%)"); ax2.legend(fontsize=8); ax2.set_title("Same PCA coloured by cohort", fontsize=10)
ax.set_title(f"Diplotype PCA on pgr-tk principal-bundle presence in the class II region ({B.shape[1]} bundles, {len(B)} haplotypes)\nlines join the two haplotypes of an individual; labels = DRB1 allele group", fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig11_classII_bundle_pca.png", dpi=160)
print("bundle matrix", B.shape)
