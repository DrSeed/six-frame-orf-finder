# Self-contained demo: plant genes in random DNA, then see which ORFs are real.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from orf_finder import find_orfs, revcomp, CODON_TABLE

SENSE = [c for c, aa in CODON_TABLE.items() if aa != "*" and c != "ATG"]


def random_dna(n, gc, rng):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choice(list("AGCT"), size=n, p=p))


def plant_genes(genome, rng, n_genes=8):
    genome = list(genome)
    truth = []
    slots = np.linspace(1000, len(genome) - 4000, n_genes).astype(int)
    for pos in slots:
        aa = int(rng.integers(150, 600))
        gene = "ATG" + "".join(rng.choice(SENSE, size=aa)) + "TAA"
        strand = "+" if rng.random() < 0.5 else "-"
        if strand == "-":
            gene = revcomp(gene)
        genome[pos:pos + len(gene)] = list(gene)
        truth.append({"start": int(pos), "end": int(pos) + len(gene), "strand": strand, "length_aa": aa})
    return "".join(genome), truth


def longest_spurious(gc, rng, n=30_000, reps=6):
    # Average over replicates: the single longest ORF is a noisy extreme value.
    return float(np.mean([max(o["length_aa"] for o in find_orfs(random_dna(n, gc, rng), min_aa=1))
                          for _ in range(reps)]))


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(7)
    genome, truth = plant_genes(random_dna(60_000, 0.5, rng), rng)
    orfs = pd.DataFrame(find_orfs(genome, min_aa=1))
    # Match on the stop codon: an upstream in-frame ATG can legitimately extend the ORF.
    stops = {("+", t["end"]) for t in truth if t["strand"] == "+"} | {("-", t["start"]) for t in truth if t["strand"] == "-"}
    orfs["planted"] = [(st, e if st == "+" else s) in stops
                       for st, s, e in zip(orfs["strand"], orfs["start"], orfs["end"])]

    rows = []
    for t in [30, 50, 75, 100, 150]:
        hits = orfs[orfs["length_aa"] >= t]
        rows.append({"min_aa": t, "orfs_called": len(hits), "planted_found": int(hits["planted"].sum()),
                     "spurious": int((~hits["planted"]).sum())})
    table = pd.DataFrame(rows)
    gc_rows = [{"gc": gc, "longest_spurious_orf_aa": longest_spurious(gc, rng)} for gc in (0.35, 0.5, 0.65, 0.72)]
    table.to_csv("results/summary.csv", index=False)
    pd.DataFrame(gc_rows).to_csv("results/gc_vs_spurious.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    ax = axes[0]
    bins = np.logspace(0, np.log10(orfs["length_aa"].max() + 1), 40)
    ax.hist(orfs.loc[~orfs["planted"], "length_aa"], bins=bins, alpha=0.7, label="spurious ORFs")
    ax.hist(orfs.loc[orfs["planted"], "length_aa"], bins=bins, alpha=0.9, label="planted genes")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("ORF length (aa)")
    ax.set_ylabel("count")
    ax.set_title("Random DNA is full of short ORFs")
    ax.legend()

    ax = axes[1]
    long = orfs[orfs["length_aa"] >= 50]
    for _, o in long.iterrows():
        ax.plot([o["start"], o["end"]], [o["frame"], o["frame"]], lw=6 if o["planted"] else 2,
                color="#2a7f62" if o["planted"] else "#bbbbbb", solid_capstyle="butt")
    ax.set_yticks([-3, -2, -1, 1, 2, 3])
    ax.set_xlabel("genome position (bp)")
    ax.set_ylabel("reading frame")
    ax.set_title("ORFs of at least 50 aa (green = planted)")

    ax = axes[2]
    g = pd.DataFrame(gc_rows)
    ax.bar([f"{int(x * 100)}%" for x in g["gc"]], g["longest_spurious_orf_aa"], color="#c0504d")
    ax.set_xlabel("GC content of random DNA")
    ax.set_ylabel("mean longest spurious ORF (aa)")
    ax.set_title("Stop codons are AT-rich, so GC hides them")
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))
    print(pd.DataFrame(gc_rows).to_string(index=False))


if __name__ == "__main__":
    main()
