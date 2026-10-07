import sys, gzip, collections, yaml
S = sys.argv[1]
mind = yaml.safe_load(open("config/settings.yaml"))["thresholds"]["min_depth"]
d = f"results/{S}/qc/{S}.cov"
reg = [l.rstrip("\n").split("\t") for l in gzip.open(d + ".regions.bed.gz", "rt")]
thr_lines = [l.rstrip("\n").split("\t") for l in gzip.open(d + ".thresholds.bed.gz", "rt")]
hdr = thr_lines[0]; thr = thr_lines[1:]
col = {h: i for i, h in enumerate(hdr)}          # e.g. '100X' -> column index
TS = ["100X", "250X", "500X", "1000X"]
for t in TS:
    assert t in col, f"threshold {t} missing from mosdepth output; rerun with --thresholds 100,250,500,1000"
G = collections.defaultdict(lambda: {"n": 0, "dm": 0.0, "t": {k: 0 for k in TS}, "reg": 0, "low": 0})
with open(f"results/{S}/report/{S}.region_coverage.tsv", "w") as fo:
    fo.write("gene\tchrom\tstart\tend\tmean_depth\tbelow_min_depth\n")
    for r, t in zip(reg, thr):
        n = int(r[2]) - int(r[1]); m = float(r[4]); g = r[3]; x = G[g]
        x["n"] += n; x["dm"] += m * n; x["reg"] += 1; x["low"] += (m < mind)
        for k in TS: x["t"][k] += int(t[col[k]])
        fo.write(f"{g}\t{r[0]}\t{r[1]}\t{r[2]}\t{m:.0f}\t{'YES' if m < mind else ''}\n")
with open(f"results/{S}/report/{S}.gene_coverage.tsv", "w") as fo:
    fo.write("gene\tmean_depth\tpct_100x\tpct_250x\tpct_500x\tpct_1000x\tregions\tregions_below_min_depth\n")
    for g, x in sorted(G.items()):
        p = [f"{100*x['t'][k]/x['n']:.1f}" for k in TS]
        fo.write(f"{g}\t{x['dm']/x['n']:.0f}\t{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\t{x['reg']}\t{x['low']}\n")
