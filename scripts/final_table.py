import sys, gzip, re, csv, yaml
from collections import Counter
vcf, bed, settings, mane_f, out = sys.argv[1:6]
cfg = yaml.safe_load(open(settings))["thresholds"]
LCV, LCA = cfg.get("low_conf_vaf", 0.10), cfg.get("low_conf_alt_reads", 10)

panel = set()
for l in open(bed):
    for g in l.split("\t")[3].strip().split(","):
        panel.add(g.replace("_promoter", ""))
panel.discard("GAPDH")

mane = {}
for l in gzip.open(mane_f, "rt"):
    if l.startswith("#"): continue
    c = l.rstrip("\n").split("\t")
    if c[9] == "MANE Select": mane[c[3]] = c[5].split(".")[0]

aa = dict(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",His="H",Ile="I",
          Leu="L",Lys="K",Met="M",Phe="F",Pro="P",Ser="S",Thr="T",Trp="W",Tyr="Y",Val="V",Ter="*")
def p1(p):
    m = re.match(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*)$", p)
    return f"p.({aa.get(m[1],m[1])}{m[2]}{aa.get(m[3],m[3])})" if m else p

REP = ("missense","stop_gained","stop_lost","start_lost","frameshift","splice_donor","splice_acceptor","inframe")
rows = []
for line in open(vcf):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
    info = dict(x.split("=", 1) if "=" in x else (x, "1") for x in f[7].split(";"))
    fmt = dict(zip(f[8].split(":"), f[9].split(":")))
    altn = int(fmt["AD"].split(",")[1]); dp = int(fmt["DP"]); vaf = float(fmt["AF"])
    ents = [e.split("|") for e in info.get("ANN", "").split(",") if e]
    pe = [e for e in ents if e[3] in panel]
    if not pe: continue
    gene = pe[0][3]
    same = [e for e in pe if e[3] == gene]
    pref = [e for e in same if e[6].split(".")[0] == mane.get(gene)]
    e, tflag = (pref[0], "MANE") if pref else (same[0], "NOT_MANE")
    cons = e[1]
    popaf = float(info["POPAF"].split(",")[0]) if "POPAF" in info else None
    cs = info.get("CLNSIG", "-")
    benign = cs.lower().startswith(("benign", "likely_benign"))
    common = popaf is not None and popaf < 2.0
    indel = len(ref) != len(alt)
    lod = cfg["lod_indel_vaf"] if indel else cfg["lod_snv_vaf"]
    flags = []
    if vaf < LCV or altn < LCA: flags.append("low_confidence")
    if (ref, alt) in (("C","T"),("G","A")) and vaf < LCV: flags.append("possible_FFPE_deamination")
    if abs(len(alt) - len(ref)) >= 20: flags.append("large_indel_check_IGV")
    if tflag == "NOT_MANE": flags.append("no_MANE_transcript_match")
    if dp < cfg["min_depth"] or altn < cfg["min_alt_reads"] or vaf < lod: cat = "1_excluded_low_support"
    elif not any(k in cons for k in REP): cat = "2_not_reported_consequence"
    elif common or benign: cat = "3_likely_benign_or_common"
    else: cat = "0_REVIEW"
    loc = (("Intron " if "intron" in cons else "Exon ") + e[8]) if e[8] else ""
    vcv = f"VCV{int(vid):09d}" if vid.isdigit() else ""
    g = f"{chrom}:g.{pos}{ref}>{alt}" if not indel else f"{chrom}:{pos} {ref[:10]}>{alt[:10]} (indel, HGVS g. manual)"
    rows.append([cat, gene, e[6], tflag, (e[9] if len(e[9]) <= 40 else e[9][:37] + "...(long)"),
                 p1(e[10]) if e[10] else "", cons, e[2], loc, g, chrom, pos, ref, alt, dp, altn, f"{vaf:.3f}",
                 popaf if popaf is not None else "NA", cs, vcv, info.get("ALLELEID", ""), ";".join(flags)])

rows.sort(key=lambda r: (r[0], -float(r[16])))
hdr = ["category","gene","transcript","transcript_source","HGVS_c","HGVS_p","consequence","impact","location",
       "genomic","chrom","pos","ref","alt","depth","alt_reads","VAF","POPAF","ClinVar","ClinVar_VCV","ClinVar_AlleleID","flags"]
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr); w.writerows(rows)
print(Counter(r[0] for r in rows))
