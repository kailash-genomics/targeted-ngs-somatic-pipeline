#!/usr/bin/env python3
"""Fill the Clinical Somatic NGS worksheet for one sample.
Usage (from the project root):  python scripts/fill_report.py SAMPLE_ID
Yellow cells = ENTER MANUALLY, or SUGGESTED (pipeline suggestion a person must confirm)."""
import sys, os, re, gzip, json, csv, copy, subprocess, datetime, hashlib
from math import comb, log10
from collections import defaultdict
import yaml, docx
from docx.enum.text import WD_COLOR_INDEX

S = sys.argv[1]
os.chdir(os.environ.get("PIPE_ROOT", "."))
CFG = yaml.safe_load(open("config/settings.yaml"))
TH = CFG["thresholds"]
R = f"results/{S}"
TODAY = datetime.date.today().isoformat()
NA = "Not available"

class Manual(str): pass     # human must enter
class Sugg(str): pass       # pipeline suggestion, human confirms

# ---------------------------------------------------------------- helpers
def rd(p): return open(p).read() if os.path.exists(p) else ""
def jload(p): return json.load(open(p)) if os.path.exists(p) else {}
def num(n): return f"{int(n):,}"
def f1(x, nd=1): return NA if x is None else f"{x:.{nd}f}"
def run(cmd):
    try: return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=180)
    except Exception: return None
def ver(cmd, rx):
    r = run(cmd)
    if not r: return None
    m = re.search(rx, (r.stdout or "") + (r.stderr or ""))
    return m.group(0) if m else None
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""): h.update(b)
    return h.hexdigest()
def kv(p):
    d = {}
    for l in rd(p).splitlines():
        if "=" in l: k, v = l.split("=", 1); d[k] = v
    return d
def tsv(p):
    if not os.path.exists(p): return []
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

# ---------------------------------------------------------------- inputs
meta = {}
for r in csv.DictReader(open("config/metadata.csv")):
    if r["sample_id"] == S: meta = {k: (v or "").strip() for k, v in r.items()}
def M(k, hint=None):
    v = meta.get(k, "")
    return v if v else Manual(hint or "")

samples = {r["sample_id"]: r for r in csv.DictReader(open("config/samples.csv"))}
FQ1, FQ2 = samples[S]["r1"], samples[S]["r2"]

fj = jload(f"{R}/qc/{S}.fastp.json")
def rs(which, side):
    d = fj.get(f"read{side}_{which}_filtering", {})
    tb, tr, q30 = d.get("total_bases", 0), d.get("total_reads", 0), d.get("q30_bases", 0)
    qc = d.get("quality_curves", {}).get("mean", []); cc = d.get("content_curves", {})
    gc = 100 * (sum(cc["G"]) / len(cc["G"]) + sum(cc["C"]) / len(cc["C"])) if cc.get("G") and cc.get("C") else None
    nn = 100 * sum(cc["N"]) / len(cc["N"]) if cc.get("N") else None
    return dict(reads=tr, bases=tb, q30=100 * q30 / tb if tb else None,
                mq=sum(qc) / len(qc) if qc else None, gc=gc, n=nn, cycles=d.get("total_cycles"))
B1, B2, A1, A2 = rs("before", 1), rs("before", 2), rs("after", 1), rs("after", 2)
raw_reads = B1["reads"] + B2["reads"]; trim_reads = A1["reads"] + A2["reads"]

fs = {}
for l in rd(f"{R}/qc/{S}.flagstat.txt").splitlines():
    m = re.match(r"(\d+) \+ (\d+) (.*)", l)
    if m:
        nm = re.sub(r"\s*\(.*\)$", "", m.group(3)).strip()
        pc = re.search(r"\(([\d.]+)%", l)
        fs[nm] = (int(m.group(1)), float(pc.group(1)) if pc else None)
prim = fs.get("primary", (0,))[0]; pmap = fs.get("primary mapped", (0,))[0]

def pic(p, first):
    L = rd(p).splitlines()
    for i, l in enumerate(L):
        if l.startswith(first): return dict(zip(l.split("\t"), L[i + 1].split("\t")))
    return {}
HS = pic(f"{R}/qc/{S}.hsmetrics.txt", "BAIT_SET")
DUP = pic(f"{R}/qc/{S}.dupmetrics.txt", "LIBRARY")
def hs(k):
    v = HS.get(k)
    try: return float(v)
    except (TypeError, ValueError): return None

BS = kv(f"{R}/qc/{S}.bamstats.txt")
stx = rd(f"{R}/qc/{S}.samtools_stats.txt").splitlines()
def sn(key):
    for l in stx:
        if l.startswith("SN\t" + key): return float(l.split("\t")[2])
isz = [(int(l.split("\t")[1]), int(l.split("\t")[2])) for l in stx if l.startswith("IS\t")]
def med_is():
    tot = sum(c for _, c in isz); acc = 0
    for s, c in isz:
        acc += c
        if tot and acc >= tot / 2: return s

GC_ = {r["gene"]: r for r in tsv(f"{R}/report/{S}.gene_coverage.tsv")}
RC_ = tsv(f"{R}/report/{S}.region_coverage.tsv")
VT = tsv(f"{R}/report/{S}.variants.tsv")
REVIEW = [v for v in VT if v["category"] == "0_REVIEW"]
DETAIL = [v for v in REVIEW if not ({"low_confidence", "large_indel_check_IGV"} & set(v["flags"].split(";")))]

def parse_vcf(p):
    out = {}
    if not os.path.exists(p): return out
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        for l in fh:
            if l.startswith("#"): continue
            f = l.rstrip("\n").split("\t")
            info = dict(x.split("=", 1) if "=" in x else (x, "1") for x in f[7].split(";"))
            out[(f[0], f[1], f[3], f[4])] = dict(id=f[2], qual=f[5], filt=f[6], info=info, fmt=dict(zip(f[8].split(":"), f[9].split(":"))))
    return out
ANN = parse_vcf(f"{R}/variants/{S}.annotated.vcf")
def count_vcf(p):
    n = pas = 0
    if os.path.exists(p):
        with gzip.open(p, "rt") as fh:
            for l in fh:
                if l.startswith("#"): continue
                n += 1; pas += l.split("\t")[6] == "PASS"
    return n, pas
N_CAND, _ = count_vcf(f"{R}/variants/{S}.mutect2.vcf.gz")
_, N_PASS = count_vcf(f"{R}/variants/{S}.filtered.vcf.gz")

# ---------------------------------------------------------------- CIViC (best effort, exact matches only)
CV = tsv(CFG.get("civic_variants", "")); CE = tsv(CFG.get("civic_evidence", ""))
civ_coord = {(r.get("chromosome"), r.get("start"), r.get("reference_bases"), r.get("variant_bases")): r for r in CV}
civ_name = {((r.get("gene") or "").upper(), (r.get("variant") or "").upper()): r for r in CV}
def civic(gene, pshort, chrom, pos, ref, alt):
    hit = civ_coord.get((chrom.replace("chr", ""), pos, ref, alt)) or civ_name.get((gene.upper(), pshort.upper()))
    if not hit: return None, []
    mpid = hit.get("molecular_profile_id") or hit.get("variant_id")
    prof = f"{gene} {hit.get('variant', '')}".upper()
    ev = [e for e in CE if (mpid and e.get("molecular_profile_id") == mpid) or (e.get("molecular_profile") or "").upper() == prof]
    return hit, ev
def ev_line(e):
    return " | ".join(x for x in [e.get("evidence_level"), e.get("disease"), e.get("therapies"), e.get("clinical_significance")] if x)
def civ_summary(hit, ev):
    if hit is None: return dict(all="No CIViC record matched (exact coordinate or gene + protein change only)", ther=None, dxpx=None, res=None, dis=None)
    if not ev: return dict(all=f"CIViC variant record found (variant_id {hit.get('variant_id', '?')}); no evidence items matched", ther=None, dxpx=None, res=None, dis=None)
    def sel(f): return "; ".join(ev_line(e) for e in ev if f(e))[:600] or None
    t = lambda e: (e.get("evidence_type") or "").lower()
    cnt = defaultdict(int)
    for e in ev: cnt[t(e) or "?"] += 1
    return dict(all=f"CIViC: {len(ev)} evidence items (" + ", ".join(f"{k} {v}" for k, v in cnt.items()) + ")",
                ther=sel(lambda e: t(e) == "predictive"), dxpx=sel(lambda e: t(e) in ("diagnostic", "prognostic")),
                res=sel(lambda e: t(e) == "predictive" and "RESISTANCE" in (e.get("clinical_significance") or "").upper()),
                dis=", ".join(sorted({e.get("disease") for e in ev if e.get("disease")}))[:300] or None)

# ---------------------------------------------------------------- small stats
def fisher_p(a, b, c, d):
    n = a + b + c + d; r1 = a + b; c1 = a + c
    if n == 0: return 1.0
    pm = lambda x: comb(r1, x) * comb(n - r1, c1 - x) / comb(n, c1)
    p0 = pm(a)
    return min(1.0, sum(pm(x) for x in range(max(0, c1 - (n - r1)), min(r1, c1) + 1) if pm(x) <= p0 * (1 + 1e-7)))
def gnomad_af(popaf):
    if popaf in ("NA", None, ""): return NA
    p = float(popaf)
    return "No entry in gnomAD AF-only resource (Mutect2 default AF 5e-8 applied)" if p >= 7.29 else f"{10 ** -p:.2e}"
STARS = [("practice_guideline", 4), ("reviewed_by_expert_panel", 3), ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
         ("criteria_provided,_single_submitter", 1), ("criteria_provided,_conflicting", 1)]
def stars(rs_):
    for k, v in STARS:
        if rs_.startswith(k): return v
    return 0
def clean(s): return s.replace("_", " ") if s else s
def vtype(ref, alt):
    if len(ref) == len(alt): return "SNV" if len(ref) == 1 else "MNV"
    return "Insertion" if len(alt) > len(ref) else "Deletion"
def short(s, n=24): return s if len(s) <= n else s[:n] + "…"

# ---------------------------------------------------------------- versions / provenance
V = {
 "FastQC": ver("fastqc --version 2>&1", r"FastQC v[\d.]+"), "fastp": ver("fastp --version 2>&1", r"fastp [\d.]+"),
 "BWA": ver("bwa 2>&1", r"Version: [\w.\-]+"), "samtools": ver("samtools --version", r"samtools [\d.]+"),
 "GATK": ver("gatk --version 2>&1", r"\(GATK\) v[\d.]+"), "bcftools": ver("bcftools --version", r"bcftools [\d.]+"),
 "SnpEff": ver("conda run -n snpeff_env snpEff -version 2>&1", r"SnpEff\s+[\w.]+\s+[\d-]*")}
def fdate(p): return datetime.date.fromtimestamp(os.path.getmtime(p)).isoformat() if os.path.exists(p) else "?"
def clinvar_date():
    p = CFG["clinvar"]
    if os.path.exists(p):
        with gzip.open(p, "rt") as fh:
            for l in fh:
                if l.startswith("##fileDate"): return l.strip().split("=")[1]
                if not l.startswith("##"): break
    return "date not found in header"
snak = md5("Snakefile") if os.path.exists("Snakefile") else "?"
PIPE_VER = f"targeted NGS somatic pipeline v1.0 (Snakefile md5 {snak[:8]})"
def bed_ver():
    n = sum(1 for _ in open(CFG["bed"]))
    return f"{os.path.basename(CFG['bed'])} (md5 {md5(CFG['bed'])[:8]}, {n} regions)"

prov = dict(sample=S, generated=datetime.datetime.now().isoformat(timespec="seconds"), pipeline=PIPE_VER, versions=V,
            settings=CFG, snakefile_md5=snak, fastq_md5={os.path.basename(FQ1): md5(FQ1), os.path.basename(FQ2): md5(FQ2)},
            database_checksums=rd("db/CHECKSUMS.md5").splitlines(), clinvar_file_date=clinvar_date(),
            civic_download=fdate(CFG.get("civic_variants", "")), gnomad_resource_file_date=fdate(CFG["germline_resource"]),
            bed=bed_ver(), preferred_transcripts="MANE Select v1.5 (fallback: first SnpEff transcript, flagged)")
json.dump(prov, open(f"{R}/report/{S}.provenance.json", "w"), indent=2)

# ---------------------------------------------------------------- FASTQ header
hdr = ""
with gzip.open(FQ1, "rt") as fh: hdr = fh.readline().strip()
hp = hdr.lstrip("@").split(" ")[0].split(":"); ip = hdr.split(" ")[-1].split(":")[-1] if " " in hdr else ""
INST, RUNNO, FCELL, LANE = (hp + [""] * 4)[:4]

# ---------------------------------------------------------------- QC evaluation
mean_dep, pct100 = hs("MEAN_TARGET_COVERAGE"), (hs("PCT_TARGET_BASES_100X") or 0) * 100
dup_pct = float(DUP["PERCENT_DUPLICATION"]) * 100 if DUP.get("PERCENT_DUPLICATION") else None
q30min = min(x for x in (B1["q30"], B2["q30"]) if x is not None) if B1["q30"] is not None else None
on_t = int(BS["on_target_reads"]) if "on_target_reads" in BS else None
QC = []   # (metric, observed, criterion, ok, impact, action)
def add(m, obs, crit, ok, imp, act): QC.append((m, obs, crit, ok, imp, act))
if mean_dep is not None: add("Mean target depth (dedup)", f"{mean_dep:.0f}×", f"≥{TH['min_depth']}× (PLACEHOLDER)", mean_dep >= TH["min_depth"], "Sensitivity at low VAF reduced if below", "Compare with validated criterion")
add("% target bases ≥100×", f"{pct100:.1f}%", f"≥{TH['min_pct_100x']}% (PLACEHOLDER)", pct100 >= TH["min_pct_100x"], "Regions below depth cannot exclude variants", "See gene/exon coverage table")
if q30min is not None: add("Raw FASTQ % Q30 (lowest of R1/R2)", f"{q30min:.1f}%", f"≥{TH['min_q30_pct']}% (PLACEHOLDER)", q30min >= TH["min_q30_pct"], "Base-call accuracy", "—")
if dup_pct is not None: add("Duplicate rate (Picard)", f"{dup_pct:.1f}%", f"≤{TH['max_duplicate_pct']}% (PLACEHOLDER)", dup_pct <= TH["max_duplicate_pct"], "Reduces unique depth; FFPE/low input typical", "Review DNA input / library complexity")
if on_t is not None and pmap: add("On-target rate (reads overlapping BED)", f"{100 * on_t / pmap:.1f}%", Manual("validated criterion"), None, "Reads outside BED (e.g. intronic/fusion probes) lower exon depth", "Confirm capture design")
try:
    _c = rd(f"{R}/variants/{S}.contamination.table").splitlines()[1].split("\t")
    add("Cross-sample contamination (CalculateContamination)", f"{float(_c[1]) * 100:.2f}% (± {float(_c[2]) * 100:.2f}%)", Manual("validated criterion"), None,
        "Estimate is less reliable on a small panel; it feeds FilterMutectCalls", "Interpret with caution")
except Exception:
    pass
low_genes = [g for g, r in GC_.items() if float(r["pct_100x"]) < TH["min_pct_100x"]]
sev_genes = sorted([g for g, r in GC_.items() if float(r["mean_depth"]) < 30], key=lambda g: float(GC_[g]["mean_depth"]))
add(f"Genes with <{TH['min_pct_100x']}% bases at ≥{TH['min_depth']}×", f"{len(low_genes)} of {len(GC_)} genes", "PLACEHOLDER", len(low_genes) == 0,
    "‘No variant found’ is not reliable in these genes", "Listed: " + ", ".join(sorted(low_genes))[:500])
for g in sev_genes[:15]:
    r = GC_[g]; add(f"Gene {g}: very low coverage", f"mean {r['mean_depth']}×, {r['pct_100x']}% ≥100×", f"mean ≥{TH['min_depth']}×", False, "Gene effectively not assessed", "Report as inadequately covered")
add("CNV analysis", "Not performed", "Assay-validated CNV method", None, "Copy-number changes not assessed", "No normal reference pool available")
add("Gene fusions", "Not assessed", "Assay-validated fusion method", None, "Fusions not assessed by this pipeline", "Requires RNA / intronic design")
add("MSI / TMB", "Not performed", "Validated method", None, "Not reported", "—")
nflag = sum(1 for v in REVIEW if {"low_confidence", "large_indel_check_IGV"} & set(v["flags"].split(";")))
if nflag: add("Review variants with technical warnings", f"{nflag} of {len(REVIEW)}", "—", None, "Possible artefacts (low support / large indels)", "Inspect in IGV before use")
any_fail = any(q[3] is False for q in QC[:5])

# ---------------------------------------------------------------- docx helpers
d = docx.Document("config/template.docx")
T = d.tables
assert len(T) == 22 and T[8].rows[0].cells[0].text.strip() == "Gene" and T[18].rows[0].cells[0].text.strip() == "Gene", "template layout differs from expected"
def norm(s): return " ".join(s.split())
def put(cell, text, hl=False):
    p = cell.paragraphs[0]
    for e in cell.paragraphs[1:]: e._p.getparent().remove(e._p)
    run = p.runs[-1] if p.runs else p.add_run()
    for r in p.runs[:-1]: r.text = ""
    run.text = str(text)
    run.font.highlight_color = WD_COLOR_INDEX.YELLOW if hl else None
def fill(tbl, vals):
    for row in tbl.rows[1:]:
        label = norm(row.cells[0].text); vc = row.cells[1]; hint = norm(vc.text)
        v = vals.get(label, Manual(""))
        if v is None: v = Manual("")
        if isinstance(v, Manual): put(vc, "ENTER MANUALLY" + (f" — {v}" if v else (f" ({hint})" if hint else "")), True)
        elif isinstance(v, Sugg): put(vc, f"SUGGESTED: {v} — confirm", True)
        else: put(vc, v)
def rows_table(tbl, data):
    body = tbl.rows[1:]
    while len(tbl.rows) - 1 < len(data): tbl._tbl.append(copy.deepcopy(tbl.rows[-1]._tr))
    while len(tbl.rows) - 1 > max(len(data), 1): tbl._tbl.remove(tbl.rows[-1]._tr)
    for row, rec in zip(tbl.rows[1:], data):
        for cell, v in zip(row.cells, rec):
            if isinstance(v, Manual): put(cell, "ENTER MANUALLY" + (f" — {v}" if v else ""), True)
            elif isinstance(v, Sugg): put(cell, f"SUGGESTED: {v}", True)
            else: put(cell, v)
def status(ok, text=""):
    return Sugg({True: "Pass", False: "Review (below placeholder threshold)", None: "Review"}[ok] + text)

# ---------------------------------------------------------------- tables 0-7
fill(T[0], {"Patient / Sample ID": meta.get("patient_id") or S, "Case / Accession No.": M("case_no"), "Sample Type": M("sample_type", "FFPE / Blood / cfDNA / Other"),
  "Tumor Site / Specimen": M("tumor_site"), "Collection Date": M("collection_date"), "Extraction Date": M("extraction_date"), "Analysis Date": TODAY,
  "Panel": CFG.get("panel_name", "Targeted panel (set panel_name in config/settings.yaml)"), "Reference Genome": "hg19 / GRCh37", "Panel BED / Version": bed_ver(), "Pipeline Version": PIPE_VER,
  "Analyst": M("analyst"), "Reviewer": M("reviewer")})
fill(T[1], {"Sample ID / FASTQ ID": f"{S} / {os.path.basename(FQ1)}", "Sample Type": M("sample_type"), "Tumor Percentage / Estimated Tumor Purity": M("tumor_pct"),
  "Specimen Adequacy": M("specimen_adequacy", "Pass / Fail / N/A"), "DNA Input": M("dna_input_ng"), "DNA Concentration": M("dna_conc_ng_ul"), "DNA Volume": M("dna_volume_ul"),
  "DNA QC Method": M("dna_qc_method", "Qubit / Other"), "DIN / DNA Integrity": M("din"), "Library Input": M("library_input_ng"), "Library Preparation Kit / Lot": M("library_kit_lot"),
  "UMI Used?": M("umi_used", "Yes / No"), "Enrichment Method": M("enrichment_method")})
fill(T[2], {"Sequencer": meta.get("sequencer") or f"Instrument ID {INST} (model not stated in FASTQ header)", "Run ID": meta.get("run_id") or f"Run number {RUNNO} (from FASTQ header)",
  "Flow Cell": meta.get("flow_cell") or f"{FCELL} (lane {LANE})", "Read Configuration": meta.get("read_config") or f"{B1['cycles']} / {B2['cycles']} cycles (R1/R2, from FASTQ)",
  "Index / UMI Configuration": meta.get("index_config") or f"Index {ip} (from FASTQ header); UMI not detected by pipeline",
  "Total Clusters / Reads": f"{num(raw_reads)} (R1+R2 in FASTQ; clusters = {num(B1['reads'])})", "PF Reads": f"{num(raw_reads)} (FASTQ contains PF reads only)",
  "% PF": Manual("needs run-level report"), "Q30 Read 1": f1(B1["q30"]), "Q30 Read 2": f1(B2["q30"]), "Mean / Median Q Score": f"Mean {f1((B1['mq'] + B2['mq']) / 2 if B1['mq'] else None)} (fastp)",
  "Total Yield": f"{(B1['bases'] + B2['bases']) / 1e9:.2f} Gb", "Demultiplexing Status": M("demux_status", "Pass / Fail"), "Undetermined Reads": M("undetermined_reads", "run-level report")})
q3ok = q30min is not None and q30min >= TH["min_q30_pct"]
fill(T[3], {"Raw R1 Reads": num(B1["reads"]), "Raw R2 Reads": num(B2["reads"]), "Raw Total Bases": num(B1["bases"] + B2["bases"]), "R1 Mean Read Quality": f1(B1["mq"]), "R2 Mean Read Quality": f1(B2["mq"]),
  "R1 % Q30": f1(B1["q30"]), "R2 % Q30": f1(B2["q30"]), "R1 GC Content": f1(B1["gc"]), "R2 GC Content": f1(B2["gc"]),
  "Adapter Content": f"{100 * fj.get('adapter_cutting', {}).get('adapter_trimmed_reads', 0) / raw_reads:.2f} (reads with adapter trimmed)" if raw_reads else NA,
  "N Content": f1(B1["n"], 3), "Duplication Level": f1(100 * fj.get("duplication", {}).get("rate", 0)) + " (fastp estimate)", "Raw FASTQ QC Status": status(q3ok)})
cmd = fj.get("command", ""); mlen = re.search(r"--length_required\s+(\d+)", cmd)
fill(T[4], {"Tool / Version": f"fastp {fj.get('summary', {}).get('fastp_version', '')}".strip(), "Trimmed R1 Reads": num(A1["reads"]), "Trimmed R2 Reads": num(A2["reads"]),
  "Reads Retained": f"{100 * trim_reads / raw_reads:.2f}", "Reads Removed": f"{100 * (1 - trim_reads / raw_reads):.2f}",
  "Adapter Removal": f"{fj.get('adapter_cutting', {}).get('adapter_trimmed_reads', 0):,} reads", "Low-Quality Bases Removed": f"{num(B1['bases'] + B2['bases'] - A1['bases'] - A2['bases'])} bp (all trimming/filtering)",
  "Minimum Read Length": mlen.group(1) if mlen else "50 (pipeline setting)", "Post-trim R1 Q30": f1(A1["q30"]), "Post-trim R2 Q30": f1(A2["q30"]),
  "Post-trim GC": f1((A1["gc"] + A2["gc"]) / 2 if A1["gc"] else None), "Post-trim QC Status": status(A1["q30"] is not None and min(A1["q30"], A2["q30"]) >= TH["min_q30_pct"])})
uniq_pct = None
if RC_:
    tot = sum(int(r["end"]) - int(r["start"]) for r in RC_); mn = sum(float(r["mean_depth"]) * (int(r["end"]) - int(r["start"])) for r in RC_) / tot
    uniq_pct = 100 * sum(int(r["end"]) - int(r["start"]) for r in RC_ if float(r["mean_depth"]) >= 0.2 * mn) / tot
ins_mean = sn("insert size average:")
fill(T[5], {"Aligner / Version": f"BWA-MEM ({V['BWA'] or 'version n/a'})", "Reference Genome / Version": "hg19 / GRCh37 (" + CFG["reference"] + ")", "Total Reads Input": num(prim) + " (primary)",
  "Mapped Reads": num(pmap) + " (primary)", "Mapping Rate": f"{100 * pmap / prim:.2f}" if prim else NA, "Properly Paired Reads": f1(fs.get("properly paired", (0, None))[1], 2),
  "Unmapped Reads": f"{100 * (1 - pmap / prim):.2f}" if prim else NA, "Primary Alignments": num(prim), "Secondary Alignments": num(fs.get("secondary", (0,))[0]),
  "Supplementary Alignments": num(fs.get("supplementary", (0,))[0]), "Duplicate Reads": num(fs.get("primary duplicates", (0,))[0]),
  "Duplicate Rate": f1(dup_pct) + " (Picard, pair-based)", "Mean MAPQ": BS.get("mapq_mean", NA),
  "MAPQ ≥20": f"{100 * int(BS['mapq20']) / int(BS['n_primary_mapped']):.1f}" if BS else NA, "MAPQ ≥30": f"{100 * int(BS['mapq30']) / int(BS['n_primary_mapped']):.1f}" if BS else NA,
  "Insert Size Mean": f1(ins_mean), "Insert Size Median": str(med_is() or NA), "BAM File Size": f"{os.path.getsize(f'{R}/align/{S}.dedup.bam') / 1e9:.2f}",
  "BAM Index Present": "Yes" if os.path.exists(f"{R}/align/{S}.dedup.bam.bai") else "No", "BAM QC Status": status(None)})
fill(T[6], {"Panel BED Version": bed_ver(), "Number of Target Regions": str(sum(1 for _ in open(CFG["bed"]))), "Target Bases": num(hs("TARGET_TERRITORY") or 0),
  "Total Mapped Reads": num(pmap), "On-target Reads": f"{num(on_t)} (reads overlapping BED)" if on_t is not None else NA,
  "On-target Rate": f"{100 * on_t / pmap:.1f} (reads); Picard base-level {100 * (hs('PCT_SELECTED_BASES') or 0):.1f}" if on_t is not None else NA,
  "Off-target Reads": num(pmap - on_t) if on_t is not None else NA, "Off-target Rate": f"{100 * (1 - on_t / pmap):.1f}" if on_t is not None else NA,
  "Mean Target Depth": f1(mean_dep, 0) + " (duplicates excluded)", "Median Target Depth": f1(hs("MEDIAN_TARGET_COVERAGE"), 0), "Minimum Target Depth": f1(hs("MIN_TARGET_COVERAGE"), 0),
  "Maximum Target Depth": f1(hs("MAX_TARGET_COVERAGE"), 0), "Uniformity": (f"{uniq_pct:.1f} (% of target bases in regions ≥0.2× mean; region-level approximation)" if uniq_pct is not None else NA),
  "Fold-80 Base Penalty": HS.get("FOLD_80_BASE_PENALTY") if HS.get("FOLD_80_BASE_PENALTY") not in (None, "?", "") else "Not defined (Picard returned ‘?’)", "Target QC Status": status(None)})
cv = {f"% Target Bases ≥{k}×": (f1(hs(f"PCT_TARGET_BASES_{k}X") * 100) if hs(f"PCT_TARGET_BASES_{k}X") is not None else "Not in Picard output") for k in (1, 10, 20, 30, 50, 100, 250, 500, 1000)}
mind = TH["min_depth"]
below = (100 - (hs(f"PCT_TARGET_BASES_{mind}X") or 0) * 100) if hs(f"PCT_TARGET_BASES_{mind}X") is not None else None
tb = hs("TARGET_TERRITORY") or 0
cv["Target Bases Below Minimum Depth"] = f"{below:.1f}% ({num(below * tb / 100)} bp) below {mind}× (PLACEHOLDER minimum)" if below is not None else Manual(f"minimum depth {mind}× has no Picard column")
nocov = 100 - (hs("PCT_TARGET_BASES_1X") or 0) * 100
cv["No-Coverage Target Bases"] = f"{nocov:.2f}% ({num(nocov * tb / 100)} bp)"
cv["Overall Coverage QC"] = status(pct100 >= TH["min_pct_100x"])
fill(T[7], cv)

# ---------------------------------------------------------------- table 8: genes
low_reg = {}
for r in RC_:
    m = float(r["mean_depth"])
    if r["gene"] not in low_reg or m < low_reg[r["gene"]][0]: low_reg[r["gene"]] = (m, f"{r['chrom']}:{r['start']}-{r['end']} ({m:.0f}×)")
gdata = []
for g in sorted(GC_):
    r = GC_[g]; ok = float(r["pct_100x"]) >= TH["min_pct_100x"]
    gdata.append([g, r["mean_depth"], r["pct_100x"], r["pct_250x"], low_reg.get(g, (0, NA))[1], "PASS" if ok else f"FAIL (<{TH['min_pct_100x']}% ≥{mind}×, placeholder)"])
rows_table(T[8], gdata)

# ---------------------------------------------------------------- table 9, 13
fill(T[9], {"Caller / Version": f"GATK Mutect2 ({V['GATK'] or 'version n/a'})", "Tumor BAM": f"{R}/align/{S}.dedup.bam", "Normal BAM": "Not used (tumor-only)",
  "Germline Resource": f"{os.path.basename(CFG['germline_resource'])} (gnomAD AF-only b37, renamed to hg19)", "Panel of Normals (PoN)": "None used (no normal pool available)",
  "Intervals / BED": os.path.basename(CFG["bed"]), "Reference Genome": "hg19 / GRCh37", "Mutect2 Candidate Variants": num(N_CAND), "PASS Variants": num(N_PASS),
  "Filtered Variants": num(N_CAND - N_PASS), "F1R2 Orientation Model Used": "Yes", "FilterMutectCalls Used": "Yes", "Variant Calling QC Status": status(None)})
cnt = lambda c: sum(1 for v in VT if v["category"] == c)
lod_i, lod_s = TH["lod_indel_vaf"], TH["lod_snv_vaf"]
below_lod = sum(1 for v in VT if float(v["VAF"]) < (lod_i if len(v["ref"]) != len(v["alt"]) else lod_s))
fill(T[13], {"Initial Candidate Variants": num(N_CAND), "PASS After Caller Filtering": num(N_PASS), "Minimum DP Criterion": f"DP ≥{mind} and ALT reads ≥{TH['min_alt_reads']} (PLACEHOLDER)",
  "Minimum VAF Criterion": f"SNV ≥{lod_s * 100:.0f}%, indel ≥{lod_i * 100:.0f}%", "LOD of Assay": "Configurable threshold; validate for the intended assay",
  "Below-LOD Variants": f"{below_lod} (among PASS calls in panel genes)", "Population AF Filter": f"POPAF<2 (gnomAD AF>1%) or ClinVar Benign/Likely benign: {cnt('3_likely_benign_or_common')} variants set aside",
  "Germline Likelihood Filter": "FilterMutectCalls germline filter; tumor-only cannot separate germline from somatic — see VAF context",
  "Artifact Filter": f"FilterMutectCalls (orientation, strand bias, haplotype, clustered events, etc.); {nflag} review variants carry technical warnings",
  "Panel / Region Filter": "Panel BED genes only (GAPDH control excluded)", "Clinical Relevance Filter": "Missense, nonsense, frameshift, splice donor/acceptor, in-frame only",
  "Final Variants for Interpretation": str(len(REVIEW))})

# ---------------------------------------------------------------- per-variant blocks (10, 11, 12, 14, 15)
tpl_ids = [10, 11, 12, 14, 15]
purity = None
try: purity = float(meta.get("tumor_pct", "").replace("%", ""))
except ValueError: pass

def rec_for(v): return ANN.get((v["chrom"], v["pos"], v["ref"], v["alt"]), dict(id="", qual=".", filt="", info={}, fmt={}))
def block10(v):
    r = rec_for(v); f, i = r["fmt"], r["info"]
    sb = [int(x) for x in f.get("SB", "0,0,0,0").split(",")]
    fs_ = -10 * log10(max(fisher_p(sb[0], sb[1], sb[2], sb[3]), 1e-300)) if sum(sb) else None
    mbq = i.get("MBQ", "").split(","); mmq = i.get("MMQ", "").split(",")
    f1r2 = f.get("F1R2", "").split(","); f2r1 = f.get("F2R1", "").split(",")
    flags = v["flags"].replace(";", ", ")
    ok = r["filt"] == "PASS" and not flags
    return {"Variant ID": f"{v['gene']} {v['chrom']}:{v['pos']} {short(v['ref'])}>{short(v['alt'])}", "CHROM": v["chrom"], "POS": v["pos"], "REF": short(v["ref"], 30), "ALT": short(v["alt"], 30),
        "FILTER": r["filt"] or NA, "DP": v["depth"], "AD REF": f.get("AD", "0,0").split(",")[0], "AD ALT": v["alt_reads"], "VAF / AF": f"{100 * float(v['VAF']):.1f}",
        "QUAL": "Not reported by Mutect2 (‘.’); TLOD " + i.get("TLOD", NA), "Base Quality of ALT Reads": f"median {mbq[1]}" if len(mbq) > 1 else NA,
        "MAPQ of ALT Reads": f"median {mmq[1]}" if len(mmq) > 1 else NA, "Forward ALT Reads": str(sb[2]), "Reverse ALT Reads": str(sb[3]),
        "Strand Bias / FS": f"FS {fs_:.2f} (Phred-scaled Fisher, computed from SB counts {','.join(map(str, sb))})" if fs_ is not None else NA,
        "Read Position Bias": f"median distance from read end {i['MPOS']} bp (MPOS)" if "MPOS" in i else NA,
        "F1R2 / Orientation Bias": f"ALT F1R2 {f1r2[1] if len(f1r2) > 1 else '?'} / ALT F2R1 {f2r1[1] if len(f2r1) > 1 else '?'}; orientation filter: {'not triggered' if 'orientation' not in r['filt'] else 'TRIGGERED'}",
        "Read End Bias": Manual("not computed; see MPOS above or inspect in IGV"), "Artifact Flag / Evidence": f"FILTER={r['filt']}" + (f"; pipeline flags: {flags}" if flags else "; no pipeline flags"),
        "Technical Review": Sugg("Pass" if ok else "Review")}
def block11(v):
    r = rec_for(v); nm, _, ver_ = v["transcript"].partition(".")
    cons = v["consequence"]
    spl = "Splice donor/acceptor site" if ("splice_donor" in cons or "splice_acceptor" in cons) else ("Splice region" if "splice_region" in cons else "None")
    return {"CHROM:POS": f"{v['chrom']}:{v['pos']}", "Gene": v["gene"], "Transcript ID": nm, "Transcript Version": ver_ or NA,
        "Transcript Selection Rule": "MANE Select (v1.5)" if v["transcript_source"] == "MANE" else "NOT MANE — first SnpEff transcript (verify)", "Exon": v["location"] or NA,
        "cDNA HGVS (c.)": v["HGVS_c"], "Protein HGVS (p.)": v["HGVS_p"] or "—", "HGVS Genomic (g.)": v["genomic"], "Consequence": cons, "Codon Change": "Not extracted (SnpEff ANN gives HGVS c. only)",
        "Amino Acid Change": v["HGVS_p"] or "—", "Splice Region / Site": spl, "SnpEff / VEP Effect": f"SnpEff: {cons}", "SnpEff / VEP Impact": f"SnpEff: {v['impact']}"}
def civ_for(v):
    m = re.match(r"p\.\(([A-Z*])(\d+)([A-Z*])\)", v["HGVS_p"] or "")
    ps = f"{m[1]}{m[2]}{m[3]}" if m else ""
    hit, ev = civic(v["gene"], ps, v["chrom"], v["pos"], v["ref"], v["alt"])
    return civ_summary(hit, ev)
def block12(v):
    r = rec_for(v); rsx = r["info"].get("CLNREVSTAT", ""); cs = v["ClinVar"]
    c = civ_for(v)
    return {"ClinVar Variation ID": (r["id"] if r["id"].isdigit() else "Not in ClinVar") + (f" ({v['ClinVar_VCV']}; version suffix not in VCF)" if v["ClinVar_VCV"] else ""),
        "ClinVar Allele ID": v["ClinVar_AlleleID"] or ("Not in ClinVar" if not r["id"].isdigit() else "Not extracted"),
        "ClinVar Clinical Significance": clean(cs) if cs != "-" else "Not in ClinVar", "ClinVar Review Status / Stars": f"{clean(rsx)} ({stars(rsx)}★)" if rsx else "—",
        "ClinVar Conflicting?": "Yes" if "onflicting" in cs else "No", "gnomAD AF": gnomad_af(v["POPAF"]),
        "gnomAD Popmax AF": "Not available (AF-only resource; panel subset not downloaded)", "gnomAD Allele Count / Number": "Not available (AF-only resource)",
        "COSMIC ID": "Not queried (license required)", "COSMIC Occurrence / Cancer Type": "Not queried (license required)", "OncoKB Evidence / Level": "Not queried (license required)",
        "CIViC Evidence": c["all"], "Other Database IDs": "dbSNP rs: not queried (dbSNP not downloaded)"}
def block14(v):
    c = civ_for(v); cs = v["ClinVar"]
    return {"Gene / Variant": f"{v['gene']} {v['HGVS_c']} {v['HGVS_p']}".strip(), "Variant Type": vtype(v["ref"], v["alt"]),
        "Oncogenicity / Pathogenicity Assessment": Manual(f"ClinVar: {clean(cs)}" if cs != "-" else "not in ClinVar"), "Evidence Sources": f"ClinVar ({clinvar_date()}); CIViC (nightly {fdate(CFG.get('civic_variants', ''))}); MANE v1.5",
        "Cancer Type Relevance": c["dis"] or Manual("CIViC diseases: none matched"), "Therapeutic Evidence": c["ther"] or Manual("none matched in CIViC"),
        "Diagnostic / Prognostic Evidence": c["dxpx"] or Manual("none matched in CIViC"), "Resistance Evidence": c["res"] or Manual("none matched in CIViC"),
        "Clinical Trial / Other Evidence": Manual(), "Somatic Tier": Manual("I / II / III / IV — not assigned by pipeline"), "Interpretive Comment": Manual(), "Reportable?": Sugg("Review")}
def block15(v):
    vaf = 100 * float(v["VAF"]); exp = purity / 2 if purity else None
    if exp is None: vi = Manual("tumor purity not in metadata.csv")
    elif vaf > 1.5 * exp + 10: vi = Sugg(f"VAF {vaf:.0f}% far exceeds expected ~{exp:.0f}% for a heterozygous clonal variant at {purity:.0f}% purity: consider LOH/copy-neutral LOH, germline origin, or underestimated purity")
    else: vi = Sugg("VAF consistent with expectation for purity")
    germ = Sugg("possible (VAF ≥ 90%)") if vaf >= 90 else (Sugg("possible (VAF near 50%)") if 40 <= vaf <= 60 else Sugg("less likely from VAF alone"))
    return {"Estimated Tumor Purity": M("tumor_pct"), "Variant VAF": f"{vaf:.1f}", "Expected VAF Based on Context": f"{exp:.1f} (heterozygous, diploid, clonal)" if exp else Manual("needs purity"),
        "Copy Number Status": "Unknown (CNV not assessed)", "LOH Evidence": "Unknown", "Possible Germline Contribution": germ, "Possible Subclonality": Sugg("VAF <30% may indicate subclonal or artefact") if vaf < 30 else Sugg("no indication from VAF"), "VAF Interpretation": vi}
BLOCKS = {10: block10, 11: block11, 12: block12, 14: block14, 15: block15}

def new_label(text):
    p = d.add_paragraph(); r = p.add_run(text); r.bold = True
    el = p._p; el.getparent().remove(el); return el
if DETAIL:
    for idx in tpl_ids:
        tbl = T[idx]; base = tbl._tbl; made = []
        for k in range(1, len(DETAIL)):
            made.append(docx.table.Table(copy.deepcopy(base), tbl._parent))
        tables = [tbl] + made
        prev = base
        for k, t in enumerate(tables):
            if k > 0: prev.addnext(t._tbl); prev = t._tbl
        for k, (t, v) in enumerate(zip(tables, DETAIL)):
            fill(t, BLOCKS[idx](v))
            t._tbl.addprevious(new_label(f"Variant {k + 1} of {len(DETAIL)}: {v['gene']} {v['HGVS_c']} {v['HGVS_p']}".strip()))
else:
    for idx in tpl_ids: fill(T[idx], {"Variant ID": Manual("no candidate variants")})

# ---------------------------------------------------------------- 16, 17
fill(T[16], {"CNV Method / Software": "Not performed", "Reference / Control Samples": "None (no normal reference pool)", "Normalization Method": "N/A", "Genome-wide CNV QC": "N/A",
  "Gene-level CNV Calls": "Not assessed", "Amplifications Detected": "Not assessed", "Deletions Detected": "Not assessed", "Estimated Copy Number": "Not assessed", "CNV Confidence / Score": "N/A",
  "CNV LOD / Validated Threshold": "N/A — CNV analysis not implemented or validated in this workflow", "CNV Reportable?": "No — not assessed"})
fill(T[17], {"MSI Status": "N/A", "MSI Score / Metric": "N/A", "TMB": "Not calculated (panel too small / not validated)", "TMB Classification / Cutoff": "N/A", "Other Biomarker": "None",
  "Method / Software Version": "N/A", "QC Status": "N/A"})

# ---------------------------------------------------------------- 18, 19
sm = [[v["gene"], v["HGVS_c"], v["HGVS_p"] or "—", v["consequence"], f"{100 * float(v['VAF']):.1f}%", v["depth"],
       ("ClinVar: " + clean(v["ClinVar"])) if v["ClinVar"] != "-" else "ClinVar: not listed", Manual("tier"), Sugg("Review")] for v in REVIEW]
rows_table(T[18], sm or [["No variants for interpretation"] + [""] * 8])
qrows = []
for m_, obs, crit, ok, imp, act in QC:
    qrows.append([m_, obs, crit, Sugg({True: "Pass", False: "Fail"}[ok]) if ok is not None else "N/A / Review", imp, act])
rows_table(T[19], qrows)

# ---------------------------------------------------------------- 21, 22
fill(T[20], {"Overall Sample QC": Sugg("LIMITED" if any_fail else "PASS"), "Variant Interpretation QC": Sugg("REVIEW"), "Final Report Status": Sugg("Hold"),
  "Primary Analyst": M("analyst"), "Analyst Date": Manual(), "Reviewer": M("reviewer"), "Reviewer Date": Manual(),
  "Comments": "Auto-generated by the targeted NGS somatic pipeline. All thresholds are unvalidated PLACEHOLDERS. CNV, fusion, MSI/TMB not assessed. Low-coverage genes are listed in sections 8 and 19. Human review required before any release."})
cvd = clinvar_date()
fill(T[21], {"FastQC": V["FastQC"] or Manual(), "fastp": V["fastp"] or Manual(), "BWA": V["BWA"] or Manual(), "samtools": V["samtools"] or Manual(), "GATK": (V["GATK"] or "") + " (includes Picard)",
  "bcftools": V["bcftools"] or Manual(), "SnpEff": (V["SnpEff"] or "").strip() + f" — database {CFG['snpeff_db']} (RefSeq)", "SnpSift": "Not used (bcftools annotate used for ClinVar)",
  "ClinVar Release / Date": f"ClinVar GRCh37 VCF, file date {cvd}", "gnomAD Release / Date": f"gnomAD AF-only resource (GATK best-practices b37), file date {fdate(CFG['germline_resource'])}",
  "COSMIC Release / Date": "Not used (license)", "OncoKB Version / Date": "Not used (license)",
  "Other Reference / Database": f"MANE Select v1.5 (transcripts); CIViC nightly downloaded {fdate(CFG.get('civic_variants', ''))}; hg19 reference; provenance file {S}.provenance.json"})

# ---------------------------------------------------------------- legend, save
leg = d.add_paragraph()
leg.add_run(f"Auto-filled by the targeted NGS somatic pipeline on {TODAY} for sample {S}. ").bold = True
r = leg.add_run("Yellow cells = ENTER MANUALLY, or SUGGESTED (pipeline suggestion that a person must confirm). Pass/Fail thresholds are unvalidated placeholders from config/settings.yaml. Checklist (section 20) and sign-off are intentionally left for people.")
r.font.highlight_color = WD_COLOR_INDEX.YELLOW
sub = d.paragraphs[1]._p; sub.addnext(leg._p)
out = f"{R}/report/{S}.QC_worksheet.docx"
d.save(out); print("saved", out, "| detail variants:", len(DETAIL), "| review variants:", len(REVIEW))
