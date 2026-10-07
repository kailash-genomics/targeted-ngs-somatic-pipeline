#!/usr/bin/env python3
"""Pre-flight checks. Usage: preflight.py [R1 R2] [--full]   (run inside the 'somatic' environment)"""
import os, sys, shutil, subprocess, re, gzip
import yaml
root = os.environ.get("PIPE_ROOT", os.getcwd()); os.chdir(root)
full = "--full" in sys.argv
args = [a for a in sys.argv[1:] if not a.startswith("--")]
R1, R2 = (args + [None, None])[:2]
fails = warns = 0
def _hook(t, v, tb):
    print(f"  FAIL  checker error ({t.__name__}: {v}) - treat as a failed pre-flight and tell Claude")
    sys.exit(1)
sys.excepthook = _hook
def ok(m): print(f"  OK    {m}")
def bad(m):
    global fails; fails += 1; print(f"  FAIL  {m}")
def warn(m):
    global warns; warns += 1; print(f"  WARN  {m}")
def have(p, what=None):
    if os.path.exists(p) and os.path.getsize(p) > 0: ok(what or p)
    else: bad(f"missing or empty: {p}")
def sh(c, t=300):
    try: return subprocess.run(c, shell=True, capture_output=True, text=True, timeout=t)
    except Exception as e: return None

print("[1] Settings"); 
try:
    cfg = yaml.safe_load(open("config/settings.yaml")); ok("config/settings.yaml readable")
except Exception as e:
    print("  FAIL  cannot read config/settings.yaml:", e); sys.exit(1)
for k in ("reference", "bed", "snpeff_db", "snpeff_datadir", "clinvar", "germline_resource", "civic_variants", "civic_evidence", "mane", "tmp_dir", "threads", "java_mem_gb", "thresholds"):
    ok(f"setting {k}") if k in cfg else bad(f"setting missing: {k}")
th = cfg.get("thresholds", {})
for k in ("lod_snv_vaf", "lod_indel_vaf", "min_depth", "min_alt_reads", "min_pct_100x", "min_q30_pct", "max_duplicate_pct", "low_conf_vaf", "low_conf_alt_reads"):
    ok(f"threshold {k}={th[k]}") if isinstance(th.get(k), (int, float)) else bad(f"threshold missing/not numeric: {k}")
if cfg.get("threads", 0) > (os.cpu_count() or 1): warn(f"threads={cfg['threads']} exceeds CPU count {os.cpu_count()}")

if fails:
    print(f"\nPRE-FLIGHT: {fails} failed. Fix config/settings.yaml first (later checks need these settings).")
    sys.exit(1)
print("[2] Reference and databases")
ref = cfg["reference"]
for ext in ("", ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"): have(ref + ext)
have(ref.replace(".fa", ".dict"))
have(cfg["bed"])
for f in (cfg["clinvar"], cfg["germline_resource"]): have(f); have(f + ".tbi")
have(os.path.join(cfg["snpeff_datadir"], cfg["snpeff_db"], "snpEffectPredictor.bin"))
for k in ("mane", "civic_variants", "civic_evidence"): have(cfg[k])
fai = [l.split("\t")[0] for l in open(ref + ".fai")] if os.path.exists(ref + ".fai") else []
bedchr = {l.split("\t")[0] for l in open(cfg["bed"])} if os.path.exists(cfg["bed"]) else set()
(ok("BED chromosome names match reference") if bedchr and bedchr <= set(fai) else bad(f"BED chromosomes not in reference: {sorted(bedchr - set(fai))[:5]}"))
if os.path.exists(cfg["bed"]):
    rows = [l.rstrip("\n").split("\t") for l in open(cfg["bed"]) if l.strip() and not l.startswith(("#", "track", "browser"))]
    nb = sum(1 for r in rows if len(r) < 3 or not r[1].isdigit() or not r[2].isdigit() or int(r[2]) <= int(r[1]))
    (bad(f"BED has {nb} invalid lines (need chrom, start, end, with end > start)") if nb else ok(f"BED: {len(rows)} regions, all valid"))
    tb = sum(int(r[2]) - int(r[1]) for r in rows if len(r) >= 3 and r[1].isdigit() and r[2].isdigit() and int(r[2]) > int(r[1]))
    ok(f"BED target size {tb / 1e6:.2f} Mb")
    rs = sh(f"sort -k1,1 -k2,2n -c {cfg['bed']}")
    if rs is not None and rs.returncode != 0: bad("BED is not sorted (run: sort -k1,1 -k2,2n)")
    else:
        ok("BED is sorted")
        rm = sh(f"bedtools merge -i {cfg['bed']}")
        if rm is not None and rm.returncode == 0 and rm.stdout.strip():
            nm = len(rm.stdout.strip().splitlines())
            if nm < len(rows):
                warn(f"BED has overlapping/touching regions ({len(rows)} lines merge to {nm}): coverage may be double counted; merge the BED first")
            else:
                ok("BED has no overlapping regions")
    labs = [r[3] for r in rows if len(r) >= 4]
    if not labs: bad("BED has no 4th column; gene names are required")
    else:
        names = {g.replace("_promoter", "") for x in labs for g in x.split(",")}
        odd = sorted(n for n in names if not re.match(r"^[A-Za-z0-9._-]+$", n))
        if odd: bad(f"BED 4th column is not plain gene symbols (e.g. {odd[:3]}); variants would be silently dropped")
        elif os.path.exists(cfg["mane"]):
            with gzip.open(cfg["mane"], "rt") as fh: syms = {l.split("\t")[3] for l in fh if not l.startswith("#")}
            fr = sum(n in syms for n in names) / len(names)
            msg = f"{len(names)} gene labels, {100 * fr:.0f}% found in MANE symbols"
            (ok(msg) if fr >= 0.9 else (warn(msg + " (old/alias symbols? check)") if fr >= 0.6 else bad(msg + ": 4th column does not look like HGNC gene symbols")))
for f in (cfg["clinvar"], cfg["germline_resource"]):
    r = sh(f"tabix -l {f} | head -3")
    (ok(f"{os.path.basename(f)} uses chr-style contigs") if r and r.stdout.startswith("chr") else bad(f"{os.path.basename(f)} contigs not chr-style (rename step missing?)"))
if full and os.path.exists("db/CHECKSUMS.md5"):
    r = sh("cd db && md5sum -c CHECKSUMS.md5", t=1800)
    (ok("database checksums match") if r and r.returncode == 0 else bad("database checksum mismatch:\n" + (r.stdout if r else "")))

print("[3] Workflow files")
for f in ("Snakefile", "config/samples.csv", "config/metadata.csv", "config/template.docx", "scripts/final_table.py", "scripts/coverage_tables.py", "scripts/fill_report.py", "scripts/bamstats.sh"): have(f)
try:
    import docx
    n = len(docx.Document("config/template.docx").tables)
    ok("worksheet template has 22 tables") if n == 22 else bad(f"template has {n} tables, expected 22")
except Exception as e: bad(f"python-docx/template problem: {e}")
if os.path.exists("Snakefile"):
    sf = open("Snakefile").read(); names = re.findall(r"^rule (\w+):", sf, re.M)
    dup = [n for n in set(names) if names.count(n) > 1]
    bad(f"duplicate rules in Snakefile: {dup}") if dup else ok(f"Snakefile has {len(names)} unique rules")

print("[4] Tools")
for t in ("fastp", "fastqc", "bwa", "samtools", "bcftools", "tabix", "bgzip", "gatk", "mosdepth", "bedtools", "snakemake", "java"):
    ok(f"{t} -> {shutil.which(t)}") if shutil.which(t) else bad(f"{t} not found (is the 'somatic' environment active?)")
cb = (sh("conda info --base") or type("x", (), {"stdout": ""})).stdout.strip()
(ok("SnpEff environment present") if cb and os.path.exists(f"{cb}/envs/snpeff_env/bin/snpEff") else bad("snpeff_env/bin/snpEff not found"))
for m in ("yaml", "pandas", "docx", "pysam"):
    try: __import__(m); ok(f"python module {m}")
    except Exception: bad(f"python module {m} missing")

print("[5] Machine")
mem = {l.split(":")[0]: int(l.split()[1]) for l in open("/proc/meminfo")}
tot, av = mem["MemTotal"] / 1e6, mem["MemAvailable"] / 1e6
(ok(f"memory total {tot:.1f} GB") if tot >= 9 else bad(f"memory total {tot:.1f} GB (<9 GB: BWA may be killed)"))
(ok(f"memory available {av:.1f} GB") if av >= 6 else warn(f"only {av:.1f} GB memory available now; close other programs"))
d = shutil.disk_usage("/mnt/c" if os.path.exists("/mnt/c") else ".")
free = d.free / 1e9
(ok(f"free disk {free:.0f} GB") if free >= 60 else (warn(f"free disk {free:.0f} GB (<60)") if free >= 25 else bad(f"free disk {free:.0f} GB (<25): not enough")))
try: _r = subprocess.run(["pgrep", "-f", "bin/snakemake"], capture_output=True, text=True).stdout.split()
except Exception: _r = []
bad("another Snakemake run is active (pids " + ",".join(_r) + ")") if _r else ok("no other pipeline run active")

if R1 and R2:
    print("[6] FASTQ files")
    for f in (R1, R2):
        if not os.path.exists(f): bad(f"not found: {f}"); continue
        (ok(f"{os.path.basename(f)} ({os.path.getsize(f) / 1e6:.0f} MB)") if os.path.getsize(f) > 1e6 else bad(f"{f} is smaller than 1 MB"))
        (ok("name ends in .fastq.gz") if re.search(r"\.(fastq|fq)\.gz$", f) and not f.endswith(".gz.gz") else bad(f"name must end in .fastq.gz (no .gz.gz): {f}"))
    if os.path.exists(R1) and os.path.exists(R2):
        def head(f):
            with gzip.open(f, "rt") as fh: name = fh.readline().split()[0]
            return name[:-2] if name.endswith(("/1", "/2")) else name
        try:
            h1, h2 = head(R1), head(R2)
            (ok("first read names of R1 and R2 match") if h1 == h2 else bad(f"R1/R2 read names differ: {h1} vs {h2} (not a pair?)"))
        except Exception as e: bad(f"cannot read FASTQ header: {e}")
        procs = [subprocess.Popen(["bash", "-c", f"set -o pipefail; gzip -dc '{f}' | wc -l"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) for f in (R1, R2)]
        counts = []
        for f, p in zip((R1, R2), procs):
            out, _ = p.communicate()
            if p.returncode != 0: bad(f"corrupt gzip: {os.path.basename(f)}"); counts.append(None)
            else:
                n = int(out.strip()); counts.append(n)
                (ok(f"{os.path.basename(f)} intact, {n // 4:,} reads") if n % 4 == 0 else bad(f"{os.path.basename(f)}: line count {n} not divisible by 4 (truncated?)"))
        if None not in counts:
            (ok("R1 and R2 have the same number of reads") if counts[0] == counts[1] else bad(f"read counts differ: R1 {counts[0] // 4:,} vs R2 {counts[1] // 4:,}"))

print(f"\nPRE-FLIGHT: {fails} failed, {warns} warnings")
sys.exit(1 if fails else 0)
