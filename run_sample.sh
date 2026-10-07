#!/bin/bash
# One command for the whole pipeline (targeted NGS somatic pipeline v1.0).
#   ./run_sample.sh SAMPLE_ID R1.fastq.gz R2.fastq.gz   new sample: FASTQ -> Word worksheet
#   ./run_sample.sh SAMPLE_ID                           refresh an existing sample (e.g. after editing config/metadata.csv)
#   ./run_sample.sh status [SAMPLE_ID]                  progress and time estimate
#   ./run_sample.sh check [--full]                      pre-flight checks only (--full also verifies database checksums)
set -euo pipefail
ROOT="${PIPE_ROOT:-$(cd "$(dirname "$0")" && pwd)}"
cd "$ROOT"
CB="$(conda info --base)"
usage() { sed -n '2,6p' "$0"; exit 1; }
[ $# -ge 1 ] || usage

# ---------------------------------------------------------------- status
if [ "$1" = "status" ]; then
  ID="${2:-$(tail -n 1 config/samples.csv | cut -d, -f1)}"
  LOG="results/$ID.run.log"
  echo "Sample: $ID"
  if pgrep -f "bin/snakemake" >/dev/null; then echo "Pipeline: RUNNING"; else echo "Pipeline: not running"; fi
  [ -f "$LOG" ] || { echo "No log yet ($LOG)"; exit 0; }
  grep -a -E "steps \(" "$LOG" | tail -n 1 || echo "No step finished yet"
  grep -a -E "^(local)?rule " "$LOG" | tail -n 1 || true
  if pgrep -f "bwa mem" >/dev/null; then
    python3 - "$ID" <<'PY'
import json, re, subprocess, sys
s = sys.argv[1]; d = f"results/{s}"
tot = json.load(open(f"{d}/qc/{s}.fastp.json"))["summary"]["after_filtering"]["total_reads"]
done = sum(int(x) for x in re.findall(r"Processed (\d+) reads", open(f"{d}/logs/bwa.log", errors="ignore").read()))
pid = subprocess.run(["pgrep", "-f", "bwa mem"], capture_output=True, text=True).stdout.split()[0]
el = int(subprocess.run(["ps", "-o", "etimes=", "-p", pid], capture_output=True, text=True).stdout)
print(f"Alignment: {done:,} of {tot:,} reads ({100*done/tot:.1f}%), running {el/60:.0f} min", end="")
print(f", about {el*(tot-done)/done/60:.0f} min left" if done else "")
PY
  fi
  grep -a -q "ALL_DONE" "$LOG" && echo "FINISHED. Report: final_reports/$ID/$ID.QC_worksheet.docx"
  grep -a -i -E "error|exception" "$LOG" | tail -n 3 || true
  exit 0
fi

source "$CB/etc/profile.d/conda.sh"; conda activate somatic

# ---------------------------------------------------------------- check only
if [ "$1" = "check" ]; then
  python3 scripts/preflight.py "${@:2}"; exit $?
fi

# ---------------------------------------------------------------- new / refresh
ID="$1"; R1="${2:-}"; R2="${3:-}"
[[ "$ID" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "Sample ID may only contain letters, digits, _ and -"; exit 1; }
if pgrep -f "bin/snakemake" >/dev/null; then
  echo "Another pipeline run is active. Wait for it (./run_sample.sh status) and try again."; exit 1
fi

if [ -n "$R1" ]; then
  [ -n "$R2" ] || usage
  for f in "$R1" "$R2"; do [ -s "$f" ] || { echo "File not found: $f"; exit 1; }; done
  for f in "$R1" "$R2"; do case "$f" in *.fastq.gz|*.fq.gz) ;; *) echo "FASTQ must end in .fastq.gz (check for .gz.gz): $f"; exit 1;; esac; done
  case "$R1$R2" in *.gz.gz*) echo "Double .gz extension found"; exit 1;; esac
  mkdir -p samples
  B1="samples/$(basename "$R1")"; B2="samples/$(basename "$R2")"
  for p in "$R1:$B1" "$R2:$B2"; do
    src="${p%%:*}"; dst="${p#*:}"
    [ "$(readlink -f "$src")" = "$(readlink -f "$dst" 2>/dev/null || true)" ] || { echo "Copying $(basename "$src") into the fast Ubuntu folder..."; cp "$src" "$dst"; }
  done
  echo "Running pre-flight checks (includes a full read of both FASTQ files; a few minutes)..."
  python3 scripts/preflight.py "$B1" "$B2" || { echo "Pre-flight FAILED. Nothing started. Fix the FAIL lines above and run again."; exit 1; }
  python3 - "$ID" "$B1" "$B2" <<'PY'
import csv, sys
sid, r1, r2 = sys.argv[1:4]
rows = list(csv.DictReader(open("config/samples.csv")))
hit = [r for r in rows if r["sample_id"] == sid]
if hit and (hit[0]["r1"], hit[0]["r2"]) != (r1, r2):
    sys.exit(f"Sample {sid} already exists with different FASTQ files. Choose a new ID.")
if not hit:
    rows.append(dict(sample_id=sid, r1=r1, r2=r2))
    with open("config/samples.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample_id", "r1", "r2"], lineterminator="\n"); w.writeheader(); w.writerows(rows)
mrows = list(csv.DictReader(open("config/metadata.csv")))
hdr = list(mrows[0].keys()) if mrows else open("config/metadata.csv").readline().strip().split(",")
if not any(r["sample_id"] == sid for r in mrows):
    mrows.append({k: (sid if k == "sample_id" else "") for k in hdr})
    with open("config/metadata.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr, lineterminator="\n"); w.writeheader(); w.writerows(mrows)
print("Registered", sid)
PY
else
  grep -q "^$ID," config/samples.csv || { echo "Unknown sample $ID. Give R1 and R2 for a new sample."; exit 1; }
  python3 scripts/preflight.py || { echo "Pre-flight FAILED. Nothing started."; exit 1; }
fi

TARGET="results/$ID/report/$ID.QC_worksheet.docx"
mkdir -p results tmp
snakemake -n --cores 6 "$TARGET" --rerun-triggers mtime > "results/$ID.dryrun.log" 2>&1 || { tail -n 15 "results/$ID.dryrun.log"; echo "Dry run failed, nothing started."; exit 1; }
if [ -z "$R1" ] && grep -a -q -E "^(align|fastp) +[0-9]" "results/$ID.dryrun.log"; then
  echo "Refresh would redo alignment for $ID (inputs look changed). Not starting. See results/$ID.dryrun.log"; exit 1
fi
echo "Planned jobs:"; grep -a -E "^[a-z_0-9]+ +[0-9]+$" "results/$ID.dryrun.log" | sort -u | tr '\n' ';'; echo

LOCAL="$ROOT/final_reports/$ID"
nohup bash -c "
  source '$CB/etc/profile.d/conda.sh'; conda activate somatic; cd '$ROOT'
  snakemake --cores 6 '$TARGET' --rerun-triggers mtime &&
  mkdir -p '$LOCAL' && for f in $ID.QC_worksheet.docx $ID.variants.tsv $ID.gene_coverage.tsv $ID.region_coverage.tsv $ID.provenance.json; do cp results/$ID/report/\$f '$LOCAL/'; done &&
  echo ALL_DONE
" > "results/$ID.run.log" 2>&1 &
echo
echo "Started for $ID (it keeps running if you close this window)."
echo "  Progress:  ./run_sample.sh status $ID"
echo "  Report:    $LOCAL/$ID.QC_worksheet.docx"
echo "  Fill config/metadata.csv for $ID while it runs; the report step reads it last."
