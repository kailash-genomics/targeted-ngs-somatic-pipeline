# targeted NGS somatic pipeline v1.0 -- tumor-only targeted panel, hg19. One file, no includes.
import os
import pandas as pd
configfile: "config/settings.yaml"

S = (
    pd.read_csv("config/samples.csv").set_index("sample_id")
    if os.path.exists("config/samples.csv")
    else pd.DataFrame(columns=["r1", "r2"]).rename_axis("sample_id")
)
SAMPLES = list(S.index)
T = config["threads"]
MEM = config["java_mem_gb"]
REF = config["reference"]
BED = config["bed"]
DICT = REF.replace(".fa", ".dict")
TMP = config["tmp_dir"]

wildcard_constraints:
    s="[A-Za-z0-9_-]+"

rule all:
    input: expand("results/{s}/report/{s}.QC_worksheet.docx", s=SAMPLES)

rule worksheets:
    input: expand("results/{s}/report/{s}.QC_worksheet.docx", s=SAMPLES)

# ------------------------------------------------------------------ QC, trimming, alignment
rule check_inputs:
    input:
        r1=lambda w: S.loc[w.s, "r1"],
        r2=lambda w: S.loc[w.s, "r2"],
        ref=REF, bed=BED
    output: "results/{s}/qc/{s}.inputs_ok.txt"
    shell:
        """
        mkdir -p results/{wildcards.s}/qc
        gzip -t {input.r1} && gzip -t {input.r2}
        test -s {input.ref}.bwt && test -s {input.ref}.fai && test -s {DICT}
        echo OK > {output}
        """

rule fastp:
    input:
        ok="results/{s}/qc/{s}.inputs_ok.txt",
        r1=lambda w: S.loc[w.s, "r1"],
        r2=lambda w: S.loc[w.s, "r2"]
    output:
        r1="results/{s}/trim/{s}_R1.trim.fastq.gz",
        r2="results/{s}/trim/{s}_R2.trim.fastq.gz",
        js="results/{s}/qc/{s}.fastp.json",
        html="results/{s}/qc/{s}.fastp.html"
    threads: T
    log: "results/{s}/logs/fastp.log"
    shell:
        """
        mkdir -p results/{wildcards.s}/trim results/{wildcards.s}/logs
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
          --detect_adapter_for_pe --trim_poly_g --cut_tail --cut_tail_mean_quality 20 \
          --length_required 50 -w {threads} -j {output.js} -h {output.html} 2> {log}
        """

rule fastqc:
    input: r1="results/{s}/trim/{s}_R1.trim.fastq.gz"
    output: "results/{s}/qc/{s}_R1_fastqc.html"
    threads: 2
    log: "results/{s}/logs/fastqc.log"
    shell:
        """
        d=results/{wildcards.s}/qc
        ln -sf "$(readlink -f {input.r1})" $d/{wildcards.s}_R1.fastq.gz
        fastqc -t {threads} -o $d $d/{wildcards.s}_R1.fastq.gz > {log} 2>&1
        rm -f $d/{wildcards.s}_R1.fastq.gz
        """

rule align:
    input:
        r1="results/{s}/trim/{s}_R1.trim.fastq.gz",
        r2="results/{s}/trim/{s}_R2.trim.fastq.gz"
    output: temp("results/{s}/align/{s}.sorted.bam")
    threads: T
    log: "results/{s}/logs/bwa.log"
    params: tmp=TMP
    shell:
        """
        mkdir -p {params.tmp} results/{wildcards.s}/logs
        RG="@RG\\tID:{wildcards.s}\\tSM:{wildcards.s}\\tPL:ILLUMINA\\tLB:{wildcards.s}_lib"
        bwa mem -t {threads} -R "$RG" {REF} {input.r1} {input.r2} 2> {log} \
          | samtools sort -@ 2 -m 1G -T {params.tmp}/{wildcards.s} -o {output} -
        """

rule markdup:
    input: "results/{s}/align/{s}.sorted.bam"
    output:
        bam="results/{s}/align/{s}.dedup.bam",
        met="results/{s}/qc/{s}.dupmetrics.txt"
    log: "results/{s}/logs/markdup.log"
    params: tmp=TMP
    shell:
        """
        mkdir -p {params.tmp}
        gatk --java-options "-Xmx{MEM}g" MarkDuplicates -I {input} -O {output.bam} \
          -M {output.met} --TMP_DIR {params.tmp} --CREATE_INDEX false > {log} 2>&1
        """

rule index:
    input: "results/{s}/align/{s}.dedup.bam"
    output: "results/{s}/align/{s}.dedup.bam.bai"
    shell: "samtools index {input}"

rule flagstat:
    input: bam="results/{s}/align/{s}.dedup.bam", bai="results/{s}/align/{s}.dedup.bam.bai"
    output: "results/{s}/qc/{s}.flagstat.txt"
    shell: "samtools flagstat {input.bam} > {output}"

rule intervals:
    input: bed=BED, ref=REF
    output: "results/{s}/qc/{s}.targets.interval_list"
    shell: "gatk BedToIntervalList -I {input.bed} -O {output} -SD {DICT}"

rule hsmetrics:
    input:
        bam="results/{s}/align/{s}.dedup.bam",
        bai="results/{s}/align/{s}.dedup.bam.bai",
        iv="results/{s}/qc/{s}.targets.interval_list"
    output: "results/{s}/qc/{s}.hsmetrics.txt"
    log: "results/{s}/logs/hsmetrics.log"
    shell:
        """
        gatk --java-options "-Xmx{MEM}g" CollectHsMetrics -I {input.bam} -O {output} \
          -R {REF} --BAIT_INTERVALS {input.iv} --TARGET_INTERVALS {input.iv} > {log} 2>&1
        """

rule bamstats:
    input:
        bam="results/{s}/align/{s}.dedup.bam",
        bai="results/{s}/align/{s}.dedup.bam.bai",
        bed=BED
    output:
        txt="results/{s}/qc/{s}.bamstats.txt",
        st="results/{s}/qc/{s}.samtools_stats.txt"
    threads: 4
    shell: "bash scripts/bamstats.sh {input.bam} {input.bed} {threads} {output.txt} {output.st}"

# ------------------------------------------------------------------ variant calling (Mutect2, tumor-only)
rule mutect2:
    input:
        bam="results/{s}/align/{s}.dedup.bam",
        bai="results/{s}/align/{s}.dedup.bam.bai"
    output:
        vcf="results/{s}/variants/{s}.mutect2.vcf.gz",
        f1r2="results/{s}/variants/{s}.f1r2.tar.gz"
    threads: 4
    log: "results/{s}/logs/mutect2.log"
    shell:
        """
        mkdir -p results/{wildcards.s}/variants
        gatk --java-options "-Xmx{MEM}g" Mutect2 -R {REF} -I {input.bam} \
          -L {BED} --germline-resource {config[germline_resource]} \
          --f1r2-tar-gz {output.f1r2} --native-pair-hmm-threads {threads} \
          -O {output.vcf} > {log} 2>&1
        """

rule orientation:
    input: "results/{s}/variants/{s}.f1r2.tar.gz"
    output: "results/{s}/variants/{s}.orientation.tar.gz"
    log: "results/{s}/logs/orientation.log"
    shell: 'gatk --java-options "-Xmx{MEM}g" LearnReadOrientationModel -I {input} -O {output} > {log} 2>&1'

rule pileup:
    input:
        bam="results/{s}/align/{s}.dedup.bam",
        bai="results/{s}/align/{s}.dedup.bam.bai"
    output: "results/{s}/variants/{s}.pileup.table"
    log: "results/{s}/logs/pileup.log"
    shell:
        """
        mkdir -p results/{wildcards.s}/variants
        gatk --java-options "-Xmx{MEM}g" GetPileupSummaries -I {input.bam} \
          -V {config[germline_resource]} -L {BED} -O {output} > {log} 2>&1
        """

rule contamination:
    input: "results/{s}/variants/{s}.pileup.table"
    output: "results/{s}/variants/{s}.contamination.table"
    log: "results/{s}/logs/contamination.log"
    shell: 'gatk --java-options "-Xmx{MEM}g" CalculateContamination -I {input} -O {output} > {log} 2>&1'

rule filter:
    input:
        vcf="results/{s}/variants/{s}.mutect2.vcf.gz",
        ori="results/{s}/variants/{s}.orientation.tar.gz",
        con="results/{s}/variants/{s}.contamination.table"
    output: "results/{s}/variants/{s}.filtered.vcf.gz"
    log: "results/{s}/logs/filter.log"
    shell:
        """
        gatk --java-options "-Xmx{MEM}g" FilterMutectCalls -R {REF} -V {input.vcf} \
          --ob-priors {input.ori} --contamination-table {input.con} -O {output} > {log} 2>&1
        """

rule normalize:
    input: "results/{s}/variants/{s}.filtered.vcf.gz"
    output:
        vcf="results/{s}/variants/{s}.norm.vcf.gz",
        tbi="results/{s}/variants/{s}.norm.vcf.gz.tbi"
    shell:
        """
        bcftools norm -m -any -f {REF} {input} -Oz -o {output.vcf}
        tabix -p vcf -f {output.vcf}
        """

# ------------------------------------------------------------------ annotation, tables, worksheet
rule annotate:
    input:
        vcf="results/{s}/variants/{s}.norm.vcf.gz",
        tbi="results/{s}/variants/{s}.norm.vcf.gz.tbi"
    output: "results/{s}/variants/{s}.annotated.vcf"
    log: "results/{s}/logs/annotate.log"
    shell:
        """
        P=results/{wildcards.s}/variants/{wildcards.s}
        bcftools view -f PASS {input.vcf} -Oz -o $P.pass.vcf.gz
        tabix -p vcf -f $P.pass.vcf.gz
        bcftools annotate -a {config[clinvar]} -c ID,INFO/ALLELEID,INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN $P.pass.vcf.gz -Oz -o $P.pass.clinvar.vcf.gz
        source $(conda info --base)/etc/profile.d/conda.sh
        set +u; conda activate snpeff_env; set -u
        snpEff -Xmx4g -dataDir $PWD/{config[snpeff_datadir]} -noStats {config[snpeff_db]} $P.pass.clinvar.vcf.gz > {output} 2> {log}
        """

rule variant_table:
    input:
        vcf="results/{s}/variants/{s}.annotated.vcf",
        bed=BED,
        cfg="config/settings.yaml",
        scr="scripts/final_table.py"
    output: "results/{s}/report/{s}.variants.tsv"
    shell:
        """
        mkdir -p results/{wildcards.s}/report
        python scripts/final_table.py {input.vcf} {input.bed} {input.cfg} {config[mane]} {output}
        """

rule coverage_tables:
    input:
        bam="results/{s}/align/{s}.dedup.bam",
        bai="results/{s}/align/{s}.dedup.bam.bai",
        bed=BED,
        cfg="config/settings.yaml",
        scr="scripts/coverage_tables.py"
    output:
        gene="results/{s}/report/{s}.gene_coverage.tsv",
        region="results/{s}/report/{s}.region_coverage.tsv"
    threads: 4
    shell:
        """
        mkdir -p results/{wildcards.s}/report
        mosdepth -t {threads} -n --by {input.bed} --thresholds 100,250,500,1000 results/{wildcards.s}/qc/{wildcards.s}.cov {input.bam}
        python scripts/coverage_tables.py {wildcards.s}
        """

rule fill_worksheet:
    input:
        tpl="config/template.docx",
        meta="config/metadata.csv",
        cfg="config/settings.yaml",
        scr="scripts/fill_report.py",
        var="results/{s}/report/{s}.variants.tsv",
        gcov="results/{s}/report/{s}.gene_coverage.tsv",
        rcov="results/{s}/report/{s}.region_coverage.tsv",
        bst="results/{s}/qc/{s}.bamstats.txt",
        st="results/{s}/qc/{s}.samtools_stats.txt",
        hs="results/{s}/qc/{s}.hsmetrics.txt",
        dup="results/{s}/qc/{s}.dupmetrics.txt",
        fl="results/{s}/qc/{s}.flagstat.txt",
        fj="results/{s}/qc/{s}.fastp.json",
        ann="results/{s}/variants/{s}.annotated.vcf",
        mut="results/{s}/variants/{s}.mutect2.vcf.gz",
        flt="results/{s}/variants/{s}.filtered.vcf.gz",
        contam="results/{s}/variants/{s}.contamination.table"
    output:
        docx="results/{s}/report/{s}.QC_worksheet.docx",
        prov="results/{s}/report/{s}.provenance.json"
    shell: "python scripts/fill_report.py {wildcards.s}"
