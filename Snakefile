
# Full imputation-prep Snakefile with liftover/updatebuild/none modes
# and a QC manifest with before/after counts at each step.

configfile: "config.yaml"

import os

# ---------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------
plinkIn = str(config.get("plink_genotype_file", "")).strip()
if not plinkIn:
    raise ValueError("config.yaml must define plink_genotype_file")

reqScripts = str(config.get("req_scripts", ".")).rstrip("/")
chainStrandIn = str(config.get("strand_chain_file", "")).strip()
ref_file = str(config.get("ref_file", "")).strip()
MODE = str(config.get("mode", "none")).lower()

MIND = float(config.get("mind", 0.05))
GENO = float(config.get("geno", 0.05))
MAF = float(config.get("maf", 0.01))
HWE = float(config.get("hwe", 1e-6))

update_build_script = str(config.get("update_build_script", f"{reqScripts}/update_build.sh")).strip()

# ---------------------------------------------------------------------
# Chromosomes
# ---------------------------------------------------------------------
CHROMS = [str(i) for i in range(1, 23)]

if config.get("include_chrX", False):
    CHROMS.append("X")

if config.get("include_mt", False):
    CHROMS.append("MT")

extra = config.get("extra_chroms", [])
if extra:
    CHROMS += [str(x) for x in extra]

_autosomes = [str(i) for i in range(1, 23)]
_non_auto = [c for c in CHROMS if c not in _autosomes]
CHROMS = _autosomes + _non_auto

# ---------------------------------------------------------------------
# Mode / build prefix
# ---------------------------------------------------------------------
if MODE == "liftover":
    build_prefix = "build38/subjects.lifted.hg38.final"
elif MODE == "updatebuild":
    build_prefix = "build38/subjects"
else:
    build_prefix = "build38/subjects"

print("Running pipeline MODE =", MODE)
print("Processing chromosomes:", CHROMS)
print("Final build prefix:", build_prefix)

# ---------------------------------------------------------------------
# Helpers for QC text files
# ---------------------------------------------------------------------
def count_lines(path):
    with open(path, "r") as handle:
        return sum(1 for _ in handle)

def pct_removed(before, after):
    if before == 0:
        return 0.0
    return ((before - after) / before) * 100.0

def format_report(title, before_samples, after_samples, before_snps, after_snps):
    return (
        f"{title}\n"
        f"Samples before: {before_samples}\n"
        f"Samples after:  {after_samples}\n"
        f"Samples removed: {before_samples - after_samples}\n"
        f"Samples percent removed: {pct_removed(before_samples, after_samples):.2f}\n"
        f"SNPs before: {before_snps}\n"
        f"SNPs after:  {after_snps}\n"
        f"SNPs removed: {before_snps - after_snps}\n"
        f"SNPs percent removed: {pct_removed(before_snps, after_snps):.2f}\n\n"
    )

# ---------------------------------------------------------------------
# Final targets
# ---------------------------------------------------------------------
rule all:
    input:
        "run_manifest.txt",
        "qc_stats/liftover.txt",
        "qc_stats/mind.txt",
        "qc_stats/geno.txt",
        "qc_stats/maf.txt",
        "qc_stats/hwe.txt",
        "qc_stats/hrc.txt",
        "qc_stats/final.txt",
        expand("vcf_MIS/subjects-updated-chr{chrom}.vcf.gz", chrom=CHROMS),
        expand("vcf_MIS/subjects-updated-chr{chrom}.vcf.gz.tbi", chrom=CHROMS)

# ---------------------------------------------------------------------
# Build conversion
# ---------------------------------------------------------------------
if MODE == "liftover":

    rule triple_liftover_hg19_to_hg38:
        input:
            bed = f"{plinkIn}.bed",
            bim = f"{plinkIn}.bim",
            fam = f"{plinkIn}.fam"
        output:
            bed = build_prefix + ".bed",
            bim = build_prefix + ".bim",
            fam = build_prefix + ".fam"
        params:
            reqScripts = reqScripts,
            raw_prefix = "build38/subjects",
            liftover_prefix = "build38/subjects_liftover_1",
            lifted_prefix = "build38/subjects.lifted",
            updated_prefix = "build38/subjects.lifted.hg38",
            final_prefix = "build38/subjects.lifted.hg38.final"
        conda:
            "envs/buildtools.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p build38

            # Normalize chromosome naming for PLINK
            plink --bfile {plinkIn} --output-chr MT --make-bed --out {params.raw_prefix}

            # Triple-liftOver
            perl {params.reqScripts}/triple-liftOver-main/tripleliftover_v133.pl --bim {params.raw_prefix}.bim --base hg19 --target hg38 --outprefix {params.liftover_prefix}
            awk 'NR>1 {{print $1}}' {params.liftover_prefix}.invr.txt > build38/variant.name.in.liftOver.output
            awk 'NR>1 {{print $1 "\t" $3}}' {params.liftover_prefix}.invr.txt > build38/variant.name.hg38.position
            awk 'NR>1 && $5 == "inverted" {{print $1}}' {params.liftover_prefix}.invr.txt > build38/inverted.SNV.name

            plink --bfile {params.raw_prefix} --extract build38/variant.name.in.liftOver.output --make-bed --out {params.lifted_prefix}
            plink --bfile {params.lifted_prefix} --update-map build38/variant.name.hg38.position --make-bed --out {params.updated_prefix}
            plink --bfile {params.updated_prefix} --flip build38/inverted.SNV.name --make-bed --out {params.final_prefix}
            """

elif MODE == "updatebuild":

    rule build:
        input:
            plinkIn + ".bed",
            plinkIn + ".bim",
            plinkIn + ".fam"
        params:
            inp = plinkIn,
            stInp = chainStrandIn,
            updateBuild = update_build_script,
            out = "build38/subjects"
        output:
            "build38/subjects.bed",
            "build38/subjects.bim",
            "build38/subjects.fam"
        conda:
            "envs/buildtools.yaml"
        run:
            shell("mkdir -p build38")
            shell("bash {params.updateBuild} {params.inp} {params.stInp} {params.out}")

else:

    rule pass_through_build:
        input:
            bed = f"{plinkIn}.bed",
            bim = f"{plinkIn}.bim",
            fam = f"{plinkIn}.fam"
        output:
            bed = build_prefix + ".bed",
            bim = build_prefix + ".bim",
            fam = build_prefix + ".fam"
        conda:
            "envs/buildtools.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p build38
            plink --bfile {plinkIn} --make-bed --out {build_prefix}
            """

# ---------------------------------------------------------------------
# QC pipeline
# ---------------------------------------------------------------------
rule plink_sub_call_rate:
    input:
        build_prefix + ".bed",
        build_prefix + ".bim",
        build_prefix + ".fam"
    output:
        temp("plink_sub_cr/subjects.bed"),
        temp("plink_sub_cr/subjects.bim"),
        temp("plink_sub_cr/subjects.fam")
    params:
        mind = MIND
    conda:
        "envs/imputation.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p plink_sub_cr
        plink --bfile {build_prefix} --mind {params.mind} --make-bed --out plink_sub_cr/subjects
        """

rule plink_geno:
    input:
        "plink_sub_cr/subjects.bed",
        "plink_sub_cr/subjects.bim",
        "plink_sub_cr/subjects.fam"
    output:
        temp("plink_geno/subjects.bed"),
        temp("plink_geno/subjects.bim"),
        temp("plink_geno/subjects.fam")
    params:
        geno = GENO
    conda:
        "envs/imputation.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p plink_geno
        plink --bfile plink_sub_cr/subjects --geno {params.geno} --make-bed --out plink_geno/subjects
        """

rule plink_maf:
    input:
        "plink_geno/subjects.bed",
        "plink_geno/subjects.bim",
        "plink_geno/subjects.fam"
    output:
        temp("plink_maf/subjects.bed"),
        temp("plink_maf/subjects.bim"),
        temp("plink_maf/subjects.fam")
    params:
        maf = MAF
    conda:
        "envs/imputation.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p plink_maf
        plink --bfile plink_geno/subjects --maf {params.maf} --make-bed --out plink_maf/subjects
        """

rule plink_hwe:
    input:
        "plink_maf/subjects.bed",
        "plink_maf/subjects.bim",
        "plink_maf/subjects.fam"
    output:
        temp("plink_snp_cr/subjects.bed"),
        temp("plink_snp_cr/subjects.bim"),
        temp("plink_snp_cr/subjects.fam")
    params:
        hwe = HWE
    conda:
        "envs/imputation.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p plink_snp_cr
        plink --bfile plink_maf/subjects --hwe {params.hwe} --make-bed --out plink_snp_cr/subjects
        """

rule plink_afterQC_freq:
    input:
        "plink_snp_cr/subjects.bed",
        "plink_snp_cr/subjects.bim",
        "plink_snp_cr/subjects.fam"
    output:
        temp("plink_afterQC_freq/subjects.frq")
    conda:
        "envs/imputation.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p plink_afterQC_freq
        plink --bfile plink_snp_cr/subjects --freq --out plink_afterQC_freq/subjects
        """

rule HRC_check:
    input:
        bed = "plink_snp_cr/subjects.bed",
        bim = "plink_snp_cr/subjects.bim",
        fam = "plink_snp_cr/subjects.fam",
        freq = "plink_afterQC_freq/subjects.frq",
        ref = ref_file,
        hrc1000G = reqScripts + "/HRC-1000G-check-bim.pl"
    output:
        "plink_snp_cr/subjects-updated.bed",
        "plink_snp_cr/subjects-updated.bim",
        "plink_snp_cr/subjects-updated.fam"
    conda:
        "envs/buildtools.yaml"
    shell:
        r"""
        set -euo pipefail
        perl {input.hrc1000G} -b {input.bim} -f {input.freq} -r {input.ref} -h -v -n 
        chmod u+x plink_snp_cr/Run-plink.sh
        sed -i '/--recode vcf/d' plink_snp_cr/Run-plink.sh
        sed -i 's|rm TEMP\*|rm -f plink_snp_cr/TEMP*|' plink_snp_cr/Run-plink.sh
        # Remove per-chromosome output generation
        sed -i '/--real-ref-alleles --make-bed --chr /d' plink_snp_cr/Run-plink.sh
        sed -i '/--real-ref-alleles --recode vcf --chr /d' plink_snp_cr/Run-plink.sh
        sed -i '/--recode vcf/d' plink_snp_cr/Run-plink.sh
        bash plink_snp_cr/Run-plink.sh
        """


rule make_vcf:
    input:
        "plink_snp_cr/subjects-updated.bed",
        "plink_snp_cr/subjects-updated.bim",
        "plink_snp_cr/subjects-updated.fam"
    output:
        vcf = "vcf_MIS/subjects-updated-chr{chrom}.vcf.gz",
        tbi = "vcf_MIS/subjects-updated-chr{chrom}.vcf.gz.tbi"
    params:
        chrom = lambda wildcards: wildcards.chrom
    conda:
        "envs/imputation.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p vcf_MIS

        plink --bfile plink_snp_cr/subjects-updated --real-ref-alleles --output-chr chrMT --recode vcf-iid --chr {params.chrom}  --out vcf_MIS/tmp_chr{params.chrom}
        bcftools view -Oz -o {output.vcf} vcf_MIS/tmp_chr{params.chrom}.vcf
        bcftools index -t {output.vcf}
        rm -f vcf_MIS/tmp_chr{params.chrom}.vcf
        """

# ---------------------------------------------------------------------
# QC reports
# ---------------------------------------------------------------------
rule qc_liftover:
    input:
        before_fam = f"{plinkIn}.fam",
        before_bim = f"{plinkIn}.bim",
        after_fam = build_prefix + ".fam",
        after_bim = build_prefix + ".bim"
    output:
        "qc_stats/liftover.txt"
    run:
        os.makedirs("qc_stats", exist_ok=True)
        before_samples = count_lines(input.before_fam)
        after_samples = count_lines(input.after_fam)
        before_snps = count_lines(input.before_bim)
        after_snps = count_lines(input.after_bim)

        with open(output[0], "w") as fh:
            if MODE == "liftover":
                fh.write(format_report("LIFTOVER (hg19 -> hg38)", before_samples, after_samples, before_snps, after_snps))
            else:
                fh.write(f"LIFTOVER: not applied (mode={MODE})\n\n")

rule qc_mind:
    input:
        before_fam = f"{plinkIn}.fam",
        after_fam = "plink_sub_cr/subjects.fam",
        before_bim = f"{plinkIn}.bim",
        after_bim = "plink_sub_cr/subjects.bim"
    output:
        "qc_stats/mind.txt"
    run:
        os.makedirs("qc_stats", exist_ok=True)
        before_samples = count_lines(input.before_fam)
        after_samples = count_lines(input.after_fam)
        before_snps = count_lines(input.before_bim)
        after_snps = count_lines(input.after_bim)

        with open(output[0], "w") as fh:
            fh.write(format_report("MIND FILTER", before_samples, after_samples, before_snps, after_snps))

rule qc_geno:
    input:
        before_fam = "plink_sub_cr/subjects.fam",
        after_fam = "plink_geno/subjects.fam",
        before_bim = "plink_sub_cr/subjects.bim",
        after_bim = "plink_geno/subjects.bim"
    output:
        "qc_stats/geno.txt"
    run:
        os.makedirs("qc_stats", exist_ok=True)
        before_samples = count_lines(input.before_fam)
        after_samples = count_lines(input.after_fam)
        before_snps = count_lines(input.before_bim)
        after_snps = count_lines(input.after_bim)

        with open(output[0], "w") as fh:
            fh.write(format_report("GENO FILTER", before_samples, after_samples, before_snps, after_snps))

rule qc_maf:
    input:
        before_fam = "plink_geno/subjects.fam",
        after_fam = "plink_maf/subjects.fam",
        before_bim = "plink_geno/subjects.bim",
        after_bim = "plink_maf/subjects.bim"
    output:
        "qc_stats/maf.txt"
    run:
        os.makedirs("qc_stats", exist_ok=True)
        before_samples = count_lines(input.before_fam)
        after_samples = count_lines(input.after_fam)
        before_snps = count_lines(input.before_bim)
        after_snps = count_lines(input.after_bim)

        with open(output[0], "w") as fh:
            fh.write(format_report("MAF FILTER", before_samples, after_samples, before_snps, after_snps))

rule qc_hwe:
    input:
        before_fam = "plink_maf/subjects.fam",
        after_fam = "plink_snp_cr/subjects.fam",
        before_bim = "plink_maf/subjects.bim",
        after_bim = "plink_snp_cr/subjects.bim"
    output:
        "qc_stats/hwe.txt"
    run:
        os.makedirs("qc_stats", exist_ok=True)
        before_samples = count_lines(input.before_fam)
        after_samples = count_lines(input.after_fam)
        before_snps = count_lines(input.before_bim)
        after_snps = count_lines(input.after_bim)

        with open(output[0], "w") as fh:
            fh.write(format_report("HWE FILTER", before_samples, after_samples, before_snps, after_snps))

rule qc_hrc:
    input:
        before_fam = "plink_snp_cr/subjects.fam",
        after_fam = "plink_snp_cr/subjects-updated.fam",
        before_bim = "plink_snp_cr/subjects.bim",
        after_bim = "plink_snp_cr/subjects-updated.bim"
    output:
        "qc_stats/hrc.txt"
    run:
        os.makedirs("qc_stats", exist_ok=True)
        before_samples = count_lines(input.before_fam)
        after_samples = count_lines(input.after_fam)
        before_snps = count_lines(input.before_bim)
        after_snps = count_lines(input.after_bim)

        with open(output[0], "w") as fh:
            fh.write(format_report("HRC FILTER", before_samples, after_samples, before_snps, after_snps))

rule qc_final:
    input:
        fam = "plink_snp_cr/subjects-updated.fam",
        bim = "plink_snp_cr/subjects-updated.bim"
    output:
        "qc_stats/final.txt"
    run:
        os.makedirs("qc_stats", exist_ok=True)
        samples = count_lines(input.fam)
        snps = count_lines(input.bim)

        with open(output[0], "w") as fh:
            fh.write("FINAL DATASET (post-HRC)\n")
            fh.write(f"Samples: {samples}\n")
            fh.write(f"SNPs: {snps}\n")
            fh.write("\n")

# ---------------------------------------------------------------------
# Final manifest
# ---------------------------------------------------------------------
rule write_manifest:
    input:
        initial_fam = f"{plinkIn}.fam",
        initial_bim = f"{plinkIn}.bim",
        liftover = "qc_stats/liftover.txt",
        mind = "qc_stats/mind.txt",
        geno = "qc_stats/geno.txt",
        maf = "qc_stats/maf.txt",
        hwe = "qc_stats/hwe.txt",
        hrc = "qc_stats/hrc.txt",
        final = "qc_stats/final.txt"
    output:
        "run_manifest.txt"
    shell:
        r"""
        set -euo pipefail
        {{
            echo "### RUN MANIFEST ###"
            echo ""

            echo "Date (UTC): $(date -u)"
            echo ""

            echo "User: $USER"
            echo ""

            echo "Working Dir: $(pwd)"
            echo ""

            echo "MODE: {MODE}"
            echo ""

            echo "Processing CHROMS: {CHROMS}"
            echo ""

            echo "PLINK input prefix: {plinkIn}"
            echo ""

            echo "Strand/Chain file: {chainStrandIn}"
            echo ""

            echo "Reference file: {ref_file}"
            echo ""

            echo "QC thresholds:"
            echo "  mind={MIND}"
            echo "  geno={GENO}"
            echo "  maf={MAF}"
            echo "  hwe={HWE}"
            echo ""

            echo "### INITIAL DATA ###"
            echo ""

            echo "Samples: $(wc -l < {input.initial_fam})"
            echo "SNPs: $(wc -l < {input.initial_bim})"
            echo ""

            echo "### QC SUMMARY ###"
            echo ""

            cat {input.liftover}
            cat {input.mind}
            cat {input.geno}
            cat {input.maf}
            cat {input.hwe}
            cat {input.hrc}
            cat {input.final}
        }} > {output}
        """
