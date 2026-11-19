from pathlib import Path

from ....engine.params import ParamSpec
from ....services.settings_service import get_setting
from ....services.tool_service import find_gatk_java_home, gatk_env, resolve
from .base import (CallCommand, CallJob, FilterRule, VariantCaller, common_filter_params,
                   common_filter_rules, register_caller, tool_defaults_param)


@register_caller
class GatkHaplotypeCaller(VariantCaller):
    key = "gatk_hc"
    name = "GATK HaplotypeCaller"
    tool = "gatk"
    requires_read_group = True
    description = ("Local de-novo re-assembly of haplotypes in active regions; the GATK Best Practices "
                   "germline caller. Strong for indels. Output here is single-sample, unfiltered.")
    params = [
        tool_defaults_param("min call confidence 30, read filter MAPQ ≥ 20, base quality ≥ 10, "
                            "PCR indel model CONSERVATIVE."),
        ParamSpec("stand_call_conf", "Min call confidence (QUAL)", "float", 30.0,
                  "--standard-min-confidence-threshold-for-calling", min=0, max=1000),
        ParamSpec("min_MQ", "Min mapping quality", "int", 20, "--minimum-mapping-quality", min=0, max=60),
        ParamSpec("min_BQ", "Min base quality", "int", 10, "--min-base-quality-score", min=0, max=60),
        ParamSpec("pcr_indel_model", "PCR indel model", "choice", "CONSERVATIVE",
                  "--pcr-indel-model: use NONE for PCR-free libraries.",
                  choices=[("CONSERVATIVE", "CONSERVATIVE (PCR libraries, e.g. exome)"),
                           ("NONE", "NONE (PCR-free)"), ("AGGRESSIVE", "AGGRESSIVE"), ("HOSTILE", "HOSTILE")]),
        ParamSpec("dragen_mode", "DRAGEN-GATK mode", "bool", False,
                  "--dragen-mode: DRAGEN's genotyping model inside GATK (base-quality dropout, foreign-read "
                  "detection, DRAGEN MAPQ transform, posterior-based QUAL, soft-clipping of low-quality ends, "
                  "dynamic read disqualification, 150 bp indel padding). Values set explicitly on this form "
                  "(min MAPQ, min call confidence, base quality, PCR model) still win over the mode's own."),
        ParamSpec("dragstr_calibrate", "DRAGEN mode: calibrate DragSTR on this BAM", "bool", True,
                  "With DRAGEN mode: models indel errors in short tandem repeats from this BAM "
                  "(ComposeSTRTableFile once per reference, cached outside the reference folder, then "
                  "CalibrateDragstrModel per run -> --dragstr-params-path). Off = GATK's default STR "
                  "parameters. Adds ~10-30 min (+ a one-off table build per reference)."),
        ParamSpec("java_xmx_gb", "Java heap (GB)", "int", 8, "", min=2, max=48, advanced=True),
    ] + common_filter_params(min_qual=0) + [
        ParamSpec("filter_gatk_hard", "Filter: GATK Best-Practices hard filters", "bool", True,
                  "SNPs: QD<2, FS>60, MQ<40, SOR>3, MQRankSum<-12.5, ReadPosRankSum<-8. "
                  "Indels: QD<2, FS>200, SOR>10, ReadPosRankSum<-20. (Used when VQSR is not possible, "
                  "e.g. a single exome.)", group="Filtering"),
    ]

    def filter_rules(self, params):
        rules = common_filter_rules(params)
        if params.get("filter_gatk_hard"):
            rules.append(FilterRule(
                "HardFilterSNP",
                'TYPE="snp" && (INFO/QD<2 || INFO/FS>60 || INFO/MQ<40 || INFO/SOR>3 || '
                'INFO/MQRankSum<-12.5 || INFO/ReadPosRankSum<-8)',
                "GATK SNP hard filter (QD<2|FS>60|MQ<40|SOR>3|MQRankSum<-12.5|ReadPosRankSum<-8)"))
            rules.append(FilterRule(
                "HardFilterIndel",
                'TYPE="indel" && (INFO/QD<2 || INFO/FS>200 || INFO/SOR>10 || INFO/ReadPosRankSum<-20)',
                "GATK indel hard filter (QD<2|FS>200|SOR>10|ReadPosRankSum<-20)"))
        return rules

    def available(self):
        gatk = resolve("gatk", get_setting("gatk_path"))
        if not gatk:
            return False, "GATK launcher not found (Settings → GATK launcher path)."
        if not find_gatk_java_home(get_setting("gatk_java_home", "") or None):
            return False, "No Java 17-21 found for GATK."
        return True, ""

    def build_commands(self, job: CallJob):
        p = job.params
        gatk, env = get_setting("gatk_path"), gatk_env(get_setting("gatk_java_home", "") or None)
        java = ["--java-options", f"-Xmx{p['java_xmx_gb']}g"]
        regions = ["-L", job.regions_bed, "--interval-padding", "0"] if job.regions_bed else []
        cmds: list[CallCommand] = []
        argv = [gatk, *java, "HaplotypeCaller", "-R", job.reference, "-I", job.bam, "-O", job.out_vcf,
                "--native-pair-hmm-threads", str(min(job.threads, 8))]
        if not p.get("tool_defaults"):
            if p.get("dragen_mode"):
                argv += ["--dragen-mode"]          # first, so the explicit values below win
                if p.get("dragstr_calibrate", True):
                    table = dragstr_table_path(job.reference)
                    if not table.exists():          # one-off per reference; built beside, then renamed
                        table.parent.mkdir(parents=True, exist_ok=True)
                        tmp = table.with_name(table.name + ".building.zip")
                        cmds.append(CallCommand([gatk, *java, "ComposeSTRTableFile", "-R", job.reference,
                                                 "-O", tmp, "--tmp-dir", job.work_dir], env=env))
                        cmds.append(CallCommand(["mv", "-n", tmp, table]))
                    dparams = Path(job.work_dir) / "dragstr_params.txt"
                    cmds.append(CallCommand([gatk, *java, "CalibrateDragstrModel", "-R", job.reference,
                                             "-I", job.bam, "-str", table, "-O", dparams,
                                             "--tmp-dir", job.work_dir, *regions], env=env))
                    argv += ["--dragstr-params-path", dparams]
            argv += ["--standard-min-confidence-threshold-for-calling", str(p["stand_call_conf"]),
                     "--minimum-mapping-quality", str(p["min_MQ"]),
                     "--min-base-quality-score", str(p["min_BQ"])]
            if not (p.get("dragen_mode") and p.get("dragstr_calibrate", True)):
                argv += ["--pcr-indel-model", p["pcr_indel_model"]]   # DragSTR replaces the PCR indel model
        argv += ["--tmp-dir", job.work_dir, *regions]   # padding already applied to the regions
        cmds.append(CallCommand(argv, env=env))
        return cmds


def dragstr_table_path(reference) -> Path:
    """STR table of a reference for DragSTR, cached OUTSIDE the reference folder (never modified)."""
    import hashlib

    from flask import current_app
    ref = Path(reference)
    base = Path(get_setting("dragstr_dir", "") or Path(current_app.instance_path) / "dragstr")
    key = hashlib.sha1(str(ref.resolve()).encode()).hexdigest()[:8]
    return base / f"{ref.stem}_{key}.str_table.zip"
