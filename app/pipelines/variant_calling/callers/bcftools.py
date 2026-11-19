from ....engine.params import ParamSpec
from .base import CallCommand, CallJob, VariantCaller, common_filter_params, register_caller, tool_defaults_param


@register_caller
class BcftoolsCaller(VariantCaller):
    key = "bcftools"
    name = "BCFtools mpileup/call"
    tool = "bcftools"
    description = ("Pileup-based caller with the multiallelic model (-m). Fast and transparent; "
                   "weaker than haplotype-based callers for indels in repeats.")
    params = [
        tool_defaults_param("mpileup MAPQ ≥ 0, base quality ≥ 13, max depth 250; call -m -v (multiallelic "
                            "caller, variant sites only — 'call' needs a model). AD/DP/GQ are still "
                            "requested: they are annotations only and do not change the calls."),
        ParamSpec("min_MQ", "Min mapping quality", "int", 20, "mpileup -q: ignore reads with MAPQ below this.",
                  min=0, max=60),
        ParamSpec("min_BQ", "Min base quality", "int", 20, "mpileup -Q: ignore bases below this Phred quality.",
                  min=0, max=60),
        ParamSpec("max_depth", "Max depth per file", "int", 1000, "mpileup -d: reads per position kept.",
                  min=50, max=100000),
        ParamSpec("variants_only", "Variant sites only", "bool", True, "call -v: do not output ref-only sites."),
    ] + common_filter_params(min_qual=20)

    def build_commands(self, job: CallJob):
        p = job.params
        raw = bool(p.get("tool_defaults"))
        mpileup = ["bcftools", "mpileup", "--threads", str(max(1, job.threads // 2)), "-f", job.reference]
        if not raw:
            mpileup += ["-q", str(p["min_MQ"]), "-Q", str(p["min_BQ"]), "-d", str(p["max_depth"])]
        mpileup += ["-a", "FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD", "-Ou"]
        if job.regions_bed:
            mpileup += ["-R", job.regions_bed]
        mpileup.append(job.bam)
        call = ["bcftools", "call", "--threads", str(max(1, job.threads // 2)), "-m", "-a", "FORMAT/GQ", "-Oz", "-o", job.out_vcf]
        if raw or p["variants_only"]:
            call.insert(4, "-v")
        return [CallCommand([mpileup, call])]
