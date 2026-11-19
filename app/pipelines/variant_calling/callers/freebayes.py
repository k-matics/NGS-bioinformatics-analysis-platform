from ....engine.params import ParamSpec
from .base import (CallCommand, CallJob, FilterRule, VariantCaller, common_filter_params,
                   common_filter_rules, register_caller, tool_defaults_param)


@register_caller
class FreeBayesCaller(VariantCaller):
    key = "freebayes"
    name = "FreeBayes"
    tool = "freebayes"
    description = ("Bayesian haplotype-based caller (variants in short windows). Single-threaded here, "
                   "so slower on a full exome. Reports complex/MNP alleles that are decomposed later.")
    params = [
        tool_defaults_param("MAPQ ≥ 1, base quality ≥ 0, min ALT fraction 0.05, min ALT reads 2, "
                            "min coverage 0."),
        ParamSpec("min_MQ", "Min mapping quality", "int", 20, "-m / --min-mapping-quality", min=0, max=60),
        ParamSpec("min_BQ", "Min base quality", "int", 20, "-q / --min-base-quality", min=0, max=60),
        ParamSpec("min_alt_fraction", "Min alternate allele fraction", "float", 0.2,
                  "-F: an allele must reach this fraction in a sample to be considered.", min=0.01, max=1),
        ParamSpec("min_alt_count", "Min alternate read count", "int", 2, "-C", min=1, max=1000),
        ParamSpec("min_coverage", "Min coverage", "int", 5, "--min-coverage", min=0, max=10000),
    ] + common_filter_params(min_qual=20, min_gq=None) + [
        ParamSpec("filter_strand", "Filter: ALT seen on one strand only", "bool", True,
                  "FILTER=StrandBias when all ALT reads are on one strand (INFO/SAF=0 or SAR=0). "
                  "Common FreeBayes artefact filter.", group="Filtering"),
    ]

    def filter_rules(self, params):
        rules = common_filter_rules(params)
        if params.get("filter_strand"):
            rules.append(FilterRule("StrandBias", "INFO/SAF=0 || INFO/SAR=0",
                                    "All ALT-supporting reads on one strand"))
        return rules

    def build_commands(self, job: CallJob):
        p = job.params
        fb = ["freebayes", "-f", job.reference]
        if not p.get("tool_defaults"):
            fb += ["-m", str(p["min_MQ"]), "-q", str(p["min_BQ"]),
                   "-F", str(p["min_alt_fraction"]), "-C", str(p["min_alt_count"]),
                   "--min-coverage", str(p["min_coverage"])]
        if job.regions_bed:
            fb += ["-t", job.regions_bed]
        fb.append(job.bam)
        return [CallCommand([fb, ["bgzip", "-c", "-@", str(min(job.threads, 4))]], stdout_path=job.out_vcf)]
