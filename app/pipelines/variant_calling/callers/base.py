"""Variant-caller plug-in interface.

Every caller is a subclass of ``VariantCaller`` and registers itself with
``@register_caller``. A caller declares:

* ``params``          - its OWN parameter set (ParamSpec). Parameters are not
                        shared between callers on purpose: "minimum base
                        quality" or "minimum allele fraction" are implemented
                        differently by each tool, so they are set per caller.
* ``build_commands``  - the argv list(s) that produce a bgzipped raw VCF.
* ``filter_rules``    - soft filters applied after normalisation (FILTER column is
                        set, records are never deleted). Each caller defines its own
                        rules because callers annotate different INFO/FORMAT fields.
* ``mode``            - "germline" now; "somatic" callers (e.g. Mutect2) will be
                        separate plug-ins with tumour/normal inputs - germline and
                        somatic logic are never mixed.

Adding a caller = one new file in this package + an import in __init__.py.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ....engine.params import ParamSpec


@dataclass
class CallJob:
    bam: Path
    reference: Path
    regions_bed: Path | None     # padded calling regions (0-based BED)
    out_vcf: Path                # bgzipped raw VCF to produce
    work_dir: Path
    threads: int
    sample_id: str
    params: dict                 # this caller's validated params
    label: str = ""              # unique label (run id) e.g. for container names


@dataclass
class CallCommand:
    commands: list               # argv or list of argv (a pipe)
    stdout_path: Path | None = None
    env: dict | None = None
    cleanup: list | None = None  # argv run (errors ignored) if the command fails or is cancelled


@dataclass
class FilterRule:
    """One soft filter: records matching ``expr`` (bcftools -e syntax) get FILTER=name."""
    name: str
    expr: str
    description: str


def tool_defaults_param(defaults_text: str) -> ParamSpec:
    """'Raw' switch: run the tool with its OWN defaults. Only what is needed to run is passed
    (reference, BAM, output, the calling regions) plus options that change speed/memory but
    never the calls (threads, temp folder, Java heap). The other caller options on the form
    are then not used; the filter options still are (set them to 0/off for an unfiltered VCF)."""
    return ParamSpec("tool_defaults", "Raw command: use the tool's own defaults", "bool", False,
                     "Ignore this caller's options above and run the bare command. Tool defaults: "
                     + defaults_text + " The exact command is in the run's command log.", group="Mode")


def common_filter_params(min_qual: float, min_depth: int = 10, min_vaf: float = 0.15,
                         min_gq: int | None = 20) -> list[ParamSpec]:
    """Filter parameters shared in meaning (not in value) by all germline callers."""
    specs = [
        ParamSpec("filter_min_qual", "Filter: min QUAL", "float", min_qual,
                  "Records with QUAL below this get FILTER=LowQual (0 = off). QUAL scales differ "
                  "between callers, so the default differs too.", min=0, max=10000, group="Filtering"),
        ParamSpec("filter_min_depth", "Filter: min read depth (FORMAT/DP)", "int", min_depth,
                  "FILTER=LowDP below this depth (0 = off).", min=0, max=100000, group="Filtering"),
        ParamSpec("filter_min_vaf", "Filter: min variant allele fraction", "float", min_vaf,
                  "FILTER=LowVAF when ALT reads / (REF+ALT reads) is below this (0 = off). "
                  "Germline hets are expected near 0.5.", min=0, max=1, group="Filtering"),
    ]
    if min_gq is not None:
        specs.append(ParamSpec("filter_min_gq", "Filter: min genotype quality (GQ)", "int", min_gq,
                               "FILTER=LowGQ below this GQ (0 = off).", min=0, max=99, group="Filtering"))
    return specs


def common_filter_rules(p: dict) -> list[FilterRule]:
    rules = [FilterRule("NoAltGT", 'GT="RR" || GT="mis"',
                        "Genotype is homozygous reference or missing (no variant in this sample)")]
    if p.get("filter_min_qual"):
        rules.append(FilterRule("LowQual", f"QUAL<{p['filter_min_qual']}",
                                f"QUAL < {p['filter_min_qual']}"))
    if p.get("filter_min_depth"):
        rules.append(FilterRule("LowDP", f"FMT/DP<{p['filter_min_depth']}",
                                f"FORMAT/DP < {p['filter_min_depth']}"))
    if p.get("filter_min_gq"):
        rules.append(FilterRule("LowGQ", f"FMT/GQ<{p['filter_min_gq']}", f"GQ < {p['filter_min_gq']}"))
    if p.get("filter_min_vaf"):
        v = p["filter_min_vaf"]
        rules.append(FilterRule(
            "LowVAF", f"(FMT/AD[0:0]+FMT/AD[0:1])>0 && FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])<{v}",
            f"ALT allele fraction (AD) < {v}"))
    return rules


class VariantCaller:
    key: str = ""
    name: str = ""
    tool: str = ""               # executable recorded for provenance/version
    mode: str = "germline"
    requires_read_group: bool = False
    description: str = ""
    params: list[ParamSpec] = []

    def available(self) -> tuple[bool, str]:
        """Can this caller run on this server right now? (tool installed, image present...)"""
        return True, ""

    def build_commands(self, job: CallJob) -> list[CallCommand]:
        raise NotImplementedError

    def filter_rules(self, params: dict) -> list[FilterRule]:
        """Soft filters for this caller (applied to the normalised VCF)."""
        return common_filter_rules(params)

    def version_info(self) -> dict:
        """Extra provenance (e.g. container image digest)."""
        return {}


_CALLERS: dict[str, VariantCaller] = {}


def register_caller(cls):
    _CALLERS[cls.key] = cls()
    return cls


def get_caller(key: str) -> VariantCaller:
    if key not in _CALLERS:
        raise KeyError(f"Unknown variant caller '{key}'")
    return _CALLERS[key]


def all_callers(mode: str | None = None) -> list[VariantCaller]:
    return [c for c in _CALLERS.values() if mode is None or c.mode == mode]
