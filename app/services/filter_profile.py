"""Extract the filtering profile of a VCF and turn it into a parameter preset.

"Profile" = everything the file itself says about how it was filtered, from two sources:

1. **Header** – the variant caller and its command line (##DRAGENCommandLine, ##GATKCommandLine,
   ##bcftools_*Command, ##commandline, ##source), and every ##FILTER definition. Many
   definitions contain the rule itself (DRAGEN: "Set if true:QUAL < 0.4"); these are parsed.
2. **Records** – which FILTER labels were actually applied and how often (a header often
   defines filters that were never used, e.g. somatic filters in a germline VCF), and the
   *observed* lowest QUAL / DP / GQ / VAF among the PASS calls. The observed values are the
   effective thresholds, whatever the header claims.

The mapping onto our per-caller soft filters is explicit and documented per value in the
returned ``notes``: each is "mapped" (same meaning), "approximate" (similar, different scale)
or "not mappable". Nothing is guessed silently.
"""
from __future__ import annotations

import re
import shlex
from collections import Counter
from pathlib import Path

import pysam

RULE_RE = re.compile(r"\b(QUAL|DP|GQ|AF|VAF|MQ|QD|FS|SOR|VQL|JGQ|MQRankSum|ReadPosRankSum)\s*(<=|>=|<|>|==)\s*"
                     r"(-?\d+(?:\.\d+)?)")
MAX_RECORDS = 5_000_000


def _caller_from_header(lines: list[str]) -> str:
    text = "\n".join(lines)
    if "##DRAGENCommandLine" in text or "##source=DRAGEN" in text:
        return "DRAGEN"
    if "DeepVariant" in text:
        return "DeepVariant"
    if "##GATKCommandLine" in text and "HaplotypeCaller" in text:
        return "GATK HaplotypeCaller"
    if "##bcftools_callCommand" in text:
        return "bcftools"
    if "freebayes" in text.lower():
        return "FreeBayes"
    return "unknown"


def _command_lines(header: pysam.VariantHeader) -> list[dict]:
    out = []
    for rec in header.records:
        key = rec.key or ""
        if "command" not in key.lower() and key not in ("source", "commandline"):
            continue
        if rec.type == "GENERIC":
            out.append({"key": key, "id": None, "command": str(rec.value)})
        else:
            items = dict(rec.items())
            out.append({"key": key, "id": items.get("ID"), "version": items.get("Version"),
                        "command": items.get("CommandLineOptions") or items.get("CommandLine") or str(items)})
    return out


def _options(cmd: str) -> dict:
    """--option value pairs from a command line (best effort, for display and a few checks)."""
    try:
        toks = shlex.split(cmd)
    except ValueError:
        toks = cmd.split()
    opts = {}
    for i, t in enumerate(toks):
        if t.startswith("--"):
            nxt = toks[i + 1] if i + 1 < len(toks) and not toks[i + 1].startswith("--") else True
            opts[t] = nxt
    return opts


def _gatk_variantfiltration(cmds: list[dict]) -> list[dict]:
    rules = []
    for c in cmds:
        if c.get("id") != "VariantFiltration":
            continue
        toks = shlex.split(c["command"]) if c["command"] else []
        exprs = [toks[i + 1] for i, t in enumerate(toks) if t in ("--filter-expression", "-filter") and i + 1 < len(toks)]
        names = [toks[i + 1] for i, t in enumerate(toks) if t in ("--filter-name",) and i + 1 < len(toks)]
        for n, e in zip(names, exprs):
            rules.append({"filter": n, "expression": e, "source": "GATK VariantFiltration"})
    return rules


def _pct(values: list[float], q: float):
    if not values:
        return None
    values = sorted(values)
    return round(values[min(len(values) - 1, int(q * len(values)))], 4)


def extract_profile(vcf: Path, max_records: int = MAX_RECORDS) -> dict:
    vcf = Path(vcf)
    with pysam.VariantFile(str(vcf)) as fh:
        header = fh.header
        header_lines = str(header).splitlines()
        cmds = _command_lines(header)
        filters = []
        for fid, f in header.filters.items():
            desc = f.description or ""
            parsed = [{"field": m.group(1), "op": m.group(2), "value": float(m.group(3))}
                      for m in RULE_RE.finditer(desc)]
            filters.append({"id": fid, "description": desc, "rules": parsed})
        samples = list(header.samples)

        counts: Counter = Counter()
        combos: Counter = Counter()
        obs = {"QUAL": [], "DP": [], "GQ": [], "VAF": []}
        n = 0
        for rec in fh:
            n += 1
            if n > max_records:
                raise ValueError(f"{vcf.name}: more than {max_records:,} records; profile extraction stopped "
                                 "(no silent sampling).")
            fl = [k for k in rec.filter.keys()] or ["."]
            combos[";".join(fl)] += 1
            for k in fl:
                counts[k] += 1
            if fl != ["PASS"]:
                continue
            if rec.qual is not None:
                obs["QUAL"].append(rec.qual)
            if samples:
                s = rec.samples[0]
                dp = s.get("DP") if "DP" in rec.format else None
                if dp is None and "DP" in rec.info:
                    dp = rec.info["DP"]
                if dp is not None:
                    obs["DP"].append(dp)
                if "GQ" in rec.format and s.get("GQ") is not None:
                    obs["GQ"].append(s.get("GQ"))
                ad = s.get("AD") if "AD" in rec.format else None
                if ad and len(ad) >= 2 and None not in ad[:2] and sum(ad) > 0:
                    obs["VAF"].append(max(ad[1:]) / sum(ad))

    used = {k for k in counts if k not in ("PASS", ".")}
    for f in filters:
        f["records"] = counts.get(f["id"], 0)
        f["used"] = f["id"] in used or f["id"] == "PASS"
    observed = {k: {"n": len(v), "min": _pct(v, 0), "p01": _pct(v, 0.01), "p05": _pct(v, 0.05),
                    "median": _pct(v, 0.5)} for k, v in obs.items()}
    caller = _caller_from_header(header_lines)
    main_cmd = next((c for c in cmds if c["key"] == "DRAGENCommandLine" and c.get("id") == "dragen"), None)
    opts = _options(main_cmd["command"]) if main_cmd else {}
    facts = []
    if caller == "DRAGEN":
        facts.append({"key": "--vc-target-bed",
                      "value": opts.get("--vc-target-bed", "not set"),
                      "meaning": ("Variant calling was restricted to target regions."
                                  if "--vc-target-bed" in opts else
                                  "Not set: DRAGEN called the whole genome, including off-target reads. "
                                  "Compare only inside the target regions.")})
        for k in sorted(opts):
            if k.startswith("--vc-") and k != "--vc-target-bed":
                facts.append({"key": k, "value": opts[k], "meaning": ""})
    return {"vcf": str(vcf), "caller": caller, "samples": samples, "records": n,
            "filter_counts": dict(counts.most_common()), "filter_combinations": dict(combos.most_common(20)),
            "filters": filters, "gatk_filters": _gatk_variantfiltration(cmds),
            "command_lines": [{**c, "command": c["command"][:4000]} for c in cmds],
            "facts": facts, "observed_pass": observed}


# ---------------------------------------------------------------------------------- mapping
def _rule(profile: dict, field: str) -> tuple[dict, dict] | None:
    """A germline filter with a rule on `field`: one that was applied to records first, else one that is
    defined in the header but never triggered (e.g. DRAGEN's QUAL < 0.4 hard filter)."""
    found = []
    for f in profile["filters"]:
        if f["id"].startswith(("Mosaic", "Targeted", "Recombinant")):
            continue
        for r in f["rules"]:
            if r["field"] == field and r["op"] in ("<", "<="):
                found.append((not f["used"], f, r))
    if not found:
        return None
    found.sort(key=lambda x: x[0])
    _, f, r = found[0]
    return f, r


def _desc(f: dict) -> str:
    used = f"{f['records']:,} records" if f["used"] else "defined in the header, never triggered"
    return f"{f['id']}: '{f['description']}' ({used})"


def profile_to_params(profile: dict) -> tuple[dict, list[dict]]:
    """Suggested preset parameters + notes [{param, value, status, reason}]."""
    notes: list[dict] = []
    obs = profile["observed_pass"]
    cp: dict = {"bcftools": {}, "gatk_hc": {}, "freebayes": {}, "deepvariant": {}}

    def note(param, value, status, reason):
        notes.append({"param": param, "value": value, "status": status, "reason": reason})

    # --- depth
    r = _rule(profile, "DP")
    if r:
        f, rule = r
        v = int(rule["value"]) + (1 if rule["op"] == "<=" else 0)
        note("filter_min_depth", v, "mapped", f"{_desc(f)} → our LowDP keeps DP ≥ {v}.")
    else:
        v = 0
        note("filter_min_depth", v, "mapped", f"The VCF applies no depth filter (lowest PASS DP observed: "
             f"{obs['DP']['min']}) → LowDP off.")
    for c in cp:
        cp[c]["filter_min_depth"] = v

    # --- genotype quality
    r = _rule(profile, "GQ")
    gq = int(r[1]["value"]) if r else 0
    note("filter_min_gq", gq, "mapped",
         (_desc(r[0]) if r else
          f"No germline GQ filter is used (lowest PASS GQ observed: {obs['GQ']['min']}). "
          "DRAGEN's TargetedLowGQ applies only to its targeted-repeat caller.") + f" → LowGQ {'≥ ' + str(gq) if gq else 'off'}.")
    for c in ("bcftools", "gatk_hc", "deepvariant"):
        cp[c]["filter_min_gq"] = gq

    # --- allele fraction
    r = _rule(profile, "AF") or _rule(profile, "VAF")
    vaf = float(r[1]["value"]) if r else 0.0
    mosaic = next((f for f in profile["filters"] if f["id"] == "MosaicLowAF" and f["used"]), None)
    note("filter_min_vaf", vaf, "mapped",
         (_desc(r[0]) if r else
          f"No germline allele-fraction filter (lowest PASS VAF observed: {obs['VAF']['min']})")
         + (f". MosaicLowAF ({mosaic['records']:,} records) applies to DRAGEN's mosaic calls only and is not mapped"
            if mosaic else "") + f" → LowVAF {'≥ ' + str(vaf) if vaf else 'off'}.")
    for c in cp:
        cp[c]["filter_min_vaf"] = vaf

    # --- QUAL (different scale per caller!)
    r = _rule(profile, "QUAL")
    q = float(r[1]["value"]) if r else 0.0
    why = (_desc(r[0]) if r else "No QUAL filter")
    note("filter_min_qual (bcftools, FreeBayes)", q, "approximate",
         f"{why}. QUAL is computed differently by every caller, so the same number is not the same confidence; "
         "used as the most permissive equivalent.")
    for c in ("bcftools", "freebayes"):
        cp[c]["filter_min_qual"] = q
    stand = max(q, 1.0)
    cp["gatk_hc"]["stand_call_conf"] = stand
    note("gatk_hc.stand_call_conf", stand, "approximate",
         f"GATK only emits calls above this QUAL (in-house 30). Lowered to {stand} so GATK also reports the "
         "low-confidence calls the source VCF passes; filter them afterwards if needed.")

    # --- caller-specific hard filters
    has_gatk_vf = bool(profile.get("gatk_filters"))
    cp["gatk_hc"]["filter_gatk_hard"] = has_gatk_vf
    note("gatk_hc.filter_gatk_hard", has_gatk_vf, "mapped",
         "The source VCF used GATK VariantFiltration." if has_gatk_vf else
         "The source VCF applies no QD/FS/MQ/SOR hard filters → off.")
    strand = any(f["used"] and any(r["field"] in ("FS", "SOR") for r in f["rules"]) for f in profile["filters"])
    cp["freebayes"]["filter_strand"] = strand
    note("freebayes.filter_strand", strand, "mapped",
         "A strand-bias filter is used in the source." if strand else "No strand-bias filter in the source → off.")

    # --- not mappable
    covered = {"QUAL", "DP", "GQ", "AF", "VAF"}
    for f in profile["filters"]:
        if not f["used"] or f["id"] == "PASS":
            continue
        fields = {x["field"] for x in f["rules"]}
        if f["id"].startswith(("Mosaic", "Targeted", "Recombinant")):
            note(f["id"], None, "not mapped",
                 f"'{f['description']}' ({f['records']:,} records): applies to DRAGEN's mosaic/targeted callers, "
                 "which have no equivalent here.")
        elif fields and fields <= covered:
            continue   # same field as a mapped rule above
        else:
            note(f["id"], None, "not mappable",
                 f"'{f['description']}' ({f['records']:,} records): no equivalent in our filters.")
    note("min_MQ / min_BQ / caller emission", None, "not mappable",
         "Read-level settings of the source caller are not recorded in the VCF; in-house values are kept.")
    return {"caller_params": cp}, notes
