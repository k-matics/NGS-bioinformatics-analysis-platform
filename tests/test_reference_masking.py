"""DRAGEN-style ALT masking: the reference check must verify the FASTA, not trust the sidecar."""
import json
import random

import pysam

from app.services.genome_build import classify_contigs
from app.services.reference_service import inspect_reference, verify_masking


def _write(tmp_path, mask_really: bool):
    rnd = random.Random(1)
    seqs = {"chr1": "".join(rnd.choice("ACGT") for _ in range(9000)),
            "chr1_KI270706v1_alt": "".join(rnd.choice("ACGT") for _ in range(6000))}
    bed = tmp_path / "mask.bed"
    bed.write_text("chr1_KI270706v1_alt\t1000\t5000\n")
    if mask_really:
        s = seqs["chr1_KI270706v1_alt"]
        seqs["chr1_KI270706v1_alt"] = s[:1000] + "N" * 4000 + s[5000:]
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{c}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n"
                             for c, s in seqs.items()))
    pysam.faidx(str(fasta))
    (tmp_path / "ref.fa.masking.json").write_text(json.dumps(
        {"method": "test", "mask_bed": str(bed), "mask_bed_copy": str(bed), "intervals": 1, "masked_bp": 4000}))
    return fasta


def test_masking_verified(tmp_path):
    fasta = _write(tmp_path, mask_really=True)
    m = verify_masking(fasta, {"chr1": 9000, "chr1_KI270706v1_alt": 6000})
    assert m["verified"] and m["checked"] == 1
    rep = inspect_reference(fasta)
    assert rep["alt_masking"]["verified"]
    assert not any("NOT ALT-aware" in w for w in rep["warnings"])


def test_sidecar_without_masking_is_an_error(tmp_path):
    fasta = _write(tmp_path, mask_really=False)
    rep = inspect_reference(fasta)
    assert rep["alt_masking"]["verified"] is False
    assert rep["status"] == "error"


def test_no_sidecar_keeps_alt_warning(tmp_path):
    fasta = _write(tmp_path, mask_really=True)
    (tmp_path / "ref.fa.masking.json").unlink()
    rep = inspect_reference(fasta)
    assert "alt_masking" not in rep
    assert any("NOT ALT-aware" in w for w in rep["warnings"])


def test_1000g_decoys_classified_as_decoy():
    c = classify_contigs({"chr1": 1, "chrUn_JTFH01000001v1_decoy": 1, "chrUn_KI270302v1": 1, "chrEBV": 1,
                          "HLA-A*01:01:01:01": 1})
    assert c["decoy"] == 2 and c["unplaced"] == 1 and c["hla"] == 1
