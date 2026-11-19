"""Sample intake validation workflow.

Purpose: before any processing, confirm that every input file is intact,
correctly formatted, and on the same genome build as our reference.
This is the first gate in a clinical pipeline: garbage-in must be caught
here, not discovered after a 6-hour alignment.
"""
from ...engine.registry import WorkflowDefinition, register
from .steps import (
    BamValidationStep,
    ChecksumStep,
    FastqIntegrityStep,
    IntakeReportStep,
    ReferenceCheckStep,
    VcfValidationStep,
)


def build_steps(params: dict):
    return [
        ReferenceCheckStep(),
        ChecksumStep(),
        FastqIntegrityStep(),
        BamValidationStep(),
        VcfValidationStep(),
        IntakeReportStep(),
    ]


def check_inputs(sample, files):
    return [] if files else ["No input files registered for this sample."]


register(WorkflowDefinition(
    key="intake_validation",
    name="Sample intake validation",
    version="1.0",
    description=("Reference check, SHA-256 checksums, FASTQ gzip/format/pairing checks, "
                 "BAM/VCF header and genome-build compatibility checks."),
    build_steps=build_steps,
    default_params={"fastq_check_records": 100_000},
    check_inputs=check_inputs,
    estimate_disk_gb=lambda sample, files: 0.05,
))
