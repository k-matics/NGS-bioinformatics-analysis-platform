"""Variant-caller plug-ins. Importing this package registers all callers."""
from . import bcftools, deepvariant, freebayes, gatk  # noqa: F401
from .base import CallCommand, CallJob, VariantCaller, all_callers, get_caller  # noqa: F401
