"""
Shared setup for the fast unit tests.

The workflow scripts live in hybracter/workflow/scripts/ and are normally run by
Snakemake's `script:` directive (which injects a `snakemake` object). Each script
guards its top-level call with `if "snakemake" in globals():`, so importing them
here is side-effect free and we can unit-test their functions directly.
"""

import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "hybracter" / "workflow" / "scripts"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
