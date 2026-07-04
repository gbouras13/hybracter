"""
Shared setup for the fast unit tests.

The workflow scripts live in hybracter/workflow/scripts/ and are normally run by
Snakemake's `script:` directive (which injects a `snakemake` object). Each script
guards its top-level call with `if "snakemake" in globals():`, so importing them
here is side-effect free and we can unit-test their functions directly.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "hybracter" / "workflow" / "scripts"

# scripts dir: import the workflow scripts as top-level modules
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

# repo root: import the installed-style `hybracter` package (e.g. for the CLI
# tests) without needing `pip install .` in the fast unit lane
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
