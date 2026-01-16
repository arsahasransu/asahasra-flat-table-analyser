
#!/usr/bin/env python3
"""
check_env.py — Verify that the current Python environment satisfies required dependencies.

- Validates package presence and minimum versions.
- Emits a structured report (human readable + JSON).
- Exits non-zero on failure (useful for CI).
"""

from __future__ import annotations
import json
import os
import sys
import platform
import importlib
import traceback
from dataclasses import dataclass
from typing import Optional, Dict, Any


try:
    from packaging.version import Version, InvalidVersion
except Exception:
    print("ERROR: 'packaging' is required (`pip install packaging`).", file=sys.stderr)
    sys.exit(2)

# Prefer Python 3.8+: importlib.metadata is stdlib (backport is 'importlib_metadata')
try:
    from importlib.metadata import version as get_dist_version, PackageNotFoundError
except Exception:
    try:
        from importlib_metadata import version as get_dist_version, PackageNotFoundError  # type: ignore
    except Exception:
        print("ERROR: 'importlib_metadata' is required for Python<3.8.", file=sys.stderr)
        sys.exit(2)

# ---------------------------------------------------------------------------
# Customize this block for your project:
# - key: import name (module), value: minimum acceptable version (or None for "presence only")
# - You can pin exact versions e.g. "==2.1.0" by encoding the logic in `meets_min_version`.
# ---------------------------------------------------------------------------

REQUIRED: Dict[str, Optional[str]] = {
    # numerics
    "numpy": "2.4.0",

    # ML / DL
    # "torch": "2.1.0",            # PyTorch

    # HEP / ROOT
    "ROOT": None,                # presence-only (PyROOT often isn't semantically versioned the same way)
    "awkward": "2.8.0",          # Awkward Array v2
    "pyarrow": "22.0.0",        # Apache Arrow / Parquet - Dependency of awkward

    "yaml": "6.0.0",
    "IPython": "9.7.0"
}

# ---------------------------------------------------------------------------

@dataclass
class CheckResult:
    name: str
    required: Optional[str]
    found: Optional[str]
    ok: bool
    note: Optional[str] = None
    extra: Dict[str, Any] = None

def import_module_safely(modname: str):
    try:
        return importlib.import_module(modname), None
    except Exception as e:
        return None, e

def get_version_from_metadata(dist_name: str) -> Optional[str]:
    try:
        return get_dist_version(dist_name)
    except PackageNotFoundError:
        # Some packages have different dist name vs import name (e.g., "Pillow" vs "PIL").
        return None
    except Exception:
        return None

def get_module_version(mod) -> Optional[str]:
    # try __version__ as a fallback
    return getattr(mod, "__version__", None)

def meets_min_version(found: Optional[str], required: Optional[str]) -> bool:
    if required is None:
        return found is not None  # presence-only
    if found is None:
        return False
    try:
        # Support basic ">=x.y.z" semantics
        req = Version(required)
        cur = Version(found)
        return cur >= req
    except InvalidVersion:
        # If version strings are non-standard, fallback to equality or presence
        return found == required

def check_one(name: str, min_version: Optional[str]) -> CheckResult:
    mod, err = import_module_safely(name)
    if mod is None:
        return CheckResult(name, min_version, None, False, note=f"ImportError: {err}")
    # Try metadata, then __version__
    found_ver = get_version_from_metadata(name) or get_module_version(mod)
    ok = meets_min_version(found_ver, min_version)
    note = None
    if not ok and min_version:
        note = f"Found {found_ver!r}, require >= {min_version}"
    return CheckResult(name, min_version, found_ver, ok, note=note)

def main() -> int:
    print("=" * 70)
    print("Environment check")
    print("=" * 70)
    print(f"Python  : {sys.version.split()[0]}  ({sys.executable})")
    print(f"Platform: {platform.system()} {platform.release()}  [{platform.machine()}]")
    print(f"Conda   : {os.environ.get('CONDA_DEFAULT_ENV', '<none>')}")
    print("-" * 70)

    results: list[CheckResult] = []
    for modname, minver in REQUIRED.items():
        r = check_one(modname, minver)
        status = "OK" if r.ok else "FAIL"
        req_str = f">={minver}" if minver else "present"
        print(f"{status:4}  {modname:12}  required: {req_str:10}  found: {str(r.found):10}  {r.note or ''}")
        results.append(r)

    ok_all = all(r.ok for r in results)
    print("-" * 70)
    print("SUMMARY:", "OK ✅" if ok_all else "FAILED ❌")

    # Emit a JSON artifact for CI or downstream tooling
    artifact = {
        "python": sys.version,
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "conda_env": os.environ.get("CONDA_DEFAULT_ENV"),
        "requirements": [
            {
                "name": r.name,
                "required": r.required,
                "found": r.found,
                "ok": r.ok,
                "note": r.note,
                "extra": r.extra,
            }
            for r in results
        ],
        "ok": ok_all,
    }
    with open("env_report.json", "w", encoding="utf-8") as f:
        json.dump(artifact, f, indent=2)
    print("Wrote env_report.json")

    return 0 if ok_all else 1

if __name__ == "__main__":
    sys.exit(main())
