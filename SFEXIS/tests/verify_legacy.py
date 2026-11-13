"""Verify every archived original byte for byte (legacy files use -text)."""
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
manifest = json.loads((ROOT / "SFEXIS/docs/legacy-manifest.json").read_text())
archive = ROOT / manifest["destination_directory"]
expected = {entry["path"]: entry["git_blob_sha1"] for entry in manifest["files"]}
actual = {p.relative_to(archive).as_posix() for p in archive.rglob("*") if p.is_file()}
if actual != set(expected):
    raise SystemExit(f"Archive file list changed: missing={set(expected)-actual}, extra={actual-set(expected)}")
for relative, sha in expected.items():
    path = (archive / relative).relative_to(ROOT).as_posix()
    found = subprocess.check_output(
        ["git", "hash-object", "--no-filters", path], cwd=ROOT, text=True
    ).strip()
    if found != sha:
        raise SystemExit(f"Original content changed: {path}: {found} != {sha}")
print(f"PASS: all {len(expected)} original Git blobs preserved from {manifest['source_commit']}.")
