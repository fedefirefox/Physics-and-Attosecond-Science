"""Validate CFF 1.2.0 metadata; optional first argument is a local schema file."""
import json
from pathlib import Path
import sys
from urllib.request import urlopen

import jsonschema
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCHEMA_URL = "https://raw.githubusercontent.com/citation-file-format/citation-file-format/1.2.0/schema.json"
if len(sys.argv) > 1:
    schema = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
else:
    with urlopen(SCHEMA_URL, timeout=30) as response:
        schema = json.load(response)

contents = [(ROOT / name).read_text(encoding="utf-8")
            for name in ("CITATION.cff", "SFEXIS/CITATION.cff")]
if contents[0] != contents[1]:
    raise SystemExit("Root and SFEXIS citation metadata differ.")
data = yaml.safe_load(contents[0])
jsonschema.Draft7Validator(schema, format_checker=jsonschema.FormatChecker()).validate(data)
paper = data["preferred-citation"]
expected_authors = ["Vismarra", "Borrego-Varillas", "Wu", "Mocci", "Nisoli", "Lucchini"]
if ([a["family-names"] for a in paper["authors"]] != expected_authors
        or paper["doi"] != "10.1088/2515-7647/ac7991"
        or str(paper["start"]) != "034006" or str(paper["volume"]) != "4"
        or paper["year"] != 2022 or "doi" in data):
    raise SystemExit("Unexpected article metadata or a paper DOI assigned to software.")
print("PASS: both CFF files match and validate against CFF 1.2.0; article metadata verified.")
