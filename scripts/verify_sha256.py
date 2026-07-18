#!/usr/bin/env python3
import csv, hashlib, sys
from pathlib import Path

index = Path(sys.argv[1])
base = index.parent.parent
failed = 0
with index.open(encoding="utf-8-sig", newline="") as handle:
    for row in csv.DictReader(handle):
        value = row.get("relative_path") or row.get("file_name") or ""
        expected = row.get("sha256", "").upper()
        path = base / value
        if not path.exists():
            print(f"MISSING {value}")
            failed += 1
            continue
        h = hashlib.sha256(path.read_bytes()).hexdigest().upper()
        status = "PASS" if h == expected else "FAIL"
        print(f"{status} {value}")
        failed += status == "FAIL"
raise SystemExit(1 if failed else 0)
