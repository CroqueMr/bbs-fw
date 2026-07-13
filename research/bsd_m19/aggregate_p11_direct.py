#!/usr/bin/env python3
from __future__ import annotations
import glob, hashlib, json, os, sys

root, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
files = sorted(glob.glob(os.path.join(root, "**", "p11_direct_shard_*.json"), recursive=True))
if len(files) != 128:
    raise SystemExit(f"expected 128 shard files, found {len(files)}")
rows = [json.load(open(f, encoding="utf-8")) for f in files]
rows.sort(key=lambda r: r["shard"])

expected = list(range(128))
if [r["shard"] for r in rows] != expected:
    raise SystemExit("shard indices are not exactly 0..127")
keys = ["case", "engine", "nshards", "p", "q", "ell1", "ell2", "eta1", "eta2", "n"]
for key in keys:
    vals = {json.dumps(r[key], sort_keys=True) for r in rows}
    if len(vals) != 1:
        raise SystemExit(f"metadata mismatch for {key}: {vals}")
for i, r in enumerate(rows):
    exp_lo = 1 + (r["n"] - 1) * i // 128
    exp_hi = (r["n"] - 1) * (i + 1) // 128
    if (r["lo"], r["hi"]) != (exp_lo, exp_hi):
        raise SystemExit(f"interval mismatch at shard {i}")
    if i and rows[i-1]["hi"] + 1 != r["lo"]:
        raise SystemExit(f"gap or overlap before shard {i}")

residue = sum(r["residue"] for r in rows) % rows[0]["p"]
manifest = []
for f in files:
    data = open(f, "rb").read()
    manifest.append({"file": os.path.basename(f), "sha256": hashlib.sha256(data).hexdigest()})
out = {
    "case": rows[0]["case"],
    "engine": rows[0]["engine"],
    "independent_of_m17_table_numba_engine": True,
    "shards": 128,
    "coverage": [rows[0]["lo"], rows[-1]["hi"]],
    "expected_coverage": [1, rows[0]["n"] - 1],
    "p": rows[0]["p"],
    "q": rows[0]["q"],
    "ell1": rows[0]["ell1"],
    "ell2": rows[0]["ell2"],
    "eta1": rows[0]["eta1"],
    "eta2": rows[0]["eta2"],
    "n": rows[0]["n"],
    "delta_n_mod_p": residue,
    "nonzero": residue != 0,
    "total_nonzero_weights": sum(r["nonzero_weights"] for r in rows),
    "total_cpu_s": sum(r["elapsed_s"] for r in rows),
    "pari_versions": sorted({r["pari_version"] for r in rows}),
    "sage_versions": sorted({r["sage_version"] for r in rows}),
    "shard_manifest": manifest,
}
raw = json.dumps(out, indent=2, sort_keys=True).encode()
out["certificate_sha256_without_self"] = hashlib.sha256(raw).hexdigest()
path = os.path.join(outdir, "p11_direct_reproduction.json")
open(path, "w", encoding="utf-8").write(json.dumps(out, indent=2, sort_keys=True))
print(json.dumps(out, sort_keys=True))
