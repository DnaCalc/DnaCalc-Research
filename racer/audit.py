"""Reconcile both backend runs, scoped identities, first kills and reachability.

All Excel input/answer/provenance fields come from the existing witness trace.
The union of model predictions is a reachability constraint, never a pass rate.
"""
from pathlib import Path
import collections
import hashlib
import json
import sys


def read(path):
    return [json.loads(line) for line in path.open()]


def collected(directory, pattern):
    records = {}
    for path in sorted(directory.glob(pattern)):
        for record in read(path):
            assert record["order"] not in records
            records[record["order"]] = record
    return records


def audit(directory):
    root = Path(directory)
    rows = read(root / "frozen/scoped.jsonl")
    native = root / "native"
    software = root / "software"
    ns = json.loads((native / "summary.json").read_text())
    ss = json.loads((software / "summary.json").read_text())
    assert {k: v for k, v in ns.items() if k not in ["backend", "reachability_seconds"]} == {
        k: v for k, v in ss.items() if k not in ["backend", "reachability_seconds"]}
    n = collected(native, "reachability-*.jsonl")
    s = collected(software, "reachability-*.jsonl")
    assert n == s
    assert set(n) == set(range(len(rows)))
    nr = collected(native, "unreachable-roots-*.jsonl")
    sr = collected(software, "unreachable-roots-*.jsonl")
    assert nr == sr
    assert len(nr) == ns["whole_family_unreachable"]
    for file in ["first-kills.jsonl", "programs.jsonl", "survivors.json", "witness-kills.jsonl"]:
        assert (native / file).read_bytes() == (software / file).read_bytes(), file
    kills = read(native / "first-kills.jsonl")
    assert len(kills) == len({k["candidate_id"] for k in kills}) == 39366
    assert {k["candidate_id"] for k in kills} == set(range(39366))
    counts = collections.Counter()
    by_source = collections.defaultdict(collections.Counter)
    for k in kills:
        r = rows[k["order"]]
        assert k["witness"] == r["raw"]
        assert k["trace_index"] == r["trace_index"]
        assert k["expected"] == r["raw"]["expected"] != k["got"]
        counts[(k["order"], r["raw"]["id"])] += 1
    statistics = read(native / "witness-kills.jsonl")
    assert len(statistics) == len(rows)
    remaining = 39366
    with (root / "unreachable-witnesses.jsonl").open("w") as gaps, (root / "reachability-ordered.jsonl").open("w") as reachable:
        for order, r in enumerate(rows):
            record = n[order]
            assert record["trace_index"] == r["trace_index"]
            assert record["id"] == r["raw"]["id"] and record["source"] == r["raw"]["source"]
            assert record["expected"] == r["raw"]["expected"]
            assert record["reachable"] == (record["expected"] in record["outputs"])
            stats = statistics[order]
            assert stats["killed"] == counts[(order, r["raw"]["id"])]
            remaining -= stats["killed"]
            assert stats["remaining"] == remaining
            assert stats["id"] == r["raw"]["id"] and stats["trace_index"] == r["trace_index"]
            assert stats["reading"] == r["raw"]["reading"]
            by_source[r["raw"]["source"]]["scoped"] += 1
            by_source[r["raw"]["source"]]["unreachable"] += not record["reachable"]
            reachable.write(json.dumps(record, separators=(",", ":")) + "\n")
            if not record["reachable"]:
                assert order in nr
                gaps.write(json.dumps({"order": order, "witness": r["raw"], "frozen_leaves": r["leaves"],
                    "prior_gap": r["prior_gap"], "family_outputs": record["outputs"]}, separators=(",", ":")) + "\n")
    assert remaining == 0
    summary = {"native_software_prediction_differences": 0, "native_software_raw_root_differences": 0,
               "scoped_rows_reconciled": len(rows), "programs_killed_once": len(kills),
               "zero_kill_witnesses": len(rows) - len(counts),
               "nonzero_kill_witnesses": [{"order": o, "id": ident, "killed": count} for (o, ident), count in counts.items() if count],
               "unreachable_by_source": dict(by_source),
               "unreachable_witnesses_sha256": hashlib.sha256((root / "unreachable-witnesses.jsonl").read_bytes()).hexdigest(),
               "no_predictions_scored_as_excel_answers": True}
    # Accesses to Counter during reconciliation add zero entries; recount directly.
    summary["zero_kill_witnesses"] = sum(s["killed"] == 0 for s in statistics)
    (root / "audit-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    audit(sys.argv[1])
