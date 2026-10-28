"""Regrade saved model answers with the current grader (after a grader fix) and write
results/REGRADED.md. Only the grading changes; the stored answers are not touched otherwise."""
import glob, json, os
from grade import grade

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    items = {}
    for fn in ("limitbench.jsonl", "limitbench_gaps.jsonl"):
        for l in open(os.path.join(HERE, fn)):
            it = json.loads(l); items[it["id"]] = it
    out = ["# Model results, regraded with the current grader", ""]
    for p in sorted(glob.glob(os.path.join(HERE, "results", "*anthropic_*.jsonl"))):
        rows = [json.loads(l) for l in open(p)]
        for r in rows: r.update(grade(items[r["id"]], r["answer_tail"]))
        with open(p, "w") as fh:
            for r in rows: fh.write(json.dumps(r) + "\n")
        pct = lambda s, k: f"{100 * sum(r[k] for r in s) / len(s):.1f}%"
        out += [f"### {os.path.basename(p)[:-6]}  ({len(rows)} items)", "",
                "| subset | n | verdict right | exact answer | unparsed |", "|---|---|---|---|---|"]
        subs = [("all", rows)] + [(f"tier {t}", [r for r in rows if r["tier"] == t]) for t in range(1, 7)] + \
               [(c, [r for r in rows if r["category"] == c]) for c in sorted({r["category"] for r in rows}) if c.startswith(("gap", "decoy"))] + \
               [("path-test traps", [r for r in rows if r["trap"]])]
        for nm, s in subs:
            if s: out.append(f"| {nm} | {len(s)} | {pct(s, 'verdict')} | {pct(s, 'exact')} | {sum(not r['parsed'] for r in s)} |")
        wrong = [r["id"] for r in rows if not r["exact"]]
        out += ["", f"Not exact: {', '.join(wrong) or 'none'}", ""]
    text = "\n".join(out) + "\n"
    open(os.path.join(HERE, "results", "REGRADED.md"), "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
