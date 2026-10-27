"""Evaluate a model on LimitBench, and the built-in reference methods.

  python evaluate.py --method paths                 the common "limits along 8 paths" method
  python evaluate.py --method fmsl                  FMSL itself (sanity check: must be 100%)
  python evaluate.py --method anthropic --model M   a Claude model via the Anthropic API
                                                    (needs ANTHROPIC_API_KEY in the environment)
Options: --limit N (N items spread evenly over the set), --workers K (parallel API calls, default 4).
Writes results/<method>[_<model>].jsonl and prints a summary table (also appended to results/SUMMARY.md).
"""
import argparse, json, math, os, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
from grade import prompt_for, grade


def load(limit=None):
    items = [json.loads(l) for l in open(os.path.join(HERE, "limitbench.jsonl"))]
    if not limit or limit >= len(items): return items
    step = len(items) / limit                         # an even spread over all tiers, deterministic
    return [items[int(k * step)] for k in range(limit)]


def answer_paths(item):
    if item["path_test"] is None: return None
    vals = item["path_test"]["values"]
    if len(vals) == 1 and vals[0] not in ("+oo", "-oo"): return f"ANSWER: LIMIT {vals[0]}"
    return "ANSWER: SET " + " u ".join("{" + v + "}" for v in vals)


def answer_fmsl(item):
    if len(item["variables"]) == 2:
        from fmsl.multivar import profile
        st = profile(item["f"])["cluster"]["signed_text"]
    else:
        from fmsl.space import profile3
        st = profile3(item["f"])["signed_text"]
    st = st.replace(" (algebraic)", "")
    if st.startswith("{") and st.count("{") == 1 and " u " not in st and "oo" not in st:
        return "ANSWER: LIMIT " + st[1:-1]
    return "ANSWER: SET " + st


import threading
USAGE = {"input_tokens": 0, "output_tokens": 0}
_LOCK = threading.Lock()


def answer_anthropic(item, model, max_tokens):
    import anthropic
    client = anthropic.Anthropic()
    for attempt in range(5):
        try:
            msg = client.messages.create(model=model, max_tokens=max_tokens,
                                         messages=[{"role": "user", "content": prompt_for(item)}])
            with _LOCK:
                USAGE["input_tokens"] += getattr(msg.usage, "input_tokens", 0) or 0
                USAGE["output_tokens"] += getattr(msg.usage, "output_tokens", 0) or 0
            return "".join(b.text for b in msg.content if getattr(b, "type", "") == "text")
        except Exception as exc:                      # rate limits, overload: back off and retry
            if attempt == 4: return f"ERROR: {type(exc).__name__}: {exc}"
            time.sleep(10 * (attempt + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True, choices=["paths", "fmsl", "anthropic"])
    ap.add_argument("--model", default="claude-sonnet-5-5")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-tokens", type=int, default=8000)
    a = ap.parse_args()
    items = load(a.limit)
    if a.method == "paths": f = answer_paths
    elif a.method == "fmsl": f = answer_fmsl
    else: f = lambda it: answer_anthropic(it, a.model, a.max_tokens)
    with ThreadPoolExecutor(max_workers=a.workers if a.method == "anthropic" else 1) as ex:
        answers = list(ex.map(f, items))
    errors = [a for a in answers if isinstance(a, str) and a.startswith("ERROR:")]
    if errors:
        # failed API calls are not answers: never score them; stop if the run did not really happen
        print(f"{len(errors)} of {len(answers)} calls failed, e.g. {errors[0][:300]}")
        if len(errors) > len(answers) // 10:
            sys.exit("Too many failed calls: nothing recorded. Fix the cause (e.g. API credits) and rerun.")
    rows = []
    for it, ans in zip(items, answers):
        if ans is None or (isinstance(ans, str) and ans.startswith("ERROR:")): continue
        g = grade(it, ans)
        rows.append({"id": it["id"], "tier": it["tier"], "category": it["category"], "trap": bool(it["path_test"] and it["path_test"]["trap"]),
                     **g, "answer_tail": ans[-400:]})
    name = a.method + (f"_{a.model}" if a.method == "anthropic" else "")
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", name + ".jsonl"), "w") as fh:
        for r in rows: fh.write(json.dumps(r) + "\n")
    def pct(sub, k): return f"{100 * sum(r[k] for r in sub) / len(sub):.1f}%" if sub else "-"
    lines = [f"### {name}  ({len(rows)} items)", "",
             "| subset | n | verdict right | exact answer | unparsed |", "|---|---|---|---|---|"]
    subsets = [("all", rows)] + [(f"tier {t}", [r for r in rows if r["tier"] == t]) for t in (1, 2, 3, 4)] + \
              [("path-test traps", [r for r in rows if r["trap"]])]
    for nm, sub in subsets:
        if sub:
            lines.append(f"| {nm} | {len(sub)} | {pct(sub, 'verdict')} | {pct(sub, 'exact')} | "
                         f"{sum(1 for r in sub if not r['parsed'])} |")
    if a.method == "anthropic":
        lines += ["", f"Tokens used: {USAGE['input_tokens']:,} input, {USAGE['output_tokens']:,} output."]
    text = "\n".join(lines) + "\n"
    print(text)
    with open(os.path.join(HERE, "results", "SUMMARY.md"), "a") as fh: fh.write(text + "\n")


if __name__ == "__main__":
    main()
