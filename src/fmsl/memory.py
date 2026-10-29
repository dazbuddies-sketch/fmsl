"""Memory and learning: FMSL starts each session where the last one left off.

Two parts, both off by default (FMSL's answers never depend on them):

1. Proven results. Every proven answer is saved with the FMSL version that proved it; asking the
   same question again returns it at once. A new version starts a fresh memory, so a fix in the
   engine is never hidden by an old answer.

2. Strategy learning (a bandit, the simplest form of reinforcement learning). Where FMSL has a
   choice that does not affect the answer but can change the time a proof takes by a large factor
   (for example which leading direction the squaring map uses for a distance in three variables),
   it records, for each kind of problem, how each choice did: reward 1/(1 + seconds) when the proof
   succeeds within its time slice, 0 when it does not. Choices are then tried best first (UCB1:
   the mean reward plus an exploration bonus, untried choices first). Learning only orders the
   attempts; every answer is still proven by the engine, so a bad choice costs time, never
   correctness.

    import fmsl.memory as memory
    memory.enable()                      # ~/.cache/fmsl/memory.json, or memory.enable("file.json")
    ...                                  # use fmsl as usual
    memory.stats()                       # what has been learned

Setting the environment variable FMSL_MEMORY to a file path enables it at import.
"""
import json, math, os, threading

__all__ = ["enable", "disable", "enabled", "stats", "lookup", "store", "order", "reward"]

_LOCK = threading.Lock()
_STATE = {"path": None, "data": None}


def _version():
    from . import __version__
    return __version__


def enable(path=None):
    path = path or os.path.join(os.path.expanduser("~"), ".cache", "fmsl", "memory.json")
    data = {"version": _version(), "results": {}, "strategy": {}}
    try:
        with open(path) as fh:
            old = json.load(fh)
        if old.get("version") == data["version"]:
            data = old
        else:                                       # a new engine: keep what was learned about speed only
            data["strategy"] = old.get("strategy", {})
    except (OSError, ValueError):
        pass
    _STATE.update(path=path, data=data)
    return path


def disable():
    _STATE.update(path=None, data=None)


def enabled():
    return _STATE["data"] is not None


def _save():
    path, data = _STATE["path"], _STATE["data"]
    if not path or data is None: return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh: json.dump(data, fh)
    os.replace(tmp, path)


def _key(kind, expr, variables, point=None, side=None):
    return json.dumps([kind, expr.replace(" ", ""), list(variables or []), point, side], sort_keys=True, default=str)


def lookup(kind, expr, variables, point=None, side=None):
    """a saved proven result (a copy), or None"""
    if not enabled(): return None
    with _LOCK:
        r = _STATE["data"]["results"].get(_key(kind, expr, variables, point, side))
    return json.loads(json.dumps(r)) if r is not None else None


def store(kind, expr, variables, result, point=None, side=None):
    """save a result if it is proven"""
    if not enabled() or not result or not result.get("proven"): return
    with _LOCK:
        _STATE["data"]["results"][_key(kind, expr, variables, point, side)] = result
        _save()


def order(bucket, n_actions, c=0.3):
    """the actions 0..n-1, best first, for this kind of problem (UCB1; untried first, in order)"""
    if not enabled() or n_actions <= 1: return list(range(n_actions))
    with _LOCK:
        st = _STATE["data"]["strategy"].get(bucket, {})
    total = sum(v[0] for v in st.values()) or 1
    def score(a):
        n, s = st.get(str(a), (0, 0.0))
        if n == 0: return (1, -a)                   # untried: first, in the default order
        return (0, s / n + c * math.sqrt(math.log(total + 1) / n))
    return sorted(range(n_actions), key=score, reverse=True)


def reward(bucket, action, value):
    if not enabled(): return
    with _LOCK:
        st = _STATE["data"]["strategy"].setdefault(bucket, {})
        n, s = st.get(str(action), (0, 0.0))
        st[str(action)] = (n + 1, s + float(value))
        _save()


def stats():
    if not enabled(): return {"enabled": False}
    d = _STATE["data"]
    return {"enabled": True, "path": _STATE["path"], "version": d["version"], "results": len(d["results"]),
            "strategy": {b: {a: {"tries": n, "mean_reward": round(s / n, 3)} for a, (n, s) in st.items()}
                         for b, st in d["strategy"].items()}}


if os.environ.get("FMSL_MEMORY"):
    enable(os.environ["FMSL_MEMORY"])
