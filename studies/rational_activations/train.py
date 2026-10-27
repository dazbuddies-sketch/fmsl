"""Train small MLPs with the rational activation of Boulle, Nakatsukasa & Townsend (2020) and log,
at every step, everything the pre-registered analysis needs (see PREREGISTRATION.md).

R(z) = P(z)/Q(z), P = a0 + a1 z + a2 z^2 + a3 z^3, Q = b0 + b1 z + b2 z^2, one R per hidden layer,
initialised at the paper's Table 1 coefficients.  Exact gradients in NumPy (checked against
finite differences by `python train.py --gradcheck`).

Run:  python train.py            (all 40 runs, writes runs/lr{lr}_seed{s}.npz)
"""
import os, sys
import numpy as np

A0 = np.array([1.1915, 1.5957, 0.5000, 0.0218])     # numerator, Table 1 of arXiv:2004.01902
B0 = np.array([2.3830, 0.0000, 1.0000])             # denominator
WIDTH, STEPS, BATCH = 32, 3000, 128
LRS = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
SEEDS = list(range(8))
OUT = os.path.join(os.path.dirname(__file__), "runs")


def data():
    rng = np.random.default_rng(12345)
    x = rng.uniform(-2, 2, size=(2000, 1))
    y = np.sin(3 * x) + 0.3 * x ** 2
    return x, y


def init(seed):
    rng = np.random.default_rng(seed)
    def glorot(n_in, n_out):
        lim = np.sqrt(6 / (n_in + n_out)); return rng.uniform(-lim, lim, size=(n_in, n_out))
    return {"W1": glorot(1, WIDTH), "c1": np.zeros(WIDTH), "W2": glorot(WIDTH, WIDTH), "c2": np.zeros(WIDTH),
            "W3": glorot(WIDTH, 1), "c3": np.zeros(1),
            "a1": A0.copy(), "b1": B0.copy(), "a2": A0.copy(), "b2": B0.copy()}


def rational(z, a, b):
    P = a[0] + z * (a[1] + z * (a[2] + z * a[3]))
    Q = b[0] + z * (b[1] + z * b[2])
    return P, Q


def forward(p, x):
    z1 = x @ p["W1"] + p["c1"]; P1, Q1 = rational(z1, p["a1"], p["b1"]); h1 = P1 / Q1
    z2 = h1 @ p["W2"] + p["c2"]; P2, Q2 = rational(z2, p["a2"], p["b2"]); h2 = P2 / Q2
    out = h2 @ p["W3"] + p["c3"]
    return out, (x, z1, P1, Q1, h1, z2, P2, Q2, h2)


def loss_and_grads(p, x, y):
    out, (x, z1, P1, Q1, h1, z2, P2, Q2, h2) = forward(p, x)
    n = x.shape[0]
    err = out - y
    loss = float(np.mean(err ** 2))
    d_out = 2 * err / n
    g = {"W3": h2.T @ d_out, "c3": d_out.sum(0)}
    d_h2 = d_out @ p["W3"].T
    def act_back(d_h, z, P, Q, a, b):
        dP = a[1] + z * (2 * a[2] + z * 3 * a[3]); dQ = b[1] + 2 * b[2] * z
        dz = d_h * (dP * Q - P * dQ) / Q ** 2
        ga = np.array([np.sum(d_h * z ** i / Q) for i in range(4)])
        gb = np.array([np.sum(-d_h * P * z ** j / Q ** 2) for j in range(3)])
        return dz, ga, gb
    d_z2, g["a2"], g["b2"] = act_back(d_h2, z2, P2, Q2, p["a2"], p["b2"])
    g["W2"] = h1.T @ d_z2; g["c2"] = d_z2.sum(0)
    d_h1 = d_z2 @ p["W2"].T
    d_z1, g["a1"], g["b1"] = act_back(d_h1, z1, P1, Q1, p["a1"], p["b1"])
    g["W1"] = x.T @ d_z1; g["c1"] = d_z1.sum(0)
    return loss, g


def gradcheck():
    x, y = data(); x, y = x[:64], y[:64]
    p = init(0)
    _, g = loss_and_grads(p, x, y)
    worst = 0.0
    rng = np.random.default_rng(1)
    for k in p:
        for _ in range(5):
            idx = tuple(rng.integers(0, s) for s in p[k].shape)
            old = p[k][idx]; h = 1e-6
            p[k][idx] = old + h; lp, _ = loss_and_grads(p, x, y)
            p[k][idx] = old - h; lm, _ = loss_and_grads(p, x, y)
            p[k][idx] = old
            num = (lp - lm) / (2 * h)
            rel = abs(num - g[k][idx]) / max(1e-8, abs(num) + abs(g[k][idx]))
            worst = max(worst, rel)
    print(f"gradient check: worst relative difference {worst:.2e}")
    return worst


def run(lr, seed):
    x, y = data()
    p = init(seed)
    m = {k: np.zeros_like(v) for k, v in p.items()}; v = {k: np.zeros_like(w) for k, w in p.items()}
    rng = np.random.default_rng(1000 + seed)
    log = {k: [] for k in ("loss", "gnorm", "a1", "b1", "a2", "b2", "r1", "r2")}
    with np.errstate(all="ignore"):
        for t in range(1, STEPS + 1):
            # full-set forward: loss and pre-activation ranges (what the analysis uses)
            out, (_, z1, _, _, _, z2, _, _, _) = forward(p, x)
            log["loss"].append(float(np.mean((out - y) ** 2)))
            log["r1"].append((float(np.min(z1)), float(np.max(z1))))
            log["r2"].append((float(np.min(z2)), float(np.max(z2))))
            for k in ("a1", "b1", "a2", "b2"): log[k].append(p[k].copy())
            idx = rng.integers(0, x.shape[0], size=BATCH)
            _, g = loss_and_grads(p, x[idx], y[idx])
            gn = float(np.sqrt(sum(np.sum(gg ** 2) for gg in g.values())))
            log["gnorm"].append(gn)
            if not np.isfinite(gn): break                 # parameters would become NaN
            for k in p:                                   # Adam
                m[k] = 0.9 * m[k] + 0.1 * g[k]; v[k] = 0.999 * v[k] + 0.001 * g[k] ** 2
                mh = m[k] / (1 - 0.9 ** t); vh = v[k] / (1 - 0.999 ** t)
                p[k] = p[k] - lr * mh / (np.sqrt(vh) + 1e-8)
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(os.path.join(OUT, f"lr{lr}_seed{seed}.npz"),
                        **{k: np.array(v) for k, v in log.items()}, lr=lr, seed=seed)
    return log


if __name__ == "__main__":
    if "--gradcheck" in sys.argv:
        sys.exit(0 if gradcheck() < 1e-5 else 1)
    which = sys.argv[1:] or None
    for lr in LRS:
        for s in SEEDS:
            log = run(lr, s)
            L = log["loss"]
            print(f"lr={lr:<6} seed={s}  steps={len(L)}  final loss={L[-1]:.4g}  max loss={max(L):.4g}", flush=True)
