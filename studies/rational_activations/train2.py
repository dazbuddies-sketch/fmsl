"""Experiment 2 (PREREGISTRATION_2.md): unsafe Pade activation units, degrees (5, 4),
R(z) = (a0 + ... + a5 z^5) / (1 + b1 z + ... + b4 z^4), initialised at Table 6 of
arXiv:1907.06732 (LeakyReLU, slope 0.2).  Same network, task and logging as train.py; the log
stores Q's full coefficient list (1, b1, ..., b4) so analyze.py's exact test applies unchanged.
Run:  python train2.py --gradcheck ;  python train2.py
"""
import os, sys
import numpy as np
import train as T
_base_init = T.init

A0 = np.array([0.02557776, 0.66182815, 1.58182975, 2.94478759, 0.95287794, 0.23319681])
B0 = np.array([0.50962605, 4.18376890, 0.37832090, 0.32407314])          # b1..b4
LRS = [0.002, 0.01, 0.03, 0.1, 0.3]
OUT = os.path.join(os.path.dirname(__file__), "runs2")


def rational(z, a, b):
    P = np.zeros_like(z); 
    for c in a[::-1]: P = P * z + c
    S = np.zeros_like(z)
    for c in b[::-1]: S = (S + c) * z                                     # b1 z + ... + b4 z^4
    return P, 1 + S


def act_back(d_h, z, P, Q, a, b):
    dP = np.zeros_like(z)
    for i in range(len(a) - 1, 0, -1): dP = dP * z + i * a[i]
    dQ = np.zeros_like(z)
    for k in range(len(b), 0, -1): dQ = dQ * z + k * b[k - 1]
    dz = d_h * (dP * Q - P * dQ) / Q ** 2
    ga = np.array([np.sum(d_h * z ** i / Q) for i in range(len(a))])
    gb = np.array([np.sum(-d_h * P * z ** k / Q ** 2) for k in range(1, len(b) + 1)])
    return dz, ga, gb


def init(seed):
    p = _base_init(seed)
    for L in ("1", "2"): p["a" + L] = A0.copy(); p["b" + L] = B0.copy()
    return p


def forward(p, x):
    z1 = x @ p["W1"] + p["c1"]; P1, Q1 = rational(z1, p["a1"], p["b1"]); h1 = P1 / Q1
    z2 = h1 @ p["W2"] + p["c2"]; P2, Q2 = rational(z2, p["a2"], p["b2"]); h2 = P2 / Q2
    return h2 @ p["W3"] + p["c3"], (x, z1, P1, Q1, h1, z2, P2, Q2, h2)


def loss_and_grads(p, x, y):
    out, (x, z1, P1, Q1, h1, z2, P2, Q2, h2) = forward(p, x)
    err = out - y; loss = float(np.mean(err ** 2)); d_out = 2 * err / x.shape[0]
    g = {"W3": h2.T @ d_out, "c3": d_out.sum(0)}
    d_z2, g["a2"], g["b2"] = act_back(d_out @ p["W3"].T, z2, P2, Q2, p["a2"], p["b2"])
    g["W2"] = h1.T @ d_z2; g["c2"] = d_z2.sum(0)
    d_z1, g["a1"], g["b1"] = act_back(d_z2 @ p["W2"].T, z1, P1, Q1, p["a1"], p["b1"])
    g["W1"] = x.T @ d_z1; g["c1"] = d_z1.sum(0)
    return loss, g


# reuse the training loop of train.py with this unit
T.rational, T.init, T.forward, T.loss_and_grads, T.OUT = rational, init, forward, loss_and_grads, OUT
_orig_run = T.run
def run(lr, seed):
    log = _orig_run(lr, seed)
    fn = os.path.join(OUT, f"lr{lr}_seed{seed}.npz")
    d = dict(np.load(fn))
    for L in ("1", "2"):                                   # store Q's full coefficients 1, b1..b4
        d["b" + L] = np.concatenate([np.ones((len(d["b" + L]), 1)), d["b" + L]], axis=1)
    np.savez_compressed(fn, **d)
    return log


if __name__ == "__main__":
    if "--gradcheck" in sys.argv:
        sys.exit(0 if T.gradcheck() < 1e-5 else 1)
    for lr in LRS:
        for s in T.SEEDS:
            log = run(lr, s); L = log["loss"]
            print(f"lr={lr:<6} seed={s}  steps={len(L)}  final loss={L[-1]:.4g}  max loss after 100={max(L[100:]) if len(L) > 100 else float('nan'):.4g}", flush=True)
