"""
fmsl.matrix -- safe matrix functions for quantum dynamics (requires numpy).

Built on FMSL's scalar results:
  phi_k(z) = sum_{j>=0} z^j / (j+k)!      (the exact Laurent coefficients FMSL derives for
                                           (e^z - 1 - z - ... - z^(k-1)/(k-1)!) / z^k)
  phi_0 = exp, phi_1(z) = (e^z - 1)/z, phi_2(z) = (e^z - 1 - z)/z^2, ...

For a Hermitian Hamiltonian H = V diag(E) V^dagger and time t, with a = -i E t:
  propagator(H, t)          U = exp(-iHt)            = V diag(e^a) V^dagger
  phi(k, H, t)              phi_k(-iHt)              = V diag(phi_k(a)) V^dagger
  frechet(H, dH, t)         d/de exp(-i(H + e dH)t)  = V (Gamma o (V^dagger (-i t dH) V)) V^dagger
                            Gamma_ij = (e^a_i - e^a_j)/(a_i - a_j) = e^a_j * phi_1(a_i - a_j)
The last form is the safe one: it never divides by a small eigenvalue gap, so exact and
near-degenerate energy levels need no special-casing.

Scalar kernel: on |z| <= 1 an exact-coefficient series with a PROVEN error bound
(proven_phi_bound). Elsewhere both the series (|z| <= 8) and the recurrence
phi_k = (phi_{k-1} - 1/(k-1)!)/z are evaluated with running error analysis (Higham)
and the more accurate result is kept per element.
"""
from fractions import Fraction as Fr
from math import factorial

import numpy as np

RHO_PROVEN = 1.0   # radius of the zone with a PROVEN error bound (see proven_phi_bound)
RHO_SERIES = 8.0   # series considered up to here; 60 terms leave a tail < 8^60/60! ~ 2e-28
KTERMS = 60
KPROVEN = 24       # terms used on the proven zone |z| <= 1 (tail < 1/(24+k)! < 1e-24)
_U = 2.0 ** -53
_COEF = {}


def _coef(k):
    if k not in _COEF:
        _COEF[k] = [float(Fr(1, factorial(j + k))) for j in range(KTERMS)]
    return _COEF[k]


def _series(k, z, K=KTERMS):
    c = _coef(k)[:K]
    acc = np.full(z.shape, c[-1], dtype=complex)
    for cj in reversed(c[:-1]):
        acc = acc * z + cj
    # typical-case Horner error ~ 2u * sum |c_j| |z|^j (= 2u * phi_k(|z|)); same footing as the
    # recurrence estimate so the two can be compared fairly (worst case is in proven_phi_bound)
    r = np.abs(z)
    S = np.full(z.shape, c[-1])
    for cj in reversed(c[:-1]):
        S = S * r + cj
    return acc, 2 * _U * S


def _recurrence(k, z):
    p = np.exp(z)
    err = 2 * _U * np.abs(p)                        # exp assumed accurate to ~2 ulp
    for j in range(1, k + 1):
        num = p - 1.0 / factorial(j - 1)
        err = err + _U * np.abs(num)                # subtraction rounding
        p = num / z
        err = err / np.abs(z) + 2 * _U * np.abs(p)  # division (complex) rounding
    return p, err


def phi_scalar(k, z, return_estimate=False):
    """Vectorized phi_k(z) for complex arrays z (k >= 0).
       For each element, both the exact-coefficient series (|z| <= 8) and the recurrence are
       evaluated with a running error estimate, and the more accurate one is returned.
       |z| <= 1 always uses the series, which carries a PROVEN bound (proven_phi_bound)."""
    z = np.asarray(z, dtype=complex)
    if k == 0:
        out = np.exp(z)
        return (out, 2 * _U * np.ones(z.shape)) if return_estimate else out
    out = np.empty_like(z)
    est = np.empty(z.shape)
    r = np.abs(z)
    inner = r <= RHO_PROVEN
    mid = (r > RHO_PROVEN) & (r <= RHO_SERIES)
    outer = r > RHO_SERIES
    if inner.any():
        v, e = _series(k, z[inner], KPROVEN); out[inner] = v; est[inner] = e / np.abs(v)
    if mid.any():
        zs = z[mid]
        vs, es = _series(k, zs)
        vr, er = _recurrence(k, zs)
        es, er = es / np.abs(vs), er / np.abs(vr)
        use_s = es <= er
        out[mid] = np.where(use_s, vs, vr); est[mid] = np.where(use_s, es, er)
    if outer.any():
        v, e = _recurrence(k, z[outer]); out[outer] = v; est[outer] = e / np.abs(v)
    return (out, est) if return_estimate else out


def proven_phi_bound(k, rho=Fr(1), K=KPROVEN):
    """PROVEN relative-error bound of the series branch of phi_scalar on ALL complex |z| <= rho.
       Truncation: sum_{j>=K} |z|^j/(j+k)!  <=  rho^K/(K+k)! / (1 - rho/(K+k+1)).
       Lower bound: |phi_k(z)| >= 1/k! - sum_{j>=1} rho^j/(j+k)!  (dominant term).
       Rounding: coefficient rounding + Higham's Horner bound gamma_{2K}, with a conservative
       complex unit roundoff of 4u (covers complex multiply and add without FMA).
       Returns a Fraction (rigorous), or None if the dominant-term lower bound fails."""
    if k < 1: raise ValueError("series branch is used for k >= 1")
    rho = Fr(rho)
    u = 4 * Fr(1, 2 ** 53)
    T = rho ** K / factorial(K + k) / (1 - rho / (K + k + 1))
    S = sum(rho ** j / factorial(j + k) for j in range(K))
    L = Fr(1, factorial(k)) - sum(rho ** j / factorial(j + k) for j in range(1, K)) - T
    if L <= 0: return None
    g = 2 * K * u / (1 - 2 * K * u)
    return (1 + T / L) * (1 + (u + g * (1 + u)) * S / L) - 1


# ------------------------------ matrix functions ------------------------------
def _eigh(H):
    H = np.asarray(H)
    if not np.allclose(H, H.conj().T, atol=1e-13 * max(1.0, np.abs(H).max())):
        raise ValueError("H must be Hermitian")
    return np.linalg.eigh(H)


def propagator(H, t, eig=None):
    """U = exp(-i H t) for Hermitian H."""
    E, V = eig if eig is not None else _eigh(H)
    return (V * np.exp(-1j * E * t)) @ V.conj().T


def phi(k, H, t, eig=None):
    """phi_k(-i H t) for Hermitian H (exponential-integrator matrix functions)."""
    E, V = eig if eig is not None else _eigh(H)
    return (V * phi_scalar(k, -1j * E * t)) @ V.conj().T


def gamma_matrix(E, t):
    """Divided-difference matrix Gamma_ij = e^{a_j} phi_1(a_i - a_j), a = -iEt (safe form)."""
    a = -1j * np.asarray(E) * t
    return np.exp(a)[None, :] * phi_scalar(1, a[:, None] - a[None, :])


def frechet(H, dH, t, eig=None):
    """Derivative of U = exp(-iHt) with respect to H in direction dH (Daleckii-Krein formula).
       Exact and near-degenerate eigenvalues are handled without any threshold."""
    E, V = eig if eig is not None else _eigh(H)
    Vh = V.conj().T
    X = Vh @ (-1j * t * np.asarray(dH)) @ V
    return V @ (gamma_matrix(E, t) * X) @ Vh


def frechet_many(H, dHs, t):
    """Gradients for several control directions sharing one eigen-decomposition."""
    E, V = _eigh(H)
    Vh = V.conj().T
    G = gamma_matrix(E, t)
    return [V @ (G * (Vh @ (-1j * t * np.asarray(d)) @ V)) @ Vh for d in dHs]
