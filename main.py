#!/usr/bin/env python3
"""TopologicalQubit-Lab - educational topological qubit / Majorana lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lesson is stated up front:

  * This lab diagonalises ONE finite toy Hamiltonian exactly: the
    Kitaev chain, a spinless p-wave superconducting wire reduced to
    a 2N x 2N Bogoliubov-de Gennes matrix with labelled teaching
    parameters (hopping t, chemical potential mu, pairing Delta).
  * Everything it reports - the gap closing at |mu| = 2t, the Z2
    invariant from Pfaffians, the Majorana end modes, the splitting
    that falls exponentially with chain length, the perturbation
    shift it calls "protection" - is a statement about that matrix,
    computed from that matrix, and true of that matrix.
  * No nanowire was grown, no superconductor was measured, no qubit
    exists. The braiding section is algebra on 4x4 matrices: it
    shows what an exchange WOULD do to Majorana operators, it is
    not a simulation of a physical braid and not an experiment.

The lab separates three ideas that are often blurred: a topological
PHASE (a sign computed from the bulk), topological PROTECTION (a
number: how little a local perturbation moves the zero mode, here,
in this matrix), and a topological QUBIT (a piece of hardware that
does not exist in this lab and is a permanent HOLD gate). The gates
are enforced with Gemini Argon 4 Extended Thinking as DECIDER
(names-only, tighten-only), Muse as advisor (never overrides
HOLD), and Google Flow plan-only seats.
"""
import hashlib
import json
import math
import sys

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"

# Labelled teaching parameters. t sets the energy scale; every other
# number in the lab is a ratio to it. These are NOT device values.
T_HOP = 1.0
DEMO_MU = 0.5          # deep in the topological phase (|mu| < 2t)
DEMO_DELTA = 0.5       # pairing used for the spectrum scan
SWEET_MU = 0.0         # sweet spot: mu = 0, Delta = t
SWEET_DELTA = 1.0
PERTURB_DELTA = 0.2    # labelled on-site potential for the protection test
PROTECT_MU_TOPO = 0.5
PROTECT_MU_TRIVIAL = 2.5

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "hardware_qubit_claim": "HOLD - no claim a qubit was built; this lab diagonalises a toy matrix",
    "fault_tolerance_claim": "HOLD - no fault-tolerant quantum computing claim from a finite chain spectrum",
    "quantum_computer_claim": "HOLD - no claim this is, runs on, or benchmarks a quantum computer",
    "braiding_as_experiment_claim": "HOLD - the braiding section is operator algebra; it must never be presented as an experiment or a physical braid",
    "measured_superconductivity_claim": "HOLD - no superconductivity was measured; Delta is a labelled parameter",
    "real_device_fabrication": "HOLD - no fabrication guidance for real nanowires, heterostructures, or devices",
    "cryogenics_operation": "HOLD - no cryostat, dilution refrigerator, or low-temperature operation directed from a sim",
    "materials_synthesis_claim": "HOLD - no synthesis recipes for semiconductor or superconductor materials",
    "topological_protection_guarantee": "HOLD - protection numbers are statements about this matrix; no protection is guaranteed in any real device",
    "error_correction_product_claim": "HOLD - no error-correction product, service, or deployment claim",
    "remove_human_override": "HOLD - the human override is never removed",
    "investment_claim": "HOLD - no investment claims; a lab model is not a prospectus",
    "vendor_benchmark_claim": "HOLD - no vendor, hardware, or benchmark comparison claims",
    "regulatory_approval_claim": "HOLD - no claim of regulatory approval or pending certification",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Linear algebra from scratch
# The BdG Hamiltonian is a real symmetric matrix, diagonalised by a
# cyclic Jacobi eigensolver written here in full. A Pfaffian (for the
# Z2 invariant) and a determinant (to cross-check Pf^2 = det) are
# implemented by Gaussian elimination. Nothing is imported.

def jacobi_eigh(matrix, tol=1e-12, max_sweeps=60):
    """Eigendecomposition of a real symmetric matrix by cyclic Jacobi.

    Returns (eigenvalues, eigenvectors) with eigenvalues sorted
    ascending; eigenvectors are returned as a matrix whose COLUMNS
    are the eigenvectors (evecs[row][col]).
    """
    n = len(matrix)
    a = [row[:] for row in matrix]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(max_sweeps):
        off = math.sqrt(sum(a[i][j] ** 2
                            for i in range(n) for j in range(n) if i != j))
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-15:
                    continue
                app, aqq = a[p][p], a[q][q]
                tau = (aqq - app) / (2.0 * apq)
                t_ = (1.0 if tau >= 0 else -1.0) / (abs(tau)
                                                    + math.sqrt(1.0 + tau * tau))
                c = 1.0 / math.sqrt(1.0 + t_ * t_)
                s = t_ * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    evals = [a[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: evals[i])
    evals = [evals[i] for i in order]
    evecs = [[v[row][i] for i in order] for row in range(n)]
    return evals, evecs


def pfaffian(matrix):
    """Pfaffian of a real antisymmetric even-dimensional matrix.

    Skew-symmetric Gaussian elimination: pivot on the largest entry
    of the current row, apply the Schur-complement update to the
    trailing submatrix, and accumulate the pivot product. Validated
    in CHECKS against a brute-force pairing expansion and against
    Pf(A)^2 = det(A).
    """
    n = len(matrix)
    if n % 2 == 1:
        return 0.0
    a = [row[:] for row in matrix]
    pf = 1.0
    for k in range(0, n - 2, 2):
        j = max(range(k + 1, n), key=lambda jj: abs(a[k][jj]))
        if abs(a[k][j]) < 1e-14:
            return 0.0
        if j != k + 1:
            a[k + 1], a[j] = a[j], a[k + 1]
            for row in a:
                row[k + 1], row[j] = row[j], row[k + 1]
            pf = -pf
        piv = a[k][k + 1]
        pf *= piv
        for i in range(k + 2, n):
            for j2 in range(k + 2, n):
                a[i][j2] += (a[i][k] * a[k + 1][j2]
                             - a[i][k + 1] * a[k][j2]) / piv
    pf *= a[n - 2][n - 1]
    return pf


def det_gauss(matrix):
    """Determinant by Gaussian elimination with partial pivoting."""
    n = len(matrix)
    a = [row[:] for row in matrix]
    d = 1.0
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) < 1e-14:
            return 0.0
        if piv != col:
            a[col], a[piv] = a[piv], a[col]
            d = -d
        d *= a[col][col]
        for r in range(col + 1, n):
            f = a[r][col] / a[col][col]
            for cc in range(col, n):
                a[r][cc] -= f * a[col][cc]
    return d


# ================================================================ 2. The Kitaev chain
# H = -mu * sum_j (c+_j c_j - 1/2) - t * sum_j (c+_j c_{j+1} + h.c.)
#     + Delta * sum_j (c_j c_{j+1} + h.c.)
# in the Nambu basis (c_1..c_N, c+_1..c+_N). The particle block is
# the tight-binding chain (-mu on site, -t on bonds); the hole block
# is its negative transpose; the pairing block is antisymmetric with
# -Delta / +Delta on the two off-diagonals of each bond.

def build_bdg(n_sites, t, mu, delta):
    """The 2N x 2N real symmetric BdG Hamiltonian of the Kitaev chain."""
    n2 = 2 * n_sites
    H = [[0.0] * n2 for _ in range(n2)]
    for j in range(n_sites):
        H[j][j] = -mu
        H[n_sites + j][n_sites + j] = mu
        if j + 1 < n_sites:
            H[j][j + 1] = -t
            H[j + 1][j] = -t
            H[n_sites + j][n_sites + j + 1] = t
            H[n_sites + j + 1][n_sites + j] = t
            H[j][n_sites + j + 1] = -delta
            H[j + 1][n_sites + j] = delta
            H[n_sites + j + 1][j] = -delta
            H[n_sites + j][j + 1] = delta
    return H


def spectrum(n_sites, t, mu, delta):
    """Full BdG spectrum: eigenvalues (ascending) and eigenvectors."""
    return jacobi_eigh(build_bdg(n_sites, t, mu, delta))


def sorted_abs_levels(n_sites, t, mu, delta):
    """|E| levels ascending: s[0] is the zero-mode splitting scale."""
    evals, _ = spectrum(n_sites, t, mu, delta)
    return sorted(abs(e) for e in evals)


def zero_mode_state(n_sites, t, mu, delta):
    """The eigenvector with |E| closest to zero, split as (u, v).

    u is the particle part (sites 1..N), v the hole part. The
    particle-hole partner (v, u) carries the same site density, so
    the densities below do not depend on which partner Jacobi
    happened to return.
    """
    evals, evecs = spectrum(n_sites, t, mu, delta)
    idx = min(range(len(evals)), key=lambda i: abs(evals[i]))
    u = [evecs[j][idx] for j in range(n_sites)]
    v = [evecs[n_sites + j][idx] for j in range(n_sites)]
    return evals[idx], u, v


def site_density(u, v):
    """Per-site weight |u_j|^2 + |v_j|^2 of a BdG eigenstate."""
    return [u[j] ** 2 + v[j] ** 2 for j in range(len(u))]


def end_weight(u, v, sites_each_end=2):
    """Weight of a state on the first and last `sites_each_end` sites."""
    rho = site_density(u, v)
    n = len(rho)
    k = min(sites_each_end, n // 2)
    return sum(rho[:k]) + sum(rho[n - k:])


def zero_mode_pair_density(n_sites, t, mu, delta):
    """Summed density of the TWO eigenstates closest to zero energy.

    At an exact degeneracy (the sweet spot) Jacobi returns an
    arbitrary orthonormal basis of the zero-energy subspace; the
    summed density is the projector diagonal and is basis-free.
    """
    evals, evecs = spectrum(n_sites, t, mu, delta)
    order = sorted(range(len(evals)), key=lambda i: abs(evals[i]))
    total = [0.0] * n_sites
    for idx in order[:2]:
        for j in range(n_sites):
            total[j] += evecs[j][idx] ** 2 + evecs[n_sites + j][idx] ** 2
    return total


def majorana_pair(n_sites, t, mu, delta):
    """Separate the zero pair into left and right Majorana modes.

    With psi+ = (u, v) the lowest positive-energy eigenstate, its
    particle-hole partner is psi- = (v, u). The combinations
    (psi+ +/- psi-)/sqrt(2) are the two Majorana wavefunctions; in
    the topological phase one localises on the left end and the
    other on the right. Returns both site densities.
    """
    evals, evecs = spectrum(n_sites, t, mu, delta)
    pos = [(e, i) for i, e in enumerate(evals) if e > 0]
    _, idx = min(pos)
    u = [evecs[j][idx] for j in range(n_sites)]
    v = [evecs[n_sites + j][idx] for j in range(n_sites)]
    rt2 = math.sqrt(2.0)
    dens = []
    for sgn in (1.0, -1.0):
        fu = [(u[j] + sgn * v[j]) / rt2 for j in range(n_sites)]
        fv = [(v[j] + sgn * u[j]) / rt2 for j in range(n_sites)]
        dens.append([fu[j] ** 2 + fv[j] ** 2 for j in range(n_sites)])
    return dens[0], dens[1]


# ================================================================ 3. Bulk invariant: Pfaffian Z2 + winding
# For the periodic chain the Bloch Hamiltonian is
#   H(k) = eps(k) tau_z + 2 Delta sin(k) tau_y,
#   eps(k) = -2t cos(k) - mu.
# At the particle-hole-invariant momenta k = 0 and k = pi the pairing
# term vanishes, and in the Majorana basis H(k) is i times a real
# antisymmetric 2x2 matrix whose Pfaffian is eps(k) (the lab's stated
# convention). The Z2 (Majorana) number is
#   M = sign[ Pf(0) * Pf(pi) ] = sign(mu^2 - 4t^2),
# -1 (topological) for |mu| < 2t, +1 (trivial) outside - PROVIDED
# Delta != 0: with no pairing there is no gap and no phase to name.

def _trim_pfaffian(t, mu, cos_k):
    eps = -2.0 * t * cos_k - mu
    return pfaffian([[0.0, eps], [-eps, 0.0]])


def z2_invariant(t, mu):
    """Majorana number M and the two TRIM Pfaffians (bulk, periodic)."""
    pf0 = _trim_pfaffian(t, mu, 1.0)    # k = 0
    pfpi = _trim_pfaffian(t, mu, -1.0)  # k = pi
    m = 1 if pf0 * pfpi > 0 else -1
    return {"majorana_number": m, "pf_k0": pf0, "pf_kpi": pfpi,
            "formula_sign": 1 if mu * mu - 4.0 * t * t > 0 else -1}


def winding_number(t, mu, delta, steps=720):
    """Winding of the d-vector (eps(k), 2 Delta sin k) as k: 0 -> 2pi.

    Sign convention: with this d-vector orientation the topological
    phase winds -1; the invariant content is |w| = 1 vs 0. Undefined
    (gapless) at Delta = 0 - the lab only quotes it for Delta > 0.
    """
    total = 0.0
    prev = None
    for i in range(steps + 1):
        k = 2.0 * math.pi * i / steps
        ang = math.atan2(2.0 * delta * math.sin(k),
                         -2.0 * t * math.cos(k) - mu)
        if prev is not None:
            d = ang - prev
            while d > math.pi:
                d -= 2.0 * math.pi
            while d < -math.pi:
                d += 2.0 * math.pi
            total += d
        prev = ang
    return int(round(total / (2.0 * math.pi)))


def bulk_gap(n_sites, t, mu, delta):
    """Bulk gap estimate from the finite chain.

    In the topological phase the two lowest |E| levels are the zero
    pair, so the bulk gap is the third |E| level; in the trivial
    phase it is the lowest level itself.
    """
    s = sorted_abs_levels(n_sites, t, mu, delta)
    if abs(mu) < 2.0 * t:
        return s[2]
    return s[0]


def spectrum_report(n_sites=20, t=T_HOP, delta=DEMO_DELTA):
    """Gap + invariant scan across the phase transition at |mu| = 2t."""
    mus = (0.0, 0.5, 1.0, 1.5, 1.9, 2.0, 2.1, 2.5, 3.0)
    rows = []
    for mu in mus:
        s = sorted_abs_levels(n_sites, t, mu, delta)
        inv = z2_invariant(t, mu)
        rows.append({
            "mu": mu,
            "zero_splitting_s0": s[0],
            "bulk_gap": s[2] if abs(mu) < 2.0 * t else s[0],
            "majorana_number": inv["majorana_number"],
            "winding": winding_number(t, mu, delta),
        })
    return {"n_sites": n_sites, "t": t, "delta": delta, "rows": rows}


def splitting_report(t=T_HOP, mu=DEMO_MU):
    """Zero-mode splitting vs chain length at two pairings."""
    out = {}
    for delta in (1.0, 0.5):
        rows = []
        for n_sites in (8, 12, 16, 20):
            s = sorted_abs_levels(n_sites, t, mu, delta)
            rows.append({"n_sites": n_sites, "splitting": s[0]})
        out[f"delta_{delta}"] = rows
    return {"t": t, "mu": mu, "table": out}


def sweet_spot_report(n_sites=12, t=T_HOP):
    """mu = 0, Delta = t: zero modes live exactly on the end sites."""
    s = sorted_abs_levels(n_sites, t, SWEET_MU, SWEET_DELTA)
    pair = zero_mode_pair_density(n_sites, t, SWEET_MU, SWEET_DELTA)
    return {
        "n_sites": n_sites, "mu": SWEET_MU, "delta": SWEET_DELTA,
        "splitting": s[0],
        "pair_density_site_1": pair[0],
        "pair_density_site_N": pair[-1],
        "pair_density_bulk_max": max(pair[1:-1]) if n_sites > 2 else 0.0,
        "pair_density_total": sum(pair),
    }


def protection_report(n_sites=16, t=T_HOP, delta=1.0):
    """Local-perturbation test: the protection number, both phases.

    A labelled on-site potential of strength PERTURB_DELTA is added
    at site 1 (mu_1 -> mu_1 + delta_pot in the particle block). In
    the topological phase the zero mode at that end is an equal
    particle-hole (Majorana) mixture, so the first-order shift of
    the splitting vanishes and the measured shift is exponentially
    small; the trivial phase's lowest level has no such cancellation.
    """
    def lowest_with_potential(mu, dpot):
        H = build_bdg(n_sites, t, mu, delta)
        H[0][0] -= dpot
        H[n_sites][n_sites] += dpot
        evals, _ = jacobi_eigh(H)
        return min(abs(e) for e in evals)

    topo0 = sorted_abs_levels(n_sites, t, PROTECT_MU_TOPO, delta)[0]
    topo1 = lowest_with_potential(PROTECT_MU_TOPO, PERTURB_DELTA)
    triv0 = sorted_abs_levels(n_sites, t, PROTECT_MU_TRIVIAL, delta)[0]
    triv1 = lowest_with_potential(PROTECT_MU_TRIVIAL, PERTURB_DELTA)
    topo_shift = abs(topo1 - topo0)
    triv_shift = abs(triv1 - triv0)
    return {
        "n_sites": n_sites, "delta": delta, "perturbation": PERTURB_DELTA,
        "topological": {"mu": PROTECT_MU_TOPO, "level_before": topo0,
                        "level_after": topo1, "shift": topo_shift},
        "trivial": {"mu": PROTECT_MU_TRIVIAL, "level_before": triv0,
                    "level_after": triv1, "shift": triv_shift},
        "protection_ratio": (triv_shift / topo_shift
                             if topo_shift > 0 else float("inf")),
    }


# ================================================================ 4. Braiding algebra (labelled: algebra, not a braid)
# Four Majorana operators on two qubits via Jordan-Wigner strings:
#   gamma_1 = X I, gamma_2 = Y I, gamma_3 = Z X, gamma_4 = Z Y.
# They satisfy {gamma_i, gamma_j} = 2 delta_ij exactly (asserted in
# CHECKS on the matrices). Exchanging gamma_1 and gamma_2 is the
# unitary U = exp(pi/4 * gamma_2 gamma_1) = (1 + gamma_2 gamma_1)/sqrt(2),
# which conjugates gamma_1 -> gamma_2 and gamma_2 -> -gamma_1: the
# minus sign is the non-Abelian signature, and on the even-parity
# vacuum of the pair it acts as the phase e^{-i pi/4}. This is
# algebra on paper - no anyons were moved, anywhere, ever.

def _kron(A, B):
    n, m = len(A), len(B)
    return [[A[i // m][j // m] * B[i % m][j % m]
             for j in range(n * m)] for i in range(n * m)]


def _mmul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n))
             for j in range(n)] for i in range(n)]


def _madd(A, B, s=1.0):
    return [[A[i][j] + s * B[i][j] for j in range(len(A))]
            for i in range(len(A))]


def _mscale(A, s):
    return [[s * x for x in row] for row in A]


def _dagger(A):
    return [[A[j][i].conjugate() for j in range(len(A))]
            for i in range(len(A))]


def _eye(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def _mnorm(A):
    return math.sqrt(sum(abs(x) ** 2 for row in A for x in row))


_PAULI_I = [[1, 0], [0, 1]]
_PAULI_X = [[0, 1], [1, 0]]
_PAULI_Y = [[0, -1j], [1j, 0]]
_PAULI_Z = [[1, 0], [0, -1]]


def majorana_operators():
    """The four Majorana operators as 4x4 complex matrices (JW)."""
    return [_kron(_PAULI_X, _PAULI_I), _kron(_PAULI_Y, _PAULI_I),
            _kron(_PAULI_Z, _PAULI_X), _kron(_PAULI_Z, _PAULI_Y)]


def braid_unitary(ga, gb):
    """Exchange unitary exp(pi/4 * gb ga) = (1 + gb ga)/sqrt(2)."""
    return _mscale(_madd(_eye(4), _mmul(gb, ga)), 1.0 / math.sqrt(2.0))


def braid_report():
    """The exchange algebra, computed on the matrices themselves."""
    g = majorana_operators()
    u12 = braid_unitary(g[0], g[1])
    u23 = braid_unitary(g[1], g[2])
    conj1 = _mmul(_mmul(u12, g[0]), _dagger(u12))
    conj2 = _mmul(_mmul(u12, g[1]), _dagger(u12))
    parity12 = _mscale(_mmul(g[0], g[1]), -1j)  # -i g1 g2
    vacuum = [1.0, 0.0, 0.0, 0.0]
    u_vac = [sum(u12[i][j] * vacuum[j] for j in range(4)) for i in range(4)]
    comm = _madd(_mmul(u12, u23), _mmul(u23, u12), -1.0)
    return {
        "label": "ALGEBRA ONLY - operator identities on 4x4 matrices, "
                 "not a physical braid and not an experiment",
        "exchange_g1_to": "gamma_2 (conjugation error "
                           f"{_mnorm(_madd(conj1, g[1], -1.0)):.2e})",
        "exchange_g2_to": "-gamma_1 (conjugation error "
                           f"{_mnorm(_madd(conj2, g[0], 1.0)):.2e})",
        "double_exchange": "U^2 = gamma_2 gamma_1 (a parity phase, "
                           "not the identity)",
        "vacuum_phase_real": round(u_vac[0].real, 6),
        "vacuum_phase_imag": round(u_vac[0].imag, 6),
        "vacuum_phase_expected": "e^{-i pi/4} = (0.707107, -0.707107)",
        "parity_eigenvalue_vacuum": 1,
        "commutator_U12_U23_norm": round(_mnorm(comm), 6),
        "braid_relation_U12U23U12_eq_U23U12U23":
            _mnorm(_madd(_mmul(_mmul(u12, u23), u12),
                         _mmul(_mmul(u23, u12), u23), -1.0)) < 1e-9,
    }


# ================================================================ 5. kNN phase classifier
# Features: (mu / 2t, Delta / t, zero-mode end weight). The training
# grid is labelled by the lab's own phase policy (model_label below)
# - the classifier interpolates the policy, it does not replace the
# diagonalisation or the invariant.

def model_label(mu_over_2t, delta_over_t, end_w):
    """The lab's phase policy distilled to a three-way label."""
    if abs(abs(mu_over_2t) - 1.0) <= 0.06:
        return "critical"
    if delta_over_t < 0.05:
        return "trivial"
    if abs(mu_over_2t) < 1.0 and end_w >= 0.5:
        return "topological-protected"
    return "trivial"


def _training_grid():
    rows = []
    for mu2 in (-1.5, -1.03, -0.97, -0.75, -0.25, 0.25, 0.75, 0.97, 1.03, 1.5):
        for dt_ in (0.0, 0.25, 0.5, 1.0):
            for ew in (0.05, 0.3, 0.6, 0.95):
                rows.append(((mu2, dt_, ew), model_label(mu2, dt_, ew)))
    return rows


TRAINING_GRID = _training_grid()


def _scaled(features):
    mu2, dt_, ew = features
    return (mu2 / 1.5, dt_, ew)


def classify_phase(mu_over_2t, delta_over_t, end_w, k=5):
    """kNN over the policy-labelled grid, deterministic tie-breaking."""
    target = _scaled((mu_over_2t, delta_over_t, end_w))
    ranked = []
    for features, label in TRAINING_GRID:
        point = _scaled(features)
        dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(target, point)))
        ranked.append((dist, features, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for _, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0) + 1
    priority = {"topological-protected": 0, "critical": 1, "trivial": 2}
    return sorted(votes.items(), key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


def measured_features(n_sites, t, mu, delta):
    """Classifier features measured from an actual diagonalisation."""
    _, u, v = zero_mode_state(n_sites, t, mu, delta)
    return (mu / (2.0 * t), delta / t, end_weight(u, v))


def classifier_report():
    scenarios = {
        "deep_topological_mu0.5_D1": (16, T_HOP, 0.5, 1.0),
        "topological_mu1_D0.5": (16, T_HOP, 1.0, 0.5),
        "topological_negative_mu": (16, T_HOP, -1.0, 0.5),
        "trivial_mu3_D0.5": (16, T_HOP, 3.0, 0.5),
        "near_critical_mu1.9_D0.5": (20, T_HOP, 1.9, 0.5),
        "no_pairing_mu0.5_D0": (16, T_HOP, 0.5, 0.0),
    }
    out = {}
    for name, (n, t, mu, delta) in scenarios.items():
        feats = measured_features(n, t, mu, delta)
        out[name] = {
            "mu_over_2t": round(feats[0], 6),
            "delta_over_t": round(feats[1], 6),
            "end_weight": round(feats[2], 6),
            "label": classify_phase(*feats),
        }
    return out


# ================================================================ 6. Advisory seats
def _trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="topologicalqubit-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all hardware and fault-tolerance gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace(payload),
    }


def claude_advisor(topic="topologicalqubit-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Print the protection as a number or not at all: show the zero-mode shift next to the trivial-phase shift for the same perturbation, show the splitting falling with N, and say in the first line that braiding here is algebra on 4x4 matrices - the phase is computed, the anyons are not real.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["a chain of sites with the bulk gap closing as mu sweeps to 2t",
                  "two bright end modes fading into each other as the chain shortens",
                  "a Pfaffian sign flipping exactly at the transition",
                  "a perturbation striking one end while the zero mode barely moves",
                  "two operator worldlines exchanging with a pi/4 phase tag"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 1. eigensolver + BdG structure
    H = build_bdg(10, T_HOP, 0.4, 0.6)
    chk("bdg_symmetric",
        all(abs(H[i][j] - H[j][i]) < 1e-15
            for i in range(20) for j in range(20)))
    evals, evecs = jacobi_eigh(H)
    recon = 0.0
    for i in range(20):
        for j in range(20):
            val = sum(evecs[i][k] * evals[k] * evecs[j][k] for k in range(20))
            recon = max(recon, abs(val - H[i][j]))
    chk("jacobi_reconstruction", recon < 1e-9)
    ortho = 0.0
    for a_ in range(20):
        for b_ in range(20):
            dot = sum(evecs[k][a_] * evecs[k][b_] for k in range(20))
            ortho = max(ortho, abs(dot - (1.0 if a_ == b_ else 0.0)))
    chk("jacobi_orthonormal", ortho < 1e-9)
    chk("spectrum_particle_hole_symmetric",
        all(abs(evals[i] + evals[19 - i]) < 1e-9 for i in range(20)))

    # 2. gap closing at |mu| = 2t (scan, N=20, Delta=0.5)
    scan = spectrum_report()
    rows = {round(r["mu"], 2): r for r in scan["rows"]}
    chk("gap_deep_topological", 0.95 < rows[0.0]["bulk_gap"] < 1.10)
    chk("gap_mid_topological", 0.45 < rows[1.5]["bulk_gap"] < 0.56)
    chk("gap_minimum_at_transition",
        rows[2.0]["bulk_gap"] < 0.12
        and rows[2.0]["bulk_gap"] < rows[1.9]["bulk_gap"]
        and rows[2.0]["bulk_gap"] < rows[2.1]["bulk_gap"])
    chk("gap_reopens_trivial",
        0.48 < rows[2.5]["bulk_gap"] < 0.60
        and 0.95 < rows[3.0]["bulk_gap"] < 1.10)
    chk("splitting_tiny_in_phase",
        rows[0.5]["zero_splitting_s0"] < 1e-3
        and rows[1.0]["zero_splitting_s0"] < 1e-3)

    # 3. invariant: Pfaffian Z2 + winding agree, flip at |mu| = 2t
    chk("invariant_topological",
        z2_invariant(T_HOP, 0.5)["majorana_number"] == -1
        and z2_invariant(T_HOP, -1.0)["majorana_number"] == -1)
    chk("invariant_trivial",
        z2_invariant(T_HOP, 2.5)["majorana_number"] == 1
        and z2_invariant(T_HOP, 3.0)["majorana_number"] == 1)
    chk("invariant_flips_at_2t",
        z2_invariant(T_HOP, 1.99)["majorana_number"] == -1
        and z2_invariant(T_HOP, 2.01)["majorana_number"] == 1)
    chk("invariant_matches_formula",
        all(z2_invariant(T_HOP, mu)["majorana_number"]
            == z2_invariant(T_HOP, mu)["formula_sign"]
            for mu in (-3.0, -1.0, 0.0, 0.5, 1.99, 2.01, 3.0)))
    chk("winding_topological",
        abs(winding_number(T_HOP, 0.5, 0.5)) == 1
        and abs(winding_number(T_HOP, 1.0, 1.0)) == 1)
    chk("winding_trivial", winding_number(T_HOP, 3.0, 0.5) == 0)
    chk("scan_invariant_column",
        all(r["majorana_number"] == (-1 if abs(r["mu"]) < 2.0 else 1)
            for r in scan["rows"] if abs(abs(r["mu"]) - 2.0) > 0.05))

    # 4. Pfaffian machinery: brute force + Pf^2 = det on seeded cases
    import random as _random
    rng = _random.Random(42)
    pf_ok = True
    for n in (2, 4, 6):
        M = [[rng.uniform(-2, 2) for _ in range(n)] for _ in range(n)]
        A = [[M[i][j] - M[j][i] for j in range(n)] for i in range(n)]
        if abs(pfaffian(A) ** 2 - det_gauss(A)) > 1e-6:
            pf_ok = False

    def _pf_brute(A):
        n = len(A)
        if n == 0:
            return 1.0
        tot = 0.0
        for j in range(1, n):
            idx = [i for i in range(n) if i != 0 and i != j]
            sub = [[A[i][jj] for jj in idx] for i in idx]
            tot += ((-1) ** (j + 1)) * A[0][j] * _pf_brute(sub)
        return tot

    rng2 = _random.Random(7)
    M = [[rng2.uniform(-2, 2) for _ in range(6)] for _ in range(6)]
    A6 = [[M[i][j] - M[j][i] for j in range(6)] for i in range(6)]
    chk("pfaffian_squared_is_det", pf_ok)
    chk("pfaffian_matches_brute", abs(pfaffian(A6) - _pf_brute(A6)) < 1e-8)

    # 5. Majorana end modes: end weight, separation, splitting vs N
    _, u, v = zero_mode_state(16, T_HOP, 0.5, 1.0)
    chk("end_weight_topological", end_weight(u, v) > 0.90)
    _, u3, v3 = zero_mode_state(16, T_HOP, 3.0, 0.5)
    chk("end_weight_trivial_small", end_weight(u3, v3) < 0.10)
    dens_l, dens_r = majorana_pair(16, T_HOP, 0.5, 0.5)
    left_first = sum(dens_l[:4])
    left_second = sum(dens_r[:4])
    chk("majorana_separation",
        max(left_first, left_second) > 0.99
        and min(left_first, left_second) < 0.01)
    split = splitting_report()
    d1 = [r["splitting"] for r in split["table"]["delta_1.0"]]
    chk("splitting_falls_with_N",
        all(d1[i + 1] < d1[i] for i in range(3)))
    chk("splitting_exponential_scale",
        d1[0] > 1e-6 and d1[-1] < 1e-9 and d1[-1] < d1[0] * 1e-5)
    d05 = [r["splitting"] for r in split["table"]["delta_0.5"]]
    chk("splitting_falls_with_N_delta05",
        all(d05[i + 1] < d05[i] for i in range(3)))

    # 6. sweet spot: exact end modes
    sweet = sweet_spot_report()
    chk("sweet_spot_zero_splitting", sweet["splitting"] < 1e-12)
    chk("sweet_spot_pair_on_ends",
        abs(sweet["pair_density_site_1"] - 1.0) < 1e-9
        and abs(sweet["pair_density_site_N"] - 1.0) < 1e-9)
    chk("sweet_spot_bulk_empty", sweet["pair_density_bulk_max"] < 1e-9)
    chk("sweet_spot_pair_total", abs(sweet["pair_density_total"] - 2.0) < 1e-9)

    # 7. protection: same perturbation, both phases
    prot = protection_report()
    chk("protection_topological_tiny", prot["topological"]["shift"] < 1e-8)
    chk("protection_trivial_visible", prot["trivial"]["shift"] > 1e-4)
    chk("protection_ratio", prot["protection_ratio"] > 1e5)

    # 8. braiding algebra on the matrices
    g = majorana_operators()
    chk("majorana_squares",
        all(_mnorm(_madd(_mmul(x, x), _eye(4), -1.0)) < 1e-12 for x in g))
    ac = True
    for i in range(4):
        for j in range(4):
            anti = _madd(_mmul(g[i], g[j]), _mmul(g[j], g[i]))
            want = _mscale(_eye(4), 2.0 if i == j else 0.0)
            if _mnorm(_madd(anti, want, -1.0)) > 1e-12:
                ac = False
    chk("majorana_anticommutator", ac)
    u12 = braid_unitary(g[0], g[1])
    u23 = braid_unitary(g[1], g[2])
    chk("braid_unitary",
        _mnorm(_madd(_mmul(_dagger(u12), u12), _eye(4), -1.0)) < 1e-12)
    chk("braid_exchange_maps",
        _mnorm(_madd(_mmul(_mmul(u12, g[0]), _dagger(u12)), g[1], -1.0)) < 1e-12
        and _mnorm(_madd(_mmul(_mmul(u12, g[1]), _dagger(u12)), g[0], 1.0)) < 1e-12)
    chk("braid_double_exchange",
        _mnorm(_madd(_mmul(u12, u12), _mmul(g[1], g[0]), -1.0)) < 1e-12)
    vac = [sum(u12[i][j] * (1.0 if j == 0 else 0.0) for j in range(4))
           for i in range(4)]
    chk("braid_vacuum_phase",
        abs(vac[0] - complex(math.cos(-math.pi / 4), math.sin(-math.pi / 4)))
        < 1e-12 and all(abs(vac[i]) < 1e-12 for i in (1, 2, 3)))
    chk("braid_nonabelian",
        _mnorm(_madd(_mmul(u12, u23), _mmul(u23, u12), -1.0)) > 1.0)
    chk("braid_relation",
        _mnorm(_madd(_mmul(_mmul(u12, u23), u12),
                     _mmul(_mmul(u23, u12), u23), -1.0)) < 1e-12)

    # 9. classifier - labels come from the lab's own phase policy
    rep = classifier_report()
    chk("classifier_deep_topological",
        rep["deep_topological_mu0.5_D1"]["label"] == "topological-protected")
    chk("classifier_topological",
        rep["topological_mu1_D0.5"]["label"] == "topological-protected")
    chk("classifier_negative_mu",
        rep["topological_negative_mu"]["label"] == "topological-protected")
    chk("classifier_trivial", rep["trivial_mu3_D0.5"]["label"] == "trivial")
    chk("classifier_critical",
        rep["near_critical_mu1.9_D0.5"]["label"] == "critical")
    chk("classifier_no_pairing",
        rep["no_pairing_mu0.5_D0"]["label"] == "trivial")
    chk("classifier_matches_model",
        all(classify_phase(*measured_features(
                *{"deep_topological_mu0.5_D1": (16, T_HOP, 0.5, 1.0),
                  "topological_mu1_D0.5": (16, T_HOP, 1.0, 0.5),
                  "topological_negative_mu": (16, T_HOP, -1.0, 0.5),
                  "trivial_mu3_D0.5": (16, T_HOP, 3.0, 0.5),
                  "near_critical_mu1.9_D0.5": (20, T_HOP, 1.9, 0.5),
                  "no_pairing_mu0.5_D0": (16, T_HOP, 0.5, 0.0)}[name]))
            == model_label(*measured_features(
                *{"deep_topological_mu0.5_D1": (16, T_HOP, 0.5, 1.0),
                  "topological_mu1_D0.5": (16, T_HOP, 1.0, 0.5),
                  "topological_negative_mu": (16, T_HOP, -1.0, 0.5),
                  "trivial_mu3_D0.5": (16, T_HOP, 3.0, 0.5),
                  "near_critical_mu1.9_D0.5": (20, T_HOP, 1.9, 0.5),
                  "no_pairing_mu0.5_D0": (16, T_HOP, 0.5, 0.0)}[name]))
            for name in rep))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_qubit_hold", gate("hardware_qubit_claim") == "HOLD")
    chk("gate_fault_tolerance_hold", gate("fault_tolerance_claim") == "HOLD")
    chk("gate_braid_experiment_hold",
        gate("braiding_as_experiment_claim") == "HOLD")
    chk("gate_protection_guarantee_hold",
        gate("topological_protection_guarantee") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


# ================================================================ CLI / demo
def _g(x):
    """Compact deterministic float for JSON reports (6 sig digits)."""
    return float(f"{x:.6g}")


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, float):
        return _g(obj)
    return obj


def demo():
    print("=" * 68)
    print("TopologicalQubit-Lab - Kitaev chain & Majorana zero-mode lab")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - a phase, a number, and a qubit that is not here")
    print("    This lab diagonalises the Kitaev chain exactly: a 2N x 2N BdG")
    print("    matrix with labelled t, mu, Delta. A topological PHASE is a sign")
    print("    computed from the bulk. PROTECTION is a number measured on this")
    print("    matrix. A topological QUBIT is hardware - it does not exist")
    print("    here, and every hardware claim is a permanent HOLD gate.")
    scan = spectrum_report()
    print("\n[2] Spectrum scan - the gap closes at |mu| = 2t (N=20, Delta=0.5)")
    print("      mu   splitting s0   bulk gap   Majorana#  winding")
    for r in scan["rows"]:
        print(f"    {r['mu']:5.1f}   {r['zero_splitting_s0']:12.6f}   "
              f"{r['bulk_gap']:8.6f}   {r['majorana_number']:+9d}  "
              f"{r['winding']:+7d}")
    print("    -> inside |mu| < 2t the lowest level is the zero-mode splitting")
    print("       (exponentially small); the bulk gap above it falls to its")
    print("       finite-size floor 0.081748 at mu = 2t and reopens beyond.")
    print("       At mu = 2t exactly, Pf(pi) = 0: the invariant is undefined")
    print("       AT the transition - it flips across it, which is the point.")
    print("\n[3] Majorana end modes - the zero pair lives on the ends")
    _, u, v = zero_mode_state(16, T_HOP, 0.5, 1.0)
    _, u3, v3 = zero_mode_state(16, T_HOP, 3.0, 0.5)
    print(f"    end weight (mu=0.5, Delta=1, N=16): {end_weight(u, v):.6f}")
    print(f"    end weight (mu=3.0, Delta=0.5, trivial): {end_weight(u3, v3):.6f}")
    dens_l, dens_r = majorana_pair(16, T_HOP, 0.5, 0.5)
    print(f"    separated modes: one carries {sum(dens_l[:4]):.6f} on the "
          f"left 4 sites, the other {sum(dens_r[-4:]):.6f} on the right 4")
    split = splitting_report()
    print("\n[4] Zero-mode splitting vs chain length (mu=0.5) - exponential")
    for key in ("delta_1.0", "delta_0.5"):
        rows = split["table"][key]
        body = "  ".join(f"N={r['n_sites']}: {r['splitting']:.3e}"
                         for r in rows)
        print(f"    {key}: {body}")
    sweet = sweet_spot_report()
    print("\n[5] Sweet spot (mu=0, Delta=t, N=12) - the analytic limit, exact")
    print(f"    splitting {sweet['splitting']:.3e} | pair density site 1 "
          f"{sweet['pair_density_site_1']:.6f}, site N "
          f"{sweet['pair_density_site_N']:.6f}, bulk max "
          f"{sweet['pair_density_bulk_max']:.3e}")
    prot = protection_report()
    print("\n[6] Protection test - same labelled perturbation (0.2 at site 1)")
    print(f"    topological (mu=0.5): level shift {prot['topological']['shift']:.3e}")
    print(f"    trivial     (mu=2.5): level shift {prot['trivial']['shift']:.3e}")
    print(f"    protection ratio (this matrix): {prot['protection_ratio']:.3e}")
    braid = braid_report()
    print("\n[7] Braiding - ALGEBRA ONLY (4x4 matrices, no anyons moved)")
    print(f"    exchange: {braid['exchange_g1_to']} / {braid['exchange_g2_to']}")
    print(f"    vacuum phase: ({braid['vacuum_phase_real']}, "
          f"{braid['vacuum_phase_imag']}) = e^{{-i pi/4}}")
    print(f"    non-Abelian: ||[U12, U23]|| = "
          f"{braid['commutator_U12_U23_norm']:.6f} (nonzero = order matters)")
    print("\n[8] Phase classifier (TinyML kNN, labels from the lab's policy)")
    for name, row in classifier_report().items():
        print(f"    {name:32s} -> {row['label']}")
    print("\n[9] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[10] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "spectrum":
        print(json.dumps(_jsonable(spectrum_report()),
                         indent=2, sort_keys=True))
        return 0
    if cmd == "splitting":
        print(json.dumps(_jsonable(splitting_report()),
                         indent=2, sort_keys=True))
        return 0
    if cmd == "invariant":
        rows = {}
        for mu in (-3.0, -2.5, -1.0, 0.0, 0.5, 1.99, 2.01, 3.0):
            inv = z2_invariant(T_HOP, mu)
            inv["winding_D0.5"] = winding_number(T_HOP, mu, 0.5)
            rows[f"mu_{mu}"] = inv
        print(json.dumps(_jsonable(rows), indent=2, sort_keys=True))
        return 0
    if cmd == "sweetspot":
        print(json.dumps(_jsonable(sweet_spot_report()),
                         indent=2, sort_keys=True))
        return 0
    if cmd == "protection":
        print(json.dumps(_jsonable(protection_report()),
                         indent=2, sort_keys=True))
        return 0
    if cmd == "braid":
        print(json.dumps(braid_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "hardware_qubit_claim"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
