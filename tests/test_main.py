"""Tests for TopologicalQubit-Lab: Kitaev chain, invariant, protection, braiding."""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import main as lab


# ---------------- eigensolver + BdG structure ----------------

def test_bdg_is_symmetric():
    H = lab.build_bdg(8, 1.0, 0.4, 0.6)
    assert all(abs(H[i][j] - H[j][i]) < 1e-15
               for i in range(16) for j in range(16))


def test_bdg_blocks():
    H = lab.build_bdg(6, 1.0, 0.3, 0.7)
    assert H[0][0] == -0.3 and H[6][6] == 0.3
    assert H[0][1] == -1.0 and H[6][7] == 1.0
    assert H[0][7] == -0.7 and H[1][6] == 0.7  # pairing block antisymmetric


def test_jacobi_reconstruction_and_orthonormality():
    H = lab.build_bdg(8, 1.0, 0.4, 0.6)
    evals, evecs = lab.jacobi_eigh(H)
    n = 16
    for i in range(n):
        for j in range(n):
            val = sum(evecs[i][k] * evals[k] * evecs[j][k] for k in range(n))
            assert abs(val - H[i][j]) < 1e-9
    for a in range(n):
        for b in range(n):
            dot = sum(evecs[k][a] * evecs[k][b] for k in range(n))
            assert abs(dot - (1.0 if a == b else 0.0)) < 1e-9


def test_spectrum_particle_hole_symmetric():
    evals, _ = lab.spectrum(10, 1.0, 0.7, 0.5)
    assert all(abs(evals[i] + evals[len(evals) - 1 - i]) < 1e-9
               for i in range(len(evals)))


# ---------------- gap closing + invariant ----------------

def test_bulk_gap_values():
    assert 0.95 < lab.bulk_gap(20, 1.0, 0.0, 0.5) < 1.10
    assert 0.45 < lab.bulk_gap(20, 1.0, 1.5, 0.5) < 0.56
    assert 0.48 < lab.bulk_gap(20, 1.0, 2.5, 0.5) < 0.60


def test_gap_minimum_at_transition():
    g = {mu: lab.bulk_gap(20, 1.0, mu, 0.5) for mu in (1.9, 2.0, 2.1)}
    assert g[2.0] < 0.12 and g[2.0] < g[1.9] and g[2.0] < g[2.1]


def test_invariant_phases():
    assert lab.z2_invariant(1.0, 0.5)["majorana_number"] == -1
    assert lab.z2_invariant(1.0, -1.0)["majorana_number"] == -1
    assert lab.z2_invariant(1.0, 2.5)["majorana_number"] == 1
    assert lab.z2_invariant(1.0, 3.0)["majorana_number"] == 1


def test_invariant_flips_at_2t():
    assert lab.z2_invariant(1.0, 1.99)["majorana_number"] == -1
    assert lab.z2_invariant(1.0, 2.01)["majorana_number"] == 1


def test_winding():
    assert abs(lab.winding_number(1.0, 0.5, 0.5)) == 1
    assert lab.winding_number(1.0, 3.0, 0.5) == 0


def test_pfaffian_squared_is_det():
    import random
    rng = random.Random(42)
    for n in (2, 4, 6):
        M = [[rng.uniform(-2, 2) for _ in range(n)] for _ in range(n)]
        A = [[M[i][j] - M[j][i] for j in range(n)] for i in range(n)]
        assert abs(lab.pfaffian(A) ** 2 - lab.det_gauss(A)) < 1e-6


# ---------------- Majorana end modes ----------------

def test_end_weight_contrast():
    _, u, v = lab.zero_mode_state(16, 1.0, 0.5, 1.0)
    assert lab.end_weight(u, v) > 0.90
    _, u3, v3 = lab.zero_mode_state(16, 1.0, 3.0, 0.5)
    assert lab.end_weight(u3, v3) < 0.10


def test_majorana_separation():
    dl, dr = lab.majorana_pair(16, 1.0, 0.5, 0.5)
    lefts = (sum(dl[:4]), sum(dr[:4]))
    assert max(lefts) > 0.99 and min(lefts) < 0.01


def test_splitting_falls_exponentially():
    rep = lab.splitting_report()
    d1 = [r["splitting"] for r in rep["table"]["delta_1.0"]]
    assert all(d1[i + 1] < d1[i] for i in range(3))
    assert d1[0] > 1e-6 and d1[-1] < 1e-9
    d05 = [r["splitting"] for r in rep["table"]["delta_0.5"]]
    assert all(d05[i + 1] < d05[i] for i in range(3))


def test_sweet_spot_exact():
    rep = lab.sweet_spot_report()
    assert rep["splitting"] < 1e-12
    assert abs(rep["pair_density_site_1"] - 1.0) < 1e-9
    assert abs(rep["pair_density_site_N"] - 1.0) < 1e-9
    assert rep["pair_density_bulk_max"] < 1e-9


# ---------------- protection ----------------

def test_protection_contrast():
    rep = lab.protection_report()
    assert rep["topological"]["shift"] < 1e-8
    assert rep["trivial"]["shift"] > 1e-4
    assert rep["protection_ratio"] > 1e5


# ---------------- braiding algebra ----------------

def test_majorana_clifford_algebra():
    g = lab.majorana_operators()
    for i in range(4):
        for j in range(4):
            anti = lab._madd(lab._mmul(g[i], g[j]), lab._mmul(g[j], g[i]))
            want = lab._mscale(lab._eye(4), 2.0 if i == j else 0.0)
            assert lab._mnorm(lab._madd(anti, want, -1.0)) < 1e-12


def test_braid_exchange():
    g = lab.majorana_operators()
    u = lab.braid_unitary(g[0], g[1])
    c1 = lab._mmul(lab._mmul(u, g[0]), lab._dagger(u))
    c2 = lab._mmul(lab._mmul(u, g[1]), lab._dagger(u))
    assert lab._mnorm(lab._madd(c1, g[1], -1.0)) < 1e-12
    assert lab._mnorm(lab._madd(c2, g[0], 1.0)) < 1e-12


def test_braid_vacuum_phase():
    g = lab.majorana_operators()
    u = lab.braid_unitary(g[0], g[1])
    assert abs(u[0][0] - complex(math.cos(-math.pi / 4),
                                 math.sin(-math.pi / 4))) < 1e-12


def test_braid_nonabelian_and_relation():
    g = lab.majorana_operators()
    u12 = lab.braid_unitary(g[0], g[1])
    u23 = lab.braid_unitary(g[1], g[2])
    comm = lab._madd(lab._mmul(u12, u23), lab._mmul(u23, u12), -1.0)
    assert lab._mnorm(comm) > 1.0
    lhs = lab._mmul(lab._mmul(u12, u23), u12)
    rhs = lab._mmul(lab._mmul(u23, u12), u23)
    assert lab._mnorm(lab._madd(lhs, rhs, -1.0)) < 1e-12


# ---------------- classifier + gates + seats ----------------

def test_classifier_scenarios():
    rep = lab.classifier_report()
    assert rep["deep_topological_mu0.5_D1"]["label"] == "topological-protected"
    assert rep["topological_mu1_D0.5"]["label"] == "topological-protected"
    assert rep["topological_negative_mu"]["label"] == "topological-protected"
    assert rep["trivial_mu3_D0.5"]["label"] == "trivial"
    assert rep["near_critical_mu1.9_D0.5"]["label"] == "critical"
    assert rep["no_pairing_mu0.5_D0"]["label"] == "trivial"


def test_hold_gates():
    assert len(lab.HOLD_GATES) == 14
    assert lab.gate("hardware_qubit_claim") == "HOLD"
    assert lab.gate("fault_tolerance_claim") == "HOLD"
    assert lab.gate("braiding_as_experiment_claim") == "HOLD"
    assert lab.gate("run_simulation") == "PROCEED"


def test_seats_cannot_override():
    brains = lab.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["seat"] == "Gemini Argon 4 Extended Thinking"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False and brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = lab.run_checks()
    assert all(ok for _, ok in checks), [n for n, ok in checks if not ok]


def test_hold_demo_exit_2():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent / "main.py"),
         "hold-demo"], capture_output=True, text=True)
    assert proc.returncode == 2


def test_demo_deterministic():
    root = Path(__file__).resolve().parent.parent
    a = subprocess.run([sys.executable, str(root / "main.py"), "demo"],
                       capture_output=True, text=True)
    b = subprocess.run([sys.executable, str(root / "main.py"), "demo"],
                       capture_output=True, text=True)
    assert a.returncode == 0 and a.stdout == b.stdout
