import sys; sys.path.insert(0, "..")
import numpy as np
import latticefit

PHI = latticefit.PHI

def test_exact_lattice():
    """Points exactly on the lattice should give zero residuals."""
    k = np.array([0, 4, 8, 12, 20])
    x = 1.0 * PHI ** (k / 4)
    result = latticefit.fit(x, anchor=1.0, base=PHI, denom=4)
    assert result.rms < 1e-10
    np.testing.assert_array_equal(result.labels, k)

def test_sm_masses():
    masses = [5.10999e-4, 0.105658, 1.77686,
              0.00216,    1.275,    172.76,
              0.00467,    0.0934,   4.18]
    result = latticefit.fit(masses, anchor=5.10999e-4, base=PHI, denom=4)
    assert result.rms < 0.10
    assert result.labels[0] == 0   # electron anchors at q=0

def test_anchor_shift_invariance():
    """Relative q differences must be invariant under anchor change."""
    masses = [5.10999e-4, 0.105658, 1.77686]
    r1 = latticefit.fit(masses, anchor=masses[0], base=PHI, denom=4)
    r2 = latticefit.fit(masses, anchor=masses[1], base=PHI, denom=4)
    diffs1 = np.diff(r1.labels)
    diffs2 = np.diff(r2.labels)
    np.testing.assert_array_equal(diffs1, diffs2)
