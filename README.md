# Recoverability Framework: v1.2.2 Reference Implementation

Pure-Python reference implementation of the v1.2.2 recoverability metrics and the Section 5 falsification harness.

- **Python:** 3.10+ (uses `X | Y` type hints)
- **Dependencies:** none (standard library only)
- **File:** `recoverability_v1_2_2.py`

## What it does

Every recoverability score is a **skill ratio** against a null baseline:

```
R = 1 - E_model / E_null        (unclipped; negative = worse than null)
```

Edge cases are explicit rather than silent:

| Case | Result |
|---|---|
| `E_null == 0` and `E_model == 0` | `1.0` |
| `E_null == 0` and `E_model > 0` | `"undefined"` |
| `E_null > 0` | `1 - E_model / E_null` |

### Metrics

| Function | Section | Measures | Null baseline |
|---|---|---|---|
| `skill_ratio` | 2 | Core skill ratio | n/a |
| `structural_recoverability` | 3.1 | `R_s`: graph recovery, L1 off-diagonal error | Density-matched mean off-diagonal (or supplied) |
| `directional_asymmetry` | 3.2A | `A_t`, `A_t_signed`: scale-invariant directional diagnostic | n/a (diagnostic) |
| `sequential_recoverability_hamming` | 3.2B | `R_t`: sequence recovery, Hamming distance | Multiset permutation, `1 - sum(q_e^2)` |
| `causal_recoverability` | 3.3 | `R_c`: held-out classification | Balanced `1/K`, majority, or proportional |
| `global_recoverability_mse` | 3.4 | `R_g`: out-of-sample hidden-state MSE (ridge) | Training-set mean |
| `recoverability_retention` | 4 | `eta`: retention from `t1` to `t2` | Defined only when `R(t1) > 0` |

Retention returns a `status` of `retained`, `sign_crossing`, `non_positive_baseline`, or `undefined`.

### Section 5 falsification harness

`run_section5_falsification_harness` runs the five-condition control:

1. **Matching:** relative energy and state differences within tolerance.
2. **Pre-perturbation measurement.**
3. **Identical perturbation probe.**
4. **R-only predictive power:** 95% percentile bootstrap CI of CV R² must have a lower bound > 0, otherwise `falsified`.
5. **Incremental control:** ΔR²_CV from adding R on top of Z_int must have a CI lower bound > 0 **and** a signed bootstrap median at or above the practical floor (default `0.01`), otherwise `hidden_state_dependence`.

Possible classifications:

| Classification | Meaning |
|---|---|
| `irreducible_recoverability_effect` | Passed Cond 4 and survived the Cond 5 control |
| `hidden_state_dependence` | R predicts, but adds nothing beyond Z_int |
| `falsified` | R has no out-of-sample predictive power |
| `inconclusive` | Conditions 1 to 3 failed, too few samples, zero-variance target, or a degenerate CI |

**CV design:** each bootstrap replicate draws one fold assignment over the original rows and reuses it for the R-only, Z-only, and Z+R models, so ΔR² is paired. Bootstrap rows inherit the fold of their source row, so duplicates never straddle train and validation.

## Usage

```python
from recoverability_v1_2_2 import (
    skill_ratio,
    structural_recoverability,
    run_section5_falsification_harness,
)

print(skill_ratio(0.2, 1.0))   # 0.8

result = run_section5_falsification_harness(
    delta_r_samples=r_list,          # list[list[float]]
    delta_z_int_samples=z_list,      # list[list[float]]
    delta_resp_samples=y_list,       # list[float]
    max_rel_energy_diff=0.002,
    max_rel_state_diff=0.005,
    n_bootstrap=300,
    seed=20261121,                   # fixed seed for reproducibility
)
print(result["classification"], result["reason"])
```

## Verification

```bash
python recoverability_v1_2_2.py
# All v1.2.2 default-path & Section 5 harness assertions passed.
```

The built-in assertions cover skill-ratio edge cases, uniform-graph degeneracy, scale invariance and sign antisymmetry of the asymmetry diagnostic, retention domains, and the harness classifications: hidden-state, linear and nonlinear positive paths, and a pure-noise leakage check.

## Known limitations

- **Linear baseline.** The harness fits linear ridge models. If the response depends nonlinearly on `Z_int` and `R` is a proxy for `Z_int`, the linear Z-only model underfits and `R` can show spurious incremental lift, producing a false `irreducible_recoverability_effect`. The baseline model class should be at least as flexible as the R-model; this is not yet enforced.
- **Oracle null in `R_c`.** The causal baseline is computed from held-out test labels, which is stricter than a baseline built from training labels.
- **Unstandardized features.** Ridge is applied without feature scaling, so `ridge_lambda` acts differently on columns of different scale.
- **Small samples.** With `n` near 6 and 5 splits, folds hold 1 to 2 points and CV R² is mostly noise. Treat results at that size as indicative only.
- **Pure Python.** The solver is plain Gauss-Jordan on the normal equations and is not meant for large feature counts.

## Changelog

**v1.2.2 (this release)**
- Cond 5 now uses a signed median floor and a strictly positive CI lower bound. A significantly negative ΔR² no longer passes.
- Paired, leakage-free fold assignment for the bootstrap CV.
- Ridge solver raises on a singular normal equation instead of returning invalid weights.
- `global_recoverability_mse` raises a clear error on empty `hid_train`.
- Trigger key renamed `ci95_includes_zero` to `ci95_lower_not_positive`.
- Added leakage, empty-train, and fold-pairing tests.

## Author

Nile Green, independent researcher. ORCID: [0009-0007-3629-6404](https://orcid.org/0009-0007-3629-6404)

## Citation

If you use this code, please cite the Zenodo release (DOI assigned on release):

```
Green, N. (2026). Recoverability Framework v1.2.2 reference implementation [Software]. Zenodo. https://doi.org/<DOI>
```

## License

MIT. See `LICENSE`.
