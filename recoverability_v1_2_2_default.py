import math
import random
import statistics
import warnings
from typing import Any


UNDEFINED = "undefined"
SIGN_CROSSING = "sign_crossing"
NON_POSITIVE_BASELINE = "non_positive_baseline"
COMPLETE_OBSERVATION = "complete_observation"
BOOTSTRAP_CI_METHOD = "percentile_95"


def skill_ratio(
    e_model: float,
    e_null: float,
    tol: float = 1e-12,
) -> float | str:
    """
    v1.2.2 Section 2 skill ratio:
      - e_null == 0 and e_model == 0 -> 1.0
      - e_null == 0 and e_model > 0  -> "undefined"
      - e_null > 0                   -> 1.0 - (e_model / e_null)  (unclipped)
    """
    if e_model < -tol or e_null < -tol:
        raise ValueError("Reconstruction errors E_model and E_null must be non-negative.")

    e_m = 0.0 if abs(e_model) <= tol else float(e_model)
    e_n = 0.0 if abs(e_null) <= tol else float(e_null)

    if e_n == 0.0:
        if e_m == 0.0:
            return 1.0
        return UNDEFINED
    return 1.0 - (e_m / e_n)


def structural_recoverability(
    inferred_c: list[list[float]],
    ref_c: list[list[float]],
    null_c: list[list[float]] | None = None,
    tol: float = 1e-12,
) -> dict[str, Any]:
    """
    v1.2.2 Section 3.1 R_s using L1 off-diagonal error and default density-matched
    mean off-diagonal null unless explicitly supplied.
    """
    n = len(ref_c)
    if n < 2:
        raise ValueError("Graph size N must be at least 2.")

    off_vals = [
        float(ref_c[i][j])
        for i in range(n)
        for j in range(n)
        if i != j
    ]
    is_uniform_ref = (max(off_vals) - min(off_vals)) <= tol
    mean_off = sum(off_vals) / len(off_vals)

    eff_null = (
        [[0.0 if i == j else mean_off for j in range(n)] for i in range(n)]
        if null_c is None
        else null_c
    )

    e_model = sum(
        abs(float(inferred_c[i][j]) - float(ref_c[i][j]))
        for i in range(n)
        for j in range(n)
        if i != j
    )
    e_null = sum(
        abs(float(eff_null[i][j]) - float(ref_c[i][j]))
        for i in range(n)
        for j in range(n)
        if i != j
    )

    score = skill_ratio(e_model, e_null, tol=tol)
    warning = (
        "uniform_reference_graph_degeneracy"
        if (is_uniform_ref and null_c is None)
        else None
    )
    return {
        "R_s": score,
        "E_model": e_model,
        "E_null": e_null,
        "warning": warning,
    }


def directional_asymmetry(
    counts: list[list[float]],
    orientations: list[list[int]] | None = None,
    eps_scale: float = 1e-9,
    tol: float = 1e-12,
) -> dict[str, float | str]:
    """
    v1.2.2 Section 3.2A scale-invariant directional asymmetry diagnostic (A_t, A_t_signed).
    """
    n = len(counts)
    edge_fluxes: list[tuple[int, int, float]] = []
    f_tot = 0.0

    for i in range(n):
        for j in range(i + 1, n):
            n_ij = max(0.0, float(counts[i][j]))
            n_ji = max(0.0, float(counts[j][i]))
            f_ij = n_ij + n_ji
            if f_ij > tol:
                edge_fluxes.append((i, j, f_ij))
                f_tot += f_ij

    if f_tot <= tol or not edge_fluxes:
        return {"A_t": UNDEFINED, "A_t_signed": UNDEFINED}

    f_bar = f_tot / len(edge_fluxes)
    eps = eps_scale * f_bar

    a_t = 0.0
    for i, j, f_ij in edge_fluxes:
        n_ij = max(0.0, float(counts[i][j]))
        n_ji = max(0.0, float(counts[j][i]))
        delta_ij = abs(n_ij - n_ji) / (f_ij + eps)
        a_t += (f_ij / f_tot) * delta_ij

    if orientations is None:
        return {"A_t": a_t, "A_t_signed": 0.0}

    f_orient = 0.0
    weighted_signed = 0.0
    for i, j, f_ij in edge_fluxes:
        sigma_ij = int(orientations[i][j])
        if sigma_ij not in (-1, 0, 1):
            raise ValueError("Reference orientations sigma_ij must be in {-1, 0, +1}.")
        if sigma_ij != 0:
            n_ij = max(0.0, float(counts[i][j]))
            n_ji = max(0.0, float(counts[j][i]))
            f_orient += f_ij
            weighted_signed += f_ij * sigma_ij * ((n_ij - n_ji) / (f_ij + eps))

    a_t_signed = 0.0 if f_orient <= tol else (weighted_signed / f_orient)
    return {"A_t": a_t, "A_t_signed": a_t_signed}


def sequential_recoverability_hamming(
    seq_true: list[Any],
    seq_recon: list[Any],
    tol: float = 1e-12,
) -> dict[str, Any]:
    """
    v1.2.2 Section 3.2B default sequential recoverability path (Hamming distance +
    exact multiset permutation baseline E_null = 1 - sum(q_e^2)).
    """
    if len(seq_true) != len(seq_recon) or len(seq_true) == 0:
        raise ValueError("seq_true and seq_recon must be non-empty of equal length.")

    length = len(seq_true)
    mismatches = sum(1 for a, b in zip(seq_true, seq_recon) if a != b)
    e_model = mismatches / length

    sym_counts: dict[Any, int] = {}
    for sym in seq_true:
        sym_counts[sym] = sym_counts.get(sym, 0) + 1

    e_null = 1.0 - sum((cnt / length) ** 2 for cnt in sym_counts.values())
    return {
        "R_t": skill_ratio(e_model, e_null, tol=tol),
        "E_model": e_model,
        "E_null": e_null,
    }


def causal_recoverability(
    y_true: list[int],
    y_pred: list[int],
    baseline_rule: str = "majority",
    declared_prior_sum_sq: float | None = None,
    tol: float = 1e-12,
) -> dict[str, Any]:
    """
    v1.2.2 Section 3.3 causal recoverability on held-out test set.

    NOTE: the baseline is computed from the held-out test labels (oracle null),
    which is stricter than a deployable training-set baseline.
    """
    if len(y_true) != len(y_pred) or len(y_true) == 0:
        raise ValueError("y_true and y_pred must be non-empty of equal length.")

    n_test = len(y_true)
    acc_model = sum(1 for yt, yp in zip(y_true, y_pred) if yt == yp) / n_test

    counts: dict[int, int] = {}
    for yt in y_true:
        counts[yt] = counts.get(yt, 0) + 1

    k_classes = len(counts)
    props = [cnt / n_test for cnt in counts.values()]
    is_balanced = k_classes >= 2 and all(abs(p - 1.0 / k_classes) <= tol for p in props)

    if is_balanced:
        acc_base = 1.0 / k_classes
        base_mode = "balanced_1_over_K"
    elif baseline_rule == "majority":
        acc_base = max(props)
        base_mode = "majority_class"
    elif baseline_rule == "proportional":
        acc_base = (
            sum(p * p for p in props)
            if declared_prior_sum_sq is None
            else float(declared_prior_sum_sq)
        )
        base_mode = "random_proportional"
    else:
        raise ValueError(f"Unsupported baseline_rule: {baseline_rule}")

    return {
        "R_c": skill_ratio(1.0 - acc_model, 1.0 - acc_base, tol=tol),
        "acc_model": acc_model,
        "acc_base": acc_base,
        "baseline_mode": base_mode,
    }


def _solve_ridge_1d(
    x_rows: list[list[float]],
    y_vals: list[float],
    lam: float = 1e-2,
    rank_tol: float = 1e-12,
) -> list[float]:
    """Fit weights w (intercept at index 0). Raises on a singular normal equation."""
    p = len(x_rows[0]) + 1
    xtx = [[0.0] * p for _ in range(p)]
    xty = [0.0] * p

    for row, y_val in zip(x_rows, y_vals):
        aug_x = [1.0] + [float(v) for v in row]
        for i in range(p):
            xty[i] += aug_x[i] * float(y_val)
            for j in range(p):
                xtx[i][j] += aug_x[i] * aug_x[j]

    for i in range(1, p):
        xtx[i][i] += lam
    xtx[0][0] += 1e-9

    mat = [xtx[r][:] + [xty[r]] for r in range(p)]
    for col in range(p):
        piv = max(range(col, p), key=lambda r: abs(mat[r][col]))
        mat[col], mat[piv] = mat[piv], mat[col]
        if abs(mat[col][col]) <= rank_tol:
            raise ValueError(
                f"Singular normal equation at column {col}; increase ridge_lambda or "
                "remove degenerate features."
            )
        inv = 1.0 / mat[col][col]
        for c in range(col, p + 1):
            mat[col][c] *= inv
        for r in range(p):
            if r != col and mat[r][col] != 0.0:
                fac = mat[r][col]
                for c in range(col, p + 1):
                    mat[r][c] -= fac * mat[col][c]

    return [mat[i][p] for i in range(p)]


def _predict_ridge_1d(weights: list[float], x_rows: list[list[float]]) -> list[float]:
    b0 = weights[0]
    w_rest = weights[1:]
    return [
        b0 + sum(float(row[j]) * w_rest[j] for j in range(len(w_rest)))
        for row in x_rows
    ]


def global_recoverability_mse(
    obs_train: list[list[float]],
    hid_train: list[list[float]],
    obs_test: list[list[float]],
    hid_test: list[list[float]],
    ridge_lambda: float = 1e-2,
    tol: float = 1e-12,
) -> dict[str, Any]:
    """
    v1.2.2 Section 3.4 R_g using out-of-sample MSE against training-set mean baseline.
    """
    n_test = len(hid_test)
    if n_test == 0 or not hid_test[0]:
        return {"R_g": COMPLETE_OBSERVATION, "E_model": None, "E_null": None}

    n_train = len(hid_train)
    if n_train == 0:
        raise ValueError("hid_train must be non-empty when hid_test is non-empty.")

    d_h = len(hid_test[0])
    mean_h = [
        sum(float(row[c]) for row in hid_train) / n_train
        for c in range(d_h)
    ]

    preds_by_col: list[list[float]] = []
    for c in range(d_h):
        y_tr = [float(row[c]) for row in hid_train]
        w = _solve_ridge_1d(obs_train, y_tr, lam=ridge_lambda)
        preds_by_col.append(_predict_ridge_1d(w, obs_test))

    sq_err_model = 0.0
    sq_err_null = 0.0
    for i in range(n_test):
        for c in range(d_h):
            y_true_val = float(hid_test[i][c])
            y_pred_val = preds_by_col[c][i]
            sq_err_model += (y_pred_val - y_true_val) ** 2
            sq_err_null += (mean_h[c] - y_true_val) ** 2

    e_model = sq_err_model / (n_test * d_h)
    e_null = sq_err_null / (n_test * d_h)
    return {
        "R_g": skill_ratio(e_model, e_null, tol=tol),
        "E_model": e_model,
        "E_null": e_null,
    }


def recoverability_retention(
    r_t1: float | str,
    r_t2: float | str,
    tol: float = 1e-12,
) -> dict[str, Any]:
    """
    v1.2.2 Section 4 retention operator eta_X(t1 -> t2) strictly on R_X(t1) > 0
    with explicit sign-crossing handling.
    """
    if not isinstance(r_t1, (int, float)) or math.isnan(r_t1) or r_t1 <= tol:
        return {"status": NON_POSITIVE_BASELINE, "eta": None, "ratio": None}

    if not isinstance(r_t2, (int, float)) or math.isnan(r_t2):
        return {"status": UNDEFINED, "eta": None, "ratio": None}

    r1 = float(r_t1)
    r2 = 0.0 if abs(float(r_t2)) <= tol else float(r_t2)
    ratio = r2 / r1

    if r2 < 0.0:
        return {"status": SIGN_CROSSING, "eta": None, "ratio": ratio}

    return {"status": "retained", "eta": ratio, "ratio": ratio}


def _r2_score(y_true: list[float], y_pred: list[float], eps: float = 1e-12) -> float:
    mean_y = sum(y_true) / len(y_true)
    ss_res = sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred))
    ss_tot = sum((yt - mean_y) ** 2 for yt in y_true)
    if ss_tot <= eps:
        return 0.0
    return 1.0 - (ss_res / ss_tot)


def _make_folds(n: int, k: int, rng: random.Random) -> list[int]:
    """Shuffled fold assignment: returns a fold id in [0, k) for each of n rows."""
    k = min(k, n)
    perm = list(range(n))
    rng.shuffle(perm)
    ids = [0] * n
    for pos, i in enumerate(perm):
        ids[i] = pos % k
    return ids


def _cv_r2_folds(
    x_rows: list[list[float]],
    y_vals: list[float],
    fold_ids: list[int],
    ridge_lambda: float = 1e-2,
) -> float:
    """
    Cross-validated R^2 using a caller-supplied fold assignment, so competing models
    can be compared on identical folds (paired) and bootstrap duplicates can be kept
    inside a single fold (no train/validation leakage).
    """
    if max(fold_ids) < 1:
        return 0.0

    all_true: list[float] = []
    all_pred: list[float] = []

    for fold in range(max(fold_ids) + 1):
        tr = [i for i, g in enumerate(fold_ids) if g != fold]
        va = [i for i, g in enumerate(fold_ids) if g == fold]
        if not tr or not va:
            continue
        w = _solve_ridge_1d(
            [x_rows[i] for i in tr], [y_vals[i] for i in tr], lam=ridge_lambda
        )
        all_pred.extend(_predict_ridge_1d(w, [x_rows[i] for i in va]))
        all_true.extend(y_vals[i] for i in va)

    if len(all_true) < 2:
        return 0.0
    return _r2_score(all_true, all_pred)


def _percentile_ci_95(values: list[float]) -> tuple[float, float]:
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    lo_idx = max(0, min(n - 1, int(math.floor(0.025 * (n - 1)))))
    hi_idx = max(0, min(n - 1, int(math.ceil(0.975 * (n - 1)))))
    return sorted_vals[lo_idx], sorted_vals[hi_idx]


def run_section5_falsification_harness(
    delta_r_samples: list[list[float]],
    delta_z_int_samples: list[list[float]],
    delta_resp_samples: list[float],
    max_rel_energy_diff: float,
    max_rel_state_diff: float,
    measured_pre_perturbation: bool = True,
    identical_perturbation: bool = True,
    rel_energy_tol: float = 0.01,
    rel_state_tol: float = 0.02,
    n_bootstrap: int = 400,
    n_splits: int = 5,
    ridge_lambda: float = 1e-2,
    ci_degeneracy_tol: float = 1e-9,
    practical_floor: float = 0.01,
    zero_var_tol: float = 1e-12,
    seed: int = 20261121,
) -> dict[str, Any]:
    """
    v1.2.2 Section 5 five-condition falsification & control harness.

    CV design: every replicate draws ONE fold assignment over the original rows and
    reuses it for all three models (R-only, Z-only, Z+R), so the incremental
    Delta R^2_CV is paired. Bootstrap rows inherit the fold of their source row, so
    duplicates never straddle train/validation.

    Condition 5 passes only if the 95% percentile CI of incremental Delta R^2_CV has
    a strictly positive lower bound AND the (signed) bootstrap median meets the
    practical-significance floor. A significantly negative Delta R^2 does not pass.
    """
    n = len(delta_resp_samples)
    if not (len(delta_r_samples) == n and len(delta_z_int_samples) == n):
        raise ValueError("Paired sample arrays must have equal length.")

    cond1_ok = (max_rel_energy_diff <= rel_energy_tol) and (max_rel_state_diff <= rel_state_tol)
    cond2_ok = bool(measured_pre_perturbation)
    cond3_ok = bool(identical_perturbation)

    if not (cond1_ok and cond2_ok and cond3_ok):
        return {
            "classification": "inconclusive",
            "reason": "conditions_1_to_3_not_satisfied",
            "bootstrap_ci_method": BOOTSTRAP_CI_METHOD,
            "conditions_1_to_3": {
                "cond1_matching_ok": cond1_ok,
                "cond2_pre_perturbation_ok": cond2_ok,
                "cond3_identical_probe_ok": cond3_ok,
            },
        }

    if n < 6 or n_bootstrap < 20:
        return {
            "classification": "inconclusive",
            "reason": "insufficient_samples_for_bootstrap_cv",
            "bootstrap_ci_method": BOOTSTRAP_CI_METHOD,
            "n_samples": n,
        }

    if statistics.pvariance(delta_resp_samples) <= zero_var_tol:
        return {
            "classification": "inconclusive",
            "reason": "zero_variance_target",
            "bootstrap_ci_method": BOOTSTRAP_CI_METHOD,
        }

    rng = random.Random(seed)
    full_dz_plus_dr = [delta_z_int_samples[i] + delta_r_samples[i] for i in range(n)]

    # Full-data point estimates on one shared (paired) fold assignment.
    full_folds = _make_folds(n, n_splits, rng)
    r2_r_full_cv = _cv_r2_folds(delta_r_samples, delta_resp_samples, full_folds, ridge_lambda)
    r2_z_full_cv = _cv_r2_folds(delta_z_int_samples, delta_resp_samples, full_folds, ridge_lambda)
    r2_zr_full_cv = _cv_r2_folds(full_dz_plus_dr, delta_resp_samples, full_folds, ridge_lambda)
    delta_r2_full_cv = r2_zr_full_cv - r2_z_full_cv

    r2_r_boot: list[float] = []
    delta_r2_boot: list[float] = []

    for _ in range(n_bootstrap):
        boot_idx = [rng.randrange(n) for _ in range(n)]
        orig_folds = _make_folds(n, n_splits, rng)
        b_folds = [orig_folds[i] for i in boot_idx]

        b_dr = [delta_r_samples[i] for i in boot_idx]
        b_dz = [delta_z_int_samples[i] for i in boot_idx]
        b_dz_plus_dr = [full_dz_plus_dr[i] for i in boot_idx]
        b_y = [delta_resp_samples[i] for i in boot_idx]

        r2_only_r = _cv_r2_folds(b_dr, b_y, b_folds, ridge_lambda)
        r2_base_z = _cv_r2_folds(b_dz, b_y, b_folds, ridge_lambda)
        r2_joint_zr = _cv_r2_folds(b_dz_plus_dr, b_y, b_folds, ridge_lambda)

        r2_r_boot.append(r2_only_r)
        delta_r2_boot.append(r2_joint_zr - r2_base_z)

    r2_boot_median = statistics.median(r2_r_boot)
    ci4_lo, ci4_hi = _percentile_ci_95(r2_r_boot)
    ci4_width = ci4_hi - ci4_lo

    if abs(ci4_width) <= ci_degeneracy_tol:
        return {
            "classification": "inconclusive",
            "reason": "degenerate_bootstrap_ci_condition_4",
            "bootstrap_ci_method": BOOTSTRAP_CI_METHOD,
            "cond4_r2_full_cv": r2_r_full_cv,
            "cond4_r2_boot_median": r2_boot_median,
            "cond4_ci95": (ci4_lo, ci4_hi),
        }

    if ci4_lo <= 0.0:
        return {
            "classification": "falsified",
            "reason": "cond4_ci95_lower_bound_not_positive",
            "bootstrap_ci_method": BOOTSTRAP_CI_METHOD,
            "cond4_r2_full_cv": r2_r_full_cv,
            "cond4_r2_boot_median": r2_boot_median,
            "cond4_ci95": (ci4_lo, ci4_hi),
        }

    delta_r2_boot_median = statistics.median(delta_r2_boot)
    ci5_lo, ci5_hi = _percentile_ci_95(delta_r2_boot)
    ci5_width = ci5_hi - ci5_lo

    if abs(ci5_width) <= ci_degeneracy_tol:
        return {
            "classification": "inconclusive",
            "reason": "degenerate_bootstrap_ci_condition_5",
            "bootstrap_ci_method": BOOTSTRAP_CI_METHOD,
            "cond4_r2_full_cv": r2_r_full_cv,
            "cond4_r2_boot_median": r2_boot_median,
            "cond4_ci95": (ci4_lo, ci4_hi),
            "cond5_delta_r2_full_cv": delta_r2_full_cv,
            "cond5_delta_r2_boot_median": delta_r2_boot_median,
            "cond5_ci95": (ci5_lo, ci5_hi),
        }

    ci5_lower_not_positive = ci5_lo <= 0.0
    below_practical_floor = delta_r2_boot_median < practical_floor  # signed, no abs()

    if ci5_lower_not_positive or below_practical_floor:
        return {
            "classification": "hidden_state_dependence",
            "reason": "cond5_or_triggered",
            "bootstrap_ci_method": BOOTSTRAP_CI_METHOD,
            "cond5_triggers": {
                "ci95_lower_not_positive": ci5_lower_not_positive,
                "delta_r2_boot_median_below_floor": below_practical_floor,
            },
            "cond4_r2_full_cv": r2_r_full_cv,
            "cond4_r2_boot_median": r2_boot_median,
            "cond4_ci95": (ci4_lo, ci4_hi),
            "cond5_delta_r2_full_cv": delta_r2_full_cv,
            "cond5_delta_r2_boot_median": delta_r2_boot_median,
            "cond5_ci95": (ci5_lo, ci5_hi),
        }

    return {
        "classification": "irreducible_recoverability_effect",
        "reason": "passed_cond4_and_survived_cond5_control",
        "bootstrap_ci_method": BOOTSTRAP_CI_METHOD,
        "cond4_r2_full_cv": r2_r_full_cv,
        "cond4_r2_boot_median": r2_boot_median,
        "cond4_ci95": (ci4_lo, ci4_hi),
        "cond5_delta_r2_full_cv": delta_r2_full_cv,
        "cond5_delta_r2_boot_median": delta_r2_boot_median,
        "cond5_ci95": (ci5_lo, ci5_hi),
    }


def _run_verification_assertions() -> None:
    # 1. Skill ratio edge cases
    assert skill_ratio(0.0, 0.0) == 1.0
    assert skill_ratio(0.2, 0.0) == UNDEFINED
    assert abs(skill_ratio(1.5, 1.0) - (-0.5)) < 1e-12

    # 2. Structural R_s & uniform graph degeneracy
    unif = [[0.0 if i == j else 0.5 for j in range(3)] for i in range(3)]
    rs_unif = structural_recoverability(unif, unif)
    assert rs_unif["R_s"] == 1.0 and rs_unif["warning"] == "uniform_reference_graph_degeneracy"

    # 3. Directional asymmetry A_t, A_t_signed, scale invariance, and sigma-antisymmetry
    zero_mat = [[0.0] * 3 for _ in range(3)]
    da_zero = directional_asymmetry(zero_mat)
    assert da_zero["A_t"] == UNDEFINED and da_zero["A_t_signed"] == UNDEFINED

    counts_scaled = [[0.0, 100.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    counts_10x = [[v * 10.0 for v in row] for row in counts_scaled]
    orient_pos = [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]
    orient_neg = [[0, -1, 0], [1, 0, 0], [0, 0, 0]]

    a1 = directional_asymmetry(counts_scaled, orient_pos)
    a2 = directional_asymmetry(counts_10x, orient_pos)
    a_flip = directional_asymmetry(counts_scaled, orient_neg)

    assert abs(float(a1["A_t"]) - float(a2["A_t"])) < 1e-12
    assert float(a1["A_t_signed"]) > 0.999999
    assert abs(float(a1["A_t_signed"]) + float(a_flip["A_t_signed"])) < 1e-12

    # 4. Sequential Hamming R_t
    rt_exact = sequential_recoverability_hamming(["A", "B", "C"], ["A", "B", "C"])
    assert rt_exact["R_t"] == 1.0

    # 5. Causal R_c single-class test set -> undefined when imperfect
    rc_single = causal_recoverability([1, 1, 1, 1], [1, 0, 1, 1])
    assert rc_single["R_c"] == UNDEFINED

    # 6. Global R_g complete observation (d_H == 0)
    rg_full_obs = global_recoverability_mse([[1.0]], [], [[1.0]], [])
    assert rg_full_obs["R_g"] == COMPLETE_OBSERVATION

    # 6b. Non-empty hid_test with empty hid_train must raise, not divide by zero
    try:
        global_recoverability_mse([[1.0]], [], [[1.0]], [[1.0]])
        raise AssertionError("expected ValueError for empty hid_train")
    except ValueError:
        pass

    # 7. Retention eta_X domain & sign-crossing
    assert recoverability_retention(0.0, 0.5)["status"] == NON_POSITIVE_BASELINE
    assert recoverability_retention(-0.2, 0.1)["status"] == NON_POSITIVE_BASELINE
    sc = recoverability_retention(0.6, -0.3)
    assert sc["status"] == SIGN_CROSSING and abs(sc["ratio"] - (-0.5)) < 1e-12

    # 8. Zero-variance target -> explicit inconclusive with zero_variance_target reason
    z_var = run_section5_falsification_harness(
        delta_r_samples=[[0.1 * i] for i in range(10)],
        delta_z_int_samples=[[0.05 * i] for i in range(10)],
        delta_resp_samples=[1.5] * 10,
        max_rel_energy_diff=0.001,
        max_rel_state_diff=0.001,
    )
    assert z_var["classification"] == "inconclusive" and z_var["reason"] == "zero_variance_target"

    # 9. hidden_state_dependence: z_int explains y strongly, r correlates with y
    # marginally (Condition 4 passes) but adds no lift once z_int is included.
    rng = random.Random(2026)
    n_pts = 80
    z_list = [[rng.uniform(-2.0, 2.0)] for _ in range(n_pts)]
    r_list = [[0.85 * z_list[i][0] + 0.35 * rng.uniform(-1.0, 1.0)] for i in range(n_pts)]
    y_hidden = [3.0 * z_list[i][0] + 0.02 * rng.gauss(0.0, 1.0) for i in range(n_pts)]

    res_hidden = run_section5_falsification_harness(
        delta_r_samples=r_list,
        delta_z_int_samples=z_list,
        delta_resp_samples=y_hidden,
        max_rel_energy_diff=0.002,
        max_rel_state_diff=0.005,
        n_bootstrap=300,
    )
    assert res_hidden["classification"] == "hidden_state_dependence"
    assert res_hidden["reason"] == "cond5_or_triggered"

    # 10a. Linear positive-path test
    y_lin = [0.4 * z_list[i][0] + 2.5 * r_list[i][0] + 0.05 * rng.gauss(0.0, 0.2) for i in range(n_pts)]
    res_lin = run_section5_falsification_harness(
        delta_r_samples=r_list,
        delta_z_int_samples=z_list,
        delta_resp_samples=y_lin,
        max_rel_energy_diff=0.002,
        max_rel_state_diff=0.005,
        n_bootstrap=300,
    )
    assert res_lin["classification"] == "irreducible_recoverability_effect"

    # 10b. Nonlinear positive-path test (strong cubic component in r beyond z)
    r_indep = [[rng.uniform(-2.0, 2.0)] for _ in range(n_pts)]
    y_nonlin = [
        0.5 * z_list[i][0] + 1.8 * (r_indep[i][0] ** 3) + 0.05 * rng.gauss(0.0, 0.3)
        for i in range(n_pts)
    ]
    res_nonlin = run_section5_falsification_harness(
        delta_r_samples=r_indep,
        delta_z_int_samples=z_list,
        delta_resp_samples=y_nonlin,
        max_rel_energy_diff=0.002,
        max_rel_state_diff=0.005,
        n_bootstrap=300,
    )
    assert res_nonlin["classification"] == "irreducible_recoverability_effect"

    # 11. Leakage sanity: pure noise must never pass Condition 4
    noise_r = [[rng.gauss(0.0, 1.0)] for _ in range(n_pts)]
    noise_z = [[rng.gauss(0.0, 1.0)] for _ in range(n_pts)]
    noise_y = [rng.gauss(0.0, 1.0) for _ in range(n_pts)]
    res_noise = run_section5_falsification_harness(
        delta_r_samples=noise_r,
        delta_z_int_samples=noise_z,
        delta_resp_samples=noise_y,
        max_rel_energy_diff=0.002,
        max_rel_state_diff=0.005,
        n_bootstrap=300,
    )
    assert res_noise["classification"] == "falsified", res_noise

    # 12. Fold pairing: identical inputs + identical folds -> identical CV R^2
    folds = _make_folds(n_pts, 5, random.Random(7))
    a = _cv_r2_folds(r_list, y_lin, folds)
    b = _cv_r2_folds(r_list, y_lin, folds)
    assert a == b

    print("All v1.2.2 default-path & Section 5 harness assertions passed.")


if __name__ == "__main__":
    _run_verification_assertions()


