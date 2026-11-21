# Formal Recoverability Framework — Specification v1.2.2 (Frozen)

## 1. Foundational Architecture & Formal Objects

* **System & Observable Space:** Let $\Sigma_t$ denote the relational state configuration of a system at time $t$, and $\mathcal{O}_t = \Pi(\Sigma_t)$ denote the observable representation accessible at time $t$.
* **Axiom Zero (Non-Closure Primitive):** Complete nothingness or total collapse into undifferentiated identity cannot obtain, leaving a persistent non-zero structural residual ($\Phi_{\text{tot}} > 0$).
* **Foundational Proposition:** $\text{Persistence} \not\Rightarrow \text{Recoverability}$. Persistence states that a configuration or excitation continues to exist across transformation. Recoverability quantifies how accurately a target relational structure $\mathcal{T}$ can be reconstructed from post-transformation observables $\mathcal{O}_t$ relative to an uninformed null baseline.
* **Recoverability Vector:** $\mathbf{R}$ is the four-component tuple of task-specific skill scores: $\mathbf{R} = (R_s, R_t, R_c, R_g)$.
* **Dependency Chain (Conceptual Spine):** $\text{Axiom Zero} \to \text{Persistence} \to \text{Differentiation} \to \text{Relations} \to \mathbf{R} = (R_s, R_t, R_c, R_g)$.
* **Unproven Empirical Hypotheses:** Downstream mappings from $\mathbf{R}$ to accumulated memory $M(t)$, historical response $G_{ij}(\mathbf{q}, \mathbf{R}, \mathbf{M})$, geometry, and spatial capacity are empirical hypotheses subject to Section 5, not consequences deduced from definitions alone.

---

## 2. General Normalized Skill Ratio & Edge-Case Rules

For any reconstruction task $X \in \{s, t, c, g\}$, let $E_{\text{model}} \ge 0$ be the reconstruction error of a specified estimator $\mathcal{F}$ and $E_{\text{null}} \ge 0$ be the reconstruction error of a specified uninformed baseline, evaluated on the **same held-out samples** using the **same non-negative loss function** (identical units and scale).

$$R_X = \begin{cases} 1, & \text{if } E_{\text{null}} = 0 \text{ and } E_{\text{model}} = 0 \\ \text{undefined}, & \text{if } E_{\text{null}} = 0 \text{ and } E_{\text{model}} > 0 \\ 1 - \dfrac{E_{\text{model}}}{E_{\text{null}}}, & \text{if } E_{\text{null}} > 0 \end{cases}$$

### Domain & Reporting Rules
* **Codomain:** $R_X \in (-\infty, 1] \cup \{\text{undefined}\}$.
* **No artificial clipping:** Negative values ($R_X < 0$, worse than null baseline) must be preserved in primary reporting.
* **Zero-baseline convention:** When $E_{\text{null}} = 0$ and $E_{\text{model}} = 0$, $R_X = 1$ by convention. When $E_{\text{null}} = 0$ and $E_{\text{model}} > 0$, report `undefined` (baseline-degenerate comparison).
* **Non-negativity constraint:** Every discrepancy or loss $E$ must satisfy $E \ge 0$ by construction.

---

## 3. Component Specifications

### 3.1 Structural Recoverability ($R_s$)
Measures reconstruction fidelity of a reference weighted or binary relational graph $C^{\text{ref}} \in \mathbb{R}_{\ge 0}^{N \times N}$ (zero diagonal, $C_{ii} = 0$) from inferred graph $\hat{C} \in \mathbb{R}_{\ge 0}^{N \times N}$.

* **Default Error Metric ($L_1$ entrywise distance):**
  $$E_{\text{model}}^{(s)} = \sum_{i \neq j} \bigl|\hat{C}_{ij} - C^{\text{ref}}_{ij}\bigr|, \qquad E_{\text{null}}^{(s)} = \sum_{i \neq j} \bigl|C^{\text{null}}_{ij} - C^{\text{ref}}_{ij}\bigr|$$
* **Default Null Graph ($C^{\text{null}}$):** Uniform density-matched baseline where every off-diagonal entry equals the mean off-diagonal weight of $C^{\text{ref}}$: $C^{\text{null}}_{ij} = \frac{1}{N(N-1)} \sum_{u \neq v} C^{\text{ref}}_{uv}$ for $i \neq j$. For uniform reference graphs where all off-diagonal entries of $C^{\text{ref}}$ are identical, $C^{\text{null}} = C^{\text{ref}}$ makes $R_s$ degenerate; use a declared alternative null (such as the empty graph $C^{\text{null}}_{ij} = 0$ or degree-preserving randomization) in that case.
* **Score:** Apply Section 2 to $(E_{\text{model}}^{(s)}, E_{\text{null}}^{(s)})$ to obtain $R_s$.

---

### 3.2 Temporal & Sequential Recoverability ($R_t$)

#### A. Auxiliary Diagnostic — Scale-Invariant Directional Asymmetry ($A_t$)
Let $N_{i \to j} \ge 0$ be observed directed transition counts or fluxes for $i \neq j$. Define unordered edge flux $F_{ij} = N_{i \to j} + N_{j \to i}$ for $i < j$, and total trajectory flux $F_{\text{tot}} = \sum_{i < j} F_{ij} = \sum_{u \neq v} N_{u \to v}$.

* **Zero-flux guard:** If $F_{\text{tot}} = 0$, report $A_t = \text{undefined}$ and $A_t^{\text{signed}} = \text{undefined}$.
* **Scale-relative regularization:** For $F_{\text{tot}} > 0$, let $\bar{F}$ be the mean of $F_{ij}$ over all unordered pairs $i < j$ with $F_{ij} > 0$, and set $\epsilon = 10^{-9}\,\bar{F} > 0$.
* **Unsigned flux-weighted asymmetry ($A_t \in [0, 1)$):**
  $$\Delta_{ij} = \frac{\bigl|N_{i \to j} - N_{j \to i}\bigr|}{F_{ij} + \epsilon}, \qquad A_t = \sum_{i < j,\, F_{ij} > 0} \frac{F_{ij}}{F_{\text{tot}}}\,\Delta_{ij} = \frac{1}{F_{\text{tot}}} \sum_{i < j,\, F_{ij} > 0} \frac{F_{ij}\,\bigl|N_{i \to j} - N_{j \to i}\bigr|}{F_{ij} + \epsilon}$$
  *Codomain note:* Because $\epsilon > 0$, $A_t \in [0, 1)$ strictly; $1$ is a limiting bound as $\epsilon \to 0$, not an attained finite value.
* **Signed orientation alignment ($A_t^{\text{signed}} \in (-1, 1) \subset [-1, 1]$):** Given a reference orientation matrix $\sigma_{ij} \in \{-1, 0, +1\}$ satisfying $\sigma_{ji} = -\sigma_{ij}$ (with $\sigma_{ij} = 0$ on symmetric/unoriented reference edges), let $F_{\text{orient}} = \sum_{i < j,\, \sigma_{ij} \neq 0} F_{ij}$. If $F_{\text{orient}} = 0$, define $A_t^{\text{signed}} = 0$; otherwise:
  $$A_t^{\text{signed}} = \frac{1}{F_{\text{orient}}} \sum_{i < j,\, \sigma_{ij} \neq 0} F_{ij}\,\sigma_{ij}\,\frac{N_{i \to j} - N_{j \to i}}{F_{ij} + \epsilon}$$
  *Codomain note:* For finite $\epsilon > 0$, $A_t^{\text{signed}} \in (-1, 1)$ strictly; the endpoints $+1$ and $-1$ are asymptotic limits as $\epsilon \to 0$ when all oriented edges align with or oppose $\sigma_{ij}$, not attained finite values.

#### B. Primary Metric — Sequential Recoverability ($R_t$)
Let $\tau_{\text{true}} = (e_1, e_2, \dots, e_L)$ be the true ordered sequence of length $L \ge 1$, and $\hat{\tau} = (\hat{e}_1, \dots, \hat{e}_L)$ be the sequence reconstructed by estimator $\mathcal{F}$.

* **Default Order Discrepancy ($D_{\text{ord}}$ — Normalized step-indexed Hamming distance):**
  $$D_{\text{ord}}(\tau_1, \tau_2) = \frac{1}{L} \sum_{k=1}^{L} \mathbb{I}\!\left(\tau_{1,k} \neq \tau_{2,k}\right) \in [0, 1]$$
  *(Kendall-tau or Levenshtein edit distance may replace Hamming distance if declared.)*
* **Default Null Baseline ($E_{\text{null}}^{(t)}$):** Expected discrepancy under uniform random permutations $\pi$ of the true transition multiset $\tau_{\text{true}}$:
  * *For step-indexed Hamming distance only:* $E_{\text{null}}^{(t)} = 1 - \sum_{e} q_e^2$ holds analytically, where $q_e = n_e / L$ is the empirical frequency of symbol $e$ in $\tau_{\text{true}}$.
  * *For Kendall-tau, Levenshtein, or any non-Hamming discrepancy:* The analytic collision formula is invalid; $E_{\text{null}}^{(t)} = \frac{1}{M_{\text{perm}}} \sum_{m=1}^{M_{\text{perm}}} D_{\text{ord}}(\pi_m(\tau_{\text{true}}), \tau_{\text{true}})$ must be computed via $M_{\text{perm}} \ge 1000$ fixed-seed shuffles.
* **Score:** With $E_{\text{model}}^{(t)} = D_{\text{ord}}(\hat{\tau}, \tau_{\text{true}})$, apply Section 2 to $(E_{\text{model}}^{(t)}, E_{\text{null}}^{(t)})$ to obtain $R_t$.

---

### 3.3 Causal Recoverability ($R_c$)
Measures out-of-sample identification of generating historical regime $y \in \{1, \dots, K\}$ ($K \ge 2$) on a held-out test set of size $N_{\text{test}}$ via a declared classifier $\mathcal{F}$ trained strictly on training data.

* **Error Definitions:** Let $\text{Acc}_{\text{model}} = \frac{1}{N_{\text{test}}} \sum_{n=1}^{N_{\text{test}}} \mathbb{I}(\hat{y}_n = y_n)$ be held-out accuracy and $\text{Acc}_{\text{base}}$ be held-out chance accuracy. Define $E_{\text{model}}^{(c)} = 1 - \text{Acc}_{\text{model}}$ and $E_{\text{null}}^{(c)} = 1 - \text{Acc}_{\text{base}}$, yielding:
  $$R_c = \begin{cases} 1, & \text{if } \text{Acc}_{\text{base}} = 1 \text{ and } \text{Acc}_{\text{model}} = 1 \\ \text{undefined}, & \text{if } \text{Acc}_{\text{base}} = 1 \text{ and } \text{Acc}_{\text{model}} < 1 \\ \dfrac{\text{Acc}_{\text{model}} - \text{Acc}_{\text{base}}}{1 - \text{Acc}_{\text{base}}}, & \text{if } \text{Acc}_{\text{base}} < 1 \end{cases}$$
* **Single-class test guard:** Although $K \ge 2$ candidate regimes are required by design, if a degenerate held-out split contains only one class ($\text{Acc}_{\text{base}} = 1$), $E_{\text{null}}^{(c)} = 0$ and $R_c$ evaluates to `undefined` whenever $\text{Acc}_{\text{model}} < 1$.
* **Default Chance Baseline ($\text{Acc}_{\text{base}}$):**
  * Balanced test classes ($p_k = 1/K$ for all $k$): $\text{Acc}_{\text{base}} = 1/K$.
  * Imbalanced test classes with test proportions $\{p_1, \dots, p_K\}$: **default** $\text{Acc}_{\text{base}} = \max_{k} p_k$ (majority-class prior on the test set; random-proportional prior $\sum_k p_k^2$ permitted only if explicitly declared).
* **Mandatory Disclosure:** Report classifier family $\mathcal{F}$, feature scaling/selection pipeline (fit on training folds only), and train/test split protocol.

---

### 3.4 Global Recoverability ($R_g$)
Partition observables into visible subset $O \subset \Sigma$ and non-empty hidden subset $H = \Sigma \setminus O$ ($|H| = d_H \ge 1$). An estimator $\mathcal{F}$ is fit strictly on $(O_{\text{train}}, H_{\text{train}})$ and evaluated on held-out $(O_{\text{test}}, H_{\text{test}})$ of size $N_{\text{test}}$.

* **Complete-observation guard:** If $O = \Sigma$ ($d_H = 0$), $R_g$ is not evaluated; report `complete_observation`.
* **Default Loss & Baseline (MSE with training-mean predictor):**
  $$E_{\text{model}}^{(g)} = \frac{1}{N_{\text{test}}\, d_H} \sum_{n=1}^{N_{\text{test}}} \left\|\hat{H}_n(O_{\text{test}, n}) - H_{\text{test}, n}\right\|_2^2, \qquad E_{\text{null}}^{(g)} = \frac{1}{N_{\text{test}}\, d_H} \sum_{n=1}^{N_{\text{test}}} \left\|\bar{H}_{\text{train}} - H_{\text{test}, n}\right\|_2^2$$
  where $\bar{H}_{\text{train}} = \frac{1}{N_{\text{train}}} \sum_{m=1}^{N_{\text{train}}} H_{\text{train}, m}$. *(If MAE is declared instead of MSE, $\bar{H}_{\text{train}}$ must be replaced by the training-set coordinatewise median.)*
* **Score:** Apply Section 2 to $(E_{\text{model}}^{(g)}, E_{\text{null}}^{(g)})$ to obtain $R_g$.
* **Mandatory Disclosure:** Report estimator $\mathcal{F}$, regularization hyperparameters, feature selection/standardization (fit on training data only), and split protocol.

---

## 4. Temporal Memory Operators (Distinct from $\mathbf{R}$)

1. **Instantaneous recoverability:** $\mathbf{R}(t) = \bigl(R_s(t), R_t(t), R_c(t), R_g(t)\bigr)$.
2. **Accumulated recoverability (kernel-parameterized family):**
   $$\mathbf{M}_{\text{acc}}(t; K) = \int_{t_0}^{t} K(t, \tau)\,\mathbf{R}(\tau)\,d\tau$$
   Because $\mathbf{M}_{\text{acc}}$ depends on kernel $K(t, \tau)$, it is a parameterized functional family rather than a single invariant; **default kernel:** exponential decay $K(t, \tau) = e^{-(t - \tau)/T_m}$ with declared memory timescale $T_m > 0$. Unweighted cumulative integration $\mathbf{M}_{\text{acc}}(t; T_m \to \infty) = \int_{t_0}^{t} \mathbf{R}(\tau)\,d\tau$ is the asymptotic limit as $T_m \to \infty$, not an attained finite parameter.
3. **Positive-branch recoverability retention & sign-crossing rule ($t_1 \to t_2$):**
   $$\eta_X(t_1 \to t_2) \text{ is defined only on the positive baseline domain } R_X(t_1) > 0:$$
   $$\eta_X(t_1 \to t_2) = \begin{cases} \dfrac{R_X(t_2)}{R_X(t_1)}, & \text{if } R_X(t_1) > 0 \text{ and } R_X(t_2) \ge 0 \\ \text{sign\_crossing} \;\bigl(R_X(t_2)/R_X(t_1)\bigr), & \text{if } R_X(t_1) > 0 \text{ and } R_X(t_2) < 0 \\ \text{undefined (non\_positive\_baseline)}, & \text{if } R_X(t_1) \le 0 \text{ or } R_X(t_1) \text{ is undefined} \end{cases}$$
   *Domain & aggregation rule:* When $R_X(t_1) \le 0$ (including $R_X(t_1) = 0$), retention is `undefined` because no positive baseline skill exists to retain. When $R_X(t_1) > 0$ and $R_X(t_2) < 0$, recoverability has crossed below the null baseline; report as `sign_crossing` (logging the negative ratio separately) rather than a retention fraction, and **do not** average or aggregate $\eta_X$ across sign-crossing events.

---

## 5. Recoverability Independence Hypothesis — Falsification Protocol

To test whether $\mathbf{R}(t_f)$ provides independent predictive power for post-perturbation dynamics beyond standard macroscopic and observable state variables, **all five conditions** below must hold. If any condition among 1–3 is violated, if the target divergence $\Delta_{\text{resp}}$ has zero variance, or if the bootstrap confidence interval in Condition 4 or Condition 5 is degenerate (zero width due to insufficient held-out sample size), the test is classified as **inconclusive / improperly controlled** rather than confirmatory or falsified. Bootstrap CI method (default: 95% percentile bootstrap) must be disclosed.

1. **Macroscopic & observable state matching:** Across comparison groups $A$ and $B$ at $t_f$, total dissipated energy $E_{\text{diss}}$, entropy production proxy, and observable state vector $\mathbf{q}(t_f)$ are matched within pre-declared tolerance $\delta_{\text{match}}$ ($\max_k |q_{A,k} - q_{B,k}| / \text{scale}_k \le 0.02$ and $|E_{\text{diss}, A} - E_{\text{diss}, B}| / E_{\text{diss}, A} \le 0.01$), or adjusted via out-of-sample covariate residualization.
2. **Pre-perturbation measurement:** $\mathbf{R}_A(t_f)$ and $\mathbf{R}_B(t_f)$ are evaluated strictly using data up to $t_f$, before any probe perturbation is applied.
3. **Identical probe perturbation:** Both groups undergo an identical external perturbation $\delta f(t)$ over $t \in (t_f, t_f + \tau]$ with identical step size and boundary conditions.
4. **Operational out-of-sample prediction criterion:** Let $\Delta_{\text{resp}} = \|\mathbf{y}_A(t_f + \tau) - \mathbf{y}_B(t_f + \tau)\|_2$ be post-probe trajectory divergence across held-out paired runs. $\Delta \mathbf{R}(t_f) = \mathbf{R}_A(t_f) - \mathbf{R}_B(t_f)$ is scored as predictive **if and only if** an out-of-sample regression predicting $\Delta_{\text{resp}}$ from $\Delta \mathbf{R}(t_f)$ yields a non-degenerate 95% bootstrap confidence interval for cross-validated $R^2_{\text{CV}}$ whose lower bound is strictly greater than $0$. If the CI is non-degenerate and its lower bound is $\le 0$, the independence hypothesis is **falsified** on that system.
5. **Internal-variable control criterion (disjunctive OR with practical-significance floor):** Re-fit the held-out prediction of $\Delta_{\text{resp}}$ conditioning jointly on full internal state differences $\Delta \mathbf{z}_{\text{int}}(t_f)$ and $\Delta \mathbf{R}(t_f)$. Combine the following two rules with logical **OR** (either rule alone triggers the null classification, where the practical-significance floor intentionally overrides statistical significance alone): if adding $\Delta \mathbf{R}(t_f)$ yields an incremental $\Delta R^2_{\text{CV}}$ whose 95% bootstrap confidence interval includes $0$ **OR** whose bootstrap median magnitude satisfies $|\text{median}_{\text{boot}}(\Delta R^2_{\text{CV}})| < 0.01$, the divergence is classified as **standard hidden-state dependence** rather than an irreducible recoverability effect.
