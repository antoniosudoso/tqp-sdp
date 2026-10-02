
# SDP-B&B for Ternary Quadratic Problems

MATLAB implementation of a semidefinite-programming branch-and-bound framework for three ternary quadratic optimization problems:

1. **QUTO** — unconstrained ternary quadratic optimization;
2. **TQP-Linear** — ternary quadratic optimization with the constraint $\mathbf{1}^\top x=0$;
3. **TQP-Ratio** — optimization of a ratio of two ternary quadratic functions.

The repository also contains the benchmark instances used for the three problem classes.

## Problem variants

All decision variables satisfy $x\in \lbrace{-1,0,1 \rbrace}^n$.

| Variant | Formulation | Main solver | Experiment script |
|---|---|---|---|
| QUTO | $\min x^\top Qx+c^\top x$ | `solve_ternary_mosek_unc.m` | `run_SDP_TQP_QUTO.m` |
| TQP-Linear | $\min x^\top Qx+c^\top x$ subject to $\mathbf{1}^\top x=0$ | `solve_ternary_mosek_sum.m` | `run_SDP_TQP_Linear.m` |
| TQP-Ratio | $\min (x^\top Ax+a^\top x+a_0)/(x^\top Bx+b^\top x+b_0)$ | `solve_ternary_mosek_ratio.m` | `run_SDP_TQP_Ratio.m` |


## Method

The three solvers combine:

- a problem-specific SDP relaxation;
- iterative separation of triangle, RLT, split, pentagonal, and heptagonal inequalities;
- Variable Neighborhood Search (VNS) primal heuristics;
- ternary branching, which fixes a selected variable to $-1$, $0$, or $1$;

## Repository structure

```text
.
├── B&B-SDP/
│   ├── solve_ternary_mosek_unc.m
│   ├── solve_ternary_mosek_sum.m
│   ├── solve_ternary_mosek_ratio.m
│   ├── heuristics/       # VNS heuristics and multistart procedures
│   ├── inequalities/     # cutting-plane generation and separation
│   └── utils/            # branching and selection utilities
├── instancesBW1/         # Type-1 QUTO/TQP-Linear data
├── instancesBW2/         # Type-2 QUTO/TQP-Linear data
├── instancesBW3/         # Type-3 QUTO/TQP-Linear data
├── instancesRATIO/       # TQP-Ratio data
├── run_SDP_TQP_QUTO.m
├── run_SDP_TQP_Linear.m
├── run_SDP_TQP_Ratio.m
└── compute_dinkelbach.m           # Dinkelbach helper routine
```

## Requirements

- MATLAB
- YALMIP

The implementation reported in the paper uses MATLAB R2025a, MOSEK 11.0, and GUROBI 12.0.2.

## Setup

Clone the repository and start MATLAB in its root directory. Add the source folders and solver interfaces to the MATLAB path:

```matlab
addpath(genpath('./B&B-SDP'));
addpath(genpath('<path-to-mosek>')); % location of your MOSEK installation
addpath(genpath('<path-to-yalmip>')); % location of your YALMIP installation
```

## Running the SDP-B&B solvers

Select the instances and parameter values near the beginning of the corresponding script, then run one of:

```matlab
run('run_SDP_TQP_QUTO.m');    % QUTO
run('run_SDP_TQP_Linear.m');    % TQP-Linear
run('run_SDP_TQP_Ratio.m');  % TQP-Ratio
```

Each script accepts an explicit list of triples:

```matlab
instances = [
    80 25 1;
    80 50 2;
    90 75 3
];
```

The columns identify the dimension, generator parameter, and random replicate. Set `type` to `"BW1"`, `"BW2"`, or `"BW3"` for QUTO and TQP-Linear experiments.

### Main parameters

| Parameter | Description |
|---|---|
| `time_limit` | Maximum B&B running time in seconds |
| `max_nodes` | Maximum number of B&B nodes, including the root |
| `opt_gap` | Relative optimality-gap tolerance |
| `sdp_verbose` | YALMIP/MOSEK output level |
| `n_threads` | Number of MOSEK threads |
| `cp_max_iter` | Maximum iterations in each cutting-plane phase |
| `cp_eps_ineq` | Threshold to identify a violated cut |
| `cp_max_sep_ineq` | Maximum TRI/RLT/SPLIT candidates examined |
| `cp_max_new_tri` | Maximum combined TRI/RLT/SPLIT cuts added per iteration |
| `cp_max_new_penta` | Maximum combined PENTA/EPTA cuts added per iteration |
| `cp_ntimes_penta` | Multistart iterations used for pentagonal separation |
| `cp_ntimes_epta` | Multistarts iterations for heptagonal separation |

## Instances

Files follow the convention:

```text
<n>_<parameter>_<replicate>.mat
```

The second filename component takes values `25`, `50`, or `75`; its interpretation depends on the generator. Each configuration has three replicates.

| Folder | Type | Dimensions available |
|---|---|---:|
| `instancesBW1` | Type-1 | 60, 70, ..., 120 |
| `instancesBW2` | Type-2 | 60, 70, ..., 120 |
| `instancesBW3` | Type-3 | 60, 70, ..., 120 |
| `instancesRATIO` | Ratio | 60, 70, ..., 120 |

The BW files contain:

```matlab
Q  % n-by-n quadratic coefficient matrix
c  % n-by-1 linear coefficient vector
```

They are shared by the QUTO and TQP-Linear experiments. The QUTO scripts replace the diagonal of `Q` by its absolute value after loading an instance.

The ratio files contain a scalar structure named `data` with:

```matlab
data.n
data.A, data.a, data.a0   % numerator
data.B, data.b, data.b0   % denominator
```

## Output

The principal fields returned by each SDP-B&B solver are:

```matlab
result.best_ub       % best feasible objective value
result.best_x_ub     % corresponding ternary solution
result.best_lb       % best lower bound maintained by the solver
result.gap_bb        % final relative B&B gap
result.nodes         % number of processed nodes, including the root
result.time_bb       % total B&B time in seconds
```

## Citation

If you use this code or the benchmark instances, please cite:

> F. de Meijer, V. Piccialli, R. Sotirov, and A. M. Sudoso, “Beyond Binarity: Semidefinite Programming for Ternary Quadratic Problems” 2026. arXiv preprint arXiv:2603.28979

```bibtex
@article{deMeijer2026BeyondBinarity,
  author = {Frank de Meijer and Veronica Piccialli and Renata Sotirov and Antonio M. Sudoso},
  title  = {Beyond Binarity: Semidefinite Programming for Ternary Quadratic Problems},
  year   = {2026},
  note   = {arXiv preprint arXiv:2603.28979}
}
```
