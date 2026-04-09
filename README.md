# AntsNet

**AntsNet** is a comprehensive R package for exploring the mathematical isomorphisms between ant colony decision-making and machine learning. It unifies three previously separate packages into a single toolkit accompanying the research trilogy:

| Part | arXiv | Isomorphism |
|------|-------|-------------|
| I  | [2603.20328](https://arxiv.org/abs/2603.20328) | Random Forests ≅ Ant Colonies (variance reduction via decorrelation) |
| II | [2604.00038](https://arxiv.org/abs/2604.00038) | Boosting ≅ Adaptive Recruitment (bias reduction via sequential reweighting) |
| III | *forthcoming* | Neural Networks ≅ Generational Colony Learning (gradient descent via pheromone evolution) |

## Installation

```r
# From CRAN (when available)
install.packages("AntsNet")

# Development version from GitHub
remotes::install_github("ylevental/IsomorphismSim_Full")
```

## Quick Start

### Part I: Random Forest ≅ Ant Colony

```r
library(AntsNet)

# Simulate an ant colony
sim <- simulate_ant_colony(n_ants = 50, p_explore = 0.3)
cat("Colony chose site", sim$decision, "\n")

# Run variance decomposition experiment
rf_results <- variance_decomposition_experiment(n_replicates = 20)
plot_variance_decomposition(rf_results)

# Direct isomorphism test
iso <- isomorphism_test(n_replicates = 20)
plot_correlation_decay(iso)
```

### Part II: Boosting ≅ Adaptive Recruitment

```r
# AdaBoost
data <- generate_classification_data(n = 200, p = 5, noise = 0.1)
boost_res <- adaboost(data$X, data$y, M = 50)

# Ant Colony Adaptive Recruitment
sq <- c(10, 8, 6, 4, 2)
acar_res <- acar(sq, n_ants = 30, n_waves = 50)

# Compare weight/pheromone evolution
plot_weight_pheromone(boost_res, acar_res, sq)

# Weak learnability experiment
wl <- weak_learnability_experiment(n_replicates = 50)
plot_weak_learnability(wl)
```

### Part III: Neural Network ≅ Generational Colony Learning

```r
# Compare gradient descent and generational pheromone learning
plot_isomorphism(site_qualities = c(10, 7, 5, 4, 3), n_generations = 50)

# Learning curves side by side
plot_learning_curves(n_replicates = 10, n_generations = 50)

# Neural plasticity ↔ colony adaptation
plot_plasticity(n_generations = 100, env_shift_at = 50)
```

### Interactive Exploration

```r
# Shiny app for Part I (decorrelation explorer)
launch_app("part1")

# Shiny app for Part II (boosting explorer)
launch_app("part2")
```

## Package Structure

| Module | Functions | Isomorphism |
|--------|-----------|-------------|
| `part1_*` | `simulate_ant_colony()`, `variance_decomposition_experiment()`, `isomorphism_test()`, ... | RF ≅ Colony |
| `part2_*` | `adaboost()`, `acar()`, `weak_learnability_experiment()`, ... | Boosting ≅ ACAR |
| `part3_*` | `gacl()`, `simple_neural_network()`, `plot_isomorphism()`, ... | NN ≅ GACL |

### Renamed functions (from the original separate packages)

To avoid name collisions in the unified package:

| Original (separate package) | Unified (AntsNet) | Reason |
|----|----|----|
| `generate_data()` (Part I) | `generate_regression_data()` | Different signatures |
| `generate_data()` (Part II) | `generate_classification_data()` | Different signatures |
| `plot_noise_robustness()` (Part II) | `plot_noise_robustness_boost()` | Name collision with Part III |
| `plot_noise_robustness()` (Part III) | `plot_noise_robustness_nn()` | Name collision with Part II |
| `plot_convergence()` (Part II) | `plot_convergence_boost()` | Clarity |
| `noise_experiment()` (Part II) | `noise_experiment_boost()` | Clarity |
| `convergence_experiment()` (Part II) | `convergence_experiment_boost()` | Clarity |

## The Isomorphism in One Equation

All three isomorphisms share the same underlying principle. For an ensemble of *N* units with individual variance σ² and pairwise correlation ρ:

```
Var[ensemble] = ρσ² + (1 − ρ)σ²/N
```

This holds identically for:
- **Random forests**: trees decorrelated by random feature selection (θ = m_try/p)
- **Ant colonies**: ants decorrelated by stochastic exploration (θ = p_explore)

> *The ant colony is a random forest running on biological hardware;
> the random forest is an ant colony running on silicon.*

## Citation

```bibtex
@misc{fokoue2026partI,
  title={Decorrelation, Diversity, and Emergent Intelligence: The Isomorphism
         Between Social Insect Colonies and Ensemble Machine Learning},
  author={Fokou{\'e}, Ernest and Babbitt, Gregory and Levental, Yuval},
  year={2026},
  eprint={2603.20328},
  archivePrefix={arXiv},
  primaryClass={stat.ML}
}

@misc{fokoue2026partII,
  title={Isomorphic Functionalities between Ant Colony and Ensemble Learning:
         Part II --- On the Strength of Weak Learnability and the Boosting Paradigm},
  author={Fokou{\'e}, Ernest and Babbitt, Gregory and Levental, Yuval},
  year={2026},
  eprint={2604.00038},
  archivePrefix={arXiv},
  primaryClass={stat.ML}
}
```

## Authors

- **Yuval Levental** (maintainer) — Center for Imaging Science, RIT
- **Gregory Babbitt** — Gosnell School of Life Sciences, RIT
- **Ernest Fokoué** — School of Mathematics and Statistics, RIT

## License

MIT
