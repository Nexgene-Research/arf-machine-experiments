# Arf Machine Experiments

Open-source research cases from **Nexgene AI's Arf Machine** — an autonomous research engine for life-science modelling.

## About

[**Nexgene AI**](https://nexgene.ai) builds AI systems that accelerate scientific discovery across the life sciences. **Arf Machine** is our autonomous research engine: it formulates hypotheses, writes and evaluates code, iterates through experiments, and keeps only the ideas that improve held-out metrics.

This repository publishes reproducible case studies from Arf Machine runs. Each case includes the discovered predictor code, benchmark context, results, scientific notes, and figures.

## Cases

| Case | Target | Headline metric | Link |
| --- | --- | --- | --- |
| OpenBind EV-A71 2A Affinity | EV-A71 2A protease | Spearman rho 0.57, RMSE 0.70 | [cases/openbind](cases/openbind) |
| Emerald Bay Drug Sensitivity | Cell-line × drug growth rate | Mean R² 0.44, Pearson 0.69 | [cases/emeraldbay](cases/emeraldbay) |

## Repository Layout

```text
cases/
  openbind/          # First case: OpenBind EV-A71 2A affinity prediction
  emeraldbay/        # Emerald Bay drug sensitivity prediction
```

Each case directory contains its own README, predictor code, benchmark harness, dependencies, and figures.

## License

Code and documentation are released under the [MIT License](LICENSE). Upstream datasets are not redistributed; each case README links to the original data source and its license.
