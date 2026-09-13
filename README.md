# Mathematical Model and GAMS Implementation

Companion repository to

> K. Fathollahzadeh, M. Saeedi, M. Ghasempour Anaraki, and M. Rabiee,
> "Sustainable cast iron supply chain network design: Robust multi-objective
> optimization with scenario reduction via genetic algorithm,"
> *International Journal of Production Economics*, 289, 109730 (2025).
> https://doi.org/10.1016/j.ijpe.2025.109730

and to its corrigendum.

## What this repository is for

Here we provide the corrected mathematical model, the GAMS code, and all
assumptions, modeling choices and reproducibility notes, to improve the
transparency and reproducibility of the work. The PDF sets out the complete
formulation together with the rules used to prepare the input data; the GAMS
file is a runnable implementation of that formulation.

**Note on equation numbers.** The PDF numbers its equations (1)–(41) and the
data-preparation relations (P1)–(P16). This numbering is independent of the
one used in the article; the corrigendum gives the correspondence between the
two.

## Files

### `Full Model Specification and Reproducibility Notes.pdf`

The full model specification:

- sets, indices, parameters (with units) and decision variables;
- the economic, environmental and social objective functions;
- all constraints;
- the modeling assumptions;
- Part II: input-data admissibility checks, the derived inventory and
  activation bounds, and the order in which everything is evaluated before the
  solver is called.

The relations in Part II are computed before optimization. They validate the
data or produce fixed parameter values; they are not additional constraints on
the decision variables.

### `GAMS For Fathollahzadeh, et al (2025).gms`

A self-contained GAMS model with the input data embedded, so nothing else
needs to be downloaded. It covers one scenario and three planning periods and
runs, in order: data validation, the three single-objective solves for the
payoff table, and an AUGMECON sweep to produce a Pareto set.


## Software and Tested Environment

The computational model was executed and verified using:

* **GAMS:** Version 24.1.2, build r40979
* **Solver interface:** GAMS/CPLEX
* **Operating system:** Microsoft Windows, 64-bit
* **Model class:** Mixed-Integer Programming (MIP)
* **Solver status:** Normal completion
* **Model status:** Optimal

