
# Complete Mathematical Model and GAMS Implementation

## Overview

This repository provides the complete mathematical formulation and a self-contained GAMS implementation associated with the following article:

K. Fathollahzadeh, M. Saeedi, M. Ghasempour Anaraki, and M. Rabiee,
“Sustainable cast iron supply chain network design: Robust multi-objective optimization with scenario reduction via genetic algorithm,”
*International Journal of Production Economics*, Vol. 289, Article 109730, 2025.
https://doi.org/10.1016/j.ijpe.2025.109730

The materials are provided to help readers understand the complete structure of the model and reproduce its computational implementation. Additional explanations and input-data rules are included to clarify the modeling assumptions and the order in which the numerical parameters are prepared and passed to the optimization model.

## Repository Files

### `Complete_Mathematical_Model.pdf`

This document presents the complete mathematical specification of the model, including:

* Sets and indices
* Parameters and measurement units
* Decision variables
* Economic, environmental, and social objective functions
* Complete model constraints
* Modeling assumptions
* Input-data admissibility rules
* Derived inventory and activation bounds
* Initialization and model-generation procedures

The input-data relations presented in Part II are evaluated before optimization. They are used to validate or calculate fixed parameter values and do not constitute additional constraints on the decision variables.

### `IJPE-Code.gms`

This file contains a self-contained computational implementation of the mathematical model with:

* One scenario
* Three planning periods
* Embedded numerical input data
* Economic cost minimization
* Environmental-impact minimization
* Employment maximization
* Payoff-table construction
* AUGMECON-based Pareto analysis
* Input-data validation before model generation
* Export of objective values, Pareto solutions, and selected decisions

The scenario index is retained in all scenario-dependent parameters and variables, even though the supplied computational instance contains only one scenario. This preserves consistency with the complete mathematical formulation and allows the code to be extended to multiple scenarios.

## Computational Scope

The supplied code is a compact one-scenario, three-period implementation designed for model verification, transparency, and reproducibility.

It should not be interpreted as an exact reproduction of the complete multi-scenario numerical case study reported in the published article. The numerical values included in the code form a complete and internally consistent computational instance.

All numerical inputs required for this implementation are embedded directly in the GAMS file. No external spreadsheet or data file is required.

## Important Modeling Features

The implementation follows these principal modeling conditions:

* Demand is aggregated over product types for each customer group.
* Product types remain distinguished in the forward-flow and material-composition relations.
* Each customer group has one scenario-level feasibility-margin variable.
* The feasibility-margin variables do not carry demand-point or period indices.
* Period one is treated as the initialization period.
* Customer-to-collection flows are fixed to zero in period one.
* Reverse and recovery operations begin after the initialization period.
* Product returns represent exogenous return potential based on an installed-base proxy.
* Returns are not generated from current-period fulfilled sales.
* Iron is represented as one homogeneous material.
* IoT alternatives are mutually exclusive at each candidate facility.
* At most one IoT alternative can be selected for each candidate facility.
* Inventory and inbound-flow bounds are calculated and validated before model generation.
* Input-data requirements are checked before any model equation is generated.

## Input-Data Preparation

The GAMS implementation follows the sequence below:

1. Load the embedded numerical input data.
2. Verify parameter domains and nonnegativity conditions.
3. Verify the scenario probability.
4. Verify that the return-allocation proportions sum to one.
5. Verify the product-composition identity.
6. Calculate the inventory upper bounds.
7. Calculate the facility-activation bounds.
8. Fix the first-period customer-to-collection flows to zero.
9. Generate the optimization equations.
10. Construct the payoff table.
11. Execute the AUGMECON procedure.
12. Export the resulting Pareto solutions and decisions.

If a required input-data condition is violated, the program terminates before the solver is called.

## Software and Tested Environment

The computational model was executed and verified using:

* **GAMS:** Version 24.1.2, build r40979
* **Solver interface:** GAMS/CPLEX
* **Operating system:** Microsoft Windows, 64-bit
* **Model class:** Mixed-Integer Programming (MIP)
* **Solver status:** Normal completion
* **Model status:** Optimal

A licensed GAMS installation with access to a MIP-capable solver is recommended. The supplied implementation explicitly selects CPLEX.

The reported version is the environment in which the computational results were verified. The code is expected to run with newer compatible GAMS and CPLEX installations, but these versions have not been formally tested. Solver time or alternative optimal flow patterns may differ because of solver versions, numerical tolerances, and parameter settings.

## Running the Model

Download the repository and run the following command from the repository directory:

```bash
gams IJPE-Code.gms lo=2
```

The model can also be opened and executed directly through GAMS Studio or the GAMS IDE.

## Generated Results

The implementation reports or exports:

* Model status
* Solver status
* Solution time
* Maximum infeasibility
* Individual-objective anchor solutions
* Payoff-table values
* Pareto objective values
* Selected IoT alternatives
* Facility-activation decisions
* Material and product flows
* Inventory levels
* Demand-feasibility margins

Because multiple operational solutions may have identical objective values, different compatible solver versions may return different flow patterns while preserving the same objective values and feasibility status.

## Citation

If the mathematical formulation or computational implementation is used in subsequent research, please cite the original article:

```bibtex
@article{fathollahzadeh2025castiron,
  title   = {Sustainable cast iron supply chain network design:
             Robust multi-objective optimization with scenario reduction
             via genetic algorithm},
  author  = {Fathollahzadeh, Kasra and Saeedi, Mehran and
             Ghasempour Anaraki, Matin and Rabiee, Meysam},
  journal = {International Journal of Production Economics},
  volume  = {289},
  pages   = {109730},
  year    = {2025},
  doi     = {10.1016/j.ijpe.2025.109730}
}
```

## Copyright and Use

The published article remains subject to the copyright and reuse conditions of its publisher. The files in this repository are author-prepared mathematical and computational materials supplied for transparency and reproducibility.

No separate reuse license is granted through this repository unless a dedicated `LICENSE` file is subsequently added.
