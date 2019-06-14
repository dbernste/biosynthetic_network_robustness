# biosynthetic_network_robustness
This repository contains the code used in Bernstein DB, Dewhirst FE, Segre D. (2019) Metabolic network percolation quantifies biosynthetic capabilities across the human oral microbiome. *eLife*. DOI: 10.7554/eLife.39733.

**Figures/**<br />
All data and scripts necessary to recreat the figures from the original article

**algorithm_functions/**<br />
MATLAB code used to implement out algorithm

**combinatorial_theory_functions/**<br />
Code pertaining to the combinatorial equations we developed to describe our method

**example_models/**<br />
Several metabolic network models in .mat format that are used by example_calculate_PM.m

**example_results/**<br />
A directory to hold the results calculated with example_calculate_PM.m

**metabolic_networks_mat/**<br />
All metabolic network models in .mat format used in the original article.

**example_calculate_PM.m**<br />
Example script to run the functions in algorithm_functions/ and calculate the producibility metric (PM)

