# Inferring-mechanistic-synergy-OD-data

## Table of contents
* [Introduction](#introduction)
* [Technologies](#technologies)
* [Setup](#setup)
* [Usage](#usage)
* [License](#license)

## Introduction
The repository holds the code for the manuscript Hameed _et al._ (2026), "Inferring antifungal drug synergy from _Candidozyma auris_ optical density data using Bayesian mechanistic modelling".

## Technologies
The code is written in R (v4.5.2) and Stan. Details of the packages and their versions can be found in the [renv.lock](renv.lock) file. 

## Setup
First either clone or download the repository to your machine. To download the correct versions of the R packages used in this project, use the renv package and run:

```
renv::restore()
```

in the R console when prompted. The appropriate version of renv should be automatically installed once this repository has been opened in RStudio. For more details please refer to the [Introduction to renv](https://rstudio.github.io/renv/articles/renv.html) vignette.

The scripts in [`03_validation/`](R/03_validation/) and [`s1_prior_pred/`](R/s1_prior_pred/) are intended to be run using Imperial College London's High Performance Computing (HPC) services. To run these scripts yourself on Imperial's HPC facility, a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) environment with R at v4.5.2 should be created:

```
source miniconda3/bin/activate
conda create -n r452 r-base=4.5.2 -c conda-forge
source activate r452
```

and the correct versions of the R packages can then downloaded in the same manner as detailed above.

## Usage

The code is intended to be run by using the .R scripts in the [`R/`](R) folder:
* [`01_read_data.R`](R/01_read_data.R): Takes the excel spreadsheet of data and stores the data as a tibble in a "data" folder.
* [`02_plot_data.R`](R/02_plot_data.R): Visualises the experimental data.  
* [`03_validation`](R/03_validation/): Conducts k-fold cross validation. 
* [`04_post_pred.R`](R/04_post_pred.R): Performs a posterior predictive check.
* [`05_plot`](R/05_plot/): Holds the plotting scripts used to generate each the plots in the manuscript text.
* [`s1_prior_pred`](R/s1_prior_pred/): Conducts a fake data check using 5 generated fake data sets.
* [`s2_min_ode_diff.R`](R/s2_min_ode_diff.R): Finds the value for the drug-action parameters that have a minimal difference on the models' outputs.
* [`s3_plot_gompertz_params.R`](R/s3_plot_gompertz_params.R): Plots supplementary figure of inferred drug-action parameters in the Gompertz-D-HS model.

All Stan models used in the manuscript are in the [`models/`](models) folder. 

## License
Licensed under the GPLv3 license. See [LICENSE](LICENSE) for more information.