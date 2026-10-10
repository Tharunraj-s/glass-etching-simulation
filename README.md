# Mathematical Modelling of Glass Etching

MSc research project by Tharunraj Singaravelu, University of Limerick.

This curated release contains two- and three-species MATLAB differential-equation models and an exploratory Python lattice notebook. Original academic submissions and revision history remain private.

## Run

In MATLAB, run `simulate_two_species_etching_combined.m` or `three_species_etching.m`. These scripts use `ode45` and produce surface-fraction and depth-profile plots.

For Python: `pip install -r requirements.txt`, then `jupyter notebook lattice_etching.ipynb`. Start with a small lattice, such as 5 rows, 8 columns and 10 time steps.

## Scope and limitations

The MATLAB models explore selected parameters. They are not experimentally calibrated industrial models. The Python notebook is an exploratory coursework prototype, not a validated equivalent of the ODE model. Its timing rule uses absolute timestep minus layer index; this can miss etching events for layers exposed later. Notebook outputs were cleared for this release. MATLAB execution has not been verified in the release environment.
