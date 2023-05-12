# Escape calculator for the SARS-CoV-2 RBD

This repository implements the latest version of the escape calculator originally described [Greaney, Starr, & Bloom, Virus Evolution, 8:veac021 (2022)](https://doi.org/10.1093/ve/veac021).

The interactive escape calculator itself is available at [https://jbloomlab.github.io/SARS2-RBD-escape-calc/](https://jbloomlab.github.io/SARS2-RBD-escape-calc/) and that page also contains a detailed explanation of how the calculator works as well as the papers from Yunlong Cao's group from which the deep mutational scanning data are derived.
So go to that website to understand and use the calculator.

There is also a Python module for command-line implementation of the calculator, as described [here](https://jbloomlab.github.io/SARS2-RBD-escape-calc/escapecalculator.html).

## Contents of this repo

- [./Cao_data/](Cao_data): the input data used by calculator, in the form of submodules of the GitHub repos from their papers.

- [config.yaml](config.yaml): configuration for processing the data and creating the calculator.

- [environment.yml](environment.yml): `conda` environment for processing the data and building the calculator.

- [process_data.ipynb](process_data.ipynb): Jupyter notebook that processes the deep mutational scanning data to prep it for the calculator. Those processed data are placed in [./results/](results).

- [plot_calculator.ipynb](plot_calculator.ipynb): Jupyter notebook that uses the processed data to build the actual interactive calculator using Altair.

- [format_altair_html.py](format_altair_html.py): Python script that uses the configuration in [./data_for_plot_formatting/](data_for_plot_formatting) to created a formatted version of the escape calculator (adding legend and Google analytics tag) that is placed at [./docs/index.html](docs/index.html).

- [escapecalculator.py](escapecalculator.py): the Python module for command-line non-interactive escape calculations.

- [run_pipeline.bash](run_pipeline.bash): Bash script that runs all steps.

## Building the calculator

First, build and activate the `conda` environment in [environment.yml](environment.yml).

Then process the data by running the Jupyter notebook [process_data.ipynb](process_data.ipynb).

Next, build the escape calculator by running the Jupyter notebook [plot_calculator.ipynb](plot_calculator.ipynb).

Build the documentation for the Python module by running:

    pdoc escapecalculator.py -o docs/ --no-search

This must be done before the next step to avoid overwriting the index in [./docs/](docs).

Then format the escape calculator by running:

    python format_altair_html.py \
        --chart results/escape_chart.html \
        --markdown data_for_plot_formatting/legend.md \
        --site https://jbloomlab.github.io/SARS2-RBD-escape-calc \
        --title "SARS-CoV-2 RBD antibody escape calculator" \
        --description "Calculate antibody escape by mutations to SARS-CoV-2 RBD" \
        --google_analytics_tag data_for_plot_formatting/google_analytics_tag.html \
        --output docs/index.html

You can do all of these steps automatically by just running the Bash script [run_pipeline.bash](run_pipeline.bash).

## Old version of calculator
This repo replaces that for the the original version of the calculator at [https://github.com/jbloomlab//SARS2_RBD_Ab_escape_maps](https://github.com/jbloomlab//SARS2_RBD_Ab_escape_maps).
