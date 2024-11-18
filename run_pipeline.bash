#!/bin/bash

set -euo pipefail

# run `process_data.ipynb`
jupyter nbconvert \
    --to notebook \
    --inplace \
    --execute process_data.ipynb \
    --ExecutePreprocessor.timeout=-1

# run `plot_calculator.ipynb`
jupyter nbconvert \
    --to notebook \
    --inplace \
    --execute plot_calculator.ipynb \
    --ExecutePreprocessor.timeout=-1

# build docs for command-line calculator
mkdir -p docs
pdoc escapecalculator.py -o docs/ --no-search

# run `format_altair_html.py`
python format_altair_html.py \
    --chart results/escape_chart.html \
    --markdown data_for_plot_formatting/legend.md \
    --site https://jbloomlab.github.io/SARS2-RBD-escape-calc \
    --title "SARS-CoV-2 RBD antibody escape calculator" \
    --description "Calculate antibody escape by mutations to SARS-CoV-2 RBD" \
    --google_analytics_tag data_for_plot_formatting/google_analytics_tag.html \
    --output docs/index.html
