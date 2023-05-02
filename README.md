# Escape calculator for the SARS-CoV-2 RBD

## Data from Yunlong Cao et al
The data from Yunlong Cao's group is in the following submodules:

  - Their 2022 data from [Imprinted SARS-CoV-2 humoral immunity induces convergent Omicron RBD evolution](https://www.nature.com/articles/s41586-022-05644-7) from [this repo](https://github.com/jianfcpku/convergent_RBD_evolution) are in the submodule [./Cao_data/convergent_RBD_evolution](Cao_data/convergent_RBD_evolution)

  - Their 2023 data on antibodies from re-infections in [this repo](https://github.com/jianfcpku/SARS-CoV-2-reinfection-DMS) are in the submodule [./Cao_data/SARS-CoV-2-reinfection-DMS](Cao_data/SARS-CoV-2-reinfection-DMS)


Yunlong described the second dataset as follows via e-mail:

    I am excited to share a new dataset containing the DMS profiles of 1,350 SARS-CoV-2 RBD-targeting mAbs using a DMS library based on the BA.5 RBD. These mAbs primarily originate from BA.5 and BF.7 breakthrough infection convalescents, as well as BA.5/BF.7 reinfection convalescents who previously experienced BA.1 or BA.2 breakthrough infections. Note a high proportion of these mAbs are Omicron-specific.

    Following an approach similar to that in our prior work, we performed mutation calculations based on these profiles and their neutralizing activities against XBB.1.5. I anticipate that our manuscript will be available on bioRxiv within the next one or two days.

    We have uploaded most of the data and code to our GitHub repository (https://github.com/jianfcpku/SARS-CoV-2-reinfection-DMS). While the detailed README file is still in progress, let me briefly introduce the main results here. The repository contains three .csv files, which include the integrated DMS profiles:

    > "antibody_dms_merge.csv.gz" presents the raw scores from the epistasis model, normally filtered by RBD expression.

    > "antibody_dms_merge_clean.csv.gz" contains the most significant mutations derived from "antibody_dms_merge.csv.gz", filtered according to median scores, as in our previous results.

    > "calculation/antibody_dms_merge_no_filter_clean.csv.gz" features raw scores without any expression filter, as we will apply a weight on expression when calculating mutation preferences. To minimize noise, a median-based filter is also used ("_clean"). I believe that employing this file retains mutations as many as possible.

    "antibody_info.csv" contains the neutralizing activities (IC50 in ug/mL) and the cross reactivity determined by ELISA. Notably, the dataset specifically enriches for Omicron-specific antibodies, potentially introducing bias when estimating mutation preferences. I highly recommend a weighting strategy that assigns higher weights to cross-reactive mAbs, resulting in 89% cross-reactive mAbs for BA.5/BF.7 BTI cohorts and 51% for reinfection cohorts. These ratios are determined by unbiased characterization of mAbs using ELISA. This approach is detailed in the code found in "calculation/calculate_preference.ipynb".

    I hope you maintain your interest in this subject as well as sharing these data though your calculator, and that these data will prove helpful to you and others who are also intrigued.
