## SARS-CoV-2 RBD antibody escape calculator


### Citations and acknowledgements
The calculator itself is maintained by [Jesse Bloom](https://twitter.com/jbloom_lab), and the citation for the calculator is:

 - [Greaney, Starr, & Bloom, Virus Evolution, 8:veac021 (2022)](https://doi.org/10.1093/ve/veac021)

The deep mutational scanning data used by the calculator comes from the following studies by the group of [Yunlong Cao](https://twitter.com/yunlong_cao/):

 - [Cao et al, Nature, 614:521-529 (2023)](https://www.nature.com/articles/s41586-022-05644-7)
 - [Cao et al, bioRxiv, DOI 10.1101/2023.05.01.538516 (2023)](https://doi.org/10.1101/2023.05.01.538516)

We especially acknowledge [Fanchong Jian](https://twitter.com/jianfcpku) (an author on the above papers from Yunlong Cao's group), who created the GitHub repos from which the deep mutational scanning data were extracted.

If you use this calculator, please cite the papers above.

The interactive calculator was created using [Altair](https://altair-viz.github.io/), so please consider also citing that software package:

 - [VanderPlas et al, JOSS, 3:1057](https://joss.theoj.org/papers/10.21105/joss.01057)

### Detailed explanation of how the calculator works


### Code and command-line version of calculator
All the computer code used to create this calculator is at [https://github.com/jbloomlab/SARS2-RBD-escape-calc/](https://github.com/jbloomlab/SARS2-RBD-escape-calc/)

That repository also includes a Python module that implements a version of the calculator that can be run in batch over lists of RBD mutations.
The module implementing that calculator is [here](https://github.com/jbloomlab/SARS2-RBD-escape-calc/blob/main/escapecalculator.py), and documentation for using it is at **???**

### Original version of calculator
The original version of the antibody-escape calculator was implemented in the code available at [https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps](https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps)

This new calculator supersedes that earlier version as it represents a complete re-write that streamlines the code and updates the design to incorporate the newest deep mutational scanning data.
