## SARS-CoV-2 RBD antibody escape calculator
This page enables you to interactively calculate the escape from human polyclonal serum antibodies caused by mutating sites in the SARS-CoV-2 receptor binding domain (RBD).
The calculations are based on a large set of experimental measurements of antibody escape made by [yeast-display deep mutational scanning](https://doi.org/10.1016/j.chom.2020.11.007).
This calculator is described in [this paper](https://doi.org/10.1093/ve/veac021) (although the current version is substantially improved and uses a much larger dataset than described in the paper).

#### How to use the calculator
The line plot at top shows each site in the RBD.
The line height at each site indicates how much mutating that site escapes human antibodies.
Click on a site to "mutate" it, and the blue line plot will then show how much mutating each remaining site escapes human antibodies.
Click on additional sites to make more mutations
You can mouse over points for details on sites.

The blue bar below the line plot at left shows the total fraction of the antibody neutralization that is escaped by the mutations that have been made on the line plot.
Double click on this blue bar to clear selections of mutated sites.

The bar chart below the line plot at right shows the different sources (human cohorts) that are the source of the antibodies used by the calculator, and also indicates how many of those antibodies neutralize the selected virus.
To change which antibody sources are used, click on bars to select them (in which case they are shown opaque) or de-delect them (in which case they are shown transparent).

Below the plots are a number of options:

 - Choose the virus against which neutralization is measured. This is an important option--results can look quite different depending on this choice as different antibodies neutralize different viruses.
 - Choose which study or studies to use antibodies from. Generally you will want to use antibodies from any study unless you have a good reason otherwise.
 - Use only antibodies known to bind a specific RBD. For most analyses, you can keep this option as any, as only antibodies that neutralize the virus will be shown anyway assuming escape is weighted by IC50.
 - Adjust the mutation escape strength used by the calculator. Read [the paper](https://doi.org/10.1093/ve/veac021) for an explanation of this parameter, probably you do not need to adjust it.
 - Weight antibodies by the negative log IC50. Generally you want to use this option assuming you are interested in escape from antibody neutralization. Selecting no will instead show total binding rather than neutralization.
 - Re-weight antibodies from non-representative sources. This applies the re-weighting of antibodies from sources that is recommended by [Yisimayi et al (2023)](https://doi.org/10.1101/2023.05.01.538516), and you probably want to keep this option as *yes*.

#### Citations and acknowledgements
The calculator is maintained by [Jesse Bloom](https://twitter.com/jbloom_lab) and the citation for the calculator is:

 - [Greaney, Starr, & Bloom, Virus Evolution, 8:veac021 (2022)](https://doi.org/10.1093/ve/veac021)

The deep mutational scanning data used by the calculator comes from the following studies by the group of [Yunlong Cao](https://twitter.com/yunlong_cao/):

 - [Cao et al, Nature, 614:521-529 (2023)](https://www.nature.com/articles/s41586-022-05644-7)
 - [Yisimayi et al, bioRxiv, DOI 10.1101/2023.05.01.538516 (2023)](https://doi.org/10.1101/2023.05.01.538516)

Special thanks to [Fanchong Jian](https://twitter.com/jianfcpku) (an author on the above papers from Yunlong Cao's group), who created the GitHub repos from which the deep mutational scanning data were extracted.

The interactive calculator is created using [Altair](https://altair-viz.github.io/):

 - [VanderPlas et al, JOSS, 3:1057](https://joss.theoj.org/papers/10.21105/joss.01057)

#### The calculator shows escape, **not** mutational tolerance
This calculator shows the antigenic effects of mutations to the RBD.
It does **not** incorporate information about the functional effects of those mutations, such as how they affect ACE2 affinity.
In practice, we know that the RBD mutations that actually fix during SARS-CoV-2 evolution are ones that cause antibody escape while also not excessively impairing ACE2 affinity.
Therefore, the calculator is designed to answer this question: *What is the antigenic effect of mutating a specific site?*

If you instead want to try to *predict* which sites are likely to change during natural SARS-CoV-2 evolution, you will want to integrate the calculator's estimates of the antigenic effects of mutations with other data on how those mutations affect other RBD properties that contribute to viral fitness.
Currently the best data on those other properties are available at the following sites:

 - [https://jbloomlab.github.io/SARS-CoV-2-RBD_DMS_Omicron/RBD-heatmaps/](https://jbloomlab.github.io/SARS-CoV-2-RBD_DMS_Omicron/RBD-heatmaps/), which shows mutational effects on RBD ACE2 affinity and expression as measured by yeast-display deep mutational scanning in [Starr et al, PLoS Pathogens, 2022, 18:e1010951](https://doi.org/10.1371/journal.ppat.1010951)
 - [https://jbloomlab.github.io/SARS2-mut-fitness/S.html](https://jbloomlab.github.io/SARS2-mut-fitness/S.html), which shows the estimated fitness effects of mutations to the spike as estimated in [Bloom & Neher, bioRxiv, 2023, DOI 10.1101/2023.01.30.526314](https://doi.org/10.1101/2023.01.30.526314)

In addition, the amino-acid mutations that occur in natural viral evolution are predominantly ones accessible by single nucleotide mutations, especially via the [more common types of nucleotide changes](https://doi.org/10.1093/molbev/msad085).

Note that the papers by Yunlong Cao's group already integrate these different types of data to make forecasts of likely future mutations, as for instance in [Yisimayi et al, bioRxiv, DOI 10.1101/2023.05.01.538516 (2023)](https://doi.org/10.1101/2023.05.01.538516).

### Detailed explanation of how the calculator works
For each antibody $a$, let $x_{a,r}$ be the measurement of how much mutating site $r$ escapes the antibody.
For all antibodies, these site escape measurements are normalized so that the maximum escape at any site is 1; eg, $\max_r x_{a,r} = 1$.
In the line plot, the gray lines (or the blue lines when there are no mutations) show the mean of $x_{a,r}$ overall antibodies; that is, they show $\frac{\sum_a x_{a,r}}{A}$ where $A$ is the total number of antibodies.

Let $\mathcal{M}$ be the set of mutated sites.
For each antibody $a$, we compute the binding retained as
$$
b_a\left(\mathcal{M}\right) = w_a \left(\prod\limits_{r \in \mathcal{M}} \left[1 - x_{a,r}\right]\right)^s.
$$
This equation means that if the RBD is mutated at a strong site of escape for an antibody $a$, much of the contribution of that antibody is lost (if mutated at strongest site of escape, all contribution is lost).
The $s$ variable represents how dramatically binding is lost by mutations at sites of escape that are not the strongest one: larger values mean mutations even at moderate sites of escape reduce binding a lot.
The value of $s$ is set via the mutation escape strength slider below the plot.

The $w_a$ values are the "weights" of the antibodies.
By default, antibodies are weighted proportional to their negative log IC50s divided by 10 (10 is the weakest reported IC50 in the data), so that more potent antibodies have greater weights.
Specifically, the weight is calculated as
$$w_a = \max\left(0, -\log \left[IC50 / 10\right] \right)$$
where $IC50$ is the IC50 in $\mu$g/ml.
If the re-weighting option is specified, the antibodies are also re-weighted when recommended by Cao et al (eg, as in [Yisimayi et al (2023)](https://doi.org/10.1101/2023.05.01.538516)) such that the mean reweight of all antibodies from a source is one, but the balance of cross reactive and specific antibodies is what is recommended in the Cao et al papers.

The blue lines show the escape at each site **after** making the mutations $\mathcal{M}$.
For each site $r$, this is $\frac{\sum_a x_{a,r} \times b_a\left(\mathcal{M}\right)}{A}$.

The blue bar at lower left in the chart shows the total weighted fraction of all antibodies that still bind after the mutations, $\frac{\sum_a b_a\left(\mathcal{M}\right)}{\sum_a w_a}$.

#### Code, data, and command-line version of calculator
The computer code and data used to create the calculator are at [https://github.com/jbloomlab/SARS2-RBD-escape-calc/](https://github.com/jbloomlab/SARS2-RBD-escape-calc/)

That repository also includes a Python module that implements a version of the calculator that can be run in batch over lists of RBD mutations.
The module implementing that calculator is [here](https://github.com/jbloomlab/SARS2-RBD-escape-calc/blob/main/escapecalculator.py), and documentation for using it is at [https://jbloomlab.github.io/SARS2-RBD-escape-calc/escapecalculator.html](https://jbloomlab.github.io/SARS2-RBD-escape-calc/escapecalculator.html).

#### Original version of calculator
The original version of the antibody-escape calculator was implemented in the code available at [https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps](https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps)

This new calculator supersedes that earlier version as it represents a complete re-write that streamlines the code and updates the design to incorporate the newest deep mutational scanning data.
