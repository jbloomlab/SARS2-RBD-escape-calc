"""Calculate escape after some mutations to the SARS-CoV-2 RBD.

This Python module, is available at
https://github.com/jbloomlab/SARS2-RBD-escape-calc/blob/main/escapecalculator.py

It provides a programmatic method to perform the escape calculations implemented in
the interactive calculator implemented at
https://jbloomlab.github.io/SARS2-RBD-escape-calc

By default, it downloads the same data used by the interactive calculator and uses the
same initial parameter settings as the interactive calculator. You can override these
defaults by using different arguments when initializing the :class:`EscapeCalculator`.

See https://github.com/jbloomlab/SARS2-RBD-escape-calc for all code and data. The
calculator is written by Jesse Bloom, and the citation is
`this paper <https://doi.org/10.1093/ve/veac021>`_.

To use the calculator, open a Python session, import this module, and initialize an
:class:`EscapeCalculator`:

>>> calc = EscapeCalculator()

With no mutations, there is no escape (all neutralization retained):

>>> calc.binding_retained([])
1.0

But if you mutate some key antigenic sites, there will be a dramatic reduction in
neutralization retained:

>>> calc.binding_retained([440, 505]).round(3)
0.545

If you have a whole set of sequences and have tabulated which sites are mutated,
you can apply the calculator in bulk to the data frame.
First, create such a data frame of sequences:

>>> import pandas as pd
>>> seqs = pd.DataFrame.from_records(
...     [("seq1", []), ("seq2", [498]), ("seq3", [440]), ("seq4", [440, 505])],
...     columns=["name", "mutated sites"],
... )
>>> seqs
   name mutated sites
0  seq1            []
1  seq2         [498]
2  seq3         [440]
3  seq4    [440, 505]

Now apply the calculator:

>>> seqs["neutralization retained"] = seqs["mutated sites"].map(calc.binding_retained)
>>> seqs.round(3)
   name mutated sites  neutralization retained
0  seq1            []                    1.000
1  seq2         [498]                    0.769
2  seq3         [440]                    0.790
3  seq4    [440, 505]                    0.545

You can also create new calculators that compute escape relative to different viruses,
for instance BA.2:

>>> calc_ba2 = EscapeCalculator(virus="BA.2")

Now the escape will be different because different antibodies neutralize that virus:

>>> calc_ba2.binding_retained([440, 505]).round(3)
0.786

"""

__docformat__ = 'numpy'


import requests

import numpy

import pandas as pd

import yaml


class EscapeCalculator:
    """Calculates residual polyclonal antibody binding after some mutations.

    The calculator is implemented interactively at
    `https://jbloomlab.github.io/SARS2-RBD-escape-calc <https://jbloomlab.github.io/SARS2-RBD-escape-calc/>`_

    By default, this command-line calculator will exactly implement the calculations of the
    interactive calculator with default settings. You can pass parameters to override
    the default settings.

    Parameters
    ----------
    escape : str
        Path or URL of CSV containing the escape data.
    antibody_ic50s : str
        Path or URL of CSV containing the antibody IC50s.
    antibody_binding : str
        Path or URL of CSV containing the antibody binding.
    antibody_sources : str
        Path or URL of CSV containing the antibody sources.
    antibody_reweighting : str
        Path or URL of CSV containing the antibody reweightings.
    config : str
        Path or URL of YAML file containing initial settings.
    mut_escape_strength : None or float
        If not `None`, override default `init_mutation_escape_strength` in `config`.
    weight_by_neg_log_ic50 : None or bool
        If not `None`, override default `init_weight_by_neg_log_IC50` in `config`.
    study : None or str
        If not `None`, override default `init_study` in `config`.
    binds : None or str
        If not `None`, override default `init_binds` in `config`.
    virus : None or str
        If not `None`, override default `init_virus` in `config`.
    sources : None or dict
        If not `None`, override default `init_sources` in `config`.
    reweight : None or bool
        If not `None`, override default `init_reweight` in `config`.

    Example
    -------

    Initialize the calculator:

    >>> calc = EscapeCalculator()

    Calculate escape at sites after some mutations:

    >>> sites_of_interest = [403, 440, 484, 498, 505, 510]
    >>> calc.escape_per_site([440, 505]).query("site in @sites_of_interest").round(3)
         site  original_escape  retained_escape
    69    403            0.105            0.030
    101   440            0.162            0.018
    143   484            0.043            0.032
    156   498            0.169            0.076
    163   505            0.208            0.022
    167   510            0.001            0.001

    Calculate overall neutralization retained after no mutations or some mutations:

    >>> calc.binding_retained([])
    1.0
    >>> calc.binding_retained([440, 505]).round(3)
    0.545

    Now repeat tests with some non-default options:

    >>> calc2 = EscapeCalculator(
    ...     mut_escape_strength=1,
    ...     weight_by_neg_log_ic50=False,
    ...     study="Cao et al, 2022, Nature",
    ...     binds="Wuhan-Hu-1",
    ...     virus="D614G",
    ...     sources={"include_exclude": "include", "sources": ["WT convalescents"]},
    ... )

    >>> calc2.escape_per_site([484]).query("site in @sites_of_interest").round(3)
         site  original_escape  retained_escape
    62    403            0.006            0.006
    84    440            0.008            0.007
    123   484            0.212            0.029
    134   498            0.009            0.008
    141   505            0.002            0.002
    144   510            0.000            0.000

    >>> calc2.binding_retained([484]).round(3)
    0.788

    """
    def __init__(self,
        escape="https://raw.githubusercontent.com/jbloomlab/SARS2-RBD-escape-calc/main/results/escape.csv",
        antibody_ic50s="https://raw.githubusercontent.com/jbloomlab/SARS2-RBD-escape-calc/main/results/antibody_IC50s.csv",
        antibody_binding="https://raw.githubusercontent.com/jbloomlab/SARS2-RBD-escape-calc/main/results/antibody_binding.csv",
        antibody_sources="https://raw.githubusercontent.com/jbloomlab/SARS2-RBD-escape-calc/main/results/antibody_sources.csv",
        antibody_reweighting="https://raw.githubusercontent.com/jbloomlab/SARS2-RBD-escape-calc/main/results/antibody_reweighting.csv",
        config="https://raw.githubusercontent.com/jbloomlab/SARS2-RBD-escape-calc/main/config.yaml",
        *,
        mut_escape_strength=None,
        weight_by_neg_log_ic50=None,
        study=None,
        binds=None,
        virus=None,
        sources=None,
        reweight=None,
    ):
        """See main class docstring."""
        # read input data 
        self.escape = pd.read_csv(escape)
        assert set(self.escape.columns) == {"antibody", "site", "escape"}
        assert len(self.escape) == len(self.escape.groupby(["antibody", "site", "escape"]))
        assert self.escape.notnull().all().all()
        antibodies = set(self.escape["antibody"])

        self.antibody_ic50s = pd.read_csv(antibody_ic50s)
        assert set(self.antibody_ic50s.columns) == {"antibody", "virus", "IC50"}
        assert (
            len(self.antibody_ic50s)
            == len(self.antibody_ic50s.groupby(["antibody", "virus", "IC50"]))
        )
        assert self.antibody_ic50s["IC50"].max() == 10
        assert antibodies == set(self.antibody_ic50s["antibody"])

        self.antibody_binding = pd.read_csv(antibody_binding)
        assert set(self.antibody_binding.columns) == {"antibody", "binds"}
        assert (
            len(self.antibody_binding)
            == len(self.antibody_binding.groupby(["antibody", "binds"]))
        )
        assert antibodies == set(self.antibody_binding["antibody"])

        self.antibody_sources = pd.read_csv(antibody_sources)
        assert set(self.antibody_sources.columns) == {"antibody", "source", "study"}
        assert (
            len(self.antibody_sources)
            == len(self.antibody_sources.groupby(["antibody", "source", "study"]))
            == len(antibodies)
        )
        assert antibodies == set(self.antibody_sources["antibody"])

        self.antibody_reweighting = pd.read_csv(antibody_reweighting)
        assert set(self.antibody_reweighting.columns) == {"antibody", "reweight"}
        assert antibodies.issuperset(self.antibody_reweighting["antibody"])

        self.data = (
            self.escape
            .merge(self.antibody_ic50s, on="antibody")
            .merge(self.antibody_binding, on="antibody")
            .merge(self.antibody_sources, on="antibody")
            .merge(self.antibody_reweighting, on="antibody", how="left")
            .assign(reweight=lambda x: x["reweight"].fillna(1))
        )
        assert self.data.notnull().all().all()

        # get initial config
        if config.startswith("http"):
            response = requests.get(config, allow_redirects=True)
            config_text = response.content.decode("utf-8")
        else:
            with open(config) as f:
                config_text = f.read()
        config = yaml.safe_load(config_text)
        self.sites = set(range(config["sites"]["start"], config["sites"]["end"] + 1))
        studies = config["studies"]
        studies_rev = {value: key for (key, value) in studies.items()}
        assert set(self.data["study"]) == set(studies)

        if mut_escape_strength is None:
            self.mut_escape_strength = config["init_mutation_escape_strength"]
        else:
            self.mut_escape_strength = mut_escape_strength

        if weight_by_neg_log_ic50 is None:
            self.weight_by_neg_log_ic50 = config["init_weight_by_neg_log_IC50"]
        else:
            self.weight_by_neg_log_ic50 = weight_by_neg_log_ic50
        assert isinstance(self.weight_by_neg_log_ic50, bool), self.weight_by_neg_log_ic50

        if study is None:
            self.study = config["init_study"]
        else:
            if study == "any" or study in studies:
                self.study = study
            elif study in studies_rev:
                self.study = studies_rev[study]
            else:
                raise ValueError(f"invalid {study=}")

        if binds is None:
            self.binds = config["init_binds"]
        else:
            self.binds = binds

        if virus is None:
            self.virus = config["init_virus"]
        else:
            self.virus = virus

        if sources is None:
            sources = config["init_sources"]
        include_exclude = sources["include_exclude"]
        sources_list = sources["sources"]
        self.sources = set(self.data["source"])
        assert self.sources.issuperset(sources_list), f"{self.sources=}\n{sources_list=}"
        if include_exclude == "exclude":
            self.sources = self.sources - set(sources_list)
        elif include_exclude == "include":
            self.sources = set(sources_list)
        else:
            raise ValueError(f"invalid {include_exclude=} in {sources=}")

        if reweight is None:
            self.reweight = config["init_reweight"]
        else:
            self.reweight = reweight
        assert isinstance(self.reweight, bool), self.reweight

        # filter data
        if self.study != "any":
            assert self.study in set(self.data["study"])
            self.data = self.data.query("study == @self.study").drop(columns="study")
        else:
            self.data = self.data.drop(columns="study")

        if self.binds != "any":
            assert self.binds in set(self.data["binds"])
            self.data = self.data.query("binds == @self.binds").drop(columns="binds")
        else:
            self.data = self.data.drop(columns="binds").drop_duplicates()

        assert self.virus in set(self.data["virus"])
        self.data = self.data.query("virus == @self.virus").drop(columns="virus")

        self.data = self.data.query("source in @self.sources").drop(columns="source")

        assert set(self.data.columns) == {"antibody", "site", "escape", "IC50", "reweight"}
        assert len(self.data) == len(self.data.drop_duplicates())
        self.data = (
            self.data
            .assign(neg_log_ic50=lambda x: -numpy.log(x["IC50"] / 10))
            .drop(columns="IC50")
        )

        max_escape_per_antibody = (
            self.data
            .groupby("antibody")
            .aggregate(max_escape=pd.NamedAgg("escape", "max"))
        )
        assert (max_escape_per_antibody["max_escape"] == 1).all()

    def escape_per_site(self, mutated_sites):
        """Escape at each site after mutating indicated sites.

        Parameters
        ----------
        mutated_sites : array-like of integers
            List of mutated sites.

        Returns
        -------
        pandas.DataFrame
            For each site, gives the original escape and the escape
            retained after mutations.

        """
        mutated_sites = set(mutated_sites)
        if not mutated_sites.issubset(self.sites):
            raise ValueError(f"sites {mutated_sites - self.sites} not in {self.sites}")
        df = (
            self.data
            .assign(
                mutated=lambda x: x['site'].isin(mutated_sites).astype(int),
                site_bind_retain=lambda x: 1 - x["escape"] * x["mutated"],
            )
            .groupby(["antibody", "neg_log_ic50", "reweight"], as_index=False)
            .aggregate(antibody_bind_retain=pd.NamedAgg("site_bind_retain", "prod"))
            .assign(
                antibody_bind_retain=lambda x: x["antibody_bind_retain"].pow(
                    self.mut_escape_strength
                ),
                weight=lambda x: (
                    (x["neg_log_ic50"] if self.weight_by_neg_log_ic50 else 1)
                    * (x["reweight"] if self.reweight else 1)
                ),
            )
            [["antibody", "antibody_bind_retain", "weight"]]
            .merge(self.data[["antibody", "site", "escape"]])
            .assign(
                escape=lambda x: x["escape"] * x["weight"],
                retained_escape=lambda x: x["antibody_bind_retain"] * x["escape"],
            )
            .groupby("site")
            .aggregate(
                original_escape=pd.NamedAgg("escape", "sum"),
                retained_escape=pd.NamedAgg("retained_escape", "sum"),
            )
        ) / self.data["antibody"].nunique()
        return df.reset_index()

    def binding_retained(self, mutated_sites):
        """Fraction binding or neutralization retained after mutating indicated sites.

        Parameters
        ----------
        mutated_sites : array-like of integers
            List of mutated sites, must all be in :attr:`BindingCalculator.sites`.

        Returns
        -------
        float
            The fraction binding retained after these mutations.

        """
        mutated_sites = set(mutated_sites)
        if not mutated_sites.issubset(self.sites):
            raise ValueError(f"sites {mutated_sites - self.sites} not in {self.sites}")
        retained = (
            self.data
            .assign(
                mutated=lambda x: x["site"].isin(mutated_sites).astype(int),
                site_bind_retain=lambda x: 1 - x["escape"] * x["mutated"],
            )
            .groupby(["antibody", "neg_log_ic50", "reweight"], as_index=False)
            .aggregate(antibody_bind_retain=pd.NamedAgg("site_bind_retain", "prod"))
            .assign(
                antibody_bind_retain=lambda x: x["antibody_bind_retain"].pow(
                    self.mut_escape_strength
                ),
                weight=lambda x: (
                    (x["neg_log_ic50"] if self.weight_by_neg_log_ic50 else 1)
                    * (x["reweight"] if self.reweight else 1)
                ),
                weighted_antibody_bind_retain=lambda x: (
                    x["antibody_bind_retain"] * x["weight"]
                ),
            )
            [["weight", "weighted_antibody_bind_retain"]]
            .sum(axis=0)
        )
        return retained["weighted_antibody_bind_retain"] / retained["weight"]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
