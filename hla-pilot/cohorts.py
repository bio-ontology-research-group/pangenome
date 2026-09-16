"""Shared cohort order, labels and colours for the HLA figures.
HPRC r2 is split by population (data/hprc_r2_populations.tsv): Japanese (JPT), Jewish (HG002, Ashkenazi,
a single individual), other East Asian, rest. K-PanRef = 16 Korean samples from the K-PanRef graph; CPC = Chinese Pangenome Consortium Phase 1 samples from the CPC graph."""
ORDER = ["APR", "JaSaPaGe-Saudi", "HPRC-Jewish", "JaSaPaGe-Japanese", "HPRC-Japanese",
         "KPanRef-Korean", "CPC-Chinese", "HPRC-EastAsian", "HPRC-Rest"]
SHORT = {"APR": "APR (UAE Arab)", "JaSaPaGe-Saudi": "JaSaPaGe Saudi", "HPRC-Jewish": "HPRC Ashkenazi Jewish",
         "JaSaPaGe-Japanese": "JaSaPaGe Japanese", "HPRC-Japanese": "HPRC Japanese (JPT)",
         "KPanRef-Korean": "K-PanRef Korean", "CPC-Chinese": "CPC Chinese", "HPRC-EastAsian": "HPRC other East Asian", "HPRC-Rest": "HPRC rest"}
ONE = {"APR": "APR", "JaSaPaGe-Saudi": "Saudi", "HPRC-Jewish": "Jewish", "JaSaPaGe-Japanese": "Jpn (JaSaPaGe)",
       "HPRC-Japanese": "Jpn (HPRC)", "KPanRef-Korean": "Korean", "CPC-Chinese": "Chinese (CPC)", "HPRC-EastAsian": "E.Asian", "HPRC-Rest": "HPRC rest"}
COLOR = {"APR": "#1f77b4", "JaSaPaGe-Saudi": "#d62728", "HPRC-Jewish": "#9467bd", "JaSaPaGe-Japanese": "#2ca02c",
         "HPRC-Japanese": "#98df8a", "KPanRef-Korean": "#ff7f0e", "CPC-Chinese": "#e377c2", "HPRC-EastAsian": "#bcbd22", "HPRC-Rest": "#7f7f7f"}


def present(df_cohorts):
    """cohorts of ORDER that occur in the data, in ORDER"""
    s = set(df_cohorts)
    return [c for c in ORDER if c in s]
