import pandas as pd
from math import log2
import re


def count_long_ds_regions(structure: str) -> int:
    """
    Count the number of long double-stranded regions in a structure
    """
    matches = re.findall(r"\({33}", structure)

    return len(matches)


def is_innate_immunity_safe(structure: str) -> bool:
    return "(" * 33 not in structure


def get_idx_inntate_immunity_safe_designs(structures) -> list:
    """
    Return the index of the designs that are safe from innate immunity
    """
    idx_list = []
    for i, structure in enumerate(structures):
        if not is_innate_immunity_safe(structure):
            idx_list.append(i)
    return idx_list


def calc_cai(transcript):
    codon_tab = pd.read_csv("./CAI_table_human.csv", index_col=0)
    try:
        cai = 0.0
        transcript = transcript.upper()
        for i in range(0, len(transcript), 3):
            codon = transcript[i : i + 3]  # POSIBLE ERROR: case sensitive
            w_i = codon_tab.loc[codon, "X"] / codon_tab.loc[codon, "c_max"]
            cai += log2(w_i)

        answer = 2 ** (cai / (len(transcript) / 3 - 1))
    except KeyError:
        print("ERROR: Invalid codon in transcript")
        raise
    return answer


def measure_pairing_proportion(ss: str) -> float:
    LEADING_5_PRIME_LEN = 15
    proportion = ss[:LEADING_5_PRIME_LEN].count("(") / LEADING_5_PRIME_LEN

    return proportion
