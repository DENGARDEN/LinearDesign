import pandas as pd
import pprint
from typing import Tuple

# Description: This file contains the functions for text processing


def input_preprocessing(path: str, mode=None) -> str:
    """
    Preprocess the input file to remove the unnecessary line breaks
    """
    with open(path, "r") as f:
        lines = f.read().split(">")[1:]  # Remove the first empty string
        buffer = []
        for line in lines:
            temp = line.split("\n")
            if mode is None:
                buffer.append(f'>{temp[0]}\n{"".join(temp[1:])}*\n')
            elif mode == "RNA":
                buffer.append(f'>{temp[0]}\n{"".join(temp[1:])}\n')

        preprocessed_file_path = f"{path}_preprocessed.fasta"
        with open(preprocessed_file_path, "w") as d:
            d.writelines(buffer)
        return preprocessed_file_path


def make_structured_result_from_lineardesign(result: str) -> dict:
    """
    The function takes a string as input, extracts specific information from it, and returns a
    dictionary with the extracted information.

    :param result: The input parameter "result" is a string that contains information about a biological
    sequence, including its name, mRNA sequence, mRNA structure, minimum free energy (MFE), and codon
    adaptation index (CAI). The string is formatted in a specific way, with each piece of information
    separated by specific
    :type result: str
    :return: The function `make_structured_result` takes a string `result` as input and returns a
    dictionary with structured information extracted from the input string. The dictionary contains the
    following keys: "Name", "mRNA sequence", "mRNA structure", "MFE (kcal/mol)", and "CAI". The values
    for these keys are extracted from the input string using string manipulation techniques.
    """
    try:
        records = result.split(">")[1]

        # Split each record into its components
        records = records.split("\n")
        return {
            "Name": records[0],
            "mRNA sequence": records[1].split(":")[1].strip(),
            "mRNA structure": records[2].split(":")[1].strip(),
            "MFE (kcal/mol)": float(
                records[3].split(";")[0].split(":")[1].split()[0].strip()
            ),
            "CAI": float(records[3].split(";")[1].split(":")[1].strip()),
        }
    except Exception as e:
        print(e)
        print(result)
        raise Exception("Error in parsing the result from lineardesign")


def parse_best_design(records: Tuple[str, str, str, float, float]) -> dict:
    return {
        "Name": records[0].strip(),
        "mRNA sequence": records[1].strip(),
        "mRNA structure": records[2].strip(),
        "MFE (kcal/mol)": records[3],
        "CAI": records[3],
    }


def tagging_lambda_and_create_dataframe(results: list, lambda_: int) -> pd.DataFrame:
    # Tag lambda and turn the result into a dataframe
    """
    {
        "Name": ">seq1",
        "mRNA sequence": "AUGCCUAAUACUCUUGCGUGCCCGUAG",
        "mRNA structure": "..............((((....)))).",
        "MFE (kcal/mol)": -1.7999999523162842,
        "CAI": -1.7999999523162842,
    }
    """

    # Split the result into individual records
    tagged = []
    for result in results:
        result["lambda"] = lambda_
        tagged.append(result)

    pprint.pprint(tagged)
    return pd.DataFrame.from_dict(tagged)
