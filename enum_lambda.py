from concurrent.futures import ProcessPoolExecutor
import subprocess as sp
from multiprocessing import cpu_count
import pathlib
import re
import os
import pandas as pd
import numpy as np
from text_processing import *
import gc

try:
    import RNA
except ImportError:
    from ViennaRNA import RNA
from itertools import product
from RNA_toolkit import *

DATAPATH = "./data/proteins/nuclease.fasta"
DESIGNPATH = "./designs/proteins/"
LAMBDA = [0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 1000]

IDX_INNATE_IMMUNITY_SAFE = 0
IDX_SECONDARY_STRUCTURE = 1
IDX_MFE = 2
IDX_PAIRING_PROPORTION = 3
IDX_LONG_DS_REGIONS = 4
IDX_RNA_SEQUENCE = 5
IDX_CAI = 6


def rna_measurement_wrapper(sequence: str) -> pd.Series:
    ss, mfe = RNA.fold(sequence)
    robj = pd.Series(
        [
            is_innate_immunity_safe(ss),
            ss,
            mfe,
            measure_pairing_proportion(ss),
            count_long_ds_regions(ss),
            sequence,
            calc_cai(sequence),
        ]
    )

    return robj


def pipeline_lineardesign(
    split_protein_sequence: str, cmd: list, codon_table: pd.DataFrame
):
    with open(split_protein_sequence, "r") as f:
        # Create a subprocess that runs the pipeline
        # 1. Designing the CDS region except for the 5’-end leader region
        data = f.readlines()

        # Split the data into two parts: leader and follower
        name = data[0]  # Name of the sequence
        leading_sequence = data[1][:5]  # 5' leader sequence
        following_sequence = data[1][5:]  # follower sequence
        completed_process = sp.run(
            cmd, input=f"{name}\n{following_sequence}", capture_output=True, text=True
        )

        # Output redirection from stdout to object
        output, error = completed_process.stdout, completed_process.stderr
        if error != "" or output == "":
            print(error)
            raise Exception("Error in running lineardesign")

        pattern = r"(j=\d*\n)+"  # Modified regex
        output = re.sub(pattern, "", output)[
            :-1
        ]  # Remove iteration mark and tailing \n
        structured_output = make_structured_result_from_lineardesign(output)

        # 2. Enumerate all possible subsequences in the 5’-end leader region
        pooled_codon = []
        for aa in leading_sequence:
            pooled_codon.append(codon_table[codon_table["AA"] == aa].index.tolist())
        leader_candidates = ["".join(tokens) for tokens in product(*pooled_codon)]

        # 3. Concatenate and choose the best design
        test_sequences = [
            f"{x}{structured_output['mRNA sequence']}" for x in leader_candidates
        ]

        # Parallel ss computation
        best_result = {}
        with ProcessPoolExecutor(
            max_workers=cpu_count()  # DEBUG: fixed to 1
        ) as executor:  # BUG: may consume all the memory,
            # 4. Filter the design with only < 33-bp pairing in a double-stranded region: Design constraint 2
            futures = [
                executor.submit(rna_measurement_wrapper, x) for x in test_sequences
            ]

            processing_results = [f.result() for f in futures]
            df = pd.concat(processing_results, axis=1).T.astype(
                {
                    IDX_INNATE_IMMUNITY_SAFE: bool,
                    IDX_SECONDARY_STRUCTURE: str,
                    IDX_MFE: float,
                    IDX_PAIRING_PROPORTION: float,
                    IDX_LONG_DS_REGIONS: int,
                    IDX_RNA_SEQUENCE: str,
                    IDX_CAI: float,
                }
            )
            # BUG
            # df.to_csv("./test.csv")  # DEBUG
            safe_designs = df[lambda x: x[IDX_INNATE_IMMUNITY_SAFE] == True]
            if safe_designs.empty:
                safe_designs = df[
                    df[IDX_LONG_DS_REGIONS].idxmin()
                ]  # Choose the alternative design

            best_idx = safe_designs[IDX_PAIRING_PROPORTION].idxmin()
            best_design = safe_designs.loc[best_idx].to_numpy()
            # 5. Return the best design
            best_result = parse_best_design(name, best_design)

            del df, safe_designs, best_design, best_idx, processing_results, futures
            gc.collect()

        del (
            test_sequences,
            leader_candidates,
            pooled_codon,
            following_sequence,
            leading_sequence,
            name,
            data,
        )
        gc.collect()

    return best_result


def split_directory_cleanup(path: pathlib.Path) -> None:
    # get a list of all files in the directory
    file_list = os.listdir(str(path))

    # iterate over each file in the directory and remove it
    for filename in file_list:
        file_path = os.path.join(str(path), filename)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
                # print(f"Deleted {file_path}")
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")


if __name__ == "__main__":
    path = pathlib.Path(DATAPATH)
    split_path = path.parent / "split"
    try:
        split_path.mkdir(parents=True)
    except FileExistsError:
        print(f"Directory {split_path} already exists.")
        print(f"Cleaning up {split_path}...")
        split_directory_cleanup(split_path)
    processed_fasta_path = input_preprocessing(DATAPATH)
    sp.run(["split", "-l", "2", processed_fasta_path, f"./{str(split_path)}/"])
    items = list(pathlib.Path(split_path).glob("*"))

    designs = []

    codon_tab = pd.read_csv("./CAI_table_human.csv", index_col=0)

    for lambda_ in LAMBDA:
        print(f"Running lambda = {lambda_}...")
        lambda_group = []
        for item in items:
            # Define the command in the pipeline
            cmd = ["./lineardesign", "-l", str(lambda_)]
            lambda_group.append(pipeline_lineardesign(f"{str(item)}", cmd, codon_tab))

        df = tagging_lambda_and_create_dataframe(lambda_group, lambda_)
        designs.append(df)

    # Create a directory to store the designs
    pathlib.Path(DESIGNPATH).mkdir(parents=True, exist_ok=True)
    pd.concat(designs).to_csv(f"{DESIGNPATH}/{path.name}+design.csv", index=False)
