from concurrent.futures import ProcessPoolExecutor
import subprocess as sp
from multiprocessing import cpu_count
import pathlib
import re
import os
import pandas as pd
import pprint
from text_processing import input_preprocessing
from ViennaRNA import RNA

DATAPATH = "./data/proteins/testseq"
DESIGNPATH = "./designs/proteins/"
LAMBDA = [0]


def subprocess_lineardesign(split_protein_sequence: str, cmd: list):
    with open(split_protein_sequence, "r") as f:
        # Create a subprocess that runs the pipeline
        data = f.readlines()
        # leading_sequence = data[1]
        result = sp.run(cmd, input="".join(data), capture_output=True, text=True)

        # Run the pipeline and capture the output
        output, error = result.stdout, result.stderr
        if error is not None:
            print(error)

        # Modify 5’-end leader region and choose the best design

        return output


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


def parse_result(results: list, lambda_) -> pd.DataFrame:
    # Turn this result into a dataframe
    """
    >seq2
    mRNA sequence:  AUGCUGGAUCAGGUCAACAAGCUGAAGUACCCUGAGGUUUCGUUGACCUGA
    mRNA structure: ........(((((((((((((((..((....))..))))).))))))))))
    mRNA folding free energy: -20.70 kcal/mol; mRNA CAI: 0.768
    """

    # Split the result into individual records
    parsed = []
    for result in results:
        records = result.split(">")[1:][0]

        # Split each record into its components
        records = records.split("\n")

        # Data cleaning
        parsed.append(
            {
                "Name": records[0],
                "mRNA sequence": records[1].split(":")[1].strip(),
                "mRNA structure": records[2].split(":")[1].strip(),
                "MFE (kcal/mol)": float(
                    records[3].split(";")[0].split(":")[1].split()[0].strip()
                ),
                "CAI": float(records[3].split(";")[1].split(":")[1].strip()),
                "lambda": lambda_,
            }
        )

    pprint.pprint(parsed)
    return pd.DataFrame.from_dict(parsed)


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
    with ProcessPoolExecutor(max_workers=1) as executor:  # DEBUG: fixed to 1
        # output redirection from stdout to file

        for lambda_ in LAMBDA:
            print(f"Running lambda = {lambda_}...")
            lambda_group = []
            for item in items:
                # Define the command in the pipeline
                cmd = ["./lineardesign", "-l", str(lambda_)]

                future = executor.submit(subprocess_lineardesign, f"{str(item)}", cmd)
                lambda_group.append(future)

            pattern = r"(j=\d*\n)+"  # Modified regex
            lambda_group = [future.result() for future in lambda_group]
            lambda_group = list(
                map(lambda x: re.sub(pattern, "", x), lambda_group)
            )  # Remove iteration mark
            lambda_group = list(
                map(lambda x: x[:-1], lambda_group)
            )  # Remove tailing \n

            # # Print the output
            # print(output)
            pathlib.Path(DESIGNPATH).mkdir(parents=True, exist_ok=True)
            df = parse_result(lambda_group, lambda_)
            designs.append(df)

        pd.concat(designs).to_csv(f"{DESIGNPATH}/{path.name}+design.csv", index=False)
