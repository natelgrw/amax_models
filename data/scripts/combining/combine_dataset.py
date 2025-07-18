"""
raw_combine.py

Accesses cleaned versions of ChemDataExtractor, UVVisML, SMFluo1,
and FluoDB databases, writing a new dataset combining all compound
and solvent combinations along with their absorbance maxima to
a .csv file in the "data/final/dataset" directory.

Author: @natelgrw
Created: 07/16/2025
Last Edited: 07/16/2025
"""

import csv

forbidden_solvents = [
    "FC(F)(F)C(=O)O",
    "C1COCCN1CCCOS(=O)(=O)O",
    "CN=CO",
    "N=CO",
    "CCC(C)CC",
    "OCC(O)CO",
    "PBS",
    "OCC(C)OCC(C)OCC",
    "COC(CCCCC(OC)=O)=O",
    "N=C(NCCCCCCNC(=N)NC(=N)Nc1ccc(Cl)cc1)NC(=N)Nc1ccc(Cl)cc1",
    "Cl",
    "c1ccc2c(c1)c1ccccc1n2-c1ccc(-c2ccc(-n3c4ccccc4c4ccccc43)cc2)cc1",
    "O=P(c1ccccc1)(c1ccccc1)c1ccccc1Oc1ccccc1P(=O)(c1ccccc1)c1ccccc1",
    "c1cc(-c2cccc(-n3c4ccccc4c4ccccc43)c2)cc(-n2c3ccccc3c3ccccc32)c1",
    "c1ccc([Si](c2ccccc2)(c2ccc(-n3c4ccccc4c4ccccc43)cc2)c2ccc(-n3c4ccccc4c4ccccc43)cc2)cc1",
    "CN(C)P(=O)(N(C)C)N(C)C",
    "[2H]OC",
    "[2H]OCC",
    "[2H]O[2H]",
    "CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C",
    "CCCCCc1ccc(-c2ccc(C#N)cc2)cc1",
    "CCCCCCc1ccc(-c2ccc(C#N)cc2)cc1",
    "CCCCc1ccc(/N=C/c2ccc(OC)cc2)cc1",
    "[2H]C(Cl)(Cl)Cl",
    "CC(C)(C)c1ccc(-n2c3ccc([Si](c4ccccc4)(c4ccccc4)c4ccccc4)cc3c3cc([Si](c4ccccc4)(c4ccccc4)c4ccccc4)ccc32)cc1",
    "COCOC",
    "CC1CO1",
    "O=C1NCCC1C1CCCCC1",
    "O=C1CCCCCN1CO",
    "CN(C)C(=N)N(C)C",
    "CN1CCN(C)C1=O",
    "CN1CCCN(C)C1=O",
    "O=CO",
    "CN(C)C=S",
    "CCN(CC)S(=O)(=O)N(CC)CC",
]

def combine_datasets(sources):
    """
    Combine lambda max datasets from multiple sources into a single CSV file.
    Removes duplicate (compound, solvent) pairs and tags each row with its source.
    
    Parameters:
        sources (list of str): List of dataset source names to include.
    """
    base_path = f"../../intermediate/"
    output_path = "../../final/amax1_dataset.csv"
    combinations = set()

    # Open output file and write header
    with open(output_path, "w", newline="") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["compound", "solvent", "lambda_max", "source"])

    # Loop through source list
    for name in sources:
        input_path = f"{base_path}{name}_db/{name}_clean.csv"

        try:
            with open(input_path, "r") as f_in, open(output_path, "a", newline="") as f_out:
                reader = csv.reader(f_in)
                writer = csv.writer(f_out)
                header = next(reader)

                for row in reader:
                    key = (row[0], row[1])
                    if key not in combinations and row[1] not in forbidden_solvents:
                        writer.writerow([row[0], row[1], row[2], row[3] if name == "miscpapers" else name])
                        combinations.add(key)
        except FileNotFoundError:
            print(f"[WARNING] File not found: {input_path}. Skipping.")
        except Exception as e:
            print(f"[ERROR] Failed to process {input_path}: {e}")

if __name__ == "__main__":
    combine_datasets(["fluodb", "uvvisml", "smfluo1", "chemde", "photoswitch", "miscpapers"])
        


