"""
combine_smiles.py

Accesses the full AMAX-1 dataset, writing all unique compounds and solvents
to comp.smi and solv.smi lists, stored in respective directories.

Author: @natelgrw
Created: 7/17/2025
Last Edited: 7/18/2025
"""

import csv

def write_amax_smiles():
    """
    Writes unique compounds and solvents to .smi lists
    in the "data/final/compounds" and "data/final/solvents"
    directories.
    """
    amax_path = "../../final/amax1_dataset.csv"
    comp_path = "../../final/compounds/comp.smi"
    solv_path = "../../final/solvents/solv.smi"

    comp_counter = 1
    solv_counter = 1

    comps = set()
    solvs = set()

    with open(amax_path, "r") as f_amax, open(comp_path, "w") as f_comp, open(solv_path, "w") as f_solv:
        csv_reader = csv.reader(f_amax)
        header = next(csv_reader)

        for row in csv_reader:
            if row[0] not in comps:
                f_comp.write(f"{row[0]} mol{comp_counter}\n")
                comps.add(row[0])
                comp_counter += 1
            if row[1] not in solvs:
                f_solv.write(f"{row[1]} mol{solv_counter}\n")
                solvs.add(row[1])
                solv_counter += 1

if __name__ == "__main__":
    write_amax_smiles()