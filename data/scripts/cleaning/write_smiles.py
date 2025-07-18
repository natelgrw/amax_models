"""
write_smiles.py

Accesses intermediate databases, writing a comprehensive
list of compounds and solvents to comp.smi and solv.smi in respective
intermediate database directories. This is done in preparation
for descriptor calculation and analysis.

Author: @natelgrw
Created: 07/17/2025
Last Edited: 07/17/2025
"""

import csv

def write_comps_and_solvs(db_list):
    """
    Writes a comprehensive list of compounds and solvents to
    intermediate database directories for descriptor calculation 
    and analysis.
    """
    for db in db_list:
        path_database = f"../../intermediate/{db}_db/{db}_clean.csv"
        path_comp = f"../../intermediate/{db}_db/{db}_comp.smi"
        path_solv = f"../../intermediate/{db}_db/{db}_solv.smi"

        counter_comp = 1
        counter_solv = 1

        comps = set()
        solvs = set()

        with open(path_database, "r") as f_db, open(path_comp, "w") as f_comp, open(path_solv, "w") as f_solv:
            csv_reader = csv.reader(f_db)
            header = next(csv_reader)

            for row in csv_reader:
                if row[0] not in comps:
                    f_comp.write(f"{row[0]} mol{counter_comp}\n")
                    comps.add(row[0])
                    counter_comp += 1
                if row[1] not in solvs:
                    f_solv.write(f"{row[1]} mol{counter_solv}\n")
                    solvs.add(row[1])
                    counter_solv += 1

if __name__ == "__main__":
    write_comps_and_solvs(["chemde", "fluodb", "photoswitch", "smfluo1", "uvvisml", "miscpapers"])
