"""
raw_extract_comp.py

Accesses all compound-solvent combinations in the UVVisML, SMFluo1, Deep4Chem,
and ChemDe absorbance maxima databases. Running this script writes a 
comprehensive list of compounds and solvents to intermediate .smi text files
for each database.
"""

import csv

def write_smiles(database_name):
    """
    Writes compound and solvent smiles to .smi text files stored in
    "data/intermediate" directories given an input database name.
    """
    counter_comp = 1
    counter_solv = 1
    with open(f"../../data/intermediate/{database_name}_db/{database_name}_clean.csv", "r") as csvfile_clean:
        with open(f"../../data/intermediate/{database_name}_db/compounds/comp.smi", "w") as comp_file:
            with open(f"../../data/intermediate/{database_name}_db/solvents/solv.smi", "w") as solv_file:
                comp_writer = csv.writer(comp_file)
                solv_writer = csv.writer(solv_file)
                csv_reader = csv.reader(csvfile_clean)
                header = next(csv_reader)
                for row in csv_reader:
                    comp_writer.writerow([f"{row[0]} mol{counter_comp}"])
                    solv_writer.writerow([f"{row[1]} mol{counter_solv}"])
                    counter_comp += 1
                    counter_solv += 1
            solv_file.close()
        comp_file.close()
    csvfile_clean.close()

if __name__ == "__main__":
    write_smiles("uvvisml")
    write_smiles("smfluo1")
    write_smiles("chemde")
