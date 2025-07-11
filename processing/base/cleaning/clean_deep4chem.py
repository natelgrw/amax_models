"""
clean_deep4chem.py

Accesses the raw Deep4Chem absorbance maxima database. Running this script
will write a new, reformatted .csv file named deep4chem_clean.csv to the
directory path "data/intermediate/deep4chem_db".

Author: @natelgrw
Created: 07/09/2025
Last Edited: 07/10/2025
"""

import csv

def clean_deep4chem():
    """
    Writes a reformatted version of deep4chem_raw.csv to
    specified intermediate database deep4chem_clean.csv.
    """
    with open("../../../data/raw/deep4chem_db/deep4chem_raw.csv", "r") as csvfile_raw:
        csv_reader = csv.reader(csvfile_raw)
        header = next(csv_reader)
        with open("../../../data/intermediate/deep4chem_db/deep4chem_clean.csv", "w") as csvfile_clean:
            csv_writer = csv.writer(csvfile_clean)
            csv_writer.writerow(["smiles", "solvent", "lambda_max"])
            for row in csv_reader:
                if row[3] == "NaN":
                    continue
                else:
                    csv_writer.writerow([row[1], row[2], float(row[3])])
        csvfile_clean.close()
    csvfile_raw.close()

if __name__ == "__main__":
    clean_deep4chem()