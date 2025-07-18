"""
clean_smfluo1.py

Accesses the raw SMFluo1 absorbance maxima database. Running this script
will write a new, reformatted .csv file named smfluo1_clean.csv to the
directory path "data/intermediate/smfluo1_db".

Author: @natelgrw
Created: 06/29/2025
Last Edited: 07/10/2025
"""

import csv

def clean_smfluo1():
    """
    Writes a reformatted version of smfluo1_raw.csv to
    specified intermediate database smfluo1_clean.csv.
    """
    with open("../../raw/smfluo1_db/smfluo1_raw.csv", "r") as csvfile_raw:
        csv_reader = csv.reader(csvfile_raw)
        header = next(csv_reader)
        with open("../../intermediate/smfluo1_db/smfluo1_clean.csv", "w") as csvfile_clean:
            csv_writer = csv.writer(csvfile_clean)
            csv_writer.writerow(["compound", "solvent", "lambda_max", "source"])
            for row in csv_reader:
                csv_writer.writerow([row[0], row[1], float(row[2]), "smfluo1"])
        csvfile_clean.close()
    csvfile_raw.close()

if __name__ == "__main__":
    clean_smfluo1()
