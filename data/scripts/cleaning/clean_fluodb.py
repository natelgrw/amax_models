"""
clean_smfluo1.py

Accesses the raw FluoDB fluorescence database. Running this script
will write a new, reformatted .csv file named smfluo1_clean.csv to the
directory path "data/intermediate/smfluo1_db". Removes all data entries
that do not have an absorbance maxima measurement.

Author: @natelgrw
Created: 06/29/2025
Last Edited: 07/10/2025
"""

import csv

def clean_fluodb():
    """
    Writes a reformatted version of fluodb_raw.csv to
    specified intermediate database fluodb_clean.csv.
    """
    with open("../../raw/fluodb_db/fluodb_raw.csv", "r") as csvfile_raw:
        csv_reader = csv.reader(csvfile_raw)
        header = next(csv_reader)
        with open("../../intermediate/fluodb_db/fluodb_clean.csv", "w") as csvfile_clean:
            csv_writer = csv.writer(csvfile_clean)
            csv_writer.writerow(["compound", "solvent", "lambda_max", "source"])
            for row in csv_reader:
                if not row or row[0].strip() == "":
                    continue
                else:
                    csv_writer.writerow([row[4], row[5], float(row[0]), "fluodb"])
        csvfile_clean.close()
    csvfile_raw.close()

if __name__ == "__main__":
    clean_fluodb()
