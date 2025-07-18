"""
clean_photoswitch.py

Accesses the raw Photoswitch fluorescence database. Running this script
will write a new, reformatted .csv file named smfluo1_clean.csv to the
directory path "data/intermediate/photoswitch_db". Removes all data entries
that do not have an absorbance maxima measurement.

Author: @natelgrw
Created: 07/16/2025
Last Edited: 07/16/2025
"""

import csv

SOLV_DICT = {
    "ethanol": "CCO",
    "methanol": "CO",
    "water": "O",
    "toluene": "Cc1ccccc1",
    "acetonitrile": "CC#N",
    "dimethyl sulfoxide": "CS(=O)C"
}

def clean_photoswitch():
    """
    Writes a reformatted version of photoswitch_raw.csv to
    specified intermediate database photoswitch_clean.csv.
    """
    with open("../../raw/photoswitch_db/photoswitch_raw.csv", "r") as csvfile_raw:
        csv_reader = csv.reader(csvfile_raw)
        header = next(csv_reader)
        with open("../../intermediate/photoswitch_db/photoswitch_clean.csv", "w") as csvfile_clean:
            csv_writer = csv.writer(csvfile_clean)
            csv_writer.writerow(["compound", "solvent", "lambda_max", "source"])
            for row in csv_reader:
                if not row or row[6].strip() == "" or row[3].strip() == "":
                    continue
                else:
                    csv_writer.writerow([row[2], SOLV_DICT[row[6]], float(row[3]), "photoswitch"])
        csvfile_clean.close()
    csvfile_raw.close()

if __name__ == "__main__":
    clean_photoswitch()
