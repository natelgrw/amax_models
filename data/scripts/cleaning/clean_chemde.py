"""
clean_chemde.py

Accesses the raw ChemDe absorbance maxima database. Running this script
will write a new, reformatted .csv file named chemde_clean.csv to the
directory path "data/intermediate/chemde_db".

Author: @natelgrw
Created: 07/10/2025
Last Edited: 07/10/2025
"""

import csv

SOLV_DICT = {
    "tetrachloromethane": "C(Cl)(Cl)(Cl)Cl",
    "chloroform": "ClC(Cl)Cl",
    "chcl3": "ClC(Cl)Cl",
    "dichloromethane": "ClC(Cl)Cl",
    "methylene chloride": "ClC(Cl)Cl",
    "dcm": "C(Cl)Cl",
    "ch 2 cl 2": "C(Cl)Cl",
    "ch2cl2": "C(Cl)Cl",
    "dimethylformamide": "CN(C)C=O",
    "dmf": "CN(C)C=O",
    "tetrahydrofuran": "C1CCOC1",
    "thf": "C1CCOC1",
    "mthf": "CC1CCOC1",
    "methanol": "CO",
    "meoh": "CO",
    "ethanol": "CCO",
    "etoh": "CCO",
    "benzyl alcohol": "c1ccccc1CO",
    "dimethylsulfoxide": "CS(=O)C",
    "dimethyl sulfoxide": "CS(=O)C",
    "dmso": "CS(=O)C",
    "cyclohexane": "C1CCCCC1",
    "dioxane": "C1COCCO1",
    "1,4-dioxane": "C1COCCO1",
    "ethylacetate": "CCOC(=O)C",
    "ethyl acetate": "CCOC(=O)C",
    "etoac": "CCOC(=O)C",
    "water": "O",
    "h2o": "O",
    "acetonitrile": "CC#N",
    "acn": "CC#N",
    "mecn": "CC#N",
    "ch3cn": "CC#N",
    "ch 3 cn": "CC#N",
    "toluene": "Cc1ccccc1",
    "acetone": "CC(=O)C",
    "benzene": "c1ccccc1",
    "c6h6": "c1ccccc1",
    "chlorobenzene": "c1ccc(cc1)Cl",
    "pgmea": "CC(C)OCC(=O)C",
    "tfa": "FC(F)(F)C(=O)O",
    "hexane": "CCCCCC",
    "heptane": "CCCCCCC",
    "mops": "C1COCCN1CCCOS(=O)(=O)O",
    "propylene carbonate": "CC1OC(=O)O1",
    "isopropanol": "CC(O)C",
    "phcn": "c1ccccc1C#N",
    "pyridine": "c1ccncc1",
}

def clean_chemde():
    """
    Writes a reformatted version of chemde_raw.csv to
    specified intermediate database chemde_clean.csv.
    """
    with open("../../raw/chemde_db/chemde_raw.csv", "r") as csvfile_raw:
        csv_reader = csv.reader(csvfile_raw)
        header = next(csv_reader)
        with open("../../intermediate/chemde_db/chemde_clean.csv", "w") as csvfile_clean:
            csv_writer = csv.writer(csvfile_clean)
            csv_writer.writerow(["compound", "solvent", "lambda_max", "source"])
            for row in csv_reader:
                if row[1] == "NaN":
                    continue
                elif row[-1].lower() not in SOLV_DICT.keys():
                    continue
                csv_writer.writerow([row[0], SOLV_DICT[row[-1].lower()], float(row[1]), "chemde"])
        csvfile_clean.close()
    csvfile_raw.close()

if __name__ == "__main__":
    clean_chemde()
