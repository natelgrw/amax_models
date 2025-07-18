"""
raw_calculate_rdkit.py

This script runs RDKit to calculate chemical and fingerprint descriptors
for sets of compounds and solvents in the ChemDe, SMFluo1, and UVVisML
absorbance maxima databases. A full table of calculated RDKit descriptors
can be accessed in the "data" directory README.

Returns six .csv files of descriptor data for all .smi
files stored in the "data/intermediate" directory.

Author: @natelgrw
Created: 07/15/2025
Last Edited: 07/15/2025
"""