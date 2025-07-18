#!/bin/bash

# =========================================================
# raw_calculate_padel.sh
# 
# This script runs PaDEL-Descriptor.jar to calculate
# 2D descriptors for sets of compounds and
# solvents in the ChemDe, SMFluo1, and UVVisML absorbance
# maxima databases. A full table of PaDEL descriptors can 
# be accessed in the the "data" directory README.
#
# Returns six .csv files of descriptor data for all .smi
# files stored in the "data/intermediate" directory.
# 
# This script will not work if Java is not installed on
# the computer system.
# 
# To Run:
# - navigate to directory containing raw_calculate_padel.sh
# - run command: chmod +x raw_calculate_padel.sh
# - run command: ./raw_calculate_padel.sh
#
# Author: @natelgrw
# Created: 07/12/2025
# Last Edited: 07/12/2025
# =========================================================

set -e

# === PATH NAMES ===
PADEL_JAR="../../padel/PaDEL-Descriptor.jar"
INPUT_COMP_SMFLUO1="../../intermediate/smfluo1_db/compounds/comp.smi"
INPUT_COMP_UVVISML="../../intermediate/uvvisml_db/compounds/comp.smi"
INPUT_COMP_CHEMDE="../../intermediate/chemde_db/compounds/comp.smi"
INPUT_SOLV_SMFLUO1="../../intermediate/smfluo1_db/solvents/solv.smi"
INPUT_SOLV_UVVISML="../../intermediate/uvvisml_db/solvents/solv.smi"
INPUT_SOLV_CHEMDE="../../intermediate/chemde_db/solvents/solv.smi"
OUTPUT_COMP_SMFLUO1="../../intermediate/smfluo1_db/compounds/comp_descriptors.csv"
OUTPUT_COMP_UVVISML="../../intermediate/uvvisml_db/compounds/comp_descriptors.csv"
OUTPUT_COMP_CHEMDE="../../intermediate/chemde_db/compounds/comp_descriptors.csv"
OUTPUT_SOLV_SMFLUO1="../../intermediate/smfluo1_db/solvents/solv_descriptors.csv"
OUTPUT_SOLV_UVVISML="../../intermediate/uvvisml_db/solvents/solv_descriptors.csv"
OUTPUT_SOLV_CHEMDE="../../intermediate/chemde_db/solvents/solv_descriptors.csv"

# === RUNNING PADEL FOR SMFLUO ===

mkdir -p "temp"
split -l 2 "$INPUT_COMP_SMFLUO1" "temp/smi_chunk_"
for f in temp/smi_chunk_*; do mv "$f" "$f.smi"; done

echo "Starting PaDEL compound descriptor calculation for the SMFluo1 dataset..."
java -Xmx16G -jar "$PADEL_JAR" \
    -MR 1200000 \
    -2d \
    -standardizenitro \
    -removesalt \
    -threads 8 \
    -dir "temp" \
    -file "$OUTPUT_COMP_SMFLUO1";
echo "Finished with calculation. Output saved to $OUTPUT_COMP_SMFLUO1"

rm -r "temp"

mkdir -p "temp"
split -l 10 "$INPUT_SOLV_SMFLUO1" "temp/smi_chunk_"
for f in temp/smi_chunk_*; do mv "$f" "$f.smi"; done

echo "Starting PaDEL solvent descriptor calculation for the SMFluo1 dataset..."
java -Xmx16G -jar "$PADEL_JAR" \
    -maxruntime 600000 \
    -2d \
    -standardizenitro \
    -standardizetautomers \
    -removesalt \
    -threads 8 \
    -dir "temp" \
    -file "$OUTPUT_SOLV_SMFLUO1"
echo "Finished with calculation. Output saved to $OUTPUT_SOLV_SMFLUO1"

rm -r "temp"

# === RUNNING PADEL FOR UVVISML ===

mkdir -p "temp"
split -l 50 "$INPUT_COMP_UVVISML" "temp/smi_chunk_"
for f in temp/smi_chunk_*; do mv "$f" "$f.smi"; done

echo "Starting PaDEL compound descriptor calculation for the UVVisML dataset..."
java -Xmx16G -jar "$PADEL_JAR" \
    -maxruntime 1200000 \
    -2d \
    -standardizenitro \
    -standardizetautomers \
    -removesalt \
    -threads 8 \
    -dir "temp" \
    -file "$OUTPUT_COMP_UVVISML"
echo "Finished with calculation. Output saved to $OUTPUT_COMP_UVVISML"

rm -r "temp"

mkdir -p "temp"
split -l 50 "$INPUT_SOLV_UVVISML" "temp/smi_chunk_"
for f in temp/smi_chunk_*; do mv "$f" "$f.smi"; done

echo "Starting PaDEL solvent descriptor calculation for the UVVisML dataset..."
java -Xmx16G -jar "$PADEL_JAR" \
    -maxruntime 1200000 \
    -2d \
    -standardizenitro \
    -standardizetautomers \
    -removesalt \
    -threads 8 \
    -dir "temp" \
    -file "$OUTPUT_SOLV_UVVISML"
echo "Finished with calculation. Output saved to $OUTPUT_SOLV_UVVISML"

rm -r "temp"

# === RUNNING PADEL FOR ChemDe ===

mkdir -p "temp"
split -l 50 "$INPUT_COMP_CHEMDE" "temp/smi_chunk_"
for f in temp/smi_chunk_*; do mv "$f" "$f.smi"; done

echo "Starting PaDEL compound descriptor calculation for the ChemDe dataset..."
java -Xmx16G -jar "$PADEL_JAR" \
    -maxruntime 1200000 \
    -2d \
    -standardizenitro \
    -standardizetautomers \
    -removesalt \
    -threads 8 \
    -dir "temp" \
    -file "$OUTPUT_COMP_CHEMDE"
echo "Finished with calculation. Output saved to $OUTPUT_COMP_CHEMDE"

rm -r "temp"

mkdir -p "temp"
split -l 50 "$INPUT_SOLV_CHEMDE" "temp/smi_chunk_"
for f in temp/smi_chunk_*; do mv "$f" "$f.smi"; done

echo "Starting PaDEL solvent descriptor calculation for the ChemDe dataset..."
java -Xmx16G -jar "$PADEL_JAR" \
    -maxruntime 1200000 \
    -2d \
    -standardizenitro \
    -standardizetautomers \
    -removesalt \
    -threads 8 \
    -dir "temp" \
    -file "$OUTPUT_SOLV_CHEMDE"
echo "Finished with calculation. Output saved to $OUTPUT_SOLV_CHEMDE"

rm -r "temp"