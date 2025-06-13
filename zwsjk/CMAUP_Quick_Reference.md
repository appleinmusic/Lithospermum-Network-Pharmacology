# CMAUP v2.0 Quick Reference Guide

**Last Updated:** 2025-06-01

## 1. Introduction

This document provides a quick reference guide for the CMAUP v2.0 database. Its purpose is to help users quickly understand the database structure, locate relevant data files, and identify key information for their research.

**For the most detailed and accurate descriptions of each data file and their column names, please refer to `Download_Readme.txt` located in the same directory as the data files.** The `CMAUP_v2.0_Database_Architecture_Manual.md` can serve as a supplementary architectural overview.

## 2. Data File Location

All raw CMAUP v2.0 data files (original `.txt` downloads and the `Download_Readme.txt`) are located in the following directory:

`/Users/lgmoon/Desktop/zdhky/zwsjk/`

The transcriptomic change data is located in a subdirectory:

`/Users/lgmoon/Desktop/zdhky/zwsjk/005_zip_of_transcriptomic_change/`

## 3. Core Data Files Overview & Key Information Sources

The CMAUP v2.0 database includes several core `.txt` files. **Refer to `Download_Readme.txt` for complete column lists and detailed descriptions.**

*   **`Plants.txt`**: Plant information.
*   **`Ingredients_All.txt`**: All chemical ingredients.
*   **`Ingredients_onlyActive.txt`**: Pre-filtered active ingredients. (Note: "Active" criteria details in `Download_Readme.txt` or main manual).
*   **`Targets.txt`**: Molecular target information.
*   **`Plant_Ingredient_Associations_allIngredients.txt`**: Plant-to-all-ingredient associations.
    *   **CRITICAL NOTE:** Original file **may lack a header row**. Verify and add headers if needed, using `Download_Readme.txt` for column names.
*   **`Plant_Ingredient_Associations_onlyActiveIngredients.txt`**: Plant-to-active-ingredient associations.
    *   **CRITICAL NOTE:** Original file **may lack a header row**. Verify and add headers if needed.
*   **`Ingredient_Target_Associations_ActivityValues_References.txt`**: Ingredient-target associations with activity data.
*   **`Human_Oral_Bioavailability_information_of_Ingredients_All.txt`**: ADMET-related properties.
    *   **NOTE on Column Name:** The column `Prediction_Hob&Swissadme` (and potentially others listed in `Download_Readme.txt`) may contain special characters. In R, use backticks: `` `Prediction_Hob&Swissadme` ``.
*   **`Plant_Human_Disease_Associations.txt`**: Plant-disease associations.
*   **`Plant_Clinical_Trials_Associations.txt`**: Plant-clinical trial associations.
*   **`Plant_molecular_targets_overlapping_with_DEGs.txt`**: Plant targets overlapping with disease DEGs.
*   **`Download_Readme.txt`**: **PRIMARY SOURCE for detailed file contents, column names, and data descriptions.**

## 3.1. Transcriptomic Change Data (Subdirectory: `005_zip_of_transcriptomic_change/`)

**Refer to `Download_Readme.txt` (section: "CMAUPv2.0_download_Transcriptomic_profiling_of_74_diseases") for details on this data.**

*   **`001_index_of_each_disease_analyzed_in_transcriptomic_change.txt`**: Index of diseases with transcriptomic analysis.
*   **`002_meta_data/` (Directory)**: Contains sample metadata for each disease analysis (e.g., `DS1.txt`, `DS2.txt`, ... corresponding to diseases in the index).
*   **`003_Results_produced_by_DESeq2/` (Directory)**: Contains DESeq2 differential expression results for each disease (e.g., `DS1.csv`, `DS2.csv`, ...).

## 4. Key Identifiers

(Refer to `Download_Readme.txt` for context on how these IDs are used in specific files)
*   **`Plant_ID`**: Plant identifier (e.g., NPOxxxx).
*   **`Ingredient_ID`**: Ingredient identifier (e.g., NPCxxxx).
*   **`Target_ID`**: Target identifier (e.g., NPTxxxx).
*   **`Disease_ID`** (in transcriptomic data): Disease identifier (e.g., DS1, DS2 from the index).

## 5. Common Analysis Starting Points (Always verify column names with `Download_Readme.txt`)

*   **Plant's Ingredients & Targets:**
    1.  `Plants.txt` -> `Plant_Ingredient_Associations_allIngredients.txt` -> `Ingredients_All.txt` -> `Ingredient_Target_Associations_ActivityValues_References.txt` -> `Targets.txt`.
*   **ADMET Screening:**
    1.  `Ingredients_All.txt` -> `Human_Oral_Bioavailability_information_of_Ingredients_All.txt`.
*   **Disease-related Transcriptomics & Plant Targets:**
    1.  `005_zip_of_transcriptomic_change/001_index_of_each_disease_analyzed_in_transcriptomic_change.txt` (to identify disease of interest and its `DS` code).
    2.  `005_zip_of_transcriptomic_change/003_Results_produced_by_DESeq2/DSX.csv` (for DEGs).
    3.  `Plant_molecular_targets_overlapping_with_DEGs.txt` (to link plant targets to these DEGs).

## 6. Critical Reminders

*   **Consult `Download_Readme.txt` first** for the most accurate file and column details.
*   **Missing Headers:** Files like `Plant_Ingredient_Associations_allIngredients.txt` and `Plant_Ingredient_Associations_onlyActiveIngredients.txt` may lack headers.
*   **Special Characters in Column Names:** Handle these appropriately (e.g., backticks in R).
*   **Data Types:** Inspect and convert data types as needed.

## 7. For Broader Architectural Overview

For a general architectural understanding, you can also consult `CMAUP_v2.0_Database_Architecture_Manual.md`.
