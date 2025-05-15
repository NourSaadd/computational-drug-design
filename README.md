# Topological Indices Calculator

A user-friendly Streamlit web app to compute topological indices for molecular structures.  
Developed by **Joudy Allam** and **Nour Saad**, April 2025.

---

## Overview

This web application calculates three important topological indices from molecular structure files:

- **Edge Density**
- **Wiener Index**
- **Petitjean Index**

It provides an intuitive interface for users to upload molecular data and view computed indices along with an interactive visualization of the molecule.

---

## Supported Input Formats

The app supports:

- A **single `.mol` file** containing a molecular structure.
- A **single `.sdf` file** containing one or more molecules
- A **folder** containing multiple `.mol` files (via directory path or `.zip` upload)

---

## How It Works

1. **Choose the index** (or all indices) you wish to compute.
2. **Upload** a file or folder (as zip or path).
3. The app:
   - Parses the molecule(s)
   - Converts them into a graph representation
   - Computes the requested index/indices
   - Displays the result and an interactive molecule graph
4. Optionally, you can validate results using **built-in RDKit functions** (for `.mol` and `.sdf` files).

---

## Project Folder Contents

This folder contains all the materials for the Topological Indices Calculator project:

1. app.py
   → The main Python script for the Streamlit web application.
   → Allows users to calculate Edge Density, Wiener Index, and Petitjean Index from molecular input files.
   → Supports `.mol` files, `.sdf` files, folders, and zipped folders.

2. Coding_Project_Instructions.pdf
   → The official instructions and requirements document for this project.
   → Provided as part of the course "Computational Drug Design".

3. input/
   → Folder containing sample molecular structure inputs for testing the app.
   → Includes multiple formats and folder structures to validate the different input scenarios:

     ├── Antibiotics/
     │   → Contains individual `.mol` files (tests folder input with multiple .mol files).
     ├── n-pentane/
     │   ├── 7712.mol                           → Simple `.mol` structure
     │   ├── Conformer3D_COMPOUND_CID_8003.sdf  → 3D structure
     │   └── Structure2D_COMPOUND_CID_8003.sdf  → 2D structure
     ├── Steroids/
     │   → Contains individual `.mol` files (tests folder input with multiple .mol files).
     ├── Steroids.zip
     │   → Zipped folder of `.mol` files (tests folder input as compressed archive).
     ├── 24608638.mol
     │   → Another simple `.mol` file
     ├── 2424.mol
     │   → Another simple `.mol` file
     ├── Conformer3D.sdf
     │   → A multi-molecule `.sdf` file
     └── 2DStructures.sdf
         → A multi-molecule `.sdf` file

Note:
- The input files used in this application may originate from different chemical databases.
- In the case of the provided test inputs:
	- `.mol` files were typically sourced from ChemSpider.
	- `.sdf` files were typically sourced from PubChem.
- Since source databases use different identifiers and formats, and input files may be renamed or modified, no automatic linking to PubChem or ChemSpider websites is performed in the script.



