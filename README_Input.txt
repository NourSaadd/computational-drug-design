README - Topological Indices Calculator Input Guide
====================================================

This folder contains a set of example input files used to test the Topological Indices Calculator script.
It demonstrates the script's ability to handle different formats and input types.

Folder Structure and Content:
-----------------------------

input/
├── Antibiotics/
│   └── *.mol                        # Multiple .mol files for folder-based input testing
│
├── n-pentane/
│   ├── 7712.mol                     # Standard .mol file for individual molecule input
│   ├── Conformer3D_COMPOUND_CID_8003.sdf  # Single molecule with 3D coordinates
│   └── Structure2D_COMPOUND_CID_8003.sdf  # Single molecule with 2D coordinates
│
├── Steroids/
│   └── *.mol                        # Multiple .mol files for folder-based input test
│
├── Steroids.zip                     # ZIP archive of the Steroids/ folder (tests ZIP upload mode)
│
├── 24608638.mol                     # Standalone .mol file to test single file input
│
├── 2424.mol                         # Another standalone .mol file
│
├── Conformer3D.sdf                  # Multi-molecule .sdf file (3D)
│
└── 2DStructures.sdf                 # Multi-molecule .sdf file (2D)
Supported Input Types by the Script:
------------------------------------

• A single .mol file containing a molecular structure.
• A folder containing multiple `.mol` files.
• A `.zip` file containing multiple `.mol` files.
• An `.sdf` file containing multiple chemical structures.


Database Notes:
---------------
The input files used in this application may originate from different chemical databases.

In the case of the provided test inputs:
    - `.sdf` files are typically downloaded from PubChem
    - `.mol` files are typically downloaded from ChemSpider