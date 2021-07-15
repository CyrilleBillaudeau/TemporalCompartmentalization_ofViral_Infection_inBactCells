# Temporal Compartmentalization of Viral Infection in Bacterial Cells

In the framework of our study on the "Temporal compartmentalisation of viral infection in bacterial cells" (Labarde et al. published in PNAS in 2021, https://doi.org/10.1073/pnas.2018297118), we have developed an automatic analysis protocol to detect and characterize phage DNA compartments, procapsids and virion warehouses using a Matlab-based custom algorithm (Matlab R2016B version). The analysis workflow can be decomposed in three main successive steps:
1. Cell detection and identification. This step aims to generate a binary image in which cells are individualized.
2. Localization and morphology of the phage DNA compartment.
3. Localization and morphology of procapsids and virion warehouses.

The codes developed and used for this study are available in the 'Script' folder. They will require the "Image Processing" toolbox to run. 

A dataset is also provided in the 'Data' folder to test the workflow.

The analysis protocol has been described in our paper and is also included in the header of the main functions "Cell_detection_and_identification.m" and "LocMorpho_PhageDNA_compartment_procapsids_virion_warehouses.m" ('Script' folder).
