# Steps to estimate the diameter of protein domain
## Step1 - Identify the PDB structures that contain the domain of interest - can be performed using CoDIAC [BioRxiv paper](https://www.biorxiv.org/content/10.1101/2024.07.18.604100v1)
Using CoDIAC, one can identify relevant structures, download pdb files and retrive domain boundaries information.

## Step2 - PDB files will be editted to a different format to make x, y, z co-ordinate extraction simpler and saved into new pdb files

## Step3 - Retain ATOM records information from these new editted pdb files for atoms of interest. Here, we chose to make calculations for backbone atoms and we should specify the domain boundaries that is obtained using CoDIAC. 

## Step4 - Longest distance from the list of pairwise distances calculated will be represented as the diameter of the domain sphere. 

## Step5 - Average of the diameters across multiple PDBs is measured
