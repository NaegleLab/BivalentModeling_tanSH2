# BivalentModeling
Adaptation of Errington et al. Matlab code to predict tandem SH2 domain interactions


<p align="center">
    <img src="Bivalent%20Kinetics%20nomenclature.jpg" alt="Overview of Model Parameters" width="300">
</p>

## How to use this code
There are three key tools in this repository. Two are structure-based evaluation to help parameterize tandem SH2 interactions with pTyr and the third is an adapted ODE modeling approach in Matlab. 
### Linker estimation from structure
In PDB_Linker_Estimation use cal_atomicdistances.ipynb to measure the distance between the atoms on the opposite end of the linkers defined by the SH2 domains. This works specifically only for tandem SH2 domains, but can easily be extended to cover the measurement of any structure between any pair of residues. This just happens to use our predefined tanSH2_linkers.txt file, which uses the UniProt defined boundaries between domains. 

### Domain radius estimation from structure
This code measured the pairwise distances between backbone atoms (N and C) within the defined domain boundaries and records the maximum distance as the domain's diameter. It also collects information across all PDB structures for the same domain. This code includes calculations for all the PDB structures used in Portelance et al.

### Matlab code 
We received the original code from Dr. Casim Sarkar from the [2019 PNAS publication](https://pubmed.ncbi.nlm.nih.gov/31776263/). We adapted this code to specifically apply this to tandem SH2 domains and pTyr residues, giving greater access over the linker parameters and states of the bound model and to calculate an effective equilibrium dissasociation constant. To use this code, simply update the parameters (in **parameters.m** and, if needed, **EffC_Calculator.m**) as you see fit (see all parameters we used in tandem SH2 domain predictions in [SH2_prediction_parameters.csv](Bivalent_Matlab_Code/SH2_prediction_parameters.csv)

Run **Multivalent()** which will calculate a KD and KA value for the global parameters set in **[parameters.m](Bivalent_Matlab_Code/parameters.m)**. For directly changing linkers, these are in **[EffC_Calculator.m](Bivalent_Matlab_Code/EffC_Calculator.m)** and noted in [SH2_prediction_parameters.csv](Bivalent_Matlab_Code/SH2_prediction_parameters.csv) at the bottom of the columns. 
