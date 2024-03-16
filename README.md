# Virtual tumor manuscript

This GitHub repository contains all the files needed to reproduce the main outputs of the following manuscript (currently under submission in CPT:PSP):

**QSP Modelling of tumor heterogeneity in response to BH3-mimetics using virtual tumors calibrated with cell viability assays**

Thibaud Derippe1,2,3*, Sylvain Fouliard1, Xavier Decleves2, Donald E. Mager3,4

# Structure of the repository

The repository contains the following main R scripts:

+ **0_Lindner_model_PaSM_config.R:** includes the RxODE-formated modified Lindner model, along with PaSM configuration information (eg, monotonic parameter attribution),
+ **1_1_Virtual_Cells_Generation.R** uses [PaSM algorithm](https://github.com/Thibaudpmx/PaSM) to generate Virtual Cells and compute their fate (death or survival) under all 79 different drug exposure set-ups,
+ **1_2_Rearrange_VC_files.R** performs several manipulations with previously generated VCs, especially attributing them in 'bags', and selecting a single VC per bag depending on the number of drugs to consider. These Celltheques/main outputs of the two VC's related scripts are directly provided  in the folder 'Virtual_Cells/one_per_bag_combined' (other generated files were too voluminous to commit) for you can run the next scripts without the need to complete this time-consuming VCs generation, 
+ **2_Virtual_Tumor_objects.R** contains the main function to create, calibrate and evaluate/plot the virtual tumors (coded as R6 object). This script is not made to be manipulated, but used in the next one,
+ **3_Virtual_Tumor_Calibration.R** is the main script in which Virtual tumor are created and calibrated. The outputs of this script (the calibrated virtual tumors) are provided in the folder 'calibrated_VT',
+ **4_agent_based_model_Fig6.R**, **4_cell_composition_fig5.R** and **4_other_figures_article_vt.R** are the three scripts used to generate the manuscript's figures, including the heterogeneity analyses of the calibrated virtual tumors (fig5) and minimal-ABM attempts (fig 6). These 3 files are independent (no particular order to be used) and are using the final calibrated VTs already provided in this repository.


The three additional folders are:

+ **data** contains the two datasets used in this work, namely the  [cell viability assays](https://pubmed.ncbi.nlm.nih.gov/26565405/) and the  [mice data](https://www.nature.com/articles/s41375-019-0652-0), both digitized from two articles coming from D.C. Phillips *et al* (the two links above lead to the source articles).
+ **Virtual_Cells** contains the final celltheques generated after using the scripts *1_1_Virtual_Cells_Generation.R* and *1_2_Rearrange_VC_files.R*
+ **calibrated_VT** corresponds to the final outputs after using the script *3_Virtual_Tumor_Calibration.R*

# Zenodo deposit

A Zenodo deposit was created at the time of the manuscript submission to ensure sustained availability regardless of the future of this GitHhub repository: 

# Contact

For any questions or issues with the code, feel free to contact me at Thibaud.Derippe@gmail.com. However, please note I do not work in the QSP field anymore, and therefore I do not intend to extend this work or improve the algorithms. 

Thanks !

Thibaud 
