# Sociodemographic disparities in PFAS contamination
This repository holds replication code for the following study:  

Liddie, J.M., Schaider, L.S., Sunderland, E.M. Sociodemographic Factors are Associated with Abundance of PFAS Sources and Detection in U.S. Community Water Systems. Environmental Science & Technology (2023) DOI: 10.1021/acs.est.2c07255

Files are numbered according to the order in which they were run for analysis and table/figure generation. Replication datasets and additional descriptions of the data are also available on the Harvard Dataverse [here](https://doi.org/10.7910/DVN/0C06MR) and [here](https://doi.org/10.7910/DVN/8LPLCF).

# Contents
## Modelling  

- **1a_HUC processing and modelling.R**: arrange watershed data by aggregating the community water system data; spatial error regressions with watershed data
- **1b_Sources and demographics.R**: logistic regressions estimating the association between sociodemographic factors and the presence of sources in watersheds of CWS
- **1c_PFAS-demo modelling.R**: logistic regressions estimating the association between sociodemographic factors and PFAS detections (>5 ng/L) and detections above the lowest state-level MCL in CWS
- **1d_Secondary + sensitivity analysis modelling.R**: secondary analyses, including analyses stratified by urban/rural classification 

Additional secondary analyses are also available in the sub-folder ("Additional secondary analyses"). They include analyses incorporating additional area-level measures of racial/ethnic composition and socioeconomic status.

## Tables and figures  

- **2_Final tables and figures.RMD**: all primary and supplemental figures and tables

# Authors  

- [Jahred Liddie](https://scholar.harvard.edu/jmliddie), Department of Environmental Health, Harvard T.H. Chan School of Public Health
- Laurel A. Schaider, Silent Spring Institute
- [Elsie M. Sunderland](https://bgc.seas.harvard.edu/), Department of Environmental Health, Harvard T.H. Chan School of Public Health; Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University

# Contact

Jahred Liddie
Harvard T.H. Chan School of Public Health
Department of Environmental Health
jliddie@g.harvard.edu