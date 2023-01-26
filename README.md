# Sociodemographic disparities in PFAS contamination
This repository holds replication code from the project titled, "Sociodemographic Factors are Associated with Abundance of PFAS Sources and Detections in U.S. Community Water Systems." 

Files are numbered according to the order in which they are run for analysis and table/figure generation. A description of the data and example data are available on the [Harvard Dataverse](https://doi.org/10.7910/DVN/8LPLCF).

# Contents
## Modelling  

- **1a_HUC processing and modelling.R**: arrange watershed data using the CWS data; spatial error regressions with watershed data
- **1b_Sources and demographics.R**: logistic regressions estimating the association. between sociodemographic factors and the presence of sources in watersheds of CWS
- **1c_PFAS-demo modelling.R**: logistic regressions estimating the association between sociodemographic factors and PFAS detections (>5 ng/L) and detections above the lowest state-level MCL in CWS
- **1d_Secondary analysis modelling.R**: secondary analysis stratifying by urban/rural classification 

Additional secondary analyses are also available in the sub-folder ("Additional secondary analyses") and they include analyses incorporating additional area-level measures of racial/ethnic composition and socioeconomic status.

## Tables and figures  

- **2_Final tables and figures.RMD**: all primary and supplemental figures and tables

# Authors  

- [Jahred Liddie](https://scholar.harvard.edu/jmliddie), Department of Environmental Health, Harvard T.H. Chan School of Public Health
- Laurel A. Schaider, Silent Spring Institute
- [Elsie M. Sunderland](https://bgc.seas.harvard.edu/), Department of Environmental Health, Harvard T.H. Chan School of Public Health; Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University