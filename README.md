# winratio-sas

This repository contains SAS programs which can be used to perform unstratified and stratified win ratio analyses. 

1. winratio.sas: Macro for calculating the win ratio and its variance using CDISC BDS TTE data structure, based on code originally provided by appendix of Dong et al. (Pharm Stat 2016).

> Dong G, Li D, Ballerstedt S, Vandemeulebroecke M. A generalized analytic solution to the win ratio to analyze a composite endpoint considering the clinical importance order among components. Pharm Stat. 2016 Sep;15(5):430-7. doi: 10.1002/pst.1763. Epub 2016 Aug 2. PMID: 27485522.

2. StratifiedWR.sas: Macro for stratified win ratio analysis using CDISC BDS TTE data structure, using methodology of Dong et al. (J Biopharm Stat 2018). 

> Dong G, Qiu J, Wang D, Vandemeulebroecke M. The stratified win ratio. J Biopharm Stat. 2018;28(4):778-796. doi: 10.1080/10543406.2017.1397007. Epub 2017 Nov 27. PMID: 29172988.

3. winratio_ExampleProgram.sas: Demonstrates application of %winratio and %StratifiedWR macros to one and two stratification factors.

---------------------------------------------------------------------------------------------------------------
Correspondence: Heather Ribaudo (ribaudo@sdac.harvard.edu), Emma Davies Smith (esmith@sdac.harvard.edu)
