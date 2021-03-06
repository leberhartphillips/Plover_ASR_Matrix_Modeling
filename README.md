# Overview
## Demographic causes of adult sex ratio variation and their consequences for parental cooperation
#### Luke J. Eberhart-Phillips, Clemens Küpper, María Cristina Carmona-Isunza, Orsolya Vincze, Sama Zefania, Medardo Cruz-López, András Kosztolányi, Tom E. X. Miller, Zoltán Barta, Innes C. Cuthill, Terry Burke, Tamás Székely, Joseph I. Hoffman, and Oliver Krüger

In this repository you can find all the necessary files needed to reproduce the analyses presented in our [paper](https://www.biorxiv.org/content/early/2017/12/14/223941).

*(Supplementary Material Appendix A)*

**`R/`**

  - **Eberhart_Phillips_et_al_SM_Appendix_A.pdf** and **Eberhart_Phillips_et_al_SM_Appendix_A.Rmd** contains the documented code for all analyses, which can be implemented after downloading the datasets provided in the **`data/`** and **`output/`** folders.

*(Supplementary Material Appendix B)*

**`data/`**

- **Appendix_B_juvenile_adult_mark-recapture_data.txt** contains the mark-recapture field data of juveniles and adults. Each row is a single uniquely marked individual identified by their *bird_ID*. The annual encounter history of an individual is expressed in their *ch*, where a "1" indicates that an individual was encountered and "0" indicates it was not encountered. *sex* describes the molecular sex-type of an individual with "M" for males and "F" for females. *age* describes the stage at which an individual was initially captured, where "J" indicates it was first captured as a chick, and "A" indicates it was first captured as an adult. *population* describes the population in which the individual was sampled from ("KIP" = Kittlitz's plover, "MP" = Madagascar plover, "WFP" = white-fronted plover, "KPT" = Kentish plover Tuzla, "KPM" = Kentish plover Maio, and "SP" = snowy plover).

- **Appendix_B_breeding_data.txt** contains the individual reproductive histories of all marked breeding adults in the population. Each row is a nesting attempt uniquely identified by the *family_ID*. *no_chicks* expresses the number of chicks that hatched from the nest. *clutch_size* indicates the number of eggs in the nest when it was initially discovered. *year* describes the year in which the nest was active. *male* and *female* indicates the unique identity of the father and mother, respectively, with "male_NA" and "female_NA" describing cases in which the other mate was not identified. *population* describes the population in which the individual was sampled from (same notation as above).

- **Appendix_B_hatching_sex_ratio_data.txt** contains the sex and origin of each chick included in our analysis to assess hatching sex ratio. Each row is a chick uniquely identified by their *chick_ID*. The family of origin for each chick is shown in their *family_ID*. *year* describes the year in which the chick hatched. A "1" in either the *male* or *female* column indicates the molecular sex-type of a given chick. *population* describes the population in which the chick was sampled (same notation as above).

- **Appendix_B_parental_care_data.txt** contains the behavioral observations of the care-system for each family. Each row is an observation of a unique family (*family_ID*) within a given *population* (same notation as above). *care_system* expresses the parental care recorded on a given observation (i.e. "male_care", "female_care", or "both_care". *hatch_date* indicates the date on which a given brood hatched. *date* indicates the date on which the observation was made.
  
**`output/`**

- **AIC_table_juvenile_adult_boot_out.txt** contains the bootstrap output for model selection of juvenile and adult survival based on the mark-recapture analysis run in Program MARK. Each row is a *model* fitted via maximum likelihood to the bootstrapped data sample of each iteration (*iter*). *Phi* describes the model structure for fitting annual survival. *p* describes the model structure for fitting annual encounter probability. *npar* reveals the number of parameters used in a given model. *AICc* is the Akaike Information Criteria statistic corrected for small sample size. *DeltaAICc* is the difference in AICc between a given model and the best fit model of a given iteration. *weight* describes the AIC weight of a given model. *Deviance* describes the deviance of a given model. *population* specifies the population from which the analysis of the iteration was based on (same notation as above).

- **ASR_boot_out_final.txt** contains the adult sex ratio estimates (*ASR_boot_out*) of each iteration of the bootstrap procedure. Each row represents an iteration (*iter*). *population* specifies the population from which the analysis of the iteration was based on (same notation as above).

- **survival_rates_boot_out_final.txt** contains the sex- and stage-specific survival estimates (*estimate*) of each iteration (*iter*) in the bootstrap procedure. Each row represents a given sex and stage (*sex_age*) in a given iteration. *population* specifies the population from which the analysis of the iteration was based on (same notation as above).

Last updated: February 20, 2018
