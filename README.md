# Ecological impacts of climate change will transform public health priorities for zoonotic and vector-borne disease

## Authors: David W Redding, Rory Gibb and Kate E Jones

A link to the pre-print of this manuscript (containing all methods) can be found [here](https://www.medrxiv.org/content/10.1101/2024.02.09.24302575v1).

## Paper Abstract: 
Climate change impacts on zoonotic/vector-borne diseases pose significant threats to humanity1 but these links are, in general, poorly understood2. Here, we project present and future geographical risk patterns for 141 infectious agents to understand likely climate change impacts, by integrating ecological models of infection hazard (climate-driven host/vector distributions and dispersal) with exposure (human populations) and vulnerability (poverty prevalence). Projections until 2050, under a medium climate change (Representative Concentration Pathway (RCP)), show a 9.6% mean increase in endemic area size for zoonotic/vector-borne diseases globally (n=101), with expansions common across continents and priority pathogen groups. Range shifts of host and vector animal species appear to drive higher disease risk for many areas near the poles by 2050 and beyond. Projections using lower climate change scenarios (RCP 2.6 & 4.5) indicated similar or slightly worse future population exposure trends than higher scenarios (RCP 6.0 & 8.5), possibly due to host and vector species being unable to track faster climatic changes. Socioeconomic development trajectories, Shared Socioeconomic Pathways (SSPs), mediate future risk through a combination of climate and demographic change, which will disrupt current, regional patterns of disease burden. Overall, our study suggests that climate change will likely exacerbate global animal-borne disease risk, emphasising the need to consider climate change as a health threat.


## Navigation
- navigation for each of the folders
- simple diagram of the workflow

## Analysis order:
- [01_species_gbif.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/01_species_gbif_25.R)
- [02_maxentX2.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/02_maxent_modelling.R)
- [03_project_rasters_present.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/03_project_rasters_present.R)
- [04_Combine_maxent_models_by_transmission_model.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/04__Combine_maxent_models_by_transmission_model.R)
- [05_combine_and_future_dispersion.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/05_combine_and_future_dispersion.R)


## Data: 

All of the data in this repository can be found in either the [data folder](https://github.com/BioDivHealth/new_global_maxent/tree/main/data) or in Zenodo links [# TO DO: ADD] where the files were too large to incorporate into the Github repository. 
