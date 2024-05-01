# Manuscript title: Ecological impacts of climate change will transform public health priorities for zoonotic and vector-borne disease

## Authors: David W Redding, Rory Gibb and Kate E Jones

A link to the pre-print of this manuscript (containing all methods) can be found [here](https://www.medrxiv.org/content/10.1101/2024.02.09.24302575v1).

## Paper Abstract: 
Climate change impacts on zoonotic/vector-borne diseases pose significant threats to humanity but these links are, in general, poorly understood. Here, we project present and future geographical risk patterns for 141 infectious agents to understand likely climate change impacts, by integrating ecological models of infection hazard (climate-driven host/vector distributions and dispersal) with exposure (human populations) and vulnerability (poverty prevalence). Projections until 2050, under a medium climate change (Representative Concentration Pathway (RCP)), show a 9.6% mean increase in endemic area size for zoonotic/vector-borne diseases globally (n=101), with expansions common across continents and priority pathogen groups. Range shifts of host and vector animal species appear to drive higher disease risk for many areas near the poles by 2050 and beyond. Projections using lower climate change scenarios (RCP 2.6 & 4.5) indicated similar or slightly worse future population exposure trends than higher scenarios (RCP 6.0 & 8.5), possibly due to host and vector species being unable to track faster climatic changes. Socioeconomic development trajectories, Shared Socioeconomic Pathways (SSPs), mediate future risk through a combination of climate and demographic change, which will disrupt current, regional patterns of disease burden. Overall, our study suggests that climate change will likely exacerbate global animal-borne disease risk, emphasising the need to consider climate change as a health threat.


## Navigation:

* 📁[data](https://github.com/BioDivHealth/new_global_maxent/tree/main/data) contains most of the data needed to complete the analyses. For any files that are too large, please see the **Data** section below.
* 📁[scripts](https://github.com/BioDivHealth/new_global_maxent/tree/main/scripts) contains all of the R code needed to complete the analyses.
* 📁[figures](https://github.com/BioDivHealth/new_global_maxent/tree/main/figures) contains all of the figures associated with this manuscript.

  
## Analysis order:
- 📝 [01_species_gbif.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/01_species_gbif_25.R)
- 📝 [02_maxentX2.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/02_maxent_modelling.R)
- 📝 [03_project_rasters_present.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/03_project_rasters_present.R)
- 📝 [04_Combine_maxent_models_by_transmission_model.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/04__Combine_maxent_models_by_transmission_model.R)
- 📝 [05_combine_and_future_dispersion.R](https://github.com/BioDivHealth/new_global_maxent/blob/main/scripts/05_combine_and_future_dispersion.R)


## Data: 

All of the data in this repository can be found in either the [data folder](https://github.com/BioDivHealth/new_global_maxent/tree/main/data) or in dropbox links (see below) where the files were too large to incorporate into the Github repository. 

* 📊**per_disease3** data can be downloaded [here](https://www.dropbox.com/scl/fo/gen8spncb15csjfz7thyj/ACuRzeswgWdM2DFIH8H8F7c?rlkey=wxnm9f13pv8yxavw018ubhum4&dl=0)
* 📊**livestock_future_2030_2050_2070_2080b.csv** can be downloaded [here](https://www.dropbox.com/scl/fi/7gvr5n5t4fvn99mho02bf/livestock_future_2030_2050_2070_2080b.csv?rlkey=5n6hix1ouu84mbkiruaoolprz&dl=0)
* 📊**disease_analyses2** data can be downloaded [here](https://www.dropbox.com/scl/fo/3ogm3f5bde2hjqs9oqjzk/ABD7HCw0V1fCx9IrQGNDoTk?rlkey=wsduuu4jj25m75vidh49y7tf0&dl=0)
* 📊**MODIS landcover** data can be downloaded [here](https://lpdaac.usgs.gov/products/mcd12q1v006/) and we have put a sample in the [MODIS data folder](https://github.com/BioDivHealth/new_global_maxent/tree/main/data/MODIS)
* 📊**Worldclim bioclimate** data can be downloaded [here](https://www.worldclim.org/data/bioclim.html) and we have put a sample in the [worldclim data folder](https://github.com/BioDivHealth/new_global_maxent/tree/main/data/worldclim)
