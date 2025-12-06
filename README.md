The repository of "Multiple key hosts and network structure shape viral prevalence across multispecies communities of bees" manuscript.

To recreate the analysis, install `renv` package (`install.packages("renv")`) and run `renv::restore()`. This function will install all required packages and their appropriate versions.

The data required for this analysis is available at figshare: (link coming soon).

The repository contains:
- `01_plotting_functions.R` - functions used to summarise and plot results
- `02_land_extrapolation.R` - harmonizing landuse data and extrapolating
- `03_networks.R` - calculating network characteristics and metrics (connectance, resource overlap, plant and pollinator richness)
- `04_patogen_data_processing.R` - harmonizing and quality-filtering viral screening data
- `05_main_data_file.R` - preparing the main dataset for the analysis
- `06_R0.R` - estimating true prevalence and R0,i for each virus and bee species
- `07_analysis_models` - statistical analysis
- `08_tables_figures` - visualisation of results (all tables and figures)

All the files are runnable independently by sourcing required upstream files.

