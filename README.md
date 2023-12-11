# MOdelling Sea Level And Inundation for Cyclones (MOSAIC)
This repository contains the code to reproduce the paper:
XXX

The scripts folder contains the scripts for the storm surge modelling (GTSM) and for the flood modelling (HydroMT-SFINCS).
The data folder contains input data for GTSM.

## GTSM: Storm Surge Modelling
The Global Tide and Surge Model (GTSM) is ran with the software environment OMUSE.
This repository contains the scripts to prepare and run the GTSM models using OMUSE, and the scripts to plot the results presented in the paper.

### Installation instructions
#### Dependencies
* AMUSE: https://github.com/amusecode/amuse                              
* OMUSE: https://github.com/omuse-geoscience/omuse.                                  

#### Installation

* Install amuse framework with `pip install amuse-framework`
* Install omuse:\
  go to omuse directory and run `pip install -e .`

## HydroMT-SFINCS: Flood modelling
### Installation instructions
Information on how to install HydroMT-SFINCS and how to run SFINCS can be found here:
* HydroMT-SFINCS: https://deltares.github.io/hydromt_sfincs/latest/
* SFINCS:https://sfincs.readthedocs.io/en/latest/

This repository contains the scripts to set-up and postprocess the SFINCS models using HydroMT-SFINCS, and the scripts to plot the results presented in the paper.

## Case studies
Three case studies have been selected for this paper:
* **TC Irma:** The tropical cyclone Irma has been forced with the coupling within OMUSE of ERA5 and the Holland model (tracks obtained from NHC). The meteorological forcing has been coupled directly to GTSM with OMUSE.
* **TC Haiyan:** The tropical cyclone Haiyan has been forced with the coupling within OMUSE of ERA5 and the Holland model (tracks obtained from JTWC). The meteorological forcing has been coupled directly to GTSM with OMUSE.
* **ETC Xynthia:** The extropical cyclone Xynthia has been forced with ERA5 downloaded and coupled to GTSM directly thanks to OMUSE.

## Model configurations
* **Baseline scenario (BS):** Default GTSM configuration as it is usually executed.
* **Temporal resolution Refined (TR):** Temporal resolution of the output of GTSM is enhanced from 1-hourly to 10-minute resolution.
* **Output resolution Refined (OR):** The output stations at which GTSM output is stored are enhanced from stations every ~5km to ~2km along the coast.
* **Improved bathymetry (IB):** The grid of GTSM is refined around the landing region of the TC/ETCs, and the bathymetry used in GTSM is enhanced by local and finer bathymetry, when available, or global but finer than Gebco2014 when no local bathymetry was available.
