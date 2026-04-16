# Juvenile Fin-Clamp: Open hardware designs, sensor processing scripts and CFD resources for minimally invasive biologging of juvenile sharks

This repository contains all hardware design files, sensor processing scripts, 
and computational fluid dynamics (CFD) resources associated with the following publication:

> Citation and DOI to be added upon publication.

## Repository Structure

### 🗂️ Fin-clamp production
Resources to reproduce the Juvenile Fin-Clamp hardware.
- **CAD resources** — CAD files in .step (universal) and .f3z (Fusion 360) formats
- **Instructions** — Step-by-step production instructions, PDF technical plans and assembly resources

### 📊 Sensor processing scripts
R scripts for processing LittleLeonardo accelerometry, depth and temperature data 
to reproduce the results of the publication. Includes an R script adapted from the 
MATLAB code of Cade et al. (2018) to run jiggle method calculations.

### 🌊 CFD
Resources for the computational fluid dynamics simulations.
- **3D models** — Simplified .stl models of the Juvenile Fin-Clamp and juvenile 
  smooth hammerhead shark used in SimScale simulations
- **Simulation results** — SimScale simulation outputs, summary table, and R scripts 
  for processing results and producing figures

## Authors

C. Robert Priester, Jorge Fontes, Fred Buyle, Diya Das, Bruno C. L. Macena, 
Colin Simpfendorfer, Esben Moland Olsen, Pedro Afonso

## Software Requirements

### R Scripts
The sensor processing and CFD result scripts were developed in R. We recommend 
using a recent version of R (≥ 4.0) and RStudio. Required packages are listed 
at the top of each script. The jiggle method script is adapted from the original 
MATLAB code of Cade et al. (2018) for use in R.

> Cade, D.E., Barr, K.R., Calambokidis, J., Friedlaender, A.S. and Goldbogen, J.A. 
> (2018). Determining forward speed from accelerometer jiggle in aquatic 
> environments. Journal of Experimental Biology, 221(2).

### CFD Simulations
Simulations were run using [SimScale](https://www.simscale.com/), a web-based 
CFD platform. The .stl files in the /CFD/3D models folder can be imported 
directly into SimScale to reproduce or modify the simulations.

## Data Availability

Raw tag data used and analyzed in this study are not publicly archived due to 
ongoing research but are available from the corresponding author upon reasonable 
request.

## License

Hardware design files (CAD, production files): CERN Open Hardware Licence v2 - Strongly Reciprocal (CERN-OHL-S v2) 
Documentation and instructions: Creative Commons Attribution 4.0 International (CC BY 4.0)

## Contact

**C. Robert Priester** (corresponding author)
robert.priester@posteo.de

For questions related to fin-clamp development and contributions, please use 
the [GitHub Issues](../../issues) function of this repository.
