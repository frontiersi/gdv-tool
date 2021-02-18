[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4547661.svg)](https://doi.org/10.5281/zenodo.4547661)

# GDV Tool (Workshop Preview)

## Introduction
Mine dewatering is the removal of unwanted groundwater to allow rock and mineral extraction. In some circumstances, this can affect the health of groundwater dependent vegetation (GDV), which relies on a stable water-table for its water requirements. Monitoring GDV is thus a requirement around mining leases but is poorly targeted.
Narrowing the search space using remote sensing imagery is the first step towards a time saving, cost efficient, and comprehensive monitoring program that can provide peace of mind that dewatering is being conducted appropriately. However, to date there is no standard tool to simplify this process.
To address this need, Curtin researchers (Dr. Todd Robinson (lead), Lewis Trotter (Research Fellow) and Dr. Adam Cross (Research Fellow) in consultation with mining companies (Roy Hill, BHP) have recently combined weighted seasonal imagery from the DEA Open Data Cube, Python and Jupyter Notebooks to map GDV.

## Objectives
This project aimed to develop validated remote sensing-based methods for detecting GDV and its change over time. It had the following objectives:
1. Develop GDV likelihood models using spatial multicriteria evaluation (SMCE) and Landsat 5, 7 and 8/Sentinel 2A and 2B satellite imagery.
2. Optimise model input using validation statistics obtained from ground-truthing data.
3. Develop trend and change models and demonstrate their applicability over three study areas with different histories.

The resulting Python code used to run this model has been uploaded online (this github). A Jupyter Notebook (GDV Tool 0.5.ipynb) has also been provided to provide a rough demonstration of the GDV Tool's functionality.

## Requirements
1. A Digital Earth Australia (DEA) Open Data Cube (ODC) Sandbox account is required to run the tool. Sign up at https://app.sandbox.dea.ga.gov.au/.
2. Basic Python and Jupyter Notebook experience.
3. Chrome or Edge web-browser.

## Limitations
Currently, the GDV Tool has only been tested in several locations within 100km of the township of Newman in the Pilbara region of Western Australia. Accuracy of this tool outside of this general area is currently unknown and needs statistical exploration.
Additionally, additional methods are currently being developed that have shown to improve the accuracy of this tool. These methods will be integrated into this tool at a later date. It is recommended that the code provided here is used only for experimentation and should not be used as part of any decision making processes.

## Setup information
A few steps are required to get the GDV Tool up and running.
1. Sign up for a DEA ODC Sandbox account here: https://app.sandbox.dea.ga.gov.au/.
2. Once logged in, open a new terminal window and execute the following: git clone https://github.com/frontiersi/GDV_TOOL
3. This will pull the latest version of the code from github and install it automatically into your Sandbox environment under the folder GDV_TOOL.
4. Open the GDV_TOOL folder and run the GDV Tool 0.5.ipynb file. THis will run a Jupyter Notebook that provides a demonstration of the tools functions.
5. Run through the Jupyter Notebook from top to bottom.

## Dependencies
- pyannkendall
- ruptures

## Licences
Apache License Version 2.0
