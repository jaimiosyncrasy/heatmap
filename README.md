# **Visual Tool for Assessing Stability of DER Configurations on Three-Phase Radial Networks**
## *Summary of Tool*
This repo contains code for evaluating control configurations of power-voltage control laws on three-phase radial distribution grid feeders. Three main processes for doing so are:
* Non-colocated placement process (NPP)
* Co-located Placement Process (CPP)
* Auto-colocated placement process (Auto-CPP)

Each process generally analyzes configurations of actuator-performance node pairs by evaluating the eigenvalue stability of the resulting state space model and generating heatmaps to illustrate the results. 

More details on these processes as well as the model used to evaluate stability is available in the associated paper at: https://arxiv.org/abs/2011.07232

You can view the full set of figures associated with the paper by viewing "imgs_forPaper.ipynb". To run/modify the code that generates these figures, follow the installation instructions below

## *Installation*
### Python Packages that Need to be Installed
* Anaconda Navigator with Jupitor Notebook 
* xlrd - [installation instructions](https://xlrd.readthedocs.io/en/latest/installation.html)
* networkx - [installation instructions](https://networkx.org/documentation/stable/install.html)
* pandas - [installation instructions](https://pypi.org/project/pandas/)
* matplotlib - [installation instructions](https://pypi.org/project/matplotlib/)
* graphviz & pygraphviz
  * to install, open Anaconda Navigator, select "Environments" on the left-hand toolbar, make sure the dropdown menu is set to "Not Installed," then search "graphviz" in the search field on the far right. Check the boxes next to "graphviz" and "pygraphviz."
* scipy - [installation instructions](https://pypi.org/project/scipy/)
## *Using the Tool*
### Getting Started - Reproducing our Results
1. Open Jupitor Notebook in Anaconda Navigator and navigate to the folder containing the cloned github repo
2. To reproduce the results described in the paper, open the python file "driver_code_forPaper.ipynb"
3. Run the code blocks in the order that they appear to produce the results and .png files of various heatmaps
4. To view all of the generated heatmaps, open the python file "imgs_forPaper.ipynb" and run the code blocks
### Introduction to py_modules
The main file is "driver_code_forPaper.ipynb,". The main file calls functions in the .py module files, which are collected into the folder "py_modules". The .py module files contain additional functions not called by the main file currently, but that may be incorporated in future releases.

***setup_nx.py***
  * reads line data from the impedance excel files to create the feeder object used for the other .py files
  
***my_feeder_funcs.py***
  * contains functions that initialize the feeder object and reset the feeder's graphviz graph

***my_impedance_funcs.py***
  * contains functions which allow you to find the impedance between any two nodes on a network

***my_detControlMatExistence_funcs.py***
  * contains the functions used to determine if a particular actuator configuration is stable, through eigenvalue analysis of the closed-loop system

***my_detLznRange_funcs.py***
  * contains functions that set up the linearized power flow model and use the model to solve for line losses. Not used in this version of the code

***my_heatmapSetup_funcs.py***
  * contains the functions called to run the placement processed and produce majority of the heatmaps. "runHeatMapProcess" runs the NPP, "placeMaxColocActs_stopAtInfeas" runs the OCPP, and "place_max_coloc_acts" runs the Auto-OCPP. 

***my_configVis_funcs.py***
  * contains functions primarily dedicated to dividing a distribution network into branches, and analyzing which branches are the best for placing actuators


