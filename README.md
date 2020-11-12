# **Visual Tool for Assessing Stability of DER Configurations on Three-Phase Radial Networks**
## *Summary of Tool*
This repo contains code to:
* analyze the stability of different actuator configurations on a distribution network
* determine the configuration which allows for the maximum number of colocated actuator performance nodes on a network
## *Installation*
### Downloading the Repository
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
The folder "py_modules" holds the .py files which contain the functions called by "driver_code_forPaper.ipynb," as well as additional functions which can be used to analyze and visualize actuator configurations on distribution networks.
*my_configVis_funcs.py
  * contains functions primarily dedicated to dividing a distribution network into branches, and determining which branches are the best for placing actuators.
*my_detControlMatExistence_funcs.py
  * contains the functions used to determine if a particular actuator configuration is stable.
*my_detLznRange_funcs.py
  *
*my_feeder_funcs.py
  *
*my_heatmapSetup_funcs.py
  *
*my_impedance_funcs.py 
  *
*setup_nx.py
  *
