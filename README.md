# Group Contribution Method
This code utilizes the tables and functions proposed by [Constantinou and Gani (1994)](https://doi.org/10.1002/aic.690401011) and [Constantinou, Gani and O'Connel (1995)](https://doi.org/10.1016/0378-3812(94)02593-P), with additional physical properties discussed in [Govindaraju & Ihme (2016)](https://doi.org/10.1016/j.ijheatmasstransfer.2016.06.079).  The code is based on Pavan B. Govindaraju's [Matlab implementation](https://github.com/gpavanb-old/GroupContribution) of the GCM, and has been expanded to include additional thermodynamic properties and mixture properties.

## Python Environment
The following conda environment is required to run this code:
~~~
conda create --name gcm-env --channel cantera cantera ipython matplotlib pandas scipy openpyxl
~~~

## Running the code
This repository includes an example `ex-mixtureProperties.py` that calculates a given mixture's density, viscosity and vapor pressure from GC x GC data of the mixture.  The results are plotted against data from NIST and NREL.

# Contributing
New contributions are always welcome.  If you have an idea for a new feature follow these steps:
1. Fork the main repository
2. Create a `newFeature` branch that contains your changes
3. Update the sphinx documentation in `newFeature`
4. Open a Pull Request (PR) from `newFeature` on your fork to branch `main` GCM-Python repository.

## Sphinx Documentation
This repository uses [Sphinx](https://www.sphinx-doc.org/en/master/usage/quickstart.html) to generate the documentation.  This requires the following Conda environment:
~~~
conda create --name sphinx-env
conda activate sphinx-env
conda install python
conda install anaconda::sphinx
conda install conda-forge::sphinx_rtd_theme 
~~~

To view the documentation locally, build the html using the following: 
~~~
cd GCM-Python/docs/sphinx
sphinx-build -M html . _build/
~~~
You should now be able to view the html by opening `GCM-Python/docs/sphinx_build/html/indx.html` in a web browser. 
