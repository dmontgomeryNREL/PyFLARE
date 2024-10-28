# Group Contribution Method
This code utilizes the tables and functions proposed by [Constantinou and Gani (1994)](https://doi.org/10.1002/aic.690401011) and [Constantinou, Gani and O'Connel (1995)](https://doi.org/10.1016/0378-3812(94)02593-P), with additional physical properties discussed in [Govindaraju & Ihme (2016)](https://doi.org/10.1016/j.ijheatmasstransfer.2016.06.079), to calculate the thermodynamic properties of jet fuels.  The code is based on Pavan Bharadwaj's [Matlab implementation](https://github.com/gpavanb-old/GroupContribution).

## Python Environment
The following conda environment is required to run this code:
~~~
conda create --name gcm-env --channel cantera cantera ipython matplotlib pandas scipy openpyxl
~~~

## Running the code
Two example cases are included: `ex-posf10325.py` and `ex-heptane.py`.  Run these scripts to see how the group contribution method can be used to calculate physical properties of the fuel mixture solely from GC x GC data of the fuel.