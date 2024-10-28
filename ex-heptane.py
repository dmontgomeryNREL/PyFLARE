import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import GroupContributionMethod as gcm
import fxns_mixtureProperties as fxns_mix

# -----------------------------------------------------------------------------
# Calculate mixture density and viscosity of Jet A (POSF-10325) based on the 
# properties and composition of the first order groups
# -----------------------------------------------------------------------------

# Fuel for GCM and data for validation
fuel_name = 'heptane'
fuel_data = 'heptane-NIST.csv'

# system temperature (C)
T_min = -60.0
T_max = 60.0

# system pressure
p = 1e5 # (Pa) or 1 bar

# droplet specs
drop = {} 
drop['d_0'] = 100*1e-6 	# initial droplet diameter (m)
drop['r_0'] = drop['d_0']/2.0 # initial droplet radius (m)

# Get the fuel properties based on the GCM
fuel = gcm.groupContribution(fuel_name)

# initial liquid mass fractions
Y_li = fuel.Y_l0

# number of compounds
num_compounds = fuel.num_compounds

# initial droplet radius (m)
drop['r_0'] = drop['d_0']/2.0

# Vector for temperature (convert from C to K)
T = fxns_mix.C2K(np.linspace(T_min,T_max,100))

# Vectors for density and viscosity
rho = np.zeros_like(T)
nu = np.zeros_like(T)

# Colors for plotting
colors = ['','tab:purple','tab:blue','tab:red','tab:orange','tab:green']

for i in range(0,len(T)): 
    # Mixture density (returns rho in kg/m^3)
    rho[i] = fxns_mix.calc_mixture_density(fuel,T[i])
    # Convert density to CGS (g/cm^3)
    rho[i] *= 1.0e-03 

    # Mixture viscosity (returns nu in m^2/s)
    nu[i] = fxns_mix.calc_mixture_viscosity(fuel,T[i],drop['r_0'],1)
    # Convert viscosity to mm^2/s
    nu[i] *= 1.0e+06

# Get experimental or NIST data for validation
dataPath = os.path.join(fuel.fuelData,"propertiesData")
data = pd.read_csv(os.path.join(dataPath,fuel_data),skiprows=[1])

# Convert data to kinematic viscosity in mm^2/s
data.Viscosity = data.Viscosity / (data.Density * 1000) * 1e6 

# Plot parameters
fsize = 26
ticksize = 24
line_thickness = 5
marker_size = 150

# Plot mixture viscosity vs. Temp
plt.figure(figsize=(12.3,8))
plt.plot(fxns_mix.K2C(T), nu, '-k',label='Model Prediction', linewidth=line_thickness)
plt.scatter(data.Temperature, data.Viscosity, 
            label="NIST Data", facecolors='tab:orange', s=marker_size)
plt.xlabel('Temperature ($^\circ$C)', fontsize=fsize)
plt.ylabel('Viscosity (mm$^2$/s)', fontsize=fsize)
plt.xlim([T_min, T_max])
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.legend(fontsize=fsize)

# Plot mixture density vs. Temp
plt.figure(figsize=(12.3,8))
plt.plot(fxns_mix.K2C(T), rho, '-k',label='Model Prediction', linewidth=line_thickness)
plt.scatter(data.Temperature, data.Density, 
            label="NIST Data", facecolors='tab:orange', s=marker_size)
plt.xlabel('Temperature ($^\circ$C)', fontsize=fsize)
plt.ylabel('Density (g/cm$^3$)', fontsize=fsize)
plt.xlim([T_min, T_max])
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.legend(fontsize=fsize)
plt.show()
