import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import GroupContributionMethod as gcm
import fxns_mixtureProperties as fxns_mix

# -----------------------------------------------------------------------------
# Calculate mixture properties from the group contribution properties
# -----------------------------------------------------------------------------

# Fuel for GCM and data for validation (see fuelData/propertiesData for options)
# fuel_name = 'decane', 'dodecane', 'heptane', , 'posf10325'
fuel_name = 'dodecane'

# droplet specs
drop = {} 
drop['d_0'] = 100*1e-6 	# initial droplet diameter (m), note: size doesn't matter
drop['r_0'] = drop['d_0']/2.0 # initial droplet radius (m)

# Get the fuel properties based on the GCM
fuel = gcm.groupContribution(fuel_name)

# initial liquid mass fractions
Y_li = fuel.Y_0

# Get data for validation
fuel_to_data = {
    'heptane': ('heptane-NIST.csv', 'NIST Data'),
    'decane': ('decane-NIST.csv', 'NIST Data'),
    'dodecane': ('dodecane-NIST.csv', 'NIST Data'),
    'posf10325': ('posf10325-GE-JetA.csv', 'NREL Data')
}
data_file, data_source = fuel_to_data.get(fuel_name)
dataPath = os.path.join(fuel.fuelData,"propertiesData")
data = pd.read_csv(os.path.join(dataPath,data_file),skiprows=[1])
T_nu_data = data.Temperature[data.Viscosity.notna()]
nu_data = data.Viscosity.dropna()
T_rho_data = data.Temperature[data.Density.notna()]
rho_data = data.Density.dropna()
T_pv_data = data.Temperature[data.VaporPressure.notna()]
pv_data = data.VaporPressure.dropna()

# Vectors for temperature (convert from C to K)
T_rho = fxns_mix.C2K(np.linspace(min(T_rho_data),max(T_rho_data),100))
T_nu = fxns_mix.C2K(np.linspace(min(T_nu_data),max(T_nu_data),100))
T_pv = fxns_mix.C2K(np.linspace(min(T_pv_data),max(T_pv_data),100))

# Vectors for density, viscosity and vapor pressure
rho = np.zeros_like(T_rho)
nu = np.zeros_like(T_nu)
pv = np.zeros_like(T_pv)

# Colors for plotting
colors = ['','tab:purple','tab:blue','tab:red','tab:orange','tab:green']

for i in range(0,len(T_rho)): 
    # Mixture density (returns rho in kg/m^3)
    rho[i] = fxns_mix.calc_mixture_density(fuel,T_rho[i],fuel.Y_0)
    # Convert density to CGS (g/cm^3)
    rho[i] *= 1.0e-03 

for i in range(0,len(T_nu)): 
    # Mixture viscosity (returns nu in m^2/s)
    nu[i] = fxns_mix.calc_mixture_viscosity(fuel,T_nu[i],fuel.Y_0,drop['r_0'])
    # Convert viscosity to mm^2/s
    nu[i] *= 1.0e+06

for i in range(0,len(T_pv)): 
    # Mixture vapor pressure (returns pv in Pa)
    pv[i] = fxns_mix.calc_vapor_pressure(fuel,T_pv[i],fuel.Y_0,drop['r_0'])
    # Convert vapor pressure to kPa
    pv[i] *= 1.0e-03

# Plot parameters
fsize = 26
ticksize = 24
line_thickness = 5
marker_size = 150

# Plot mixture viscosity vs. Temp
plt.figure(figsize=(12.3,8))
plt.plot(fxns_mix.K2C(T_nu), nu, '-k',label='Model Prediction', linewidth=line_thickness)
plt.scatter(T_nu_data, nu_data, 
            label=data_source, facecolors='tab:orange', s=marker_size)
plt.xlabel('Temperature ($^\circ$C)', fontsize=fsize)
plt.ylabel('Viscosity (mm$^2$/s)', fontsize=fsize)
plt.xlim([min(T_nu_data), max(T_nu_data)])
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.legend(fontsize=fsize)

# Plot mixture density vs. Temp
plt.figure(figsize=(12.3,8))
plt.plot(fxns_mix.K2C(T_rho), rho, '-k',label='Model Prediction', linewidth=line_thickness)
plt.scatter(T_rho_data, rho_data, 
            label=data_source, facecolors='tab:orange', s=marker_size)
plt.xlabel('Temperature ($^\circ$C)', fontsize=fsize)
plt.ylabel('Density (g/cm$^3$)', fontsize=fsize)
plt.xlim([min(T_rho_data), max(T_rho_data)])
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.legend(fontsize=fsize)

# Plot mixture vapor pressure vs. Temp
plt.figure(figsize=(12.3,8))
plt.plot(fxns_mix.K2C(T_pv), pv, '-k',label='Model Prediction', linewidth=line_thickness)
plt.scatter(T_pv_data, pv_data, 
            label=data_source, facecolors='tab:orange', s=marker_size)
plt.xlabel('Temperature ($^\circ$C)', fontsize=fsize)
plt.ylabel('Vapor Pressure (kPa)', fontsize=fsize)
plt.xlim([min(T_pv_data), max(T_pv_data)])
plt.ylim([0, max(pv_data)])
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.legend(fontsize=fsize)
plt.show()
