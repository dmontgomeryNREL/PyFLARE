import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import GroupContributionMethod as gcm

# -----------------------------------------------------------------------------
# Calculate mixture properties from the group contribution properties
# -----------------------------------------------------------------------------

# Fuel for GCM and data for validation (see fuelData/propertiesData for options)
# fuel_name = 'decane', 'dodecane', 'heptane', 'posf10264', 'posf10289', 'posf10325'
fuel_name = 'decane'

# droplet specs
drop = {} 
drop['d_0'] = 100*1e-6 	# initial droplet diameter (m), note: size doesn't matter
drop['r_0'] = drop['d_0']/2.0 # initial droplet radius (m)
drop['V_0'] = 4.0 / 3.0 * np.pi * drop['r_0'] ** 3 # initial droplet volume
def drop_mass(fuel,Yi,T):
    return drop['V_0'] / (fuel.molar_liquid_vol(T) @ Yi) * Yi * fuel.MW # (kg)

# Get the fuel properties based on the GCM
fuel = gcm.groupContribution(fuel_name)

# initial liquid mass fractions
Y_li = fuel.Y_0

# Get data for validation
fuel_to_data = {
    'heptane': ('heptane-NIST.csv', 'NIST Data'),
    'decane': ('decane-NIST.csv', 'NIST Data'),
    'dodecane': ('dodecane-NIST.csv', 'NIST Data'),
    'posf10264': ('posf10264.csv', 'Edwards Data'),
    'posf10289': ('posf10289.csv', 'Edwards Data'),
    'posf10325': ('posf10325.csv', 'Edwards Data'),
    'a-2_posf10325_Ed': ('a-2_posf10325_Ed.csv', 'Edwards Data')
}
data_file, data_source = fuel_to_data.get(fuel_name)
dataPath = os.path.join(fuel.fuelDataDir,"propertiesData")
data = pd.read_csv(os.path.join(dataPath,data_file),skiprows=[1])

# Seperate properties and associated temperatures from data
T_rho_data = data.Temperature[data.Density.notna()]
rho_data = data.Density.dropna()
T_pv_data = data.Temperature[data.VaporPressure.notna()]
pv_data = data.VaporPressure.dropna()
if "posf" not in fuel_name.lower():
    T_nu_data = data.Temperature[data.Viscosity.notna()]
    nu_data = data.Viscosity.dropna()

# Vectors for temperature (convert from C to K)
T_rho = gcm.C2K(np.linspace(min(T_rho_data),max(T_rho_data),100))
T_pv = gcm.C2K(np.linspace(min(T_pv_data),max(T_pv_data),100))
if "posf" not in fuel_name.lower():
    T_nu = gcm.C2K(np.linspace(min(T_nu_data),max(T_nu_data),100))

# Vectors for density, viscosity and vapor pressure
rho = np.zeros_like(T_rho)
pv = np.zeros_like(T_pv)
if "posf" not in fuel_name.lower():
    nu = np.zeros_like(T_nu)

for i in range(0,len(T_rho)): 
    # Mixture density (returns rho in kg/m^3)
    rho[i] = fuel.mixture_density(fuel.Y_0,T_rho[i])
    # Convert density to CGS (g/cm^3)
    rho[i] *= 1.0e-03 

for i in range(0,len(T_pv)): 
    # Mass of the droplet at current temp
    mass = drop_mass(fuel, Y_li, T_pv[i])
    # Mixture vapor pressure (returns pv in Pa)
    pv[i] = fuel.mixture_vapor_pressure(mass,T_pv[i])
    # Convert vapor pressure to kPa
    pv[i] *= 1.0e-03

if "posf" not in fuel_name.lower():
    for i in range(0,len(T_nu)): 
        # Mass of the droplet at current temp
        mass = drop_mass(fuel, Y_li, T_nu[i])
        # Mixture viscosity (returns nu in m^2/s)
        nu[i] = fuel.mixture_kinematic_viscosity(mass, T_nu[i])
        # Convert viscosity to mm^2/s
        nu[i] *= 1.0e+06

# Plotting parameters
fsize = 24
ticksize = 22
line_thickness = 4
marker_size = 75

# Plot mixture density vs. Temp
plt.figure(figsize=(10,7.2))
plt.plot(gcm.K2C(T_rho), rho, '-',label='Model Prediction', linewidth=line_thickness)
plt.scatter(T_rho_data, rho_data, 
            label=data_source, facecolors='black', s=marker_size)
plt.xlabel('Temperature ($^\circ$C)', fontsize=fsize)
plt.ylabel('Density (g/cm$^3$)', fontsize=fsize)
plt.xlim([min(T_rho_data), max(T_rho_data)])
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.legend(fontsize=fsize)

# Plot mixture vapor pressure vs. Temp
plt.figure(figsize=(10,7.2))
plt.plot(gcm.K2C(T_pv), pv, '-',label='Model Prediction', linewidth=line_thickness)
plt.scatter(T_pv_data, pv_data, 
            label=data_source, facecolors='black', s=marker_size)
plt.xlabel('Temperature ($^\circ$C)', fontsize=fsize)
plt.ylabel('Vapor Pressure (kPa)', fontsize=fsize)
plt.xlim([min(T_pv_data), max(T_pv_data)])
plt.ylim([0, max(pv_data)])
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.legend(fontsize=fsize)


# Plot mixture viscosity vs. Temp
if "posf" not in fuel_name.lower():
    plt.figure(figsize=(10,7.2))
    plt.plot(gcm.K2C(T_nu), nu, '-',label='Model Prediction', linewidth=line_thickness)
    plt.scatter(T_nu_data, nu_data, 
                label=data_source, facecolors='black', s=marker_size)
    plt.xlabel('Temperature (C)', fontsize=fsize)
    plt.ylabel('Viscosity (mm$^2$/s)', fontsize=fsize)
    plt.xlim([min(T_nu_data), max(T_nu_data)])
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
    plt.legend(fontsize=fsize)

plt.show()