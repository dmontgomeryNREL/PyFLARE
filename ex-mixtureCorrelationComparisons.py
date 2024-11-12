import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import GroupContributionMethod as gcm

# -----------------------------------------------------------------------------
# Calculate mixture properties from the group contribution properties
# -----------------------------------------------------------------------------

# Fuel for GCM and data for validation (see mixtureData/propertiesData for options)
# fuel_name = 'decane', 'dodecane', 'heptane', 'posf10325'
fuel_name = 'posf10325'

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
    'posf10325': ('posf10325-GE-JetA.csv', 'NREL Data')
}
data_file, data_source = fuel_to_data.get(fuel_name)
dataPath = os.path.join(fuel.mixtureData,"propertiesData")
data = pd.read_csv(os.path.join(dataPath,data_file),skiprows=[1])

# Seperate properties and associated temperatures from data
T_nu_data = data.Temperature[data.Viscosity.notna()]
nu_data = data.Viscosity.dropna()
T_pv_data = data.Temperature[data.VaporPressure.notna()]
pv_data = data.VaporPressure.dropna()

# Vectors for temperature (convert from C to K)
T_nu = gcm.C2K(np.linspace(min(T_nu_data),max(T_nu_data),100))
T_pv = gcm.C2K(np.linspace(min(T_pv_data),max(T_pv_data),100))

# Vectors for viscosity and vapor pressure
nu_KM = np.zeros_like(T_nu)
nu_Arr = np.zeros_like(T_nu)
pv_LK = np.zeros_like(T_pv)
pv_AW = np.zeros_like(T_pv)

for i in range(0,len(T_nu)): 
    # Mass of the droplet at current temp
    mass = drop_mass(fuel, Y_li, T_nu[i])
    # Mixture viscosity (returns nu in m^2/s)
    nu_KM[i] = fuel.mixture_kinematic_viscosity(mass, T_nu[i])
    nu_Arr[i] = fuel.mixture_kinematic_viscosity(mass, T_nu[i],'Arrhenius')
    # Convert viscosity to mm^2/s
    nu_KM[i] *= 1.0e+06
    nu_Arr[i] *= 1.0e+06

for i in range(0,len(T_pv)): 
    # Mass of the droplet at current temp
    mass = drop_mass(fuel, Y_li, T_pv[i])
    # Mixture vapor pressure (returns pv in Pa)
    pv_LK[i] = fuel.mixture_vapor_pressure(mass,T_pv[i])
    pv_AW[i] = fuel.mixture_vapor_pressure(mass,T_pv[i],'Ambrose-Walton')
    # Convert vapor pressure to kPa
    pv_LK[i] *= 1.0e-03
    pv_AW[i] *= 1.0e-03

# Plotting parameters
fsize = 26
ticksize = 24
line_thickness = 5
marker_size = 150

# Plot mixture viscosity vs. Temp
plt.figure(figsize=(12.3,8))
plt.plot(gcm.K2C(T_nu), nu_KM,'-',color='tab:blue',label='Kendal-Monroe', linewidth=line_thickness)
plt.plot(gcm.K2C(T_nu), nu_Arr,'-',color='tab:green',label='Arrhenius', linewidth=line_thickness)
plt.scatter(T_nu_data, nu_data, 
            label=data_source, facecolors='tab:orange', s=marker_size)
plt.xlabel('Temperature ($^\circ$C)', fontsize=fsize)
plt.ylabel('Viscosity (mm$^2$/s)', fontsize=fsize)
plt.xlim([min(T_nu_data), max(T_nu_data)])
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.legend(fontsize=fsize)

# Plot mixture vapor pressure vs. Temp
plt.figure(figsize=(12.3,8))
plt.plot(gcm.K2C(T_pv), pv_LK, '-',color='tab:blue',label='Lee-Kesler', linewidth=line_thickness)
plt.plot(gcm.K2C(T_pv), pv_AW, '-',color='tab:green',label='Ambrose-Walton', linewidth=line_thickness)
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
