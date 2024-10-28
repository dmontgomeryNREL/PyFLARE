import numpy as np
import cantera as ct
import fxns_singleDropletProperties as dropletFxns

def C2K(T):
    """
    Convert temperature from Celsius to Kelvin.

    Parameters:
    T (float or np.array): Temperature in Celsius.

    Returns:
    float or np.array: Temperature in Kelvin.
    """
    return T + 273.15

def K2C(T):
    """
    Convert temperature from Kelvin to Celsius.

    Parameters:
    T (float or np.array): Temperature in Kelvin.

    Returns:
    float or np.array: Temperature in Celsius.
    """
    return T - 273.15

def calc_mixture_density(fuel, T):
    """
    Calculate the mixture density of the fuel at a given temperature.

    Parameters:
    fuel: Object of the fuel class containing properties.
    T (float): Temperature in Kelvin.

    Returns:
    float: Mixture density in kg/m^3.
    """
    y = fuel.Y_l0  # Mass fractions of the fuel components
    MW = fuel.MW   # Molecular weights of the fuel components (kg/mol)
    Vm = fuel.molar_liquid_vol(T)  # Molar volume of the fuel (m^3/mol)

    # Calculate density (kg/m^3)
    density_drop = np.dot(MW, y) / np.dot(Vm, y)

    return density_drop 

def calc_components_viscosity(fuel, T):
    """
    Calculate the viscosity of individual components in the fuel at a given 
    temperature using Dutt's equation (4.23) in Viscosity of Liquids.

    Parameters:
    fuel: Object of the fuel class containing properties.
    T (float): Temperature in Kelvin.

    Returns:
    np.array: Viscosity of each component in m^2/s.
    """
    T_celsius = K2C(T)  # Convert temperature to Celsius

    # RHS of Dutt's equation (4.23) 
    rhs = -3.0171 + (442.78 + 1.6452 * K2C(fuel.Tb)) / (T_celsius + 239 - 0.19 * K2C(fuel.Tb))

    nu = np.exp(rhs)  # Viscosity in mm^2/s 

    # Convert to SI (m^2/s)
    nu *= 1e-6

    return nu 

def calc_mixture_viscosity(fuel, T, radius, model=1):
    """
    Calculate the viscosity of the fuel mixture at a given temperature and droplet radius.

    Parameters:
    fuel: Object of the fuel class containing properties.
    T (float): Temperature in Kelvin.
    radius (float): Radius of the droplet (m).
    model (int, optional): Mixing model to use (1 for Kendall-Monroe, others for Arrhenius). Default is 1.

    Returns:
    float: Mixture viscosity in mm^2/s.
    """
    nu_i = calc_components_viscosity(fuel, T)  # Viscosities of individual components

    MW = fuel.MW  # Molecular weights of the fuel components (kg/mol)
    Y_li = fuel.Y_l0  # Mass fractions of the fuel components
    massVec = dropletFxns.massVector(fuel, radius, Y_li, T)  # Mass vector of each compound

    # Calculate group mole fractions for each species
    x_l = dropletFxns.moleFracVec(massVec, MW)
    
    if model == 1:
        # Kendall-Monroe mixing model
        nu = np.sum(x_l * (nu_i ** (1.0 / 3.0))) ** (3.0)
    else:
        # Arrhenius mixing model
        nu = np.exp(np.sum(x_l * np.log(nu_i)))
    
    return nu
