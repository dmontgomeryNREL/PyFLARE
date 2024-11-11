import numpy as np
import fxns_singleDropletProperties as dropletFxns

def mixture_density(fuel, T, Yi):
    """
    Calculate the mixture density of the fuel at a given temperature.

    Parameters:
    fuel: Object of the groupContribution class containing properties.
    T (float): Temperature in Kelvin.
    Yi (np.ndarray): Mass fractions of each compound (shape: num_compounds,).

    Returns:
    float: Mixture density in kg/m^3.
    """
    MW = fuel.MW   # Molecular weights of the fuel components (kg/mol)
    Vm = fuel.molar_liquid_vol(T)  # Molar volume of the fuel (m^3/mol)

    # Calculate density (kg/m^3)
    density_drop = (MW @ Yi) / (Vm @ Yi)

    return density_drop 

def mixture_viscosity(fuel, T, Yi, radius, correlation='Kendall-Monroe'):
    """
    Calculate the viscosity of the fuel mixture at a given temperature and droplet radius.

    Parameters:
    fuel: Object of the groupContribution class containing properties.
    T (float): Temperature in Kelvin.
    Yi (np.ndarray): Mass fractions of each compound (shape: num_compounds,).
    radius (float): Radius of the droplet (m).
    correlation (str, optional): Mixing model "Kendall-Monroe", "Arrhenius". 

    Returns:
    float: Mixture viscosity in mm^2/s.
    """
    nu_i = fuel.viscosity_kinematic(T)  # Viscosities of individual components

    MW = fuel.MW  # Molecular weights of the fuel components (kg/mol)
    
    massVec = dropletFxns.massVector(fuel, radius, Yi, T)  # Mass vector of each compound

    # Calculate group mole fractions for each species
    Xi = dropletFxns.moleFracVec(massVec, MW)
    
    if (correlation.casefold() == 'Arrhenius'.casefold()):
        # Arrhenius mixing correlation
        nu = np.exp(np.sum(Xi * np.log(nu_i)))
    else:
        # Default: Kendall-Monroe mixing correlation
        nu = np.sum(Xi * (nu_i ** (1.0 / 3.0))) ** (3.0)
        
    
    return nu

def mixture_vapor_pressure(fuel, T, Yli, radius, correlation = 'Lee-Kesler'):
    """
    Calculate the vapor pressure the fuel mixture.

    Parameters:
    fuel: Object of the groupContribution class containing properties.
    T (float): Temperature in Kelvin.
    Yli (np.ndarray): Mass fractions of each compound in liquid phase (shape: num_compounds,).
    radius (float): Radius of the droplet (m).
    correlation (str, optional): "Ambrose-Walton" or "Lee-Kesler". 

    Returns:
    float: vapor pressure in Pa.
    """
    # molecular weights of fuel (kg/mol)
    MW = fuel.MW
    
    # Mass from mass fraction and droplet radius (kg)
    massVec = dropletFxns.massVector(fuel, radius, Yli, T) 
    
    # Group mole fraction for each compound
    Xi = dropletFxns.moleFracVec(massVec,MW)

    # Saturated vapor pressure for each compound (Pa)
    p_sati = fuel.psat(T,correlation)

    # Mixture vapor pressure via Raoult's law
    p_v = p_sati @ Xi
    
    return p_v