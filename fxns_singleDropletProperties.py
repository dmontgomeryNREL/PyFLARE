import numpy as np

def moleFracVec(massVector, MW):
    """
    Calculate the mole fractions from a mass vector.
    
    Parameters:
    massVector (np.ndarray): Mass of each compound (shape: num_compounds,).
    MW (np.ndarray): Molecular weight of each compound (shape: num_compounds,).

    Returns:
    Xi np.ndarray: Molar fractions of the compounds (shape: num_compounds,).
    """
    # Calculate the number of moles for each compound
    mole_num = np.divide(massVector, MW)
    
    # Normalize to get group mole fractions
    total_moles = np.sum(mole_num)
    if total_moles != 0:
        Xi = np.divide(mole_num, total_moles)
    
    return Xi

def droplet_radius_from_mass(fuel, massVector, T):
    """
    Compute the droplet radius from its mass and temperature.

    Parameters:
    fuel (object): An instance of the groupContribution class.
    massVector (np.ndarray): Mass of each compound (shape: num_compounds,).
    T (float): Temperature in Kelvin.

    Returns:
    float: Radius of the droplet in meters
    """
    rho = np.divide(fuel.MW, fuel.molar_liquid_vol(T))
    vol = np.sum(np.divide(massVector, rho))
        
    # Calculate and return the radius
    return (3 * vol / (4 * np.pi))**(1/3) if vol > 0 else 0.0

def drop_volume(radius):
    """
    Compute the volume of a droplet assuming it is spherical.

    Parameters:
    radius (float): Radius of the droplet in meters.

    Returns:
    float: Volume of the droplet in m^3.
    """
    return 4.0 / 3.0 * np.pi * radius ** 3

def droplet_density(massVector, radius):
    """
    Compute the density of the droplet from its mass and radius.

    Parameters:
    massVector (np.ndarray): Mass of each compound (shape: num_compounds,).
    radius (float): Radius of the droplet in meters.

    Returns:
    float: Density of the droplet in kg/m^3.
    """
    return np.sum(massVector) / drop_volume(radius)  # kg/m^3

def massVector(fuel, radius, Yi, T):
    """
    Calculate the mass of each compound in the fuel (kg).

    Parameters:
    fuel (object): An instance of the groupContribution class.
    radius (float): Radius of the droplet in meters.
    Yi (np.ndarray): Mole fractions of each compound (shape: num_compounds,).
    T (float): Temperature in Kelvin.

    Returns:
    np.ndarray: Mass of each compound in kg (shape: num_compounds,).
    """
    MW = fuel.MW
    volume = drop_volume(radius)
    
    # Calculate mass vector based on volume and mole fractions
    if volume > 0:
        return np.multiply((volume / np.dot(fuel.molar_liquid_vol(T), Yi)) * Yi, MW)
    else:
        return np.zeros_like(Yi)
