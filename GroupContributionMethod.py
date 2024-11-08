import pandas as pd
import numpy as np
import os

class groupContribution:
    """
    Class for handling group contribution calculations, including initialization 
    of fuel-specific data and properties from Constantinou and Gani's GCM table.
    """
    
    # Paths to input directories
    GCM_PATH = os.path.dirname(os.path.abspath(__file__))
    gcmTableDir = os.path.join(GCM_PATH, 'gcmTableData')
    fuelData = os.path.join(GCM_PATH, 'fuelData')
    compDescDir = os.path.join(fuelData, 'compoundDescriptions')
    initDataDir = os.path.join(fuelData, 'initData')

    # Default GCM table name
    input_table_name = 'gcmTable.xlsx'
    input_table = os.path.join(gcmTableDir, input_table_name)

    # Class-level variables to hold fuel-specific data and GCM table properties
    fuel = ''
    num_compounds = None
    fuel_comp_desc = None
    fuel_init_data = None

    # Initial composition and functional group data for fuel
    Y_0 = None  # Initial mass fraction for fuel (num_compounds,)
    Nij = None  # Compound vs. group matrix (num_compounds, num_groups)
    
    # GCM table data from Constantinou and Gani (num_groups,)
    N_g1 = 78
    N_g2 = 43
    Tck  = None
    Pck  = None
    Vck  = None
    Tbk  = None
    Tmk  = None
    hfk  = None
    gfk  = None
    hvk  = None
    wk   = None
    Vmk  = None
    cpak = None
    cpbk = None
    cpck = None
    mwk  = None
    
    # critical properties at standard temp (num_compounds,)
    MW     = None
    Tc     = None
    Pc     = None
    Vc     = None
    Tb     = None
    Tm     = None
    Hf     = None
    Gf     = None
    Hv_stp = None
    omega  = None
    Vm_stp = None
    Cp_stp = None
    Cp_B   = None
    Cp_C   = None
    Lv_stp = None
    
    def __init__(self, fuel="", W = 1):
        """
        Initializes fuel-specific data and GCM table properties for the specified 
        fuel. Reads GCM table and fuel data from files.
        
        Parameters:
        fuel (str): Name of the fuel to initialize data for.
        W (integer): Determines if first-order only (W = 0) approximation
        """

        self.fuel = fuel
        self.fuel_comp_desc = os.path.join(self.compDescDir, f"{fuel}.xlsx")
        self.fuel_init_data = os.path.join(self.initDataDir, f"{fuel}_init.xlsx")

        # Read functional group data for fuel (num_compounds,num_groups)
        df_Nij = pd.read_excel(self.fuel_comp_desc)
        self.Nij = df_Nij.iloc[:, 1:].to_numpy()  
        self.num_compounds = self.Nij.shape[0]
        self.num_groups = self.Nij.shape[1]

        # Read initial liquid composition of fuel and normalize to get massfrac
        fuel_init_data_df = pd.read_excel(self.fuel_init_data, usecols=[1])
        self.Y_0 = fuel_init_data_df.to_numpy().flatten().astype(float)
        self.Y_0 /= np.sum(self.Y_0)

        # Make sure fuel data is consistent:
        if (self.num_groups < self.N_g1):
            raise ValueError(
                f"Insufficient fuel description:\n"
                f"The number of columns in {self.fuel_comp_desc} is less than "
                f"the required number of first-order groups (N_g1 = {self.N_g1})."
                )
        if (self.Y_0.shape[0] != self.num_compounds):
            raise ValueError(
                f"Insufficient fuel description:\n"
                f"The number of compounds in {self.fuel_comp_desc} does not "
                f"equal the number of compounds in {self.fuel_init_data}."
                )
        
        # Read and store GCM table properties
        df_gcm_properties = pd.read_excel(self.input_table)
        gcm_properties = df_gcm_properties.loc[:, 
            ~df_gcm_properties.columns.isin(['Property', 'Units'])].to_numpy()
        
        # Determine if approximations include second-order contributions
        if (W == 0 or self.num_groups < self.N_g1 + self.N_g2):
            # Only retain first-order GCM properties
            gcm_properties = gcm_properties[:,0:self.N_g1]
            self.Nij = self.Nij[:,0:self.N_g1]
        
        # Table data for functional groups (num_compounds,)
        self.Tck  = gcm_properties[0]  # critical temperature (1)
        self.Pck  = gcm_properties[1]  # critical pressure (bar)
        self.Vck  = gcm_properties[2]  # critical volume (m^3/kmol)
        self.Tbk  = gcm_properties[3]  # boiling temperature (1)
        self.Tmk  = gcm_properties[4]  # melting point temperature (1)
        self.hfk  = gcm_properties[5]  # enthalpy of formation, (kJ/mol)
        self.gfk  = gcm_properties[6]  # Gibbs energy (kJ/mol)
        self.hvk  = gcm_properties[7]  # latent heat of vaporization (kJ/mol)
        self.wk   = gcm_properties[8]  # accentric factor (1)
        self.Vmk  = gcm_properties[9]  # liquid molar volume fraction (m^3/kmol)
        self.cpak = gcm_properties[10] # specific heat values (J/mol/K)
        self.cpbk = gcm_properties[11] # specific heat values (J/mol/K)
        self.cpck = gcm_properties[12] # specific heat values (J/mol/K)
        self.mwk  = gcm_properties[13] # molecular weights (g/mol) 
        
        # --- Compute critical properties at standard temp (num_compounds,)
        # Molecular weights
        self.MW = np.matmul(self.Nij,self.mwk) # g/mol
        self.MW *= 1e-3 # Convert to kg/mol
        
        # T_c (critical temperature)
        self.Tc = 181.128*np.log(np.matmul(self.Nij,self.Tck)) # K

        # p_c (critical pressure)
        self.Pc = 1.3705 + (np.matmul(self.Nij,self.Pck) + 0.10022)**(-2) # bar
        self.Pc *= 1e5 # Convert to Pa from bar

        # V_c (critical volume)
        self.Vc = -0.00435 + ( np.matmul(self.Nij,self.Vck) ) # m^3/kmol
        self.Vc *= 1e-3 # Convert to m^3/mol

        # T_b (boiling temperature)
        self.Tb = 204.359*np.log( np.matmul(self.Nij,self.Tbk)) # K

        # T_m (melting temperature)
        self.Tm = 102.425*np.log( np.matmul(self.Nij,self.Tmk)) # K

        # H_f (enthalpy of formation)
        self.Hf = 10.835 + np.matmul(self.Nij,self.hfk) # kJ/mol
        self.Hf *= 1e3 # Convert to J/mol

        # G_f (Gibbs free energy)
        self.Gf = -14.828 + np.matmul(self.Nij,self.gfk) # kJ/mol
        self.Gf *= 1e3 # Convert to J/mol

        # H_v,stp (enthalpy of vaporization at 298 K)
        self.Hv_stp = 6.829 + (np.matmul(self.Nij,self.hvk)) # kJ/mol
        self.Hv_stp *= 1e3 # Convert to J/mol

        # omega (accentric factor)
        self.omega = 0.4085 * np.log(np.matmul(self.Nij, self.wk) + 1.1507)**(1.0 / 0.5050)

        # V_m (molar liquid volume at 298 K)
        self.Vm_stp = 0.01211 + np.matmul(self.Nij,self.Vmk)  # m^3/kmol
        self.Vm_stp *= 1e-3 # Convert to m^3/mol

        # C_p,stp (specific heat at 298 K)
        self.Cp_stp = np.matmul(self.Nij, self.cpak) - 19.7779 # J/mol/K

        # Temperature corrections for C_p
        self.Cp_B = np.matmul(self.Nij, self.cpbk)
        self.Cp_C = np.matmul(self.Nij, self.cpck)

        # L_v,stp (latent heat of vaporization at 298 K)
        self.Lv_stp = self.Hv_stp / self.MW # J/kg

        # Lennard-Jones parameters for diffusion calculations (eqt. 30 in Govindaraju)
        self.epsVec = (0.7915 + 0.1693 * self.omega) * self.Tc
        self.SigmaVec = 1e-10 * (2.3551 - 0.0874 * self.omega) * \
                        ((1e5 * self.Tc / self.Pc)**(1 / 3))
    
    def dynamic_viscosity(self, T, rho):
        """
        Calculates liquid dynamic viscosity based on droplet temperature and 
        density using Dutt's Equation (4.23) in "Viscosity of Liquids".
        
        Parameters:
        T (float): Temperature in Kelvin.
        rho_drop (float): Density drop in kg/m^3.
        
        Returns:
        np.ndarray: Dynamic viscosity in Pa*s (shape: num_compounds,).
        """
        T = self.K2C(T) # Convert to celsius
        rho *= 1e-3 # Convert from kg/m^3 to g/cm^3
        rhs = -3.0171 + (442.78 + 1.6452 * (self.Tb - 273)) / \
              (T + (239 - 0.19 * (self.Tb - 273)))
        mu = np.exp(rhs) * rho 
        mu *= 1e-3 # Convert to Pa*s
        return mu 

    def Cp(self, T):
        """
        Computes specific heat capacity at a given temperature.
        
        Parameters:
        T (float): Temperature in Kelvin.
        
        Returns:
        np.ndarray: Specific heat capacity in J/mol/K (shape: num_compounds,).
        """
        theta = (T - 298) / 700
        cp = self.Cp_stp + self.Cp_B * theta + self.Cp_C * theta**2
        return cp
    
    def Cl(self, T):
        """
        Computes liquid specific heat capacity at a given temperature.
        
        Parameters:
        T (float): Temperature in Kelvin.
        
        Returns:
        np.ndarray: Specific heat capacity in J/kg/K (shape: num_compounds,).
        """
        cp = self.Cp(T)
        return cp / self.MW  # TODO: Why divide by MW?!?!?!

    def psat(self, T):
        """
        Computes the saturated vapor pressure.
        
        Parameters:
        T (float): Temperature in Kelvin.
        
        Returns:
        np.ndarray: Saturated vapor pressure in Pa.
        """
        Tr = T / self.Tc
        f0 = 5.92714 - (6.09648 / Tr) - 1.28862 * np.log(Tr) + 0.169347 * (Tr ** 6)
        f1 = 15.2518 - (15.6875 / Tr) - 13.4721 * np.log(Tr) + 0.43577 * (Tr ** 6)

        rhs = np.exp((f0 + (self.omega * f1)))
        psat = self.Pc * rhs 
        return psat
    
    def molar_liquid_vol(self, T):
        """
        Computes molar liquid volumes with temperature correction.
        
        Parameters:
        T (float): Temperature in Kelvin.
        
        Returns:
        np.ndarray: Molar liquid volume in m^3/mol (num_compounds,)
        """
        
        Tstp = 298.
        phi = np.zeros_like(self.Tc)
        for i in range(len(self.Tc)):
            if T > self.Tc[i]:
                phi[i] = -((1 - (Tstp / self.Tc[i])) ** (2.0 / 7.0))
            else:
                phi[i] = ((1 - (T / self.Tc[i])) ** (2.0 / 7.0)) \
                        - ((1 - (Tstp / self.Tc[i])) ** (2.0 / 7.0))
        z_vec = (0.29056 - 0.08775 * self.omega)
        Vmi = self.Vm_stp * np.power(z_vec,phi)
        return Vmi
    
    def latent_heat_vaporization(self, T):
        """
        Calculates the letent heat of vaporization, adjusted for temperature.

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        numpy.ndarray: Adjusted latent heat of vaporization in J/kg (num_compounds,)
        """
        Lvi = np.zeros_like(self.Tc)
        
        Tr    = T / self.Tc
        Trnbi = self.Tb / self.Tc

        for i in range(0,self.Tc.shape[0]):
            if (T > self.Tc[i]):
                Lvi[i] = 0.
            else:
                Lvi[i] =  self.Lv_stp[i] * (( (1. - Tr[i] ) / (1. - Trnbi[i]) )**0.38)
        
        return Lvi
    
    def diffusion_coeff(self, p, Tin):
        """
        Computes the diffusion coefficient using Lennard-Jones parameters. These 
        are from Appendix B in Govindaraju Equation (33).

        Parameters:
        p (float): Pressure in Pa.
        Tin (float): Temperature in Kelvin.

        Returns:
        numpy.ndarray: Diffusion coefficient.
        """		
        # Eps and Sigma for air, replace with desired gas phase
        eps_air = 78.6 
        Sigma_air = 3.711e-10 
        sigmaAB = (Sigma_air + self.SigmaVec) / 2
        evaVec = (self.epsVec * eps_air) ** 0.5
        Tstar = Tin / evaVec
        A = 1.06036
        B = 0.15610
        C = 0.193
        D = 0.47635
        E = 1.03587
        F = 1.52996
        G = 1.76474
        H = 3.89411

        omegaD = A / (Tstar ** B) + C / np.exp(D * Tstar) + E / np.exp(F * Tstar) + G / np.exp(H * Tstar)

        # Molecular weight for air, replace with desired gas phase
        MW_air = 28.97e-3

        MWa = 2e3 * (self.MW * MW_air) / (self.MW + MW_air)

        r = (1 / p) * 1e5 * ((3.03 - (0.98 / (MWa ** 0.5))) * 1e-27 * (Tin ** 1.5)) / ((MWa ** 0.5) * (sigmaAB ** 2) * omegaD)

        return r
