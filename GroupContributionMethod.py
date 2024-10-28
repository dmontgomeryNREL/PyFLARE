import pandas as pd
import numpy as np
import os

class groupContribution:
    """
    Class for handling group contribution calculations, including initialization 
    of fuel-specific data and properties from Constantinou and Gani's GCM table.
    """
    
    # Paths to input directories
    GCM_PATH = os.getcwd()
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
    
    # GCM table data from Constantinou and Gani (num_groups,)
    Tc1k   = None
    Pc1k   = None
    Vc1k   = None
    w1k    = None
    hv1k   = None
    Vliq1k = None
    cpak   = None
    cpbk   = None
    cpck   = None
    mwk    = None
    Tb1k   = None
    
    # Initial composition and functional group data for fuel
    Y_l0 = None  # Initial mass fraction for fuel (num_compounds,)
    Nij = None  # Compound vs. group matrix (num_compounds, num_groups)
    
    # critical properties (num_compounds,)
    MW     = None
    Tc     = None
    Pc     = None
    Vc     = None
    omega  = None
    Tb     = None
    Lv_stp = None
    Vm_stp = None
    
    def __init__(self, fuel=""):
        """
        Initializes fuel-specific data and GCM table properties for the specified 
        fuel. Reads GCM table and fuel data from files.
        
        Parameters:
        fuel (str): Name of the fuel to initialize data for.
        """
        self.fuel = fuel
        self.fuel_comp_desc = os.path.join(self.compDescDir, f"{fuel}.xlsx")
        self.fuel_init_data = os.path.join(self.initDataDir, f"{fuel}_init.xlsx")
        
        # Read and store GCM table properties
        df_gcm_properties = pd.read_excel(self.input_table)
        gcm_properties = df_gcm_properties.loc[:, 
            ~df_gcm_properties.columns.isin(['Property', 'Units'])].to_numpy()
        
        # Table data for functional groups (num_compounds,)
        self.Tc1k   = gcm_properties[0]  # critical temperature
        self.Pc1k   = gcm_properties[1]  # critical pressure
        self.Vc1k   = gcm_properties[2]  # critical volume
        self.Tb1k   = gcm_properties[3]  # boiling temperature
        self.Tm1k   = gcm_properties[4]  # melting point temperature
        self.hf1k   = gcm_properties[5]  # enthalpy of formation, kJ/mol
        self.gfc    = gcm_properties[6]  # Gibbs energy
        self.hv1k   = gcm_properties[7]  # latent heat of vaporization, kJ/mol
        self.w1k    = gcm_properties[8]  # accentric factor
        self.Vliq1k = gcm_properties[9]  # liquid molar volume fraction
        self.cpak   = gcm_properties[10] # specific heat values
        self.cpbk   = gcm_properties[11] # specific heat values
        self.cpck   = gcm_properties[12] # specific heat values
        self.mwk    = gcm_properties[13] # molecular weights, g/mol

        # Read functional group data for fuel (num_compounds,num_groups)
        df_Nij = pd.read_excel(self.fuel_comp_desc)
        self.Nij = df_Nij.iloc[:, 1:].to_numpy()  
        self.num_compounds = self.Nij.shape[0]

        # Read initial liquid composition of fuel and normalize to get massfrac
        fuel_init_data_df = pd.read_excel(self.fuel_init_data, usecols=[1])
        self.Y_l0 = fuel_init_data_df.to_numpy().flatten().astype(float)
        self.Y_l0 /= np.sum(self.Y_l0) 
        
        # --- compute critical properties/stp properties (num_compounds,)
        # Molecular weights
        self.MW = np.matmul(self.Nij,self.mwk) # g/mol
        self.MW *= 1e-3 # Convert to kg/mol
        
        # T_c (critical temperature)
        self.Tc = 181.128*np.log(np.matmul(self.Nij,self.Tc1k)) # K

        # p_c (critical pressure)
        self.Pc = (1.3705 + (np.matmul(self.Nij,self.Pc1k) + 0.10022)**(-2)) # bar
        self.Pc *= 1e5 # Convert to Pa from bar

        # omega (accentric factor)
        self.omega = 0.4085 * (np.log(np.matmul(self.Nij, self.w1k) + 1.1507)**(1 / 0.5050))

        # T_b (boiling temperature)
        self.Tb = 204.359*np.log( np.matmul(self.Nij,self.Tb1k)) # K

        # V_c (critical volume)
        self.Vc = 1e-3*(-0.00435 + ( np.matmul(self.Nij,self.Vc1k) )) # m^3/mol

        # Standard heat of vaporization
        self.Lv_stp = (6.829 + (np.matmul(self.Nij,self.hv1k)) ) * 1e3 # J/mol
        # Not sure why this is scaled by 1/molecular weight?
        self.Lv_stp = np.divide(self.Lv_stp, self.MW)

        # Molar liquid volume at standard temperature
        self.Vm_stp = 1e-3 * (0.01211 + ( np.matmul(self.Nij,self.Vliq1k) )) # m^3/mol
        
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
            
    def Cp_stp(self, T):
        """
        Computes specific heat capacity at a given temperature.
        
        Parameters:
        T (float): Temperature in Kelvin.
        
        Returns:
        np.ndarray: Specific heat capacity in J/kg/K (shape: num_compounds,).
        """
        theta = (T - 298) / 700
        cl = (np.matmul(self.Nij, self.cpak) - 19.7779) + \
             (np.matmul(self.Nij, self.cpbk) + 22.5981) * theta + \
             (np.matmul(self.Nij, self.cpck) - 10.7983) * (theta**2)
        return cl / self.MW

    def psat(self, T):
        """
        Computes the saturated vapor pressure.
        
        Parameters:
        T (float): Temperature in Kelvin.
        
        Returns:
        np.ndarray: Saturated vapor pressure in Pa.
        """
        Tr = np.divide(T,self.Tc)
        f0 = 5.92714 - np.divide(6.09648, Tr) - 1.28862 * np.log(Tr) \
            + 0.169347 * np.power(Tr, 6)
        f1 = 15.2518 - np.divide(15.6875, Tr) - 13.4721 * np.log(Tr) \
            + 0.43577  * np.power(Tr, 6)
        rhs = np.exp((f0 + (self.omega * f1)))
        psat = np.multiply(self.Pc, rhs) 
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
        for i in range(0,self.Tc.shape[0]):
            if (T > self.Tc[i]):
                phi[i] = - ((1 - (np.divide(Tstp,self.Tc[i]) ) )**(2.0/7.0))
            else:
                phi[i] = ((1 - (np.divide(T,self.Tc[i]) ) )**(2.0/7.0)) \
                            - ((1 - (np.divide(Tstp,self.Tc[i]) ) )**(2.0/7.0))
        z_vec = (0.29056 - 0.08775 * self.omega)
        Vmi = np.multiply(self.Vm_stp, (np.power(z_vec,phi)) )
        return Vmi
    
    def enthalpy_vaporization(self, T):
        """
        Calculates the enthalpy of vaporization, adjusted for temperature.

        Parameters:
        T (float): Temperature in Kelvin.

        Returns:
        numpy.ndarray: Adjusted enthalpy of vaporization in J/mol (num_compounds,)
        """
        Lvi = np.zeros_like(self.Tc)
        
        Tr    = np.divide(T,self.Tc)
        Trnbi = np.divide(self.Tb,self.Tc)

        for i in range(0,self.Tc.shape[0]):
            if (T > self.Tc[i]):
                Lvi[i] = 0.
            else:
                Lvi[i] =  np.multiply(self.Lv_stp[i], (np.divide( (1. - Tr[i] ), (1. - Trnbi[i]) )**0.38))
        
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
        sigmaAB = (Sigma_air + self.SigmaVec)/2
        evaVec = np.power((self.epsVec * eps_air), 0.5)
        Tstar  = np.divide(Tin, evaVec)
        A = 1.06036
        B = 0.15610
        C = 0.193
        D = 0.47635
        E = 1.03587
        F = 1.52996
        G = 1.76474
        H = 3.89411
        
        omegaD = np.divide(A,(np.power(Tstar,B))) + np.divide(C,np.exp(D * Tstar)) + np.divide(E,np.exp(F * Tstar)) + np.divide(G, np.exp(H * Tstar))
        
        # Molecular weight for air, replace with desired gas phase
        MW_air = 28.97e-3

        MWa = 2e3 * np.divide( (self.MW * MW_air), (self.MW + MW_air))

        r =  (1/p) * 1e5 * np.divide( (  (3.03 - np.divide(0.98,( np.power(MWa,0.5)) ) ) * 1e-27 * (Tin**1.5) ) , ( np.multiply(np.multiply( np.power(MWa,0.5), np.power(sigmaAB,2)  ) ,omegaD ) ) )
        return r
