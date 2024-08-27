# Calculating the Event Rate of FXTs for Einstein Probe, assuming 60 FXT detections in 1 year, takes redshift as input

import numpy as np
from scipy.integrate import quad

# Constants
H0 = 70  # Hubble constant in km/s/Mpc
c = 3e5  # Speed of light in km/s
Omega_m = 0.3  # Matter density parameter
Omega_Lambda = 0.7  # Dark energy density parameter

# Parameters for \rho_{FXT}
N_FXT = 60 
Omega_t = 4 * np.pi  # 4π steradians
T_t = 1 # years

# Function to calculate the Hubble parameter as a function of redshift
def H(z):
    return H0 * np.sqrt(Omega_m * (1 + z)**3 + Omega_Lambda)

# Integrand for the comoving distance
def integrand(z):
    return c / H(z)

# Function to calculate the comoving distance D_M
def comoving_distance(z):
    D_M, _ = quad(integrand, 0, z)
    return D_M

# Function to calculate the luminosity distance D_L
def luminosity_distance(z):
    D_M = comoving_distance(z)
    D_L = (1 + z) * D_M
    return D_L
    

# Function to calculate V_max
def V_max(D_L_max):
    return (4 * np.pi * D_L_max**3) / 3

# Function to calculate \rho_{FXT}
def rho_FXT(N_FXT, Omega_t, T_t, V_max):
    return (N_FXT * 4 * np.pi) / (Omega_t * T_t * V_max)


if __name__ == "__main__":
    # User input for redshift
    z = float(input("Enter the redshift: "))

    # Calculate the luminosity distance at the given redshift
    D_L_max = luminosity_distance(z)

    # Calculate V_max
    V_max_value = V_max(D_L_max)

    # Calculate \rho_{FXT}
    rho_FXT_value = rho_FXT(N_FXT, Omega_t, T_t, V_max_value)

    # Convert V_max to Gpc^3
    V_max_value_Gpc3 = V_max_value / 1e9  # 1 Gpc^3 = 10^9 Mpc^3

    # Convert \rho_{FXT} to Gpc^-3 yr^-1
    rho_FXT_value_Gpc = rho_FXT_value * 1e9  # 1 Gpc^-3 = 10^9 Mpc^-3

    # Print the results
    print(f"Luminosity distance at z={z}: {D_L_max:.2f} Mpc")
    print(f"V_max: {V_max_value:.2e} Mpc^3")
    print(f"V_max: {V_max_value_Gpc3:.2e} Gpc^3")
    print(f"ρ_FXT at z={z}: {rho_FXT_value:.2e} events per Mpc^3 per year")
    print(f"ρ_FXT at z={z}: {rho_FXT_value_Gpc:.2e} events per Gpc^3 per year")
