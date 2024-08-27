import lumfunc as lf
import numpy as np
from scipy.integrate import quad
from astropy import units as u
from scipy.integrate import romberg
from scipy.integrate import quad_vec
import matplotlib.pyplot as plt
from scipy import integrate

h = 0.7

def double_schechter_mass(M, phi_star1, phi_star2, M_star, alpha1, alpha2):
    """
    Double Schechter function for stellar mass function.
    
    Parameters:
    M (array-like): Stellar mass values
    phi_star1 (float): First normalization factor
    phi_star2 (float): Second normalization factor
    M_star (float): Characteristic mass (shared between both components)
    alpha1 (float): First low-mass slope
    alpha2 (float): Second low-mass slope
    
    Returns:
    array-like: Values of the double Schechter function
    """
    
    log_M = np.log10(M) if np.all(M > 0) else M
    log_M_star = np.log10(M_star) if M_star > 0 else M_star
    
    exp_term = np.exp(-10 ** (log_M - log_M_star))
    first_term = phi_star1 * (10 ** ((log_M - log_M_star) * (alpha1 + 1)))
    second_term = phi_star2 * (10 ** ((log_M - log_M_star) * (alpha2 + 1)))
    
    phi = np.log(10) * exp_term * (first_term + second_term)
    
    return phi


# Set up the mass range and parameters
M_range = np.logspace(8, 12, 1000)  # Mass range from 10^8 to 10^12 solar masses

# Parameters for different redshift ranges
# z < 0.06 
phi_star1_1 = 10**(-2.40)
phi_star2_1 = 10**(-3.10)
M_star_1 = 10**10.66
alpha1_1 = -0.35
alpha2_1 = -1.47  

# 0.25 <= z <= 0.75 
phi_star1_2 = 10**(-2.77)
phi_star2_2 = 10**(-3.26)
M_star_2 = 10**10.8
alpha1_2 = -0.61
alpha2_2 = -1.52

# 0.75 <= z <= 1.25
phi_star1_3 = 10**(-2.80)
phi_star2_3 = 10**(-3.26)
M_star_3 = 10**10.72
alpha1_3 = -0.46
alpha2_3 = -1.53

# 1.25 <= z <= 1.75
phi_star1_4 = 10**(-2.94)
phi_star2_4 = 10**(-3.54)
M_star_4 = 10**10.72
alpha1_4 = -0.55
alpha2_4 = -1.65

# 1.75 <= z <= 2.25 
phi_star1_5 = 10**(-3.18)
phi_star2_5 = 10**(-3.84)
M_star_5 = 10**10.77
alpha1_5 = -0.68
alpha2_5 = -1.73

# 2.25 <= z <= 2.75 
phi_star1_6 = 10**(-3.39)
phi_star2_6 = 10**(-3.78)
M_star_6 = 10**10.77
alpha1_6 = -0.62
alpha2_6 = -1.74

# 2.75 <= z <= 3.75
phi_star1_7 = 10**(-4.3)
phi_star2_7 = 10**(-3.94)
M_star_7 = 10**10.84
alpha1_7 = -0.01
alpha2_7 = -1.79

# Calculate the Schechter function values

# z < 0.06 
phi_values_1 = double_schechter_mass(M_range, phi_star1_1, phi_star2_1, M_star_1, alpha1_1, alpha2_1)

# 0.25 <= z <= 0.75 
phi_values_2 = double_schechter_mass(M_range, phi_star1_2, phi_star2_2, M_star_2, alpha1_2, alpha2_2)

# 0.75 <= z <= 1.25
phi_values_3 = double_schechter_mass(M_range, phi_star1_3, phi_star2_3, M_star_3, alpha1_3, alpha2_3)

# 1.25 <= z <= 1.75
phi_values_4 = double_schechter_mass(M_range, phi_star1_4, phi_star2_4, M_star_4, alpha1_4, alpha2_4)

# 1.75 <= z <= 2.25 
phi_values_5 = double_schechter_mass(M_range, phi_star1_5, phi_star2_5, M_star_5, alpha1_5, alpha2_5)

# 2.25 <= z <= 2.75 
phi_values_6 = double_schechter_mass(M_range, phi_star1_6, phi_star2_6, M_star_6, alpha1_6, alpha2_6)

# 2.75 <= z <= 3.75
phi_values_7 = double_schechter_mass(M_range, phi_star1_7, phi_star2_7, M_star_7, alpha1_7, alpha2_7)

# Integrate the function

# Function to convert Mpc^-3 to Gpc^-3
def mpc3_to_gpc3(value):
    return value * 1e9

# z < 0.06 
def integrand_1(log_M):
    M = 10**log_M
    return double_schechter_mass(M, phi_star1_1, phi_star2_1, M_star_1, alpha1_1, alpha2_1) * M * np.log(10)

# Integrate from 10^8 to 10^12 solar masses 
integral_1, error_1 = integrate.quad(integrand_1, 8, 12)

print(f"z < 0.06 integrated number density: {mpc3_to_gpc3(integral_1):.4e} Gpc^-3")
print(f"z < 0.06 integration error: {mpc3_to_gpc3(error_1):.4e} Gpc^-3")

# 0.25 <= z <= 0.75 
def integrand_2(log_M):
    M = 10**log_M
    return double_schechter_mass(M, phi_star1_2, phi_star2_2, M_star_2, alpha1_2, alpha2_2) * M * np.log(10)

integral_2, error_2 = integrate.quad(integrand_2, 8, 12)

print(f"0.25 <= z <= 0.75 integrated number density: {mpc3_to_gpc3(integral_2):.4e} Gpc^-3")
print(f"0.25 <= z <= 0.75 integration error: {mpc3_to_gpc3(error_2):.4e} Gpc^-3")

# 0.75 <= z <= 1.25
def integrand_3(log_M):
    M = 10**log_M
    return double_schechter_mass(M, phi_star1_3, phi_star2_3, M_star_3, alpha1_3, alpha2_3) * M * np.log(10)

integral_3, error_3 = integrate.quad(integrand_3, 8, 12)

print(f"0.75 <= z <= 1.25 integrated number density: {mpc3_to_gpc3(integral_3):.4e} Gpc^-3")
print(f"0.75 <= z <= 1.25 integration error: {mpc3_to_gpc3(error_3):.4e} Gpc^-3")

# 1.25 <= z <= 1.75
def integrand_4(log_M):
    M = 10**log_M
    return double_schechter_mass(M, phi_star1_4, phi_star2_4, M_star_4, alpha1_4, alpha2_4) * M * np.log(10)

integral_4, error_4 = integrate.quad(integrand_4, 8, 12)

print(f"1.25 <= z <= 1.75 integrated number density: {mpc3_to_gpc3(integral_4):.4e} Gpc^-3")
print(f"1.25 <= z <= 1.75 integration error: {mpc3_to_gpc3(error_4):.4e} Gpc^-3")

# 1.75 <= z <= 2.25 
def integrand_5(log_M):
    M = 10**log_M
    return double_schechter_mass(M, phi_star1_5, phi_star2_5, M_star_5, alpha1_5, alpha2_5) * M * np.log(10)

integral_5, error_5 = integrate.quad(integrand_5, 8, 12)

print(f"1.75 <= z <= 2.25 integrated number density: {mpc3_to_gpc3(integral_5):.4e} Gpc^-3")
print(f"1.75 <= z <= 2.25 integration error: {mpc3_to_gpc3(error_5):.4e} Gpc^-3")

# 2.25 <= z <= 2.75 
def integrand_6(log_M):
    M = 10**log_M
    return double_schechter_mass(M, phi_star1_6, phi_star2_6, M_star_6, alpha1_6, alpha2_6) * M * np.log(10)

integral_6, error_6 = integrate.quad(integrand_6, 8, 12)

print(f"2.25 <= z <= 2.75 integrated number density: {mpc3_to_gpc3(integral_6):.4e} Gpc^-3")
print(f"2.25 <= z <= 2.75 integration error: {mpc3_to_gpc3(error_6):.4e} Gpc^-3")

# 2.75 <= z <= 3.75
def integrand_7(log_M):
    M = 10**log_M
    return double_schechter_mass(M, phi_star1_7, phi_star2_7, M_star_7, alpha1_7, alpha2_7) * M * np.log(10)

integral_7, error_7 = integrate.quad(integrand_7, 8, 12)

print(f"2.75 <= z <= 3.75 integrated number density: {mpc3_to_gpc3(integral_7):.4e} Gpc^-3")
print(f"2.75 <= z <= 3.75 integration error: {mpc3_to_gpc3(error_7):.4e} Gpc^-3")

# Plotting
plt.figure(figsize=(12, 8))

plt.loglog(M_range, phi_values_1, label='z < 0.06')
plt.loglog(M_range, phi_values_2, label='0.25 ≤ z < 0.75')
plt.loglog(M_range, phi_values_3, label='0.75 ≤ z < 1.25')
plt.loglog(M_range, phi_values_4, label='1.25 ≤ z < 1.75')
plt.loglog(M_range, phi_values_5, label='1.75 ≤ z < 2.25')
plt.loglog(M_range, phi_values_6, label='2.25 ≤ z < 2.75')
plt.loglog(M_range, phi_values_7, label='2.75 ≤ z < 3.75')

plt.xlabel('Stellar Mass (M☉)', fontsize=14)
plt.ylabel('Φ (Gpc⁻³ dex⁻¹)', fontsize=14)
plt.title('Double Schechter Mass Function for Different Redshift Ranges', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, which="both", ls="-", alpha=0.2)

plt.tight_layout()
plt.show()

