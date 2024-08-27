import numpy as np
from scipy.integrate import quad

# Constants
H0 = 70  # Hubble constant in km/s/Mpc
c = 3e5  # Speed of light in km/s
Omega_m = 0.3  # Matter density parameter
Omega_Lambda = 0.7  # Dark energy density parameter

# Hubble parameter as a function of redshift
def H(z):
    return H0 * np.sqrt(Omega_m * (1 + z)**3 + Omega_Lambda)

# Integrand for the comoving distance
def integrand(z):
    return c / H(z)

# Function to calculate luminosity distance
def calculate_luminosity_distance(z):
    D_M, _ = quad(integrand, 0, z)
    D_L = (1 + z) * D_M
    return D_L

# Main execution
if __name__ == "__main__":
    # Get redshift input from user
    z = float(input("Enter the redshift: "))

    # Calculate the luminosity distance
    D_L = calculate_luminosity_distance(z)

    print(f"Luminosity distance at z={z}: {D_L:.2f} Mpc")
