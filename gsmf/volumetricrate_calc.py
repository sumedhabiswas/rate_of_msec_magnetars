import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# Constants for the double Schechter function
phi1 = 0.4e-3  # Gpc^-3
phi2 = 0.6e-3  # Gpc^-3
alpha1 = -1.3
alpha2 = -1.5
M_star = 1e10  # Solar masses

def double_schechter(M):
    """Double Schechter function."""
    M_ratio = M / M_star
    term1 = phi1 * (M_ratio ** alpha1) * np.exp(-M_ratio)
    term2 = phi2 * (M_ratio ** alpha2) * np.exp(-M_ratio)
    return np.log(10) * (term1 + term2)

def cumulative_number_density(mass_range):
    """Calculate the cumulative number density of galaxies."""
    cumulative_density = np.zeros_like(mass_range)
    for i, M in enumerate(mass_range):
        cumulative_density[i], _ = integrate.quad(double_schechter, mass_range[0], M)
    return cumulative_density

def main():
    # Step 1: Take user input for the rate
    rate_per_year = float(input("Enter the rate per year: "))
    print(f"Rate per year: {rate_per_year:.2e} yr^-1")

    # Step 2: Ask for the mass of the galaxy
    galaxy_mass = float(input("Enter the mass of the galaxy (in solar masses): "))
    print(f"Galaxy mass: {galaxy_mass:.2e} M_sun")

    # Step 3: Calculate R
    R = rate_per_year / galaxy_mass
    print(f"Rate per year per solar mass, R: {R:.2e} yr^-1 M_sun^-1")

    # Step 4 & 5: Integrate the double Schechter function over the mass range
    mass_range = np.logspace(8, 12, 1000)
    n_galaxies, _ = integrate.quad(double_schechter, 1e8, 1e12)
    print(f"Integrated number density of galaxies, n_galaxies: {n_galaxies:.2e} Gpc^-3")

    # Additional print statement for rate per year * n_galaxies
    rate_times_n_galaxies = rate_per_year * n_galaxies
    print(f"Rate per year multiplied by n_galaxies: {rate_times_n_galaxies:.2e} yr^-1 Gpc^-3")

    # Step 6: Calculate the value of the double Schechter function at the galaxy mass
    n_value = double_schechter(galaxy_mass)
    print(f"Number density at galaxy mass, n_value: {n_value:.2e} Gpc^-3")

    # Step 7: Define normalization value N
    N = n_galaxies / n_value
    print(f"Normalization value, N: {N:.2e}")

    # Step 8: Multiply R * n_galaxies * N
    result = R * n_galaxies * N
    print(f"Result of R * n_galaxies * N: {result:.2e} yr^-1 Gpc^-3")

    # Step 9: Plot the double Schechter function
    schechter_values = double_schechter(mass_range)
    
    plt.figure(figsize=(10, 6))
    plt.loglog(mass_range, schechter_values, label='Double Schechter Function', color='blue')
    plt.scatter(galaxy_mass, n_value, color='red', zorder=5, label=f'Galaxy Mass: {galaxy_mass:.2e} M_sun')
    plt.xlabel('Mass (Solar Masses)', fontsize=12)
    plt.ylabel('Number Density (Gpc$^{-3}$)', fontsize=12)
    plt.title('Double Schechter Function', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    
    # Adjusting the plot limits for better visibility
    plt.xlim(1e8, 1e12)  # x-axis from 10^8 to 10^12
    plt.ylim(1e-14, 1e-1)  # y-axis from 10^-14 to 10^-1
    plt.tight_layout()
    plt.show()

    # New Plot: Cumulative Number Density
    cumulative_density = cumulative_number_density(mass_range)
    
    plt.figure(figsize=(10, 6))
    plt.loglog(mass_range, cumulative_density, label='Cumulative Number Density', color='green')
    plt.xlabel('Mass (Solar Masses)', fontsize=12)
    plt.ylabel('Cumulative Number Density (Gpc$^{-3}$)', fontsize=12)
    plt.title('Cumulative Number Density of Galaxies', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    
    # Adjusting the plot limits for better visibility
    plt.xlim(1e8, 1e12)  # x-axis from 10^8 to 10^12
    plt.ylim(1e-14, 1e-1)  # y-axis from 10^-14 to 10^-1
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
