import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Define the double Schechter function
def double_schechter(m, phi_1, phi_2, alpha_1, alpha_2, M_star):
    term1 = phi_1 * (m / M_star)**alpha_1
    term2 = phi_2 * (m / M_star)**alpha_2
    return (term1 + term2) * np.exp(-m / M_star) / M_star

# Define the parameters from the table
schechter_params = [
    # Reference, z, log(M_star/M_sun), phi_1, phi_2, alpha_1, alpha_2
    ("2012baldry", 0.03, 10.66, 3.96e-3, 0.79e-3, -0.35, -1.47),
    ("2016weigel", 0.04, 10.79, 4.90e-4, 9.77e-3, -1.69, -0.79),
    ("2021hst", 0.5, 10.64, 2.34e-3, 7.76e-4, 0.25, -1.49),
    ("2021hst", 1.0, 10.51, 2.14e-3, 8.51e-4, 0.08, -1.49),
    ("2021hst", 1.5, 10.54, 1.48e-3, 4.79e-4, -0.07, -1.60),
    ("2021hst", 2.0, 10.56, 8.91e-4, 3.09e-4, -0.06, -1.63),
    ("2021hst", 2.5, 10.55, 5.25e-4, 3.16e-4, 0.02, -1.66),
    ("2021hst", 3.25, 10.64, 8.32e-5, 1.82e-4, 0.35, -1.76)
]

# Mass range for integration
mass_min = 10**8
mass_max = 10**12
mass_11 = 10**11

# Define event rates in yr^-1
event_rates = [
    (1e-6, "BWD merger rate (Sigurdsson & Rees 97)"),
    (3e-4, "Magnetars from BWD mergers (Levan+06)"),
    (1e-6, "NSWD merger rate (Kim+04) min"),
    (1e-5, "NSWD merger rate (Kim+04) max"),
    (1e-5, "NSWD merger rate (multiple papers) min"),
    (1e-3, "NSWD merger rate (multiple papers) max"),
    (1e-3, "Magnetars from massive stars (Beniamini+19) min"),
    (1e-2, "Magnetars from massive stars (Beniamini+19) max")
]

# Set font sizes for plots
plt.rcParams.update({'font.size': 14})

# Prepare for comparison plot
comparison_fig, comparison_ax = plt.subplots(figsize=(12, 8))

# Iterate over each event rate
for event_rate, label in event_rates:
    # Convert event rate to per solar mass
    r_per_unit_mass = event_rate / mass_11  # yr^-1 M_odot^-1
    
    # Lists to store results for plotting
    redshifts = []
    volumetric_rates = []

    # Iterate over each set of Schechter parameters
    for ref, z, log_M_star, phi_1, phi_2, alpha_1, alpha_2 in schechter_params:
        M_star = 10**log_M_star  # Convert log(M_star/M_sun) to M_sun

        # Integrate the function from 10^8 to 10^12 solar masses
        n_gal_mpc3, error = quad(double_schechter, mass_min, mass_max, args=(phi_1, phi_2, alpha_1, alpha_2, M_star))

        # Convert result to Gpc^-3
        n_gal_gpc3 = n_gal_mpc3 * 1e9

        # Calculate n_gal_11 at m = 10^11 solar masses
        n_gal_11_mpc3 = double_schechter(mass_11, phi_1, phi_2, alpha_1, alpha_2, M_star)

        # Convert result to Gpc^-3 M_odot^-1
        n_gal_11_gpc3 = n_gal_11_mpc3 * 1e9

        # Calculate normalization factor N
        N = n_gal_gpc3 / n_gal_11_gpc3

        # Calculate volumetric rate R
        R = r_per_unit_mass * n_gal_gpc3 * N

        # Store results for plotting
        redshifts.append(z)
        volumetric_rates.append(R)

        # Print results
        print(f"Event Rate: {label}")
        print(f"Reference: {ref}, z: {z}")
        print(f"Integrated number density: {n_gal_gpc3:.2e} Gpc^-3")
        print(f"Number density at 10^11 solar masses: {n_gal_11_gpc3:.2e} Gpc^-3 M_odot^-1")
        print(f"Normalization factor N: {N:.2e}")
        print(f"Volumetric rate R: {R:.2e} Gpc^-3 yr^-1")
        print("-" * 50)

    # Interpolate to make the curve smoother
    redshifts_interp = np.linspace(min(redshifts), max(redshifts), 100)
    volumetric_rates_interp = interp1d(redshifts, volumetric_rates, kind='cubic')(redshifts_interp)

    # Add to comparison plot
    comparison_ax.plot(redshifts_interp, volumetric_rates_interp, linestyle='-', label=f'{label}')

# Finalize comparison plot
comparison_ax.set_xlabel('Redshift')
comparison_ax.set_ylabel('Volumetric Rate R (Gpc$^{-3}$ yr$^{-1}$)')
comparison_ax.set_title('Comparison of Volumetric Rates vs. Redshift')
comparison_ax.legend()
comparison_ax.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()

