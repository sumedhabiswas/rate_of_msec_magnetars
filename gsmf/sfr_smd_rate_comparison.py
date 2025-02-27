import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Constants
M_MW = 1e10  # Stellar mass of the Milky Way in M_sun
SFR_MW = 2.0  # Star Formation Rate of the Milky Way in M_sun/yr
SFR_MW_error = 0.7  # Error in SFR of the MW

# Data from Table C.1 up to z ~ 4
z_bins = [(0.2, 0.5), (0.5, 0.8), (0.8, 1.1), (1.1, 1.5), 
          (1.5, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.5)]

# ρ* values with their errors (×10⁷ M_⊙ Mpc⁻³)
rho_star_data = [
    (15.91, 1.50, -1.43),
    (15.49, 1.25, -1.17),
    (17.51, 1.33, -1.27),
    (12.64, 0.86, -0.82),
    (7.62, 0.62, -0.63),
    (4.25, 0.58, -0.54),
    (3.51, 0.33, -0.33),
    (2.22, 0.24, -0.28),
    (0.87, 0.08, -0.08)
]

# Calculate bin centers for interpolation
z_centers = np.array([(z1 + z2)/2 for z1, z2 in z_bins])

# Convert ρ* values to M_⊙ Gpc⁻³
rho_star_values = np.array(rho_star_data)
rho_star_values *= 1e16  # Convert from (10⁷ M_⊙ Mpc⁻³) to (M_⊙ Gpc⁻³)

# Create interpolation functions
rho_star_interp = interp1d(z_centers, rho_star_values[:, 0], kind='linear', fill_value='extrapolate')
rho_star_upper = interp1d(z_centers, rho_star_values[:, 0] + rho_star_values[:, 1], kind='linear', fill_value='extrapolate')
rho_star_lower = interp1d(z_centers, rho_star_values[:, 0] + rho_star_values[:, 2], kind='linear', fill_value='extrapolate')

# r_MW values
r_MW_values = [1e-6, 1e-5, 1e-4, 1e-3]  # yr^-1

# Redshift range
z_values = np.linspace(0, 4, 100)

# SFR(z) function from Madau & Dickinson 2014
def SFR_z(z):
    """Returns SFR density in M_sun Mpc^-3 yr^-1"""
    return 0.015 * (1 + z)**2.7 / (1 + ((1 + z) / 2.9)**5.6)

# Initialize dictionaries for all rates and their uncertainties
R_SMD = {}
R_SMD_upper = {}
R_SMD_lower = {}
R_SFR = {}
R_SFR_upper = {}
R_SFR_lower = {}

# Calculate all rates and their uncertainties
for r_MW in r_MW_values:
    # R_SMD calculations with error bands
    R_SMD[r_MW] = (r_MW / M_MW) * rho_star_interp(z_values)
    R_SMD_upper[r_MW] = (r_MW / M_MW) * rho_star_upper(z_values)
    R_SMD_lower[r_MW] = (r_MW / M_MW) * rho_star_lower(z_values)
    
    # R_SFR calculation with error bands
    sfr_nominal = SFR_z(z_values) * 1e9  # Convert to Gpc^-3
    R_SFR[r_MW] = (r_MW / SFR_MW) * sfr_nominal
    R_SFR_upper[r_MW] = (r_MW / (SFR_MW - SFR_MW_error)) * sfr_nominal
    R_SFR_lower[r_MW] = (r_MW / (SFR_MW + SFR_MW_error)) * sfr_nominal

# Set global font sizes
plt.rcParams.update({'font.size': 20,
                    'axes.labelsize': 20,
                    'axes.titlesize': 20,
                    'xtick.labelsize': 20,
                    'ytick.labelsize': 20,
                    'legend.fontsize': 20})

# Plotting
plt.figure(figsize=(12, 8))

# Colorblind friendly palette (IBM ColorBlind Safe palette)
colors = ['#648FFF', '#785EF0', '#DC267F', '#FE6100']
legend_elements = []

# Add simple legend for line types
legend_elements.append(plt.Line2D([0], [0], color='gray', label='$\mathcal{R}_{SFR}$', linewidth=2))
legend_elements.append(plt.Line2D([0], [0], color='gray', linestyle='--', label='$\mathcal{R}_{SMD}$', linewidth=2))

for i, r_MW in enumerate(r_MW_values):
    # Plot R_SFR with error region
    plt.plot(z_values, R_SFR[r_MW], 
             color=colors[i],
             linewidth=2)
    plt.fill_between(z_values, 
                     R_SFR_lower[r_MW], 
                     R_SFR_upper[r_MW], 
                     color=colors[i],
                     alpha=0.1)
    
    # Plot R_SMD with error region
    plt.plot(z_values, R_SMD[r_MW], 
             color=colors[i],
             linestyle='--',
             linewidth=2)
    plt.fill_between(z_values, 
                     R_SMD_lower[r_MW], 
                     R_SMD_upper[r_MW], 
                     color=colors[i],
                     alpha=0.2)
    
    # Add r_MW value annotation on the left side
    y_pos = R_SFR[r_MW][0]  # Get y-value at z=0
    plt.annotate(f'$r_{{MW}}=10^{{{int(np.log10(r_MW))}}}$ yr$^{{-1}}$',
                xy=(0, y_pos),
                xytext=(0.83, y_pos),
                textcoords='data',
                color=colors[i],
                fontsize=20,
                horizontalalignment='right',
                verticalalignment='center')


plt.xlabel('Redshift $(z)$')
plt.ylabel('Volumetric Rate (Gpc$^{-3}$ yr$^{-1}$)')
# plt.title('Volumetric Rate vs Redshift (up to $z = 4$)')
plt.grid(True)
plt.yscale('log')
plt.legend(handles=legend_elements, loc='upper right')

# Adjust plot limits to accommodate left-side labels
plt.xlim(0, 4)
plt.tight_layout()
plt.savefig('sfr_smd_vs_rate.pdf',  dpi=300, bbox_inches='tight')
plt.show()


