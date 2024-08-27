import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogFormatterSciNotation

# Constants
flux_limit = 8.9e-10  # erg/s/cm²
luminosities = np.logspace(44, 47, num=100)  # Range of luminosities from 10^44 to 10^47 erg/s
N_FXT = 60  # Number of sources
T = 1  # Time in years
Omega = 4 * np.pi  # Solid angle in steradians

# Calculate distances
distances = np.sqrt(luminosities / (4 * np.pi * flux_limit))

# Convert distances from cm to Gpc (1 Gpc = 3.086e27 cm)
distances_gpc = distances / 3.086e27

# Calculate V_max in Gpc^3
v_max = (4 * np.pi * distances_gpc**3) / 3

# Calculate ρ_FXT
rho_FXT = (N_FXT * 4 * np.pi) / (v_max * Omega * T)

# Create figure
fig, ax = plt.subplots(figsize=(12, 9))

# Plot ρ_FXT
ax.loglog(luminosities, rho_FXT, 'k-', linewidth=2, alpha=0.7)

# Add reference points
reference_luminosities = [1e44, 1e45, 1e46, 1e47]
reference_rho_FXT = (N_FXT * 4 * np.pi) / ((4 * np.pi * (np.sqrt(np.array(reference_luminosities) / (4 * np.pi * flux_limit)) / 3.086e27)**3) / 3 * Omega * T)

ax.plot(reference_luminosities, reference_rho_FXT, 'o', color='#1f77b4', markersize=10, label='Reference Points')

# Formatting
ax.set_xlabel('Luminosity (erg s$^{-1}$)', fontsize=16)
ax.set_ylabel('$\\rho_{\\mathrm{FXT}}$ (Gpc$^{-3}$ yr$^{-1}$)', fontsize=16)
ax.set_title('Space Density Rate of Detectable Sources', fontsize=18)

ax.grid(True, which="both", ls="--", alpha=0.3)
ax.tick_params(axis='both', which='major', labelsize=14)

# Format x-axis to use scientific notation
ax.xaxis.set_major_formatter(LogFormatterSciNotation())
ax.xaxis.set_tick_params(which='minor', bottom=False)

# Format y-axis to use scientific notation
ax.yaxis.set_major_formatter(LogFormatterSciNotation())
ax.yaxis.set_tick_params(which='minor', left=False)

# Add legend
ax.legend(fontsize=14)

# Add text annotations with scientific notation (coefficient × 10^power)
for L, r in zip(reference_luminosities, reference_rho_FXT):
    power = int(np.log10(r))
    coeff = r / (10**power)
    ax.annotate(f'{coeff:.1f}$\\times10^{{{power}}}$', (L, r), textcoords="offset points", xytext=(0,10), 
                ha='center', va='bottom', fontsize=12, alpha=0.8)

# Adjust layout
plt.tight_layout()

# Show plot
plt.show()

# Print out the reference values
print("Reference values:")
for L, r in zip(reference_luminosities, reference_rho_FXT):
    print(f"Luminosity: {L:.2e} erg/s, ρ_FXT: {r:.2e} Gpc⁻³ yr⁻¹")
