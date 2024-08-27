import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FuncFormatter
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D

# Constants
G = 6.67430e-8  # Gravitational constant in cgs units
c = 2.99792458e10  # Speed of light in cm/s
M_sun = 1.989e33  # Solar mass in grams

# Magnetar parameters
eta = 0.001  # Efficiency factor for X-ray conversion

# Function to calculate moment of inertia
def calculate_I(M, R):
    return (2/5) * M * R**2

# Function to calculate luminosities
def calculate_luminosities(B, P, M, R):
    I = calculate_I(M, R)
    Omega = 2 * np.pi / (P * 1e-3)
    L_sd = (B**2 * R**6 * Omega**4) / (6 * c**3)
    L_X = eta * L_sd
    return L_sd, L_X

# Create arrays for B and P
B_range = np.logspace(14, 16, 100)
P_range = np.linspace(1, 2, 100)

# Create meshgrid
B_mesh, P_mesh = np.meshgrid(B_range, P_range)

# Custom formatter for scientific notation
def scientific_formatter(x, pos):
    exp = int(np.log10(x))
    return r'$10^{%d}$' % exp

# Function to create the contour plot
def create_contour_plot(M_NS, R_NS):
    L_sd, L_X = calculate_luminosities(B_mesh, P_mesh, M_NS, R_NS)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.rcParams.update({'font.size': 18, 'font.family': 'serif'})  # Increase base font size

    contour_sd = ax.contourf(B_mesh, P_mesh, np.log10(L_sd), levels=20, cmap='viridis')
    contour_lines_x = ax.contour(B_mesh, P_mesh, L_X, levels=LogLocator(numticks=6).tick_values(L_X.min(), L_X.max()),
                                 colors='white', linewidths=1.5)
    ax.clabel(contour_lines_x, inline=True, fontsize=12, fmt=lambda x: scientific_formatter(x, 0))
    ax.set_xscale('log')
    ax.set_xlabel('B (G)', fontsize=20)
    ax.set_ylabel('P (ms)', fontsize=20)
    ax.set_title(f'M = {M_NS/M_sun:.1f} M☉, R = {R_NS/1e5:.1f} km', fontsize=22)
    ax.xaxis.set_major_formatter(FuncFormatter(scientific_formatter))
    ax.tick_params(which='both', direction='in', top=True, right=True)

    cbar = fig.colorbar(contour_sd, ax=ax, pad=0.05)  # Add space between plot and colorbar
    cbar.ax.set_ylabel(r'$L_{\rm sd}$ [erg s$^{-1}$]', rotation=270, labelpad=25, fontsize=20)
    
    # Correct the colorbar ticks and labels
    cbar_ticks = cbar.get_ticks()
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([scientific_formatter(10**tick, 0) for tick in cbar_ticks])

    # Add legend for L_X
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color='white', lw=1.5, label=r'$L_{\rm X}$')]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=16)

    plt.tight_layout()
    filename = f'luminosity_contour_plot_M{M_NS/M_sun:.1f}_R{R_NS/1e5:.1f}.pdf'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved contour plot as {filename}")
    plt.show()

# Main execution
if __name__ == "__main__":
    # Original contour plots
    create_contour_plot(1.4 * M_sun, 12e5)  # Original: M = 1.4 M☉, R = 12 km
    create_contour_plot(2.0 * M_sun, 12e5)
    create_contour_plot(2.0 * M_sun, 10e5)  # New: M = 2.0 M☉, R = 10 km
    
