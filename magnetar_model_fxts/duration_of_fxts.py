import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker

# Use LaTeX for text rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 10,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 8,
})

# Constants
c = 3e10  # Speed of light in cm/s
I = 1e45  # Moment of inertia in g*cm^2
R_M = 1.2e6  # Magnetar radius in cm

# Function to calculate tau_EM
def tau_EM(B_p, P_i):
    return (3 * c**3 * I * P_i**2) / (B_p**2 * R_M**6 * (2 * np.pi)**2)

# Function to calculate L_0^EM
def L_0_EM(B_p, P_i):
    tau = tau_EM(B_p, P_i)
    Omega_i = 2 * np.pi / P_i
    return (I * Omega_i**2) / (2 * tau)

# Create arrays for B_p and P_i
B_p_range = np.logspace(14, 16, 1000)  # 10^14 to 10^16 G
P_i_range = np.linspace(1e-3, 2e-3, 1000)  # 1 ms to 2 ms

# Create meshgrid
B_p_mesh, P_i_mesh = np.meshgrid(B_p_range, P_i_range)

# Calculate tau_EM and L_0^EM for each combination of B_p and P_i
tau_EM_mesh = tau_EM(B_p_mesh, P_i_mesh) / 1000  # Convert to kiloseconds
L_0_EM_mesh = L_0_EM(B_p_mesh, P_i_mesh)

# Function to format axes
def format_axes(ax):
    ax.set_xscale('log')
    ax.set_xlabel(r'Magnetic Field Strength $B_p$ (G)')
    ax.set_ylabel(r'Initial Spin Period $P_i$ (ms)')
    ax.set_ylim(1, 2)
    ax.set_yticks(np.linspace(1, 2, 6))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: r'$10^{{{:.0f}}}$'.format(np.log10(x))))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.grid(True, which="both", ls="-", alpha=0.2)

# Create a single figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

# Plot 1: tau_EM
cs1 = ax1.contourf(B_p_mesh, P_i_mesh * 1000, tau_EM_mesh, levels=np.logspace(0, 4, 20), cmap='viridis', norm=LogNorm())
cbar1 = fig.colorbar(cs1, ax=ax1, label=r'$\tau_{\mathrm{EM}}$ (ks)')
cbar1.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: r'$10^{{{:.0f}}}$'.format(np.log10(x))))
cs1_specific = ax1.contour(B_p_mesh, P_i_mesh * 1000, tau_EM_mesh, levels=[10, 20], colors=['yellow'], linestyles=['solid', 'dashed'])
ax1.clabel(cs1_specific, inline=True, fmt='%1.0f ks', fontsize=8)
ax1.set_title(r'(a) Electromagnetic Spin-down Timescale ($\tau_{\mathrm{EM}}$)')
format_axes(ax1)

# Plot 2: L_0^EM
cs2 = ax2.contourf(B_p_mesh, P_i_mesh * 1000, np.log10(L_0_EM_mesh), levels=np.linspace(45, 50, 20), cmap='plasma')
cbar2 = fig.colorbar(cs2, ax=ax2, label=r'$\log_{10}(L_0^{\mathrm{EM}})$ (erg/s)')
cbar2.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: r'$10^{{{:.0f}}}$'.format(x)))
ax2.set_title(r'(b) Initial Spin-down Luminosity ($L_0^{\mathrm{EM}}$)')
format_axes(ax2)

# Adjust layout
plt.tight_layout()

# Save the figure
plt.savefig('Magnetar_FXT_Model.pdf', bbox_inches='tight')
plt.savefig('Magnetar_FXT_Model.png', bbox_inches='tight', dpi=300)

plt.show()
