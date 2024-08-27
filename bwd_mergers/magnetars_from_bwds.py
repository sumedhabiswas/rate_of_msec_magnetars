import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker

# Constants
R_NS = 1e6  # Neutron star radius in cm
R_WD = 7e8  # White dwarf radius in cm

# Functions for calculations
def calculate_b_ns(b_wd):
    return b_wd * (R_WD / R_NS)**2

# Generate data
p_ns_range = np.logspace(-3, 2, 200)  # Extended range
b_wd_range = np.logspace(4, 12, 200)  # Extended range

# Create meshgrid for contour plot
P_NS, B_WD = np.meshgrid(p_ns_range, b_wd_range)
B_NS = calculate_b_ns(B_WD)

# Plotting
plt.figure(figsize=(12, 10))

# B_NS vs P_NS (Contour plot)
contour = plt.contourf(P_NS, B_NS, B_WD, levels=20, cmap='viridis', norm=LogNorm())
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Neutron Star Period (s)', fontsize=14)
plt.ylabel('Neutron Star Magnetic Field (G)', fontsize=14)
plt.title('Neutron Star B-P Diagram', fontsize=16)

# Format tick labels
plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'$10^{{{int(np.log10(x))}}}$'))
plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'$10^{{{int(np.log10(x))}}}$'))

# Add colorbar
cbar = plt.colorbar(contour)
cbar.set_label('White Dwarf Magnetic Field (G)', fontsize=14)
cbar.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'$10^{{{int(np.log10(x))}}}$'))

# Add dashed lines for B_NS range
plt.axhline(y=1e14, color='r', linestyle='--', linewidth=2)
plt.axhline(y=1e16, color='r', linestyle='--', linewidth=2)
plt.text(1e-3, 1e14, '$10^{14}$ G', color='r', verticalalignment='bottom')
plt.text(1e-3, 1e16, '$10^{16}$ G', color='r', verticalalignment='top')

# Add dashed lines for P_NS range
plt.axvline(x=1e-3, color='b', linestyle='--', linewidth=2)
plt.axvline(x=2e-3, color='b', linestyle='--', linewidth=2)
plt.text(1e-3, 1e12, '1 ms', color='b', horizontalalignment='right', rotation=90)
plt.text(2e-3, 1e12, '2 ms', color='b', horizontalalignment='left', rotation=90)

plt.tight_layout()
plt.savefig('neutron_star_bp_diagram.png', dpi=300, bbox_inches='tight')
plt.show()

# Print some specific values
print(f"For B_WD = 1e6 G, B_NS = {calculate_b_ns(1e6):.2e} G")
print(f"For B_WD = 1e8 G, B_NS = {calculate_b_ns(1e8):.2e} G")
print(f"For B_WD = 1e10 G, B_NS = {calculate_b_ns(1e10):.2e} G")
