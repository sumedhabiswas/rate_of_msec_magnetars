import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FuncFormatter, LogLocator, AutoMinorLocator
import matplotlib.patheffects as path_effects
import csv

# Function to calculate B field from P and Pdot
def B_field(P, Pdot):
    return np.sqrt(3.2e19 * P * Pdot)

# Function to generate random points within ranges for each scenario
def generate_points(P_range, B_range, n=50):
    P = np.random.uniform(P_range[0], P_range[1], n)
    B = np.random.uniform(B_range[0], B_range[1], n)
    return P, B

# Function to format ticks in scientific notation
def sci_notation(x, pos):
    return f'$10^{{{int(np.log10(x))}}}$'

# Function to read columns from CSV file
def read_columns(csv_file, columns):
    try:
        with open(csv_file, 'r', newline='') as file:
            reader = csv.DictReader(file)
            column_data = {col: [] for col in columns}
            for row in reader:
                for col in columns:
                    column_data[col].append(float(row[col]))
            return column_data
    except FileNotFoundError:
        print("File not found. Please provide a valid file path.")
    except Exception as e:
        print("An error occurred:", e)

# Set up plot style
plt.style.use('seaborn-whitegrid')
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman'],
    'text.usetex': True,
    'axes.linewidth': 1.5,
    'axes.edgecolor': 'black',
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.minor.width': 1,
    'ytick.minor.width': 1,
})

# Professional color palette
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

# Scenario data (updated ranges based on the provided text)
scenarios = {
    'BWD Merger': {'P': [1e-3, 2e-3], 'B': [1e14, 1e16], 'rate': 1e5},
    'NS-WD Merger': {'P': [1e-3, 1e-2], 'B': [1e14, 1e16], 'rate': 2e2},
    'BNS Merger': {'P': [1e-3, 1e-1], 'B': [1e14, 1e16], 'rate': 5e2},
    'Massive Star': {'P': [1e-2, 2], 'B': [1e14, 1e15], 'rate': 5e4}
}

# Read magnetar data from CSV file
csv_file = 'Tab2.csv'
columns_to_read = ['Period', 'B']
data = read_columns(csv_file, columns_to_read)

if data:
    real_P = np.array(data['Period'])
    real_B = np.array(data['B'])
else:
    print("Failed to read magnetar data. Using empty arrays.")
    real_P = np.array([])
    real_B = np.array([])

# Function to create and save plot
def save_plot(fig, filename):
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close(fig)


# Function to create and format the plot
def create_plot(ax, x_data, y_data, labels, colors, x_label, y_label, is_pdot=False):
    markers = ['o', 's', 'D', '^']  # Different marker styles for each scenario
    for (label, data), color, marker in zip(labels.items(), colors, markers):
        P, B = generate_points(data['P'], data['B'])
        if is_pdot:
            Pdot = B**2 / (3.2e19 * P)
            ax.scatter(P, Pdot, c=[color], label=label, s=80, alpha=0.7, 
                       edgecolors='black', linewidth=0.5, marker=marker)
        else:
            ax.scatter(P, B, c=[color], label=label, s=80, alpha=0.7, 
                       edgecolors='black', linewidth=0.5, marker=marker)

    if is_pdot:
        real_Pdot = y_data**2 / (3.2e19 * x_data)
        ax.scatter(x_data, real_Pdot, c='black', marker='*', s=250, 
                   label='Known Magnetars', zorder=10, edgecolors='white', linewidth=1)
    else:
        ax.scatter(x_data, y_data, c='black', marker='*', s=250, 
                   label='Known Magnetars', zorder=10, edgecolors='white', linewidth=1)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(x_label, fontsize=22)
    ax.set_ylabel(y_label, fontsize=22)

    formatter = FuncFormatter(sci_notation)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=6))
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.tick_params(axis='both', which='major', labelsize=18, length=10, width=1.5)
    ax.tick_params(axis='both', which='minor', length=5, width=1)

    legend = ax.legend(fontsize=18, loc='best', framealpha=0.9, edgecolor='black')
    legend.get_frame().set_linewidth(1.5)

    ax.grid(True, which='major', linestyle='--', alpha=0.7)
    ax.grid(True, which='minor', linestyle=':', alpha=0.4)

# B-P Diagram
def create_bp_diagram():
    fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
    create_plot(ax, real_P, real_B, scenarios, colors, r'Period (s)', r'Magnetic Field (G)')
    save_plot(fig, 'B-P_Diagram.png')

# P-Pdot Diagram
def create_ppdot_diagram():
    fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
    create_plot(ax, real_P, real_B, scenarios, colors, r'Period (s)', r'Period Derivative (s/s)', is_pdot=True)

    # Add lines of constant B field
    B_lines = [1e12, 1e13, 1e14, 1e15, 1e16]
    for B in B_lines:
        P = np.logspace(-3, 2, 100)
        Pdot = B**2 / (3.2e19 * P)
        ax.plot(P, Pdot, 'k--', alpha=0.5)
        text = ax.text(1e2, B**2 / (3.2e19 * 1e2), f'$10^{{{int(np.log10(B))}}}$ G', 
                       rotation=-40, va='bottom', ha='left', alpha=0.7, fontsize=16)
        text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='white'),
                               path_effects.Normal()])

    save_plot(fig, 'P-Pdot_Diagram.png')

# Stacked Bar Chart
def create_stacked_bar_chart():
    fig, ax = plt.subplots(figsize=(10, 6))
    scenarios_df = pd.DataFrame(scenarios).T
    scenarios_df['rate'].plot(kind='bar', ax=ax, color=colors)
    ax.set_ylabel('Formation Rate (Gpc$^{-3}$ yr$^{-1}$)', fontsize=14)
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_title('Estimated Magnetar Formation Rates by Scenario', fontsize=16)
    plt.legend(title='Scenario', fontsize=10, title_fontsize=12)
    save_plot(fig, 'Stacked_Bar_Chart.png')

# Bubble Plot
def create_bubble_plot():
    fig, ax = plt.subplots(figsize=(10, 8))
    for (scenario, data), color in zip(scenarios.items(), colors):
        ax.scatter(np.mean(data['P']), np.mean(data['B']), s=data['rate']/100, 
                   c=[color], alpha=0.7, edgecolors='black', linewidth=1, label=scenario)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Period (s)', fontsize=14)
    ax.set_ylabel('Magnetic Field (G)', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_title('Magnetar Formation Scenarios', fontsize=16)
    plt.legend(fontsize=10, title='Scenario', title_fontsize=12)
    save_plot(fig, 'Bubble_Plot.png')


# Heatmap
def create_heatmap():
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    for (scenario, data), ax in zip(scenarios.items(), axes.flatten()):
        P = np.logspace(np.log10(data['P'][0]), np.log10(data['P'][1]), 100)
        B = np.logspace(np.log10(data['B'][0]), np.log10(data['B'][1]), 100)
        X, Y = np.meshgrid(P, B)
        Z = np.random.rand(100, 100)  # Replace with actual probability distribution if available
        c = ax.pcolormesh(X, Y, Z, cmap='viridis', shading='auto')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Period (s)', fontsize=12)
        ax.set_ylabel('Magnetic Field (G)', fontsize=12)
        ax.set_title(scenario, fontsize=14)
        fig.colorbar(c, ax=ax)
    fig.suptitle('Magnetar Formation Probability Heatmaps', fontsize=16)
    save_plot(fig, 'Heatmap.png')


# Box and Whisker Plot
def create_box_plot():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    # Prepare data
    periods = [data['P'] for data in scenarios.values()]
    fields = [data['B'] for data in scenarios.values()]
    
    # Plot for Periods
    sns.boxplot(data=periods, ax=ax1, palette=colors)
    ax1.set_ylabel('Period (s)', fontsize=14)
    ax1.set_yscale('log')
    ax1.set_xticklabels(scenarios.keys(), rotation=45, ha='right')
    ax1.tick_params(axis='both', which='major', labelsize=12)
    
    # Plot for Magnetic Fields
    sns.boxplot(data=fields, ax=ax2, palette=colors)
    ax2.set_ylabel('Magnetic Field (G)', fontsize=14)
    ax2.set_yscale('log')
    ax2.set_xticklabels(scenarios.keys(), rotation=45, ha='right')
    ax2.tick_params(axis='both', which='major', labelsize=12)
    
    fig.suptitle('Distribution of Periods and Magnetic Fields by Formation Scenario', fontsize=16)
    plt.tight_layout()
    save_plot(fig, 'Box_Whisker_Plot.png')

# Violin Plot
def create_violin_plot():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    # Prepare data
    periods = [data['P'] for data in scenarios.values()]
    fields = [data['B'] for data in scenarios.values()]
    
    # Plot for Periods
    sns.violinplot(data=periods, ax=ax1, palette=colors)
    ax1.set_ylabel('Period (s)', fontsize=14)
    ax1.set_yscale('log')
    ax1.set_xticklabels(scenarios.keys(), rotation=45, ha='right')
    ax1.tick_params(axis='both', which='major', labelsize=12)
    
    # Plot for Magnetic Fields
    sns.violinplot(data=fields, ax=ax2, palette=colors)
    ax2.set_ylabel('Magnetic Field (G)', fontsize=14)
    ax2.set_yscale('log')
    ax2.set_xticklabels(scenarios.keys(), rotation=45, ha='right')
    ax2.tick_params(axis='both', which='major', labelsize=12)
    
    fig.suptitle('Distribution of Periods and Magnetic Fields by Formation Scenario', fontsize=16)
    plt.tight_layout()
    save_plot(fig, 'Violin_Plot.png')

if __name__ == "__main__":
    # Create and save all plots
    create_bp_diagram()
    create_ppdot_diagram()
    create_stacked_bar_chart()
    create_bubble_plot()
    create_heatmap()	
    create_box_plot()
    create_violin_plot()
    
    print("All plots have been created, displayed, and saved.")
