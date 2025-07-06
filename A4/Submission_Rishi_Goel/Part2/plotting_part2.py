import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read the data from the text file
data = pd.read_csv('Part2/execution_times.txt', sep='\t')

# Create the plot
plt.figure(figsize=(10, 6),num=1)

# Plot execution time vs N with normal scales
plt.plot(data['N'], data['Execution Time (ms)'], 'bo-', linewidth=2, markersize=8)
plt.xlabel('N', fontsize=12)
plt.ylabel('Execution Time (ms)', fontsize=12)
plt.title('Execution Time vs N', fontsize=14)
plt.grid(True)

# Show the plot
plt.savefig('Plots/execution_time_plot.png', dpi=300)
plt.show()







# Read the data from the text file
data = pd.read_csv('Part2/ground_state.txt', sep='\t')

# Create the plot
plt.figure(figsize=(10, 6),num=1)

plt.plot(data['g'], data['Ground State Energy'], 'bo-', linewidth=2, markersize=4, label = "Numerical")
plt.plot(data['g'], data['Theory'], 'r-', linewidth=2, label = "Theory")
plt.xlabel('g', fontsize=12)
plt.ylabel('Ground State Energy', fontsize=12)
plt.title('Ground State Energy vs g', fontsize=14)
plt.grid(True)
plt.legend()

# Show the plot
plt.savefig('Plots/ground_state_energy_plot.png', dpi=300)
plt.show()





# Create the plot
plt.figure(figsize=(10, 6),num=1)

plt.plot(data['g'], data['Second Derivative'], 'bo-', linewidth=2, markersize=4, label = "Numerical")
plt.plot(data['g'], data['Theory Second Derivative'], 'r-', linewidth=2, label = "Theory")
plt.xlabel('g', fontsize=12)
plt.ylabel('Second Derivative', fontsize=12)
plt.title('Second Derivative vs g', fontsize=14)
plt.grid(True)
plt.legend()

# Show the plot
plt.savefig('Plots/ground_state_2nd_deriv.png', dpi=300)
plt.show()








import os

# Load the data
data = np.loadtxt("Part2/observables.txt")

# Extract columns
g_init_vals = data[:, 0]
g_vals = data[:, 1]
time = data[:, 2]
sz = data[:, 3]
sx = data[:, 4]
cxx = data[:, 5]

# Get unique g_init values
unique_g_init = np.unique(g_init_vals)

# Make sure output directory exists
os.makedirs('Plots', exist_ok=True)

# Loop over each g_init and generate a figure
for g_init in unique_g_init:
    fig_data_mask = g_init_vals == g_init
    g_subset = g_vals[fig_data_mask]
    time_subset = time[fig_data_mask]
    sz_subset = sz[fig_data_mask]
    sx_subset = sx[fig_data_mask]
    cxx_subset = cxx[fig_data_mask]

    unique_g = np.unique(g_subset)
    num_g = len(unique_g)

    fig, axs = plt.subplots(num_g, 1, figsize=(8, 4 * num_g), sharex=True)

    if num_g == 1:
        axs = [axs]  # ensure axs is iterable

    for i, g in enumerate(unique_g):
        mask = g_subset == g
        axs[i].plot(time_subset[mask], sz_subset[mask], label="sz")
        axs[i].plot(time_subset[mask], sx_subset[mask], label="sx")
        axs[i].plot(time_subset[mask], cxx_subset[mask], label="cxx")
        axs[i].set_title(f"g = {g}")
        axs[i].set_ylabel("observables")
        axs[i].legend()

    axs[-1].set_xlabel("Time")
    plt.suptitle(f"g_init = {g_init}", fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'Plots/observables_ginit_{g_init:.2f}.png', dpi=300)
    plt.close()
