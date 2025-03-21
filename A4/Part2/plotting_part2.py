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
plt.savefig('Part2/Plots_P2/execution_time_plot.png', dpi=300)
plt.show()







# Read the data from the text file
data = pd.read_csv('Part2/ground_state.txt', sep='\t')

# Create the plot
plt.figure(figsize=(10, 6),num=1)

# Plot execution time vs N with normal scales
plt.plot(data['g'], data['Ground State Energy'], 'bo-', linewidth=2, markersize=4, label = "Numerical")
plt.plot(data['g'], data['Theory'], 'r-', linewidth=2, label = "Theory")
plt.xlabel('g', fontsize=12)
plt.ylabel('Ground State Energy', fontsize=12)
plt.title('Ground State Energy vs g', fontsize=14)
plt.grid(True)
plt.legend()

# Show the plot
plt.savefig('Part2/Plots_P2/ground_state_energy_plot.png', dpi=300)
plt.show()





# Create the plot
plt.figure(figsize=(10, 6),num=1)

# Plot execution time vs N with normal scales
plt.plot(data['g'], data['Second Derivative'], 'bo-', linewidth=2, markersize=4, label = "Numerical")
plt.plot(data['g'], data['Theory Second Derivative'], 'r-', linewidth=2, label = "Theory")
plt.xlabel('g', fontsize=12)
plt.ylabel('Second Derivative', fontsize=12)
plt.title('Second Derivative vs g', fontsize=14)
plt.grid(True)
plt.legend()

# Show the plot
plt.savefig('Part2/Plots_P2/ground_state_2nd_deriv.png', dpi=300)
plt.show()







# Load the data
data = np.loadtxt("Part2/observables.txt") 

# Extract columns
time = data[:, 0]
col1 = data[:, 1]
col2 = data[:, 2]
col3 = data[:, 3]

# Plot
plt.plot(time, col1, label="sz")
plt.plot(time, col2, label="sx")
plt.plot(time, col3, label="cxx")

# Labels and title
plt.xlabel("Time")
plt.ylabel("observables")
plt.title("Time Series Data")
plt.legend()

# Show the plot
plt.savefig('Part2/Plots_P2/observables.png', dpi=300)
plt.show()